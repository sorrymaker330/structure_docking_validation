#!/usr/bin/env python3
import os
import sys
import json
import subprocess
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad

SKILL_DIR = os.path.dirname(os.path.abspath(__file__))

def _check_vina():
    try:
        res = subprocess.run(["vina", "--version"], capture_output=True, text=True, timeout=10)
        return res.stdout.strip() or res.stderr.strip().split("\n")[0]
    except Exception:
        raise OSError("AutoDock Vina not found in PATH.")

def run_structure_docking_validation(input_data, params_dict=None, default_params=None, output_dir=None):
    adata = input_data if isinstance(input_data, ad.AnnData) else ad.read(input_data)
    params = {**{"compound_library": "fda_approved", "docking_exhaustiveness": 8, "docking_n_poses": 10, "output_dir": "docking_results"}, **(default_params or {}), **(params_dict or {})}
    out_dir = Path(output_dir or params["output_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    vina_version = _check_vina()

    target_list = params.get("target_proteins") or []
    target_list = list(dict.fromkeys(target_list))[:20]

    ligand_tmp = out_dir / "tmp_ligands"
    ligand_tmp.mkdir(exist_ok=True)
    compound_paths = []

    if params["compound_library"] == "fda_approved":
        fda_dir = Path(SKILL_DIR) / "data" / "fda_approved"
        if fda_dir.exists():
            compound_paths = list(fda_dir.glob("*.pdbqt"))
    elif params["compound_library"] == "custom_cas":
        cas_list = [cas.strip() for cas in (params.get("custom_cas_list") or []) if cas and not cas.startswith("#")]
        if cas_list:
            cas_txt_path = ligand_tmp / "cas_list.txt"
            cas_txt_path.write_text("\n".join(cas_list))
            print("[INFO] 调用同目录下的 download_cas_pubchem.py 获取配体...")
            script_path = os.path.join(SKILL_DIR, "download_cas_pubchem.py")
            subprocess.run(["python3", script_path, str(cas_txt_path), str(ligand_tmp)], check=True)
            compound_paths = list(ligand_tmp.glob("*.pdbqt"))

    if not compound_paths: 
        raise ValueError("No compounds available for docking.")

    protein_tmp = out_dir / "tmp_proteins"
    protein_tmp.mkdir(exist_ok=True)
    prepared_receptors = {} 

    for uid in target_list:
        af_id = f"AF-{uid}-F1"
        uid_dir = protein_tmp / uid
        uid_dir.mkdir(exist_ok=True)
        
        print(f"\n[INFO] 调用 download_alphafold.py 下载 {af_id}...")
        dl_script = os.path.join(SKILL_DIR, "download_alphafold.py")
        subprocess.run(["python3", dl_script, af_id, str(uid_dir)])
        
        pdb_file = uid_dir / f"{af_id}.pdb"
        if not pdb_file.exists():
            warnings.warn(f"Failed to download PDB for {uid}")
            continue
            
        print(f"[INFO] 调用 prepare_receptor.py 处理 {af_id}...")
        prep_script = os.path.join(SKILL_DIR, "prepare_receptor.py")
        # 直接调用你目录下的脚本，不再使用命令行的 -c 绕路
        subprocess.run(["python3", prep_script, str(pdb_file), str(uid_dir)])
        
        pdbqt_file = uid_dir / f"{af_id}_prepared.pdbqt"
        grid_file = uid_dir / f"{af_id}_prepared_grid.txt"
        
        if pdbqt_file.exists() and grid_file.exists():
            grid = {}
            for line in grid_file.read_text().splitlines():
                if "=" in line:
                    k, v = line.split("=")
                    grid[k.strip()] = v.strip()
            prepared_receptors[uid] = (pdbqt_file, grid)

    if not prepared_receptors: 
        raise ValueError("No protein structures could be retrieved or prepared.")

    all_results = []
    for uid, (rec_pdbqt, grid) in prepared_receptors.items():
        for lig_pdbqt in compound_paths:
            vina_out = out_dir / f"vina_out_{uid}_{lig_pdbqt.stem}.pdbqt"
            log_file = out_dir / f"vina_log_{lig_pdbqt.stem}.txt"
            
            cmd = [
                "vina", "--receptor", str(rec_pdbqt), "--ligand", str(lig_pdbqt),
                "--center_x", grid.get("center_x", "0"), "--center_y", grid.get("center_y", "0"), "--center_z", grid.get("center_z", "0"),
                "--size_x", grid.get("size_x", "30"), "--size_y", grid.get("size_y", "30"), "--size_z", grid.get("size_z", "30"),
                "--exhaustiveness", str(params["docking_exhaustiveness"]), "--num_modes", str(params["docking_n_poses"]),
                "--out", str(vina_out), "--log", str(log_file)
            ]
            print(f"[Vina] Docking {lig_pdbqt.stem} with {uid}...")
            subprocess.run(cmd, capture_output=True)
            
            records = []
            if log_file.exists():
                for line in log_file.read_text().split("\n"):
                    if line.strip().startswith("---") or "mode" in line.lower(): continue
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            aff = float(parts[0])
                            records.append({"ligand": lig_pdbqt.stem, "affinity_kcal_mol": aff, "ki_estimate_uM": round(np.exp(aff / (0.001987 * 298)) * 1e6, 3)})
                        except ValueError: pass
            if not records: 
                records.append({"ligand": lig_pdbqt.stem, "affinity_kcal_mol": float("nan"), "ki_estimate_uM": float("nan")})
            
            df = pd.DataFrame(records)
            df["target_uniprot"] = uid
            all_results.append(df)

    if not all_results: 
        raise RuntimeError("All docking runs failed.")

    results_df = pd.concat(all_results, ignore_index=True).sort_values("affinity_kcal_mol")
    ranked = results_df.groupby("target_uniprot", group_keys=False).apply(lambda g: g.sort_values("affinity_kcal_mol").head(params["docking_n_poses"]))
    
    docking_summary = {"n_targets_docked": len(prepared_receptors), "n_compounds_screened": len(compound_paths), "top_binders": {}}
    for t, grp in ranked.groupby("target_uniprot"):
        docking_summary["top_binders"][t] = [{"compound": row["ligand"], "affinity": float(row["affinity_kcal_mol"])} for _, row in grp.head(5).iterrows()]

    best = ranked.iloc[0] if not ranked.empty else None
    docking_summary["best_binding_affinity_kcal_per_mol"] = float(best["affinity_kcal_mol"]) if best is not None else None
    docking_summary["n_poses_generated"] = len(ranked)

    ranked.to_csv(out_dir / "docking_ranked_results.tsv", sep="\t", index=False)
    (out_dir / "docking_summary.json").write_text(json.dumps(docking_summary, indent=2))
    
    adata_docked = adata.copy()
    adata_docked.uns["docking_results"] = docking_summary
    return adata_docked

if __name__ == "__main__":
    result = run_structure_docking_validation(sys.argv[1], {"target_proteins": ["P17742"], "custom_cas_list": ["93479-97-1"], "compound_library": "custom_cas"}, output_dir=sys.argv[2] if len(sys.argv)>2 else "docking_results")
    print(f"\nDocking complete. Best affinity: {result.uns['docking_results'].get('best_binding_affinity_kcal_per_mol')} kcal/mol")