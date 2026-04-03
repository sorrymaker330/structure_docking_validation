#!/usr/bin/env python3
"""
structure_docking_validation — execution script.

Entry point: run_structure_docking_validation(input_data, params_dict, default_params, output_dir)

Called by SkillExecutor via exec().
Assumes: single-cell preprocessing (QC, normalization, clustering, cell annotation)
         and CellChat-based LR analysis have been completed externally.
         Target proteins are available in adata.uns['analysis_history'] or
         provided directly via params_dict['target_proteins'].

This module bridges LR-pair candidate targets → structure retrieval →
compound preparation → AutoDock Vina global docking → ranked results.
"""
import os
import sys
import json
import subprocess
import datetime
import warnings
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad


# ─────────────────────────────────────────────────────────────────────────────
# Helper: AutoDock Vina availability check
# ─────────────────────────────────────────────────────────────────────────────
def _check_vina():
    try:
        result = subprocess.run(
            ["vina", "--version"],
            capture_output=True, text=True, timeout=10
        )
        version = result.stdout.strip()
        if not version:
            version = result.stderr.strip().split("\n")[0]
        return version
    except FileNotFoundError:
        raise OSError(
            "AutoDock Vina not found in PATH. "
            "Install: conda install -c bioconda autodock-vina"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Helper: Fetch protein structure (PDB > AlphaFold)
# ─────────────────────────────────────────────────────────────────────────────
def _fetch_protein_structure(uniprot_id: str, out_dir: Path, timeout: int = 60) -> Path:
    """
    Retrieve protein structure for a UniProt ID.
    Priority: PDB (experimental) > AlphaFold (predicted).
    Saves as <uniprot_id>.pdb in out_dir.
    Returns Path to the downloaded file.
    """
    out_path = out_dir / f"{uniprot_id}.pdb"

    # Try PDB first via RCSB REST API
    pdb_api = (
        f"https://www.rcsb.org/fasta/entry/{uniprot_id}/download"
    )
    try:
        result = subprocess.run(
            ["curl", "-s", "--max-time", str(timeout), pdb_api],
            capture_output=True, text=True, timeout=timeout + 5
        )
        if result.returncode == 0 and result.stdout.strip():
            lines = result.stdout.strip().split("\n")
            # Find PDB lines (sequence lines are >PDBID:CHAIN|...)
            pdb_lines = [l for l in lines if not l.startswith(">")]
            if pdb_lines:
                pdb_content = "\n".join(pdb_lines[:5000])  # safety cap
                out_path.write_text(pdb_content)
                return out_path
    except Exception:
        pass

    # Fall back to AlphaFold
    alphafold_url = (
        f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    )
    try:
        result = subprocess.run(
            ["curl", "-s", "--max-time", str(timeout), "-o", str(out_path),
             alphafold_url],
            capture_output=True, timeout=timeout + 5
        )
        if result.returncode == 0 and out_path.exists() and out_path.stat().st_size > 100:
            return out_path
        else:
            out_path.unlink(missing_ok=True)
    except Exception:
        pass

    raise FileNotFoundError(
        f"Could not retrieve structure for UniProt ID: {uniprot_id} "
        "(tried PDB and AlphaFold)"
    )


# ─────────────────────────────────────────────────────────────────────────────
# Helper: Fetch compound from PubChem by CAS
# ─────────────────────────────────────────────────────────────────────────────
def _fetch_compound_by_cas(cas: str, out_dir: Path, timeout: int = 30) -> Path:
    """
    Resolve CAS number to SMILES via PubChem REST, then convert to
    PDBQT using OpenBabel / RDKit.
    Returns Path to the PDBQT file.
    """
    # Step 1: find PubChem CID from CAS
    try:
        cid_result = subprocess.run(
            ["curl", "-s", "--max-time", str(timeout),
             f"https://pubchem.ncbi.nlm.nih.gov/rest/pg/tabview/json/CompoundList/CAS/{cas}/Property/MolecularFormula,MolecularWeight,CanonicalSMILES/JSON"],
            capture_output=True, text=True, timeout=timeout + 5
        )
        if cid_result.returncode == 0:
            import json as _json
            data = _json.loads(cid_result.stdout)
            if data.get("Table", {}).get("Row"):
                cid = data["Table"]["Row"][0]["Cell"][0]
                smiles = data["Table"]["Row"][0]["Cell"][-1]
            else:
                # Fallback: search by name
                cid = None
                smiles = None
        else:
            cid = None
            smiles = None
    except Exception:
        cid = None
        smiles = None

    if cid:
        try:
            smiles_result = subprocess.run(
                ["curl", "-s", "--max-time", str(timeout),
                 f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON"],
                capture_output=True, text=True, timeout=timeout + 5
            )
            if smiles_result.returncode == 0:
                import json as _json
                smiles_data = _json.loads(smiles_result.stdout)
                smiles = (
                    smiles_data
                    .get("PropertyTable", {})
                    .get("Properties", [{}])[0]
                    .get("CanonicalSMILES", "")
                )
        except Exception:
            pass

    if not smiles:
        # Try RDKit's CAS resolver if PubChem fails
        try:
            from rdkit import Chem
            from rdkit.Chem import Draw
            # RDKit cannot resolve CAS directly; this is a last-resort warning
            warnings.warn(f"Could not resolve SMILES for CAS {cas}")
        except ImportError:
            pass
        raise FileNotFoundError(f"Could not fetch compound for CAS: {cas}")

    # Step 2: write SMILES to temp SDF
    sdf_path = out_dir / f"{cas}.sdf"
    pdbqt_path = out_dir / f"{cas}.pdbqt"

    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        with open(sdf_path, "w") as f:
            f.write(smiles + "\n")
            f.write("  SciTegic111\n\n")
            f.write(f"{cas}\n")
        # Convert SDF → PDBQT via OpenBabel
        subprocess.run(
            ["obabel", "-isdf", str(sdf_path), "-opdbqt", "-O", str(pdbqt_path)],
            capture_output=True, check=True, timeout=timeout
        )
        return pdbqt_path
    except Exception as e:
        raise FileNotFoundError(
            f"Failed to prepare compound for CAS {cas}: {e}"
        ) from e


# ─────────────────────────────────────────────────────────────────────────────
# Helper: Prepare protein PDBQT using AutoDockTools
# ─────────────────────────────────────────────────────────────────────────────
def _prepare_protein(protein_pdb: Path, out_dir: Path) -> Path:
    """
    Use AutoDockTools (prepare_receptor4.py) to add hydrogens,
    assign Gasteiger charges, and convert receptor to PDBQT.
    Returns Path to the PDBQT file.
    """
    receptor_pdbqt = out_dir / (protein_pdb.stem + ".pdbqt")

    # AutoDockTools prepare_receptor4.py
    cmd = [
        "python3", "-c",
        f"""
import sys
sys.path.insert(0, '')
from AutoDockTools.Utilities import prepare_receptor4
prepare_receptor4.parse_args = lambda: None
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-r', '--receptor')
parser.add_argument('-o', '--output')
parser.add_argument('-A', '--add_output_keyword_apolar')
parser.add_argument('-v', '--verbose', action='store_true')
args = parser.parse_args(['-r', '{protein_pdb}', '-o', '{receptor_pdbqt}', '-A', 'yes'])
prepare_receptor4.run(args)
"""
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=120)
    except Exception:
        # Fallback: basic PDBQT conversion using obabel
        try:
            subprocess.run(
                ["obabel", "-ipdb", str(protein_pdb), "-opdbqt",
                 "-O", str(receptor_pdbqt)],
                capture_output=True, check=True, timeout=60
            )
        except Exception as e:
            raise OSError(f"Protein preparation failed: {e}") from e

    if not receptor_pdbqt.exists():
        raise OSError(f"Protein PDBQT not generated: {receptor_pdbqt}")
    return receptor_pdbqt


# ─────────────────────────────────────────────────────────────────────────────
# Helper: Run AutoDock Vina global docking
# ─────────────────────────────────────────────────────────────────────────────
def _run_vina_docking(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    out_dir: Path,
    exhaustiveness: int = 8,
    n_poses: int = 10
) -> pd.DataFrame:
    """
    Run AutoDock Vina global docking for one receptor–ligand pair.
    Returns DataFrame with columns: ligand, affinity_kcal_mol, ki_estimate.
    """
    vina_out = out_dir / f"vina_out_{receptor_pdbqt.stem}_{ligand_pdbqt.stem}.pdbqt"

    # Determine docking box automatically (center = protein COM, size = 30Å)
    cmd_center = [
        "python3", "-c",
        f"""
import rdkit.Chem as Chem
mol = Chem.MolFromPDBFile('{receptor_pdbqt}', removeHs=False)
if mol is None:
    print('22.5,22.5,22.5')
else:
    coords = mol.GetConformer(0).GetPositions()
    cx, cy, cz = coords.mean(axis=0)
    print(f'{{cx:.3f}},{{cy:.3f}},{{cz:.3f}}')
"""
    ]
    try:
        cr = subprocess.run(cmd_center, capture_output=True, text=True, timeout=30)
        center_str = cr.stdout.strip() if cr.stdout.strip() else "22.5,22.5,22.5"
    except Exception:
        center_str = "22.5,22.5,22.5"

    cx, cy, cz = [float(x) for x in center_str.split(",")]

    vina_cmd = [
        "vina",
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x", "30", "--size_y", "30", "--size_z", "30",
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(n_poses),
        "--out", str(vina_out),
        "--log", str(out_dir / f"vina_log_{ligand_pdbqt.stem}.txt"),
    ]

    result = subprocess.run(
        vina_cmd, capture_output=True, text=True, timeout=300
    )

    # Parse binding affinities from log
    records = []
    log_file = out_dir / f"vina_log_{ligand_pdbqt.stem}.txt"
    if log_file.exists():
        content = log_file.read_text()
        for line in content.split("\n"):
            if line.strip().startswith("---") or "mode" in line.lower():
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    aff = float(parts[0])
                    records.append({
                        "ligand": ligand_pdbqt.stem,
                        "affinity_kcal_mol": aff,
                        # rough estimate: Ki = exp(aff/(0.001987*298))
                        "ki_estimate_uM": round(
                            np.exp(aff / (0.001987 * 298)) * 1e6, 3
                        )
                    })
                except ValueError:
                    continue

    if not records:
        warnings.warn(f"No docking poses retrieved for {ligand_pdbqt.stem}")
        records.append({
            "ligand": ligand_pdbqt.stem,
            "affinity_kcal_mol": float("nan"),
            "ki_estimate_uM": float("nan")
        })

    return pd.DataFrame(records)


# ─────────────────────────────────────────────────────────────────────────────
# Helper: Load built-in FDA-approved compound library
# ─────────────────────────────────────────────────────────────────────────────
def _get_fda_library(out_dir: Path) -> list:
    """
    Returns list of PDBQT Paths for built-in FDA-approved compounds.
    If library not found locally, returns empty list.
    """
    # Default expected location within the pipeline's data dir
    lib_base = Path(__file__).parent / "data" / "fda_approved"
    if not lib_base.exists():
        # Try global install
        lib_base = Path("/mnt/c/Users/zzt/Desktop/structure_docking_validation1/data/fda_approved")
    if lib_base.exists():
        return list(lib_base.glob("*.pdbqt"))
    return []


# ─────────────────────────────────────────────────────────────────────────────
# Main entry point
# ─────────────────────────────────────────────────────────────────────────────
def run_structure_docking_validation(
    input_data,
    params_dict=None,
    default_params=None,
    output_dir=None
):
    """
    Execute structure-based molecular docking validation.

    Parameters
    ----------
    input_data : str or AnnData
        Path to AnnData (.h5ad) or AnnData object.
        Must contain completed CellChat LR analysis in adata.uns['analysis_history'].
        If providing pre-computed target proteins, pass via params_dict['target_proteins'].
    params_dict : dict, optional
        Overrides for the following keys:
        target_proteins, compound_library, custom_cas_list, custom_compound_dir,
        protein_source, docking_exhaustiveness, docking_n_poses,
        drug_info_csv, output_dir, cleanup_intermediates.
    default_params : dict, optional
        Full parameter block (normally injected by executor).
    output_dir : str or Path, optional
        Explicit output directory (overrides params_dict['output_dir']).

    Returns
    -------
    AnnData with docking results stored in obs and uns.
    """

    # ── 1. Resolve input AnnData ────────────────────────────────────────────
    adata = input_data if isinstance(input_data, ad.AnnData) else ad.read(input_data)

    # ── 2. Resolve parameters ─────────────────────────────────────────────
    default_params = default_params or {
        "target_proteins": None,
        "compound_library": "fda_approved",
        "custom_cas_list": None,
        "custom_compound_dir": None,
        "protein_source": "pdb_alphafold_hybrid",
        "docking_exhaustiveness": 8,
        "docking_n_poses": 10,
        "drug_info_csv": None,
        "output_dir": "docking_results",
        "cleanup_intermediates": False,
    }
    agent_params = params_dict or {}
    params = {**default_params, **agent_params}

    out_dir = Path(output_dir or params["output_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── 3. Check Vina availability ─────────────────────────────────────────
    vina_version = _check_vina()

    # ── 4. Determine target proteins ────────────────────────────────────────
    target_list = params["target_proteins"]
    if not target_list:
        # Auto-extract from CellChat results in analysis history
        last_cellchat = None
        for entry in reversed(adata.uns.get("analysis_history", [])):
            if "cellchat" in entry.get("skill_id", "").lower() or \
               entry.get("skill_id") == "cellchat_lr":
                last_cellchat = entry
                break
        if last_cellchat is None:
            raise ValueError(
                "No CellChat LR analysis found in adata.uns['analysis_history']. "
                "Provide target_proteins explicitly via params_dict['target_proteins']."
            )
        metrics = last_cellchat.get("metrics", {})
        diff_lr = metrics.get("differential_lr_pairs", [])
        if not diff_lr:
            raise ValueError(
                "CellChat entry found but no differential LR pairs recorded. "
                "Provide target_proteins list manually."
            )
        # Extract receptor UniProt IDs from LR pairs
        target_list = [lr.get("receptor_uniprot") or lr.get("receptor_gene_symbol")
                       for lr in diff_lr if lr]
        target_list = [t for t in target_list if t]
        if not target_list:
            raise ValueError("Could not extract target proteins from LR pairs.")

    # Deduplicate and limit
    target_list = list(dict.fromkeys(target_list))[:20]  # cap at 20 targets

    # ── 5. Prepare compound library ─────────────────────────────────────────
    library_type = params["compound_library"]
    compound_paths = []

    if library_type == "fda_approved":
        compound_paths = _get_fda_library(out_dir / "tmp_ligands")
        if not compound_paths:
            warnings.warn(
                "FDA-approved compound library not found. "
                "Using empty library — no docking will be performed."
            )
    elif library_type == "custom_cas":
        cas_list_file = params.get("custom_cas_list")
        if isinstance(cas_list_file, str):
            cas_list = Path(cas_list_file).read_text().strip().split("\n")
        else:
            cas_list = cas_list_file or []
        ligand_tmp = out_dir / "tmp_ligands"
        ligand_tmp.mkdir(exist_ok=True)
        for cas in cas_list:
            cas = cas.strip()
            if not cas or cas.startswith("#"):
                continue
            try:
                pdbqt = _fetch_compound_by_cas(cas, ligand_tmp)
                compound_paths.append(pdbqt)
            except Exception as e:
                warnings.warn(f"Skipping CAS {cas}: {e}")
    elif library_type == "custom_dir":
        custom_dir = Path(params.get("custom_compound_dir", ""))
        if custom_dir.exists():
            compound_paths = list(custom_dir.glob("*.pdbqt"))
        else:
            raise FileNotFoundError(f"Custom compound directory not found: {custom_dir}")

    if not compound_paths:
        raise ValueError("No compounds available for docking. Check library configuration.")

    # ── 6. Process drug annotations if provided ───────────────────────────
    drug_annotations = {}
    drug_csv = params.get("drug_info_csv")
    if drug_csv and Path(drug_csv).exists():
        try:
            drug_df = pd.read_csv(drug_csv)
            for _, row in drug_df.iterrows():
                drug_annotations[row.get("compound_name", "")] = row.to_dict()
        except Exception as e:
            warnings.warn(f"Could not parse drug annotations CSV: {e}")

    # ── 7. Structure retrieval & protein preparation ─────────────────────────
    protein_tmp = out_dir / "tmp_proteins"
    protein_tmp.mkdir(exist_ok=True)

    prepared_proteins = {}
    for uniport_id in target_list:
        try:
            pdb_path = _fetch_protein_structure(uniport_id, protein_tmp)
            pdbqt_path = _prepare_protein(pdb_path, protein_tmp)
            prepared_proteins[uniport_id] = pdbqt_path
        except Exception as e:
            warnings.warn(f"Skipping target {uniport_id}: {e}")

    if not prepared_proteins:
        raise ValueError(
            "No protein structures could be retrieved or prepared. "
            "Check UniProt IDs and network access."
        )

    # ── 8. Docking loop ─────────────────────────────────────────────────────
    all_results = []
    docking_summary = {
        "n_targets_docked": len(prepared_proteins),
        "n_compounds_screened": len(compound_paths),
        "best_affinity_per_target": {},
        "top_binders": {},
        " vina_version": vina_version,
    }

    ligand_tmp = out_dir / "tmp_ligands"
    ligand_tmp.mkdir(exist_ok=True)

    for uniprot_id, receptor_pdbqt in prepared_proteins.items():
        for ligand_pdbqt in compound_paths:
            try:
                df = _run_vina_docking(
                    receptor_pdbqt,
                    ligand_pdbqt,
                    out_dir,
                    exhaustiveness=params["docking_exhaustiveness"],
                    n_poses=params["docking_n_poses"]
                )
                df["target_uniprot"] = uniprot_id
                df["target_gene"] = uniprot_id  # may be refined with gene symbol lookup
                if drug_annotations:
                    lig_name = ligand_pdbqt.stem
                    df["drug_annotation"] = str(drug_annotations.get(lig_name, {}))
                all_results.append(df)
            except Exception as e:
                warnings.warn(f"Docking failed for {receptor_pdbqt.stem} × {ligand_pdbqt.stem}: {e}")

    if not all_results:
        raise RuntimeError("All docking runs failed. Check Vina installation and input files.")

    results_df = pd.concat(all_results, ignore_index=True)

    # ── 9. Rank results per target ──────────────────────────────────────────
    results_df = results_df.sort_values("affinity_kcal_mol", ascending=True)
    ranked = (
        results_df.groupby("target_uniprot", group_keys=False)
        .apply(lambda g: g.sort_values("affinity_kcal_mol").head(params["docking_n_poses"]))
    )

    # Best binders summary
    for t, grp in ranked.groupby("target_uniprot"):
        best = grp.iloc[0]
        docking_summary["best_affinity_per_target"][t] = float(best["affinity_kcal_mol"])
        docking_summary["top_binders"][t] = [
            {"compound": row["ligand"], "affinity": float(row["affinity_kcal_mol"])}
            for _, row in grp.head(5).iterrows()
        ]

    best_global = ranked.iloc[0] if not ranked.empty else None
    docking_summary["best_binding_affinity_kcal_per_mol"] = (
        float(best_global["affinity_kcal_mol"]) if best_global is not None else None
    )
    docking_summary["n_poses_generated"] = len(ranked)

    # ── 10. Save outputs ─────────────────────────────────────────────────────
    results_tsv = out_dir / "docking_ranked_results.tsv"
    ranked.to_csv(results_tsv, sep="\t", index=False)

    summary_json = out_dir / "docking_summary.json"
    summary_json.write_text(json.dumps(docking_summary, indent=2))

    # ── 11. Write back to AnnData ────────────────────────────────────────────
    adata_docked = adata.copy()
    adata_docked.uns["docking_results"] = {
        "params": params,
        " vina_version": vina_version,
        "n_targets_docked": docking_summary["n_targets_docked"],
        "n_compounds_screened": docking_summary["n_compounds_screened"],
        "best_binding_affinity_kcal_per_mol": docking_summary["best_binding_affinity_kcal_mol"],
        "ranked_results_path": str(results_tsv),
        "timestamp": datetime.datetime.now().isoformat(),
    }
    adata_docked.uns["analysis_history"].append({
        "skill_id": "structure_docking_validation",
        "params": {k: v for k, v in params.items()
                   if k not in ("custom_cas_list",)},
        "metrics": {
            "n_targets_docked": docking_summary["n_targets_docked"],
            "n_compounds_screened": docking_summary["n_compounds_screened"],
            "best_binding_affinity_kcal_mol": docking_summary["best_binding_affinity_kcal_mol"],
            "n_poses_generated": docking_summary["n_poses_generated"],
            "top_binders_per_target": docking_summary["top_binders"],
        },
        "timestamp": datetime.datetime.now().isoformat(),
    })

    # ── 12. Cleanup intermediate files ───────────────────────────────────────
    if params.get("cleanup_intermediates", False):
        shutil.rmtree(protein_tmp, ignore_errors=True)
        shutil.rmtree(ligand_tmp, ignore_errors=True)

    return adata_docked


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python run.py <input.h5ad> [output_dir]")
        sys.exit(1)
    in_file = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) > 2 else "docking_results"
    result = run_structure_docking_validation(in_file, output_dir=out_dir)
    print(f"Docking complete. Results saved to {out_dir}/")
    print(f"  Targets docked: {result.uns['docking_results']['n_targets_docked']}")
    print(f"  Best affinity:  {result.uns['docking_results']['best_binding_affinity_kcal_mol']} kcal/mol")
