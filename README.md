# structure_docking_validation

**Skill ID**: `structure_docking_validation`
**Capability**: `structure_docking_validation`
**Module type**: Execution + Critic (CellForge-Agent skill interface)

---

## Overview

This skill performs **structure-based molecular docking validation** on candidate target proteins identified from single-cell ligand-receptor (LR) analysis.

It assumes that the following upstream steps have been completed by other CellForge-Agent skills:

| Upstream Step | Expected Skill / Module |
|---|---|
| scRNA-seq QC, normalization, clustering | `scRNAseq_processing` or Seurat v5 |
| Cell type annotation | `celltypist_annotate` / `scMayoMap` |
| Cell–cell communication inference | `cellchat_lr` (CellChat v2) |
| Differential LR pair identification | `cellchat_lr` comparative mode |

`structure_docking_validation` takes the differential LR pairs → retrieves protein structures → prepares compounds → runs **AutoDock Vina** global docking → returns ranked binding affinities and poses.

---

## Input Requirements

### Required

| Input | Type | Description |
|---|---|---|
| `input_data` | `.h5ad` file path or `AnnData` object | Must contain completed CellChat LR analysis in `adata.uns['analysis_history']`. Alternatively, provide `target_proteins` explicitly in `params_dict`. |
| `params_dict['target_proteins']` | `list[str]` | UniProt IDs or gene symbols of receptor proteins to dock. Auto-extracted from CellChat LR results if not provided. |

### Optional

| Parameter | Default | Description |
|---|---|---|
| `compound_library` | `"fda_approved"` | Library source: `"fda_approved"` (built-in FDA drug set), `"custom_cas"` (resolve via PubChem), `"custom_dir"` (pre-processed PDBQT files). |
| `custom_cas_list` | `None` | Path to newline-separated CAS number list (when `compound_library="custom_cas"`). |
| `custom_compound_dir` | `None` | Path to directory with `.pdbqt` ligand files (when `compound_library="custom_dir"`). |
| `docking_exhaustiveness` | `8` | Vina exhaustiveness parameter. Higher = more thorough search (scales quadratically). |
| `docking_n_poses` | `10` | Maximum binding poses to record per ligand. |
| `drug_info_csv` | `None` | Optional BindingDB-style CSV for annotating compounds with known pharmacological evidence. |
| `output_dir` | `"docking_results"` | Directory for all output files. |
| `cleanup_intermediates` | `false` | If `true`, removes temporary protein/ligand prep files after completion. |

---

## Output

### Files (written to `output_dir/`)

| File | Description |
|---|---|
| `docking_ranked_results.tsv` | All binding poses ranked by affinity (columns: `ligand`, `target_uniprot`, `affinity_kcal_mol`, `ki_estimate_uM`). |
| `docking_summary.json` | Machine-readable summary with best binders per target and overall metrics. |
| `tmp_proteins/` | Downloaded + prepared receptor PDBQT files (deleted if `cleanup_intermediates=true`). |
| `tmp_ligands/` | Fetched + prepared ligand PDBQT files (deleted if `cleanup_intermediates=true`). |

### AnnData Updates

The returned `AnnData` object includes:

```
adata.uns['docking_results'] = {
    'n_targets_docked': int,
    'n_compounds_screened': int,
    'best_binding_affinity_kcal_mol': float,
    'ranked_results_path': str,
    'timestamp': str,
}
adata.uns['analysis_history'] += [{
    'skill_id': 'structure_docking_validation',
    'params': {...},
    'metrics': {
        'n_targets_docked': int,
        'n_compounds_screened': int,
        'best_binding_affinity_kcal_mol': float,
        'n_poses_generated': int,
        'top_binders_per_target': {uniprot_id: [{'compound', 'affinity'}, ...]},
    },
}]
```

---

## Dependencies

### Required system packages

| Package | Install command |
|---|---|
| **AutoDock Vina** | `conda install -c bioconda autodock-vina` |
| **OpenBabel** | `conda install -c conda-forge openbabel` |
| **RDKit** | `conda install -c conda-forge rdkit` |
| **curl** | (usually pre-installed on Linux/WSL) |

### Required Python packages

```
anndata>=0.8
numpy
pandas
```

### External resources accessed at runtime

- RCSB PDB REST API (`https://www.rcsb.org`)
- AlphaFold EBI (`https://alphafold.ebi.ac.uk`)
- PubChem REST API (`https://pubchem.ncbi.nlm.nih.gov`)

---

## Quality Thresholds (Critic)

The `critic.py` post-processor evaluates results using these thresholds:

| Metric | Threshold | Warning |
|---|---|---|
| `n_targets_docked` | > 0 | CRITICAL if 0 |
| `n_compounds_screened` | > 0 | CRITICAL if 0 |
| `best_binding_affinity_kcal_mol` | < 0 | UNLIKELY if > 0 |
| `best_binding_affinity_kcal_mol` | < -5 | May be weak if > -5 |
| `n_poses_generated` | ≥ 30% of expected | LOW YIELD warning |
| All top binders identical | — | DIVERSITY warning |

---

## Example Usage

### Via Python API

```python
from structure_docking_validation.run import run_structure_docking_validation
import anndata as ad

adata = ad.read("diabetic_nephropathy_cellchat.h5ad")

result = run_structure_docking_validation(
    input_data=adata,
    params_dict={
        "target_proteins": ["P62988", "Q9UNQ0"],   # UniProt IDs (PPIA, BSG from scDock paper)
        "compound_library": "fda_approved",
        "docking_exhaustiveness": 8,
        "output_dir": "docking_results",
    }
)
```

### Via CellForge-Agent executor

```json
{
  "skill_id": "structure_docking_validation",
  "input_data": "diabetic_nephropathy_cellchat.h5ad",
  "params_dict": {
    "target_proteins": null,
    "compound_library": "custom_cas",
    "custom_cas_list": "/path/to/cas_list.txt",
    "docking_exhaustiveness": 16,
    "output_dir": "docking_results"
  }
}
```

---

## Limitations

- Docking results are **computational estimates** of binding affinity and do not constitute direct evidence of functional drug–target interaction. Experimental validation is required.
- `cleanup_intermediates=true` deletes raw structure downloads — only keep if disk space is constrained and reprocessing is acceptable.
- CAS number resolution via PubChem may fail for non-standard or trade-name compounds; in such cases, use `custom_dir` mode with manually curated PDBQT files.
