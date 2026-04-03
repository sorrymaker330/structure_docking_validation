#!/usr/bin/env python3
"""
structure_docking_validation — critic / post-processing script.

Validates execution results and extracts quality metrics.
Called separately from the execution script.
"""
import warnings


def critic_post_process(adata, context=None):
    """
    Evaluate molecular docking results for quality and biological relevance.

    Returns
    -------
    dict with keys:
        metrics      : dict of extracted metrics
        warnings     : list of warning strings
        success      : bool
        context_used : str
    """
    context = context or {}
    metrics = {}
    warnings_list = []

    # ── Locate docking results in analysis history ──────────────────────────
    last_docking = None
    for entry in reversed(adata.uns.get("analysis_history", [])):
        if entry.get("skill_id") == "structure_docking_validation":
            last_docking = entry
            break

    if last_docking is None:
        warnings_list.append(
            "structure_docking_validation results not found in adata.uns['analysis_history']. "
            "Ensure run.py executed successfully."
        )
        return {
            "metrics": {},
            "warnings": warnings_list,
            "success": False,
            "context_used": "default",
        }

    dock_metrics = last_docking.get("metrics", {})
    metrics.update(dock_metrics)

    # ── Extract key values ────────────────────────────────────────────────────
    n_targets = dock_metrics.get("n_targets_docked", 0)
    n_compounds = dock_metrics.get("n_compounds_screened", 0)
    best_affinity = dock_metrics.get("best_binding_affinity_kcal_mol")
    n_poses = dock_metrics.get("n_poses_generated", 0)

    # ── Quality checks ───────────────────────────────────────────────────────

    # Check 1: at least one target successfully docked
    if n_targets == 0:
        warnings_list.append(
            "CRITICAL: No target proteins were successfully docked. "
            "Check UniProt ID validity and AlphaFold/PDB network access."
        )

    # Check 2: at least one compound in library
    if n_compounds == 0:
        warnings_list.append(
            "CRITICAL: No compounds were screened. "
            "Verify compound library (FDA-approved path, CAS list, or custom directory)."
        )

    # Check 3: reasonable binding affinity range
    if best_affinity is not None:
        if best_affinity > 0:
            warnings_list.append(
                f"UNLIKELY BEST AFFINITY: {best_affinity:.2f} kcal/mol is unfavorable. "
                "All docking poses may be unreliable — consider increasing exhaustiveness."
            )
        elif best_affinity < -15:
            warnings_list.append(
                f"VERY STRONG BINDING: {best_affinity:.2f} kcal/mol exceeds typical "
                "small-molecule binding range. May indicate a covalent or highly charged binder; "
                "verify compound identity and receptor preparation."
            )
        elif best_affinity > -5:
            warnings_list.append(
                f"WEAK BEST BINDING: {best_affinity:.2f} kcal/mol. "
                "Consider extending the docking search (higher exhaustiveness) "
                "or verifying the binding site location."
            )

    # Check 4: pose coverage
    if n_targets > 0 and n_poses == 0:
        warnings_list.append(
            "CRITICAL: No docking poses were generated. "
            "AutoDock Vina may have crashed — check Vina installation and receptor/ligand file formats."
        )

    expected_poses = n_targets * 10  # default n_poses = 10
    if n_poses < expected_poses * 0.3:
        warnings_list.append(
            f"LOW POSE YIELD: Only {n_poses} poses generated "
            f"(expected ~{expected_poses}). Some target-ligand pairs failed. "
            "Review warning messages from the docking run."
        )

    # Check 5: compound diversity
    top_binders = dock_metrics.get("top_binders_per_target", {})
    all_top_ligands = set()
    for ligands in top_binders.values():
        for item in ligands:
            all_top_ligands.add(item["compound"])
    if len(all_top_ligands) == 1 and n_targets > 1:
        warnings_list.append(
            "SAME LIGAND FOR ALL TARGETS: The top binder is identical across all targets. "
            "This may indicate a promiscuous binder or overly permissive binding site definition."
        )

    # ── Determine success ─────────────────────────────────────────────────────
    success = (
        n_targets > 0
        and n_compounds > 0
        and n_poses > 0
        and (best_affinity is not None and best_affinity < 0)
    )

    # ── Scientific context warnings ──────────────────────────────────────────
    if success and best_affinity is not None:
        if -12 < best_affinity < -8:
            metrics["significance"] = (
                "Moderate binding affinity — consistent with typical drug-target interactions. "
                "Experimental validation (e.g., surface plasmon resonance) is recommended."
            )
        elif best_affinity <= -12:
            metrics["significance"] = (
                "High binding affinity — strong candidate for experimental follow-up. "
                "Further selectivity profiling against off-target proteins is advised."
            )

    return {
        "metrics": metrics,
        "warnings": warnings_list,
        "success": success,
        "context_used": "default",
    }
