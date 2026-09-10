#!/usr/bin/env python3
"""Run locked Phase 3A-1 donor-aware broad retinal cellular architecture.

The script consumes the finalized R8 MAGMA/scDRS gene sets and the full exact
official-ID healthy reference.  Primary inference is donor x broad class;
cell-level scDRS group results and MAGMA gene-property tests are secondary or
orthogonal support.  No fine-state analysis is performed here.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import os
import shlex
import subprocess
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import cosine, spearmanr, t, shapiro

# scDRS 1.0.2 still references np.float_, which was removed in NumPy 2.0.
# Keep this compatibility shim local to the analysis process; do not alter the
# locked scDRS environment.
if not hasattr(np, "float_"):
    np.float_ = np.float64


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "phase3a1"
LOGS = OUT / "logs"
FIGURES = ROOT / "figures" / "phase3a1"
H5AD = ROOT / "data" / "single_cell" / "HRCA" / "HRCA_EXACT_MAPPED_HEALTHY_CORE.h5ad"
GS = ROOT / "results" / "phase3ar8" / "retinaplex_traits.gs"
MAGMA_SCORES = ROOT / "results" / "phase3ar8" / "magma_gene_scores.tsv"
FINE_READY = OUT / "fine_state_readiness.tsv"
LOCK = ROOT / "config" / "phase3a1_analysis_lock.yaml"
SCDRS_PYTHON = Path(os.environ.get("SCDRS_PYTHON", "python"))
MAGMA = Path(os.environ.get("MAGMA_EXECUTABLE", "magma"))
TRAITS = ["AMD", "POAG", "RE"]
CLASSES = [
    "Muller glia",
    "RGC",
    "amacrine cells",
    "astrocytes",
    "bipolar cells",
    "cones",
    "horizontal cells",
    "microglia",
    "rods",
]
CLASS_VAR = {x: x.replace(" ", "_") for x in CLASSES}
EXPECTED_OVERLAP = {"AMD": 962, "POAG": 937, "RE": 955}
BASE_SEED = 20260905
N_BOOT = 5000


def bh(values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p values, retaining NaN as NaN."""
    values = np.asarray(values, dtype=float)
    result = np.full(values.shape, np.nan, dtype=float)
    finite = np.isfinite(values)
    if not finite.any():
        return result
    p = values[finite]
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty_like(q)
    out[order] = q
    result[finite] = out
    return result


def finite_fraction(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.isfinite(values).mean()) if len(values) else math.nan


def run_command(command: list[str], log_path: Path) -> subprocess.CompletedProcess[str]:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as handle:
        handle.write("COMMAND\t" + " ".join(shlex.quote(x) for x in command) + "\n\n")
        result = subprocess.run(command, cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        handle.write(result.stdout or "")
        handle.write(f"\nEXIT_CODE\t{result.returncode}\n")
    return result


def require_lock() -> None:
    if not LOCK.exists():
        raise RuntimeError(f"locked Phase 3A-1 configuration missing: {LOCK}")
    text = LOCK.read_text(encoding="utf-8")
    for token in ("status: LOCKED", "HRCA_EXACT_MAPPED_HEALTHY_CORE", "donor_class_minimum_cells: 30", "primary_controls: 1000"):
        if token not in text:
            raise RuntimeError(f"Phase 3A-1 lock does not contain required token: {token}")
    digest = hashlib.sha256(H5AD.read_bytes()).hexdigest()
    if f"h5ad_sha256: {digest}" not in text:
        raise RuntimeError("reference H5AD checksum does not match the Phase 3A-1 lock")


def ensure_dirs() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)


def load_normalized_reference():
    import scdrs

    adata = scdrs.util.load_h5ad(str(H5AD), flag_filter_data=False, flag_raw_count=True)
    if not sparse.issparse(adata.X):
        raise RuntimeError("Phase 3A-1 requires a sparse expression matrix after scDRS loading")
    if adata.n_obs != 98834 or adata.n_vars != 36601:
        raise RuntimeError(f"unexpected Phase 3A-1 reference shape: {adata.shape}")
    required = {"donor_id", "sample", "harmonized_broad_class", "total_counts", "n_genes", "pct_counts_mt"}
    missing = required - set(adata.obs.columns)
    if missing:
        raise RuntimeError(f"reference missing technical/annotation columns: {sorted(missing)}")
    if set(adata.obs["harmonized_broad_class"].astype(str)) != set(CLASSES):
        raise RuntimeError("reference broad-class set differs from the locked nine classes")
    return adata


def build_covariates(adata) -> pd.DataFrame:
    donor = adata.obs["donor_id"].astype(str)
    cov = pd.get_dummies(donor, prefix="donor", drop_first=True, dtype=float)
    cov.index = adata.obs_names
    n_genes = pd.to_numeric(adata.obs["n_genes"], errors="coerce").to_numpy(dtype=float)
    sd = np.nanstd(n_genes)
    if not np.isfinite(sd) or sd == 0:
        raise RuntimeError("n_genes covariate has no finite variation")
    cov["n_genes"] = (n_genes - np.nanmean(n_genes)) / sd
    cov.to_csv(OUT / "scdrs_covariates.tsv", sep="\t", index_label="CELL")
    return cov


def load_gene_sets(adata) -> tuple[dict[str, tuple[list[str], np.ndarray]], dict[str, int]]:
    import scdrs

    dict_gs = scdrs.util.load_gs(str(GS), to_intersect=list(adata.var_names))
    if set(dict_gs) != set(TRAITS):
        raise RuntimeError(f"unexpected scDRS traits: {sorted(dict_gs)}")
    overlaps: dict[str, int] = {}
    audit = pd.read_csv(ROOT / "results" / "phase3ar8" / "gene_overlap_audit.tsv", sep="\t")
    for trait in TRAITS:
        genes, weights = dict_gs[trait]
        overlaps[trait] = len(genes)
        expected = EXPECTED_OVERLAP[trait]
        recorded = int(audit.loc[audit["trait"] == trait, "matched_genes"].iloc[0])
        if len(genes) != expected or recorded != expected:
            raise RuntimeError(f"{trait} gene overlap changed: loaded={len(genes)} recorded={recorded} expected={expected}")
        if len(genes) < 100:
            raise RuntimeError(f"{trait} matched scDRS gene set is too small")
    return dict_gs, overlaps


def score_traits(adata, cov: pd.DataFrame, dict_gs, overlaps: dict[str, int]):
    import scdrs

    # Reuse the completed, QC-passed score files after an interrupted run. The
    # scores are locked inputs to the downstream donor/group/MAGMA stages, so
    # this avoids recomputing 3 x 1000 control sets while preserving validation.
    score_paths = {trait: OUT / f"{trait}.score.gz" for trait in TRAITS}
    full_paths = {trait: OUT / f"{trait}.full_score.gz" for trait in TRAITS}
    existing_qc = OUT / "scdrs_score_qc.tsv"
    if existing_qc.exists() and all(path.exists() for path in [*score_paths.values(), *full_paths.values()]):
        qc = pd.read_csv(existing_qc, sep="\t")
        if set(qc.get("trait", [])) == set(TRAITS) and (qc["status"] == "PASS").all():
            donor_delta_by_trait: dict[str, dict[str, dict[str, float]]] = {}
            technical_rows: list[dict[str, object]] = []
            for trait in TRAITS:
                df = pd.read_csv(score_paths[trait], sep="\t", index_col=0)
                if len(df) != adata.n_obs or not np.isfinite(df["norm_score"].to_numpy(dtype=float)).all():
                    break
                donor_delta_by_trait[trait] = compute_donor_deltas(
                    adata.obs, df["norm_score"].to_numpy(dtype=float), trait
                )
                technical_rows.extend(
                    technical_covariate_rows(adata.obs, df["norm_score"].to_numpy(dtype=float), trait)
                )
            else:
                pd.DataFrame(technical_rows).to_csv(
                    OUT / "scdrs_technical_covariate_qc.tsv", sep="\t", index=False
                )
                # The cached score files do not carry the in-memory normalized
                # expression matrix. Reapply the locked one-time preprocessing
                # so the downstream neighbor graph sees the same input scale
                # expected by scDRS group analysis.
                with np.errstate(all="ignore"):
                    scdrs.pp.preprocess(adata, cov=cov, n_mean_bin=20, n_var_bin=20, copy=False)
                return qc, score_paths, full_paths, donor_delta_by_trait

    with warnings.catch_warnings(record=True) as caught_pre:
        warnings.simplefilter("always")
        with np.errstate(all="ignore"):
            scdrs.pp.preprocess(adata, cov=cov, n_mean_bin=20, n_var_bin=20, copy=False)
    gene_stats = adata.uns["SCDRS_PARAM"]["GENE_STATS"]
    if not np.isfinite(gene_stats[["mean", "var", "var_tech"]].to_numpy(dtype=float)).all():
        raise RuntimeError("scDRS preprocessing produced non-finite gene statistics")
    preprocess_warning_text = "; ".join(sorted(set(f"{w.category.__name__}: {w.message}" for w in caught_pre))) or "none"

    qc_rows: list[dict[str, object]] = []
    score_paths: dict[str, Path] = {}
    full_paths: dict[str, Path] = {}
    donor_delta_by_trait: dict[str, dict[str, dict[str, float]]] = {}
    technical_rows: list[dict[str, object]] = []

    for trait in TRAITS:
        gene_list, gene_weights = dict_gs[trait]
        start = time.time()
        with warnings.catch_warnings(record=True) as caught_score:
            warnings.simplefilter("always")
            with np.errstate(all="ignore"):
                df = scdrs.score_cell(
                    adata,
                    gene_list,
                    gene_weight=gene_weights,
                    ctrl_match_key="mean_var",
                    n_ctrl=1000,
                    n_genebin=200,
                    weight_opt="vs",
                    copy=False,
                    return_ctrl_raw_score=False,
                    return_ctrl_norm_score=True,
                    random_seed=0,
                    verbose=False,
                )
        df = df.astype(np.float32)
        score_path = OUT / f"{trait}.score.gz"
        full_path = OUT / f"{trait}.full_score.gz"
        df.iloc[:, :6].to_csv(score_path, sep="\t", index=True, compression="gzip")
        df.to_csv(full_path, sep="\t", index=True, compression="gzip")
        score_paths[trait] = score_path
        full_paths[trait] = full_path
        raw = pd.to_numeric(df["raw_score"], errors="coerce").to_numpy(dtype=float)
        norm = pd.to_numeric(df["norm_score"], errors="coerce").to_numpy(dtype=float)
        controls = [x for x in df.columns if str(x).startswith("ctrl_norm_score_")]
        finite_raw = np.isfinite(raw)
        finite_norm = np.isfinite(norm)
        nan_count = int(pd.isna(df[["raw_score", "norm_score"]]).to_numpy().sum())
        inf_count = int(np.isinf(df[["raw_score", "norm_score"]].to_numpy(dtype=float)).sum())
        warning_text = "; ".join(sorted(set(f"{w.category.__name__}: {w.message}" for w in caught_score))) or "none"
        qc_rows.append({
            "trait": trait,
            "input_cells": int(adata.n_obs),
            "output_cells": int(len(df)),
            "trait_genes": 1000,
            "matched_genes": overlaps[trait],
            "finite_raw_scores": int(finite_raw.sum()),
            "finite_normalized_scores": int(finite_norm.sum()),
            "finite_norm_fraction": finite_fraction(norm),
            "NaN": nan_count,
            "Inf": inf_count,
            "control_sets": len(controls),
            "runtime_seconds": round(time.time() - start, 2),
            "runtime_warnings": warning_text,
            "preprocess_warnings": preprocess_warning_text,
            "status": "PASS" if finite_fraction(norm) >= 0.99 and len(controls) == 1000 else "FAIL",
        })
        donor_delta_by_trait[trait] = compute_donor_deltas(adata.obs, norm, trait)
        technical_rows.extend(technical_covariate_rows(adata.obs, norm, trait))
        del df

    qc = pd.DataFrame(qc_rows)
    qc.to_csv(OUT / "scdrs_score_qc.tsv", sep="\t", index=False)
    pd.DataFrame(technical_rows).to_csv(OUT / "scdrs_technical_covariate_qc.tsv", sep="\t", index=False)
    if (qc.loc[qc["trait"].isin(["POAG", "RE"]), "finite_norm_fraction"] < 0.99).any():
        raise RuntimeError("POAG or RE normalized scDRS score finite fraction is below 99%; biological interpretation stopped")
    return qc, score_paths, full_paths, donor_delta_by_trait


def technical_covariate_rows(obs: pd.DataFrame, norm: np.ndarray, trait: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    score = pd.Series(norm, index=obs.index, dtype=float)
    continuous = ["total_counts", "n_genes", "pct_counts_mt"]
    for variable in continuous:
        x = pd.to_numeric(obs[variable], errors="coerce")
        valid = x.notna() & score.notna() & np.isfinite(x) & np.isfinite(score)
        if valid.sum() >= 3 and x.loc[valid].nunique() > 1:
            rho, pval = spearmanr(x.loc[valid], score.loc[valid])
        else:
            rho, pval = math.nan, math.nan
        within: list[float] = []
        for _, idx in obs.loc[valid].groupby("donor_id", observed=True).groups.items():
            xx = x.loc[idx]
            yy = score.loc[idx]
            if len(idx) >= 3 and xx.nunique() > 1 and yy.nunique() > 1:
                r, _ = spearmanr(xx, yy)
                if np.isfinite(r):
                    within.append(float(r))
        max_abs = max((abs(x) for x in within), default=math.nan)
        rows.append({
            "trait": trait,
            "technical_variable": variable,
            "cells": int(valid.sum()),
            "global_spearman_rho": float(rho) if np.isfinite(rho) else math.nan,
            "global_spearman_p": float(pval) if np.isfinite(pval) else math.nan,
            "n_donors_within_test": len(within),
            "median_within_donor_rho": float(np.median(within)) if within else math.nan,
            "max_abs_within_donor_rho": max_abs,
            "major_dependence_flag": "FLAG_MAJOR_TECHNICAL_DEPENDENCE" if np.isfinite(max_abs) and max_abs >= 0.30 else "NO_MAJOR_DEPENDENCE_DETECTED",
            "status": "PASS" if valid.sum() >= 3 else "NOT_ESTIMABLE",
            "notes": "within-donor Spearman; descriptive QC only; scoring model not changed post hoc",
        })
    # In this reference each sample belongs to one donor, so sample is
    # perfectly confounded with donor and cannot be tested within donor.
    rows.append({
        "trait": trait,
        "technical_variable": "sample",
        "cells": int(len(obs)),
        "global_spearman_rho": math.nan,
        "global_spearman_p": math.nan,
        "n_donors_within_test": 0,
        "median_within_donor_rho": math.nan,
        "max_abs_within_donor_rho": math.nan,
        "major_dependence_flag": "NOT_ESTIMABLE_SAMPLE_CONFOUNDED_WITH_DONOR",
        "status": "NOT_ESTIMABLE",
        "notes": "one sample per donor in the eligible core; sample effect is represented by donor-aware analysis and mapping sensitivities",
    })
    return rows


def eligible_pairs(obs: pd.DataFrame) -> tuple[pd.DataFrame, dict[tuple[str, str], np.ndarray]]:
    classes = obs["harmonized_broad_class"].astype(str)
    donors = obs["donor_id"].astype(str)
    count = pd.crosstab(donors, classes).reindex(columns=CLASSES, fill_value=0)
    eligible = count >= 30
    masks: dict[tuple[str, str], np.ndarray] = {}
    for donor in count.index:
        for cls in CLASSES:
            masks[(str(donor), cls)] = np.flatnonzero((donors.to_numpy() == donor) & (classes.to_numpy() == cls))
    return eligible, masks


def compute_donor_deltas(obs: pd.DataFrame, norm: np.ndarray, trait: str) -> dict[str, dict[str, float]]:
    eligible, masks = eligible_pairs(obs)
    classes = obs["harmonized_broad_class"].astype(str).to_numpy()
    donors = obs["donor_id"].astype(str).to_numpy()
    result: dict[str, dict[str, float]] = {}
    for cls in CLASSES:
        for donor in eligible.index:
            if not bool(eligible.loc[donor, cls]):
                continue
            class_idx = masks[(str(donor), cls)]
            other_classes = [x for x in CLASSES if bool(eligible.loc[donor, x]) and x != cls]
            other_idx = np.concatenate([masks[(str(donor), x)] for x in other_classes]) if other_classes else np.array([], dtype=int)
            if len(other_idx) == 0:
                continue
            class_score = norm[class_idx]
            other_score = norm[other_idx]
            if not np.isfinite(class_score).all() or not np.isfinite(other_score).all():
                continue
            result[f"{donor}::{cls}"] = {
                "trait": trait,
                "donor_id": str(donor),
                "broad_class": cls,
                "n_cells_class": int(len(class_idx)),
                "n_cells_other": int(len(other_idx)),
                "class_mean_norm_score": float(np.mean(class_score)),
                "other_mean_norm_score": float(np.mean(other_score)),
                "delta": float(np.mean(class_score) - np.mean(other_score)),
            }
    return result


def sign_flip_p(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return math.nan
    if len(values) > 15:
        rng = np.random.default_rng(BASE_SEED)
        signs = rng.choice(np.array([-1.0, 1.0]), size=(200000, len(values)))
        null = signs @ values / len(values)
        return float((np.abs(null) >= abs(np.mean(values)) - 1e-12).mean())
    signs = np.array(np.meshgrid(*([[-1.0, 1.0]] * len(values)))).T.reshape(-1, len(values))
    null = signs @ values / len(values)
    return float((np.abs(null) >= abs(np.mean(values)) - 1e-12).mean())


def bootstrap_mean(values: np.ndarray, seed: int) -> tuple[float, float, float, float]:
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    boot = rng.choice(values, size=(N_BOOT, len(values)), replace=True).mean(axis=1)
    return float(np.median(boot)), float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975)), float(np.mean(boot > 0))


def primary_effects(obs: pd.DataFrame, donor_delta_by_trait: dict[str, dict[str, dict[str, float]]]):
    rows: list[dict[str, object]] = []
    donor_long: list[dict[str, object]] = []
    loo_rows: list[dict[str, object]] = []
    for trait in TRAITS:
        deltas = donor_delta_by_trait[trait]
        for cls in CLASSES:
            values = [v for k, v in deltas.items() if k.endswith(f"::{cls}")]
            values = sorted(values, key=lambda x: x["donor_id"])
            donor_values = np.array([x["delta"] for x in values], dtype=float)
            for item in values:
                donor_long.append(item)
            n = len(donor_values)
            mean = float(np.mean(donor_values)) if n else math.nan
            sd = float(np.std(donor_values, ddof=1)) if n > 1 else math.nan
            se = sd / math.sqrt(n) if n > 1 else math.nan
            tcrit = float(t.ppf(0.975, n - 1)) if n > 1 else math.nan
            ci_low = mean - tcrit * se if np.isfinite(se) else math.nan
            ci_high = mean + tcrit * se if np.isfinite(se) else math.nan
            if n >= 5 and np.isfinite(donor_values).all():
                primary_p = float(2 * t.sf(abs(mean / se), n - 1)) if np.isfinite(se) and se > 0 else (1.0 if mean == 0 else 0.0)
                normality_p = float(shapiro(donor_values).pvalue) if n >= 3 and np.unique(donor_values).size >= 3 else math.nan
                parametric_status = "PASS_FINITE_N_GE_5"
                flip_p = sign_flip_p(donor_values)
                boot_med, boot_low, boot_high, boot_pos = bootstrap_mean(donor_values, BASE_SEED + len(rows) + 1)
            else:
                primary_p = math.nan
                normality_p = math.nan
                parametric_status = "INSUFFICIENT_DONORS"
                flip_p = math.nan
                boot_med = boot_low = boot_high = boot_pos = math.nan
            donor_pos = int((donor_values > 0).sum()) if n else 0
            donor_neg = int((donor_values < 0).sum()) if n else 0
            row = {
                "trait": trait,
                "broad_class": cls,
                "n_donors": n,
                "n_cells": int(sum(x["n_cells_class"] for x in values)),
                "n_other_cells": int(sum(x["n_cells_other"] for x in values)),
                "mean_donor_delta": mean,
                "SE": se,
                "CI95_lower": ci_low,
                "CI95_upper": ci_high,
                "median_donor_delta": float(np.median(donor_values)) if n else math.nan,
                "donors_delta_gt_0": donor_pos,
                "donors_delta_lt_0": donor_neg,
                "donor_direction_fraction_positive": donor_pos / n if n else math.nan,
                "primary_t_p": primary_p,
                "sign_flip_p": flip_p,
                "normality_shapiro_p": normality_p,
                "parametric_test_status": parametric_status,
                "bootstrap_median": boot_med,
                "bootstrap_CI95_lower": boot_low,
                "bootstrap_CI95_upper": boot_high,
                "bootstrap_fraction_gt_0": boot_pos,
                "global_broad_FDR": math.nan,
                "within_trait_FDR": math.nan,
                "status": "PASS" if n >= 5 else "INSUFFICIENT_DONORS",
            }
            rows.append(row)
            donors = [x["donor_id"] for x in values]
            for removed in donors:
                remain = donor_values[np.array([d != removed for d in donors])]
                loo_rows.append({
                    "trait": trait,
                    "broad_class": cls,
                    "removed_donor": removed,
                    "remaining_donors": len(remain),
                    "leave_one_out_mean_delta": float(np.mean(remain)) if len(remain) else math.nan,
                    "status": "PASS" if len(remain) >= 5 else "NOT_ESTIMABLE",
                })
            if not donors:
                loo_rows.append({"trait": trait, "broad_class": cls, "removed_donor": "NA", "remaining_donors": 0, "leave_one_out_mean_delta": math.nan, "status": "NOT_ESTIMABLE"})

    effects = pd.DataFrame(rows)
    # Include all 27 prespecified tests in each FDR family; untestable rows
    # contribute p=1 rather than silently disappearing from the family.
    effects["global_broad_FDR"] = bh(effects["primary_t_p"].fillna(1.0).to_numpy())
    effects["within_trait_FDR"] = np.nan
    for trait in TRAITS:
        idx = effects["trait"] == trait
        effects.loc[idx, "within_trait_FDR"] = bh(effects.loc[idx, "primary_t_p"].fillna(1.0).to_numpy())
    loo = pd.DataFrame(loo_rows)
    loo_summary = loo.groupby(["trait", "broad_class"], observed=True).agg(
        loo_positive_fraction=("leave_one_out_mean_delta", lambda x: float((pd.to_numeric(x, errors="coerce") > 0).mean()) if np.isfinite(pd.to_numeric(x, errors="coerce")).any() else math.nan),
        loo_n_estimable=("status", lambda x: int((x == "PASS").sum())),
    ).reset_index()
    effects = effects.merge(loo_summary, on=["trait", "broad_class"], how="left", validate="one_to_one")
    effects["leave_one_donor_stability"] = np.where(effects["loo_positive_fraction"] >= 0.75, "POSITIVE_IN_ALL_OR_NEARLY_ALL", "UNSTABLE_OR_NOT_ESTIMABLE")
    effects.to_csv(OUT / "donor_class_effects.tsv", sep="\t", index=False)
    pd.DataFrame(donor_long).to_csv(OUT / "donor_class_deltas.tsv", sep="\t", index=False)
    loo.to_csv(OUT / "leave_one_donor_out.tsv", sep="\t", index=False)
    return effects, pd.DataFrame(donor_long), loo


def mapping_sensitivity(obs: pd.DataFrame, donor_delta_by_trait: dict[str, dict[str, dict[str, float]]], primary: pd.DataFrame):
    sample_frac = obs.groupby(["sample", "donor_id"], observed=True)["mapping_fraction"].first().reset_index()
    sample_frac = sample_frac.sort_values(["mapping_fraction", "sample"], ascending=[False, True])
    n_top = int(math.ceil(len(sample_frac) * 0.75))
    top_donors = set(sample_frac.head(n_top)["donor_id"].astype(str))
    lowest_donor = str(sample_frac.tail(1)["donor_id"].iloc[0])
    eligible, masks = eligible_pairs(obs)
    classes = obs["harmonized_broad_class"].astype(str).to_numpy()
    donors = obs["donor_id"].astype(str).to_numpy()
    rows: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []

    for trait in TRAITS:
        score_path = OUT / f"{trait}.score.gz"
        score_df = pd.read_csv(score_path, sep="\t", index_col=0)
        norm = score_df["norm_score"].to_numpy(dtype=float)
        # score order is the reference order written by score_cell.
        donor_to_pos = {str(d): np.flatnonzero(donors == str(d)) for d in sorted(set(donors))}
        for cls in CLASSES:
            primary_row = primary.query("trait == @trait and broad_class == @cls").iloc[0]
            for sensitivity, keep_donors, use_median in [
                ("A_TOP_75_PERCENT_MAPPING", top_donors, False),
                ("B_EXCLUDE_LOWEST_MAPPING", set(donors) - {lowest_donor}, False),
                ("C_MEDIAN_DONOR_CLASS_SCORE", set(donors), True),
            ]:
                deltas: list[float] = []
                for donor in sorted(keep_donors):
                    if donor not in donor_to_pos or not bool(eligible.loc[donor, cls]):
                        continue
                    cidx = masks[(donor, cls)]
                    others = [x for x in CLASSES if x != cls and bool(eligible.loc[donor, x])]
                    oidx = np.concatenate([masks[(donor, x)] for x in others]) if others else np.array([], dtype=int)
                    if len(oidx) == 0:
                        continue
                    deltas.append(float((np.median(norm[cidx]) if use_median else np.mean(norm[cidx])) - (np.median(norm[oidx]) if use_median else np.mean(norm[oidx]))))
                n = len(deltas)
                effect = float(np.mean(deltas)) if n else math.nan
                primary_effect = float(primary_row["mean_donor_delta"]) if np.isfinite(primary_row["mean_donor_delta"]) else math.nan
                reversal = bool(np.isfinite(effect) and np.isfinite(primary_effect) and np.sign(effect) != np.sign(primary_effect) and effect != 0 and primary_effect != 0)
                status = "PASS" if n >= 5 else "NOT_ESTIMABLE"
                rows.append({
                    "trait": trait,
                    "broad_class": cls,
                    "sensitivity": sensitivity,
                    "donors": n,
                    "excluded_donor_or_rule": lowest_donor if sensitivity == "B_EXCLUDE_LOWEST_MAPPING" else ("top_75_percent_mapping" if sensitivity == "A_TOP_75_PERCENT_MAPPING" else "median_score"),
                    "effect": effect,
                    "primary_effect": primary_effect,
                    "direction_reversal": reversal,
                    "status": status,
                })
            cls_rows = [x for x in rows if x["trait"] == trait and x["broad_class"] == cls]
            estimable = [x for x in cls_rows if x["status"] == "PASS"]
            summaries.append({
                "trait": trait,
                "broad_class": cls,
                "mapping_sensitivity_min_effect": float(min(x["effect"] for x in estimable)) if estimable else math.nan,
                "mapping_sensitivity_direction_reversals": int(sum(bool(x["direction_reversal"]) for x in estimable)),
                "mapping_sensitivity_status": "PASS" if estimable and not any(bool(x["direction_reversal"]) for x in estimable) else ("NOT_ESTIMABLE" if not estimable else "SEVERE_SENSITIVITY_FAILURE"),
            })
        del score_df
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "mapping_selection_sensitivity.tsv", sep="\t", index=False)
    return result, pd.DataFrame(summaries)


def cell_balance_sensitivity(obs: pd.DataFrame, primary: pd.DataFrame):
    # Use the already generated scores, while downsampling only the biological
    # unit's contributing cells.  This is deterministic and secondary.
    rng = np.random.default_rng(BASE_SEED)
    eligible, masks = eligible_pairs(obs)
    classes = obs["harmonized_broad_class"].astype(str).to_numpy()
    donors = obs["donor_id"].astype(str).to_numpy()
    balanced: dict[tuple[str, str], np.ndarray] = {}
    for donor in sorted(set(donors)):
        for cls in CLASSES:
            idx = masks[(donor, cls)]
            balanced[(donor, cls)] = rng.choice(idx, size=min(500, len(idx)), replace=False) if len(idx) else np.array([], dtype=int)
    rows: list[dict[str, object]] = []
    for trait in TRAITS:
        score_df = pd.read_csv(OUT / f"{trait}.score.gz", sep="\t", index_col=0)
        norm = score_df["norm_score"].to_numpy(dtype=float)
        for cls in CLASSES:
            deltas: list[float] = []
            for donor in sorted(set(donors)):
                if not bool(eligible.loc[donor, cls]):
                    continue
                cidx = balanced[(donor, cls)]
                others = [x for x in CLASSES if x != cls and bool(eligible.loc[donor, x])]
                oidx = np.concatenate([balanced[(donor, x)] for x in others]) if others else np.array([], dtype=int)
                if len(cidx) == 0 or len(oidx) == 0:
                    continue
                deltas.append(float(np.mean(norm[cidx]) - np.mean(norm[oidx])))
            rows.append({
                "trait": trait,
                "broad_class": cls,
                "cap_cells_per_donor_class": 500,
                "n_donors": len(deltas),
                "balanced_mean_delta": float(np.mean(deltas)) if deltas else math.nan,
                "balanced_median_delta": float(np.median(deltas)) if deltas else math.nan,
                "primary_mean_delta": float(primary.query("trait == @trait and broad_class == @cls")["mean_donor_delta"].iloc[0]),
                "status": "PASS" if len(deltas) >= 5 else "NOT_ESTIMABLE",
            })
        del score_df
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "cell_count_balance_sensitivity.tsv", sep="\t", index=False)
    return result


def compute_neighbors(adata) -> None:
    from sklearn.decomposition import TruncatedSVD
    import scanpy as sc

    svd = TruncatedSVD(n_components=30, n_iter=3, random_state=0)
    adata.obsm["X_pca"] = svd.fit_transform(adata.X).astype(np.float32)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, use_rep="X_pca", random_state=0)


def official_group_analysis(adata, full_paths: dict[str, Path]):
    import scdrs

    # scDRS 1.0.2's reference Geary-C implementation loops over every edge in
    # Python. This is algebraically equivalent and uses sparse matrix-vector
    # products, while retaining the official downstream_group_analysis API and
    # all locked group-analysis parameters.
    def vectorized_gearys_c(group_adata, vals):
        graph = group_adata.obsp["connectivities"].tocsr()
        vals = np.asarray(vals, dtype=np.float64)
        n_obs = graph.shape[0]
        row_degree = np.asarray(graph.sum(axis=1)).ravel().astype(np.float64)
        col_degree = np.asarray(graph.sum(axis=0)).ravel().astype(np.float64)
        weight_sum = float(graph.data.astype(np.float64, copy=False).sum())
        mean = float(vals.mean())
        denom = 2.0 * weight_sum * float(np.square(vals - mean).sum())
        if denom == 0.0:
            return np.nan
        total = float(
            np.dot(row_degree + col_degree, np.square(vals))
            - 2.0 * np.dot(vals, graph.dot(vals))
        )
        if total < 0.0 and abs(total) < 1e-10:
            total = 0.0
        return float((n_obs - 1) * total / denom)

    scdrs.method.gearys_c = vectorized_gearys_c

    compute_neighbors(adata)
    rows: list[pd.DataFrame] = []
    for trait in TRAITS:
        df_full = pd.read_csv(full_paths[trait], sep="\t", index_col=0).astype(np.float32)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = scdrs.method.downstream_group_analysis(
                adata,
                df_full,
                ["harmonized_broad_class"],
                fdr_thresholds=[0.05, 0.1, 0.2],
            )["harmonized_broad_class"]
        result = result.reset_index().rename(columns={"group": "broad_class"})
        result["trait"] = trait
        result["runtime_warnings"] = "; ".join(sorted(set(f"{w.category.__name__}: {w.message}" for w in caught))) or "none"
        rows.append(result)
        del df_full
    group = pd.concat(rows, ignore_index=True)
    group["group_FDR"] = bh(group["assoc_mcp"].fillna(1.0).to_numpy())
    group.to_csv(OUT / "scdrs_official_group_analysis.tsv", sep="\t", index=False)
    return group


def donor_group_concordance(effects: pd.DataFrame, group: pd.DataFrame):
    merged = effects.merge(group[["trait", "broad_class", "assoc_mcz", "group_FDR"]], on=["trait", "broad_class"], how="left", validate="one_to_one")
    rows: list[dict[str, object]] = []
    for trait in TRAITS:
        x = merged.loc[merged["trait"] == trait].copy()
        valid = x["mean_donor_delta"].notna() & x["assoc_mcz"].notna()
        xx = x.loc[valid, "mean_donor_delta"].to_numpy(dtype=float)
        yy = x.loc[valid, "assoc_mcz"].to_numpy(dtype=float)
        if len(xx) >= 3:
            rho, pval = spearmanr(xx, yy)
            rank_x = pd.Series(xx).rank(method="average").to_numpy()
            rank_y = pd.Series(yy).rank(method="average").to_numpy()
            rank_concordance = float(1 - np.abs(rank_x - rank_y).mean() / max(len(xx) - 1, 1))
            discord = x.loc[valid].assign(rank_delta=np.abs(rank_x - rank_y)).query("rank_delta >= 3")["broad_class"].tolist()
        else:
            rho = pval = rank_concordance = math.nan
            discord = []
        rows.append({
            "trait": trait,
            "n_classes": int(len(xx)),
            "spearman_rho": float(rho) if np.isfinite(rho) else math.nan,
            "spearman_p": float(pval) if np.isfinite(pval) else math.nan,
            "rank_concordance": rank_concordance,
            "major_discordant_classes": ";".join(discord) if discord else "none",
            "notes": "comparison is donor-aware effect versus official scDRS assoc_mcz; cell-level group P values are secondary",
        })
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "donor_vs_scdrs_group_concordance.tsv", sep="\t", index=False)
    return result, merged


def make_pseudobulk_gene_covariates(obs: pd.DataFrame):
    import anndata

    raw = anndata.read_h5ad(H5AD)
    if not sparse.issparse(raw.X):
        raise RuntimeError("raw reference is not sparse for MAGMA pseudobulk baseline")
    eligible, masks = eligible_pairs(raw.obs)
    profile_by_class: dict[str, list[np.ndarray]] = {x: [] for x in CLASSES}
    donor_class_counts: list[dict[str, object]] = []
    for donor in sorted(eligible.index):
        for cls in CLASSES:
            if not bool(eligible.loc[donor, cls]):
                continue
            idx = masks[(str(donor), cls)]
            counts = np.asarray(raw.X[idx].sum(axis=0)).ravel().astype(float)
            lib = counts.sum()
            if lib <= 0:
                continue
            log_cpm = np.log1p(counts / lib * 1e6)
            profile_by_class[cls].append(log_cpm)
            donor_class_counts.append({"donor_id": str(donor), "broad_class": cls, "cells": int(len(idx)), "library_size": float(lib)})
    class_profile = {cls: np.vstack(rows).mean(axis=0) for cls, rows in profile_by_class.items() if rows}
    if set(class_profile) != set(CLASSES):
        raise RuntimeError("pseudobulk baseline lacks one or more locked broad classes")
    average = np.vstack([class_profile[x] for x in CLASSES]).mean(axis=0)
    score = pd.read_csv(MAGMA_SCORES, sep="\t", dtype={"gene_id": str})
    var_names = raw.var_names.astype(str).tolist()
    symbol_to_pos: dict[str, int] = {}
    for pos, symbol in enumerate(var_names):
        symbol_to_pos.setdefault(symbol, pos)
    rows: list[dict[str, object]] = []
    for item in score.itertuples(index=False):
        symbol = str(item.gene_symbol) if pd.notna(item.gene_symbol) else ""
        if symbol not in symbol_to_pos:
            continue
        pos = symbol_to_pos[symbol]
        row = {"GENE": str(item.gene_id), "Average": float(average[pos])}
        for cls in CLASSES:
            row[CLASS_VAR[cls]] = float(class_profile[cls][pos] - average[pos])
        rows.append(row)
    covar = pd.DataFrame(rows)
    if covar["GENE"].duplicated().any():
        covar = covar.drop_duplicates("GENE", keep="first")
    covar_path = OUT / "magma_broad_gene_covariates.tsv"
    covar.to_csv(covar_path, sep="\t", index=False, float_format="%.8g")
    pd.DataFrame(donor_class_counts).to_csv(OUT / "pseudobulk_donor_class_qc.tsv", sep="\t", index=False)
    del raw
    return covar_path, covar


def magma_gene_property(covar_path: Path):
    rows: list[dict[str, object]] = []
    for trait in TRAITS:
        prefix = OUT / "magma" / f"{trait}_broad_celltype"
        prefix.parent.mkdir(parents=True, exist_ok=True)
        raw_gene = ROOT / "results" / "phase3ar8" / "magma" / f"{trait}.genes.raw"
        command = [str(MAGMA), "--gene-results", str(raw_gene), "--gene-covar", str(covar_path), "--model", "direction=greater", "condition-hide=Average", "--out", str(prefix)]
        result = run_command(command, LOGS / f"{trait}_magma_broad_celltype.log")
        out = prefix.with_suffix(".gsa.out")
        if result.returncode != 0 or not out.exists():
            raise RuntimeError(f"MAGMA gene-property analysis failed for {trait}; see {LOGS / f'{trait}_magma_broad_celltype.log'}")
        frame = pd.read_csv(out, sep=r"\s+", engine="python", comment="#")
        frame.columns = [str(x).upper() for x in frame.columns]
        for _, item in frame.iterrows():
            var = str(item["VARIABLE"])
            if var not in set(CLASS_VAR.values()):
                continue
            cls = next(k for k, v in CLASS_VAR.items() if v == var)
            rows.append({
                "trait": trait,
                "broad_class": cls,
                "magma_variable": var,
                "NGENES": int(item["NGENES"]),
                "beta": float(item["BETA"]),
                "beta_std": float(item["BETA_STD"]),
                "SE": float(item["SE"]),
                "p": float(item["P"]),
                "method": "MAGMA gene-property direction=greater conditioned on Average",
                "status": "PASS",
            })
    result = pd.DataFrame(rows)
    if len(result) != 27:
        raise RuntimeError(f"MAGMA broad-cell property output has {len(result)} class rows, expected 27")
    result["global_broad_FDR"] = bh(result["p"].to_numpy(dtype=float))
    result.to_csv(OUT / "magma_broad_celltype_enrichment.tsv", sep="\t", index=False)
    return result


def cross_trait_similarity(donor_delta_by_trait: dict[str, dict[str, dict[str, float]]]):
    effect_map: dict[tuple[str, str], float] = {}
    donor_map: dict[tuple[str, str], dict[str, float]] = {}
    for trait in TRAITS:
        for cls in CLASSES:
            items = [v for k, v in donor_delta_by_trait[trait].items() if k.endswith(f"::{cls}")]
            effect_map[(trait, cls)] = float(np.mean([x["delta"] for x in items])) if items else math.nan
            donor_map[(trait, cls)] = {x["donor_id"]: float(x["delta"]) for x in items}
    pairs = [("AMD", "POAG"), ("POAG", "RE"), ("AMD", "RE")]
    rows: list[dict[str, object]] = []
    rng = np.random.default_rng(BASE_SEED)
    for a, b in pairs:
        x = np.array([effect_map[(a, cls)] for cls in CLASSES], dtype=float)
        y = np.array([effect_map[(b, cls)] for cls in CLASSES], dtype=float)
        valid = np.isfinite(x) & np.isfinite(y)
        zx = (x[valid] - np.mean(x[valid])) / np.std(x[valid], ddof=1) if valid.sum() >= 2 and np.std(x[valid], ddof=1) > 0 else np.full(valid.sum(), np.nan)
        zy = (y[valid] - np.mean(y[valid])) / np.std(y[valid], ddof=1) if valid.sum() >= 2 and np.std(y[valid], ddof=1) > 0 else np.full(valid.sum(), np.nan)
        rho, pval = spearmanr(zx, zy) if valid.sum() >= 3 and np.isfinite(zx).all() and np.isfinite(zy).all() else (math.nan, math.nan)
        cos = float(np.dot(zx, zy) / (np.linalg.norm(zx) * np.linalg.norm(zy))) if valid.sum() >= 2 and np.isfinite(zx).all() and np.isfinite(zy).all() else math.nan
        boot: list[float] = []
        donor_union = sorted(set().union(*(donor_map[(a, cls)] for cls in CLASSES), *(donor_map[(b, cls)] for cls in CLASSES)))
        for _ in range(N_BOOT):
            sampled = rng.choice(donor_union, size=len(donor_union), replace=True)
            vx: list[float] = []
            vy: list[float] = []
            for cls in CLASSES:
                aa = [donor_map[(a, cls)][d] for d in sampled if d in donor_map[(a, cls)]]
                bb = [donor_map[(b, cls)][d] for d in sampled if d in donor_map[(b, cls)]]
                vx.append(float(np.mean(aa)) if aa else math.nan)
                vy.append(float(np.mean(bb)) if bb else math.nan)
            vx = np.asarray(vx, dtype=float)
            vy = np.asarray(vy, dtype=float)
            ok = np.isfinite(vx) & np.isfinite(vy)
            if ok.sum() >= 3 and np.std(vx[ok]) > 0 and np.std(vy[ok]) > 0:
                boot.append(float(spearmanr(vx[ok], vy[ok]).statistic))
        rows.append({
            "trait_a": a,
            "trait_b": b,
            "n_classes": int(valid.sum()),
            "spearman_rho": float(rho) if np.isfinite(rho) else math.nan,
            "spearman_p": float(pval) if np.isfinite(pval) else math.nan,
            "cosine_similarity": cos,
            "bootstrap_replicates": len(boot),
            "bootstrap_median": float(np.median(boot)) if boot else math.nan,
            "bootstrap_CI95_lower": float(np.quantile(boot, 0.025)) if boot else math.nan,
            "bootstrap_CI95_upper": float(np.quantile(boot, 0.975)) if boot else math.nan,
            "bootstrap_fraction_gt_0": float(np.mean(np.asarray(boot) > 0)) if boot else math.nan,
            "standardization": "nine-class effects standardized within each trait before similarity; raw p values not compared",
        })
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "cross_trait_broad_profile_similarity.tsv", sep="\t", index=False)
    # Also save the standardized profile used by Figure P14.
    matrix = pd.DataFrame({trait: [effect_map[(trait, cls)] for cls in CLASSES] for trait in TRAITS}, index=CLASSES)
    z = matrix.apply(lambda x: (x - x.mean()) / x.std(ddof=1), axis=0)
    z.index.name = "broad_class"
    z.to_csv(OUT / "cross_trait_standardized_effect_profiles.tsv", sep="\t")
    return result, matrix, z


def consensus_and_state(effects: pd.DataFrame, group: pd.DataFrame, magma: pd.DataFrame, concordance: pd.DataFrame, mapping_summary: pd.DataFrame):
    merged = effects.merge(group[["trait", "broad_class", "assoc_mcz", "group_FDR"]], on=["trait", "broad_class"], how="left", validate="one_to_one")
    merged = merged.merge(magma[["trait", "broad_class", "beta", "SE", "p", "global_broad_FDR"]].rename(columns={"global_broad_FDR": "magma_FDR", "p": "magma_p"}), on=["trait", "broad_class"], how="left", validate="one_to_one")
    merged = merged.merge(mapping_summary, on=["trait", "broad_class"], how="left", validate="one_to_one")
    merged["tier"] = "NONE"
    tier1 = (merged["global_broad_FDR"] < 0.05) & (merged["mean_donor_delta"] > 0) & (merged["donor_direction_fraction_positive"] >= 0.75) & (merged["loo_positive_fraction"] >= 0.75) & (merged["mapping_sensitivity_status"] == "PASS")
    tier2 = (~tier1) & ((merged["within_trait_FDR"] < 0.05) | ((merged["primary_t_p"] < 0.05) & (merged["loo_positive_fraction"] >= 0.75) & (merged["magma_FDR"] < 0.05))) & (merged["mean_donor_delta"] > 0)
    merged.loc[tier1, "tier"] = "TIER1_BROAD_CLASS"
    merged.loc[tier2, "tier"] = "TIER2_BROAD_ASSOCIATION"
    merged["confidence"] = "NULL/UNSUPPORTED"
    high = (merged["tier"] == "TIER1_BROAD_CLASS") & (merged["group_FDR"] < 0.05) & (merged["magma_FDR"] < 0.05)
    moderate = (merged["tier"].isin(["TIER1_BROAD_CLASS", "TIER2_BROAD_ASSOCIATION"])) & ((merged["group_FDR"] < 0.05) | (merged["magma_FDR"] < 0.05))
    exploratory = (merged["primary_t_p"] < 0.05) | (merged["within_trait_FDR"] < 0.10) | (merged["magma_p"] < 0.05)
    merged.loc[exploratory, "confidence"] = "EXPLORATORY"
    merged.loc[moderate, "confidence"] = "MODERATE_CONFIDENCE"
    merged.loc[high, "confidence"] = "HIGH_CONFIDENCE"
    merged.to_csv(OUT / "BROAD_CELLULAR_CONSENSUS.tsv", sep="\t", index=False)

    class_labels: list[dict[str, object]] = []
    for cls in CLASSES:
        robust = {trait for trait in TRAITS if ((merged["trait"] == trait) & (merged["broad_class"] == cls) & merged["tier"].isin(["TIER1_BROAD_CLASS", "TIER2_BROAD_ASSOCIATION"]) & (merged["mean_donor_delta"] > 0)).any()}
        if len(robust) == 3:
            label = "PAN_TRAIT_CONVERGENT"
        elif robust == {"POAG", "RE"}:
            label = "POAG_RE_CONVERGENT"
        elif robust == {"AMD", "POAG"}:
            label = "AMD_POAG_CONVERGENT"
        elif robust == {"AMD", "RE"}:
            label = "AMD_RE_CONVERGENT"
        elif robust == {"POAG"}:
            label = "POAG_SPECIFIC"
        elif robust == {"RE"}:
            label = "RE_SPECIFIC"
        elif robust == {"AMD"}:
            label = "AMD_SUPPORTED"
        else:
            label = "NO_ROBUST_SIGNAL"
        class_labels.append({"broad_class": cls, "robust_traits": ";".join(sorted(robust)) if robust else "none", "cellular_convergence_classification": label})
    class_labels_df = pd.DataFrame(class_labels)
    merged = merged.merge(class_labels_df, on="broad_class", how="left")
    merged.to_csv(OUT / "BROAD_CELLULAR_CONSENSUS.tsv", sep="\t", index=False)

    fine = pd.read_csv(FINE_READY, sep="\t")
    ready_by_class = fine.loc[fine["readiness"] == "STATE_READY"].groupby("harmonized_broad_class", observed=True).agg(
        fine_state_donor_support=("donors", "max"),
        fine_states_ready=("cell_type", lambda x: ";".join(x.astype(str))),
    ).reset_index()
    state_rows: list[dict[str, object]] = []
    for cls in CLASSES:
        sub = merged.loc[merged["broad_class"] == cls]
        qualifying = sub.loc[sub["tier"].isin(["TIER1_BROAD_CLASS", "TIER2_BROAD_ASSOCIATION"]), "trait"].tolist()
        mag_support = sub.loc[(sub["tier"] == "TIER2_BROAD_ASSOCIATION") | (sub["tier"] == "TIER1_BROAD_CLASS"), "magma_FDR"].lt(0.05).any()
        fine_row = ready_by_class.loc[ready_by_class["harmonized_broad_class"] == cls]
        fine_support = int(fine_row["fine_state_donor_support"].iloc[0]) if not fine_row.empty else 0
        fine_states = str(fine_row["fine_states_ready"].iloc[0]) if not fine_row.empty else "none"
        eligible = bool(fine_support >= 5 and (len(qualifying) > 0 or mag_support))
        state_rows.append({
            "lineage": cls,
            "qualifying_traits": ";".join(sorted(set(qualifying))) if qualifying else "none",
            "broad_evidence": ";".join(f"{r.trait}:{r.tier}" for r in sub.itertuples() if r.tier != "NONE") or "none",
            "fine_state_donor_support": fine_support,
            "fine_states_ready": fine_states,
            "state_analysis_eligibility": "STATE_ANALYSIS_ELIGIBLE" if eligible else "NOT_STATE_ANALYSIS_ELIGIBLE",
        })
    state = pd.DataFrame(state_rows)
    state["priority_rank"] = np.nan
    candidate_mask = state["state_analysis_eligibility"] == "STATE_ANALYSIS_ELIGIBLE"
    if candidate_mask.any():
        rank_score = []
        for cls in state.loc[candidate_mask, "lineage"]:
            s = merged.loc[merged["broad_class"] == cls]
            rank_score.append((cls, float(np.nanmin(s["global_broad_FDR"].to_numpy()))))
        rank_score.sort(key=lambda x: x[1])
        rank_map = {cls: i + 1 for i, (cls, _) in enumerate(rank_score[:5])}
        state["priority_rank"] = state["lineage"].map(rank_map)
    state.to_csv(OUT / "state_analysis_eligibility.tsv", sep="\t", index=False)
    return merged, state, class_labels_df


def verdict(consensus: pd.DataFrame) -> tuple[str, str]:
    robust_by_class = {}
    for cls in CLASSES:
        robust_by_class[cls] = set(consensus.loc[(consensus["broad_class"] == cls) & consensus["tier"].isin(["TIER1_BROAD_CLASS", "TIER2_BROAD_ASSOCIATION"]) & (consensus["mean_donor_delta"] > 0), "trait"])
    convergent = [x for x, traits in robust_by_class.items() if len(traits) >= 2]
    trait_specific = [x for x, traits in robust_by_class.items() if len(traits) == 1]
    any_robust = bool(convergent or trait_specific)
    any_suggestive = bool((consensus["primary_t_p"] < 0.05).any() or (consensus["within_trait_FDR"] < 0.10).any() or (consensus["magma_p"] < 0.05).any())
    if convergent and trait_specific:
        return "BROAD_CELLULAR_MIXED_GO", "A. RUN PHASE 3A-2 CONVERGENT STATE ANALYSIS"
    if convergent:
        return "BROAD_CELLULAR_CONVERGENCE_GO", "A. RUN PHASE 3A-2 CONVERGENT STATE ANALYSIS"
    if trait_specific:
        return "TRAIT_SPECIFIC_BROAD_GO", "B. RUN PHASE 3A-2 TRAIT-SPECIFIC STATE ANALYSIS"
    if any_robust:
        return "BROAD_CELLULAR_SIGNAL_WEAK", "C. EXPAND/VALIDATE BROAD CELLULAR ANALYSIS"
    if any_suggestive:
        return "BROAD_CELLULAR_SIGNAL_WEAK", "C. EXPAND/VALIDATE BROAD CELLULAR ANALYSIS"
    return "BROAD_CELLULAR_NULL", "D. STOP CELL-STATE MAPPING"


def make_figures(effects: pd.DataFrame, consensus: pd.DataFrame, similarity: pd.DataFrame, magma: pd.DataFrame) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    plt.rcParams.update({
        "font.family": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 11,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 1.3,
        "legend.frameon": False,
        "svg.fonttype": "none",
    })
    blue = "#0F4D92"
    teal = "#42949E"
    violet = "#9A4D8E"
    red = "#B64342"
    outputs: list[Path] = []
    z = pd.read_csv(OUT / "cross_trait_standardized_effect_profiles.tsv", sep="\t", index_col=0).reindex(CLASSES)
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    vmax = max(1.5, float(np.nanmax(np.abs(z.to_numpy()))))
    im = ax.imshow(z.to_numpy(), cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax), aspect="auto")
    ax.set_xticks(range(3), TRAITS)
    ax.set_yticks(range(len(CLASSES)), CLASSES)
    for i, cls in enumerate(CLASSES):
        for j, trait in enumerate(TRAITS):
            row = consensus.query("broad_class == @cls and trait == @trait").iloc[0]
            mark = "*" if row["global_broad_FDR"] < 0.05 else ("·" if row["within_trait_FDR"] < 0.05 else "")
            ax.text(j, i, mark, ha="center", va="center", color="black", fontsize=13, fontweight="bold")
    ax.set_title("P14  Donor-aware broad cellular effects\nwithin-trait standardized effects")
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("standardized donor-aware effect")
    fig.tight_layout(pad=1.2)
    for ext in ("png", "svg"):
        path = FIGURES / f"Figure_P14_donor_aware_effect_heatmap.{ext}"
        fig.savefig(path, dpi=300, bbox_inches="tight")
        outputs.append(path)
    plt.close(fig)

    sel = consensus.loc[consensus["tier"].isin(["TIER1_BROAD_CLASS", "TIER2_BROAD_ASSOCIATION"])].copy()
    if sel.empty:
        sel = consensus.sort_values(["global_broad_FDR", "primary_t_p"], na_position="last").groupby("trait", observed=True).head(1).copy()
    sel["label"] = sel["trait"] + " | " + sel["broad_class"]
    deltas = pd.read_csv(OUT / "donor_class_deltas.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(9.5, max(3.5, 0.55 * len(sel) + 1.5)))
    y_positions = np.arange(len(sel))
    for y, (_, row) in enumerate(sel.iterrows()):
        d = deltas.query("trait == @row.trait and broad_class == @row.broad_class")
        ax.scatter(d["delta"], np.full(len(d), y), s=35, color=blue, alpha=0.75, zorder=3)
        ax.errorbar(row["mean_donor_delta"], y, xerr=[[row["mean_donor_delta"] - row["CI95_lower"]], [row["CI95_upper"] - row["mean_donor_delta"]]], fmt="D", color=red, capsize=3, markersize=5, zorder=4)
    ax.axvline(0, color="0.35", lw=1)
    ax.set_yticks(y_positions, sel["label"])
    ax.set_xlabel("norm_score class − other eligible classes\n(donor-level delta)")
    ax.set_title("P15  Donor-level effects and 95% CI")
    fig.tight_layout(pad=1.2)
    for ext in ("png", "svg"):
        path = FIGURES / f"Figure_P15_donor_level_forest_dot.{ext}"
        fig.savefig(path, dpi=300, bbox_inches="tight")
        outputs.append(path)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    colors = [blue, violet, teal]
    x = np.arange(len(similarity))
    for i, row in similarity.iterrows():
        y = row["spearman_rho"]
        lo = y - row["bootstrap_CI95_lower"] if np.isfinite(y) and np.isfinite(row["bootstrap_CI95_lower"]) else 0
        hi = row["bootstrap_CI95_upper"] - y if np.isfinite(y) and np.isfinite(row["bootstrap_CI95_upper"]) else 0
        ax.errorbar(x[i], y, yerr=[[lo], [hi]], fmt="o", color=colors[i], capsize=4, markersize=7)
    ax.axhline(0, color="0.35", lw=1)
    ax.set_xticks(x, [f"{a}–{b}" for a, b in zip(similarity["trait_a"], similarity["trait_b"])])
    ax.set_ylabel("Spearman profile similarity")
    ax.set_title("P16  Cross-trait broad cellular profiles")
    fig.tight_layout(pad=1.2)
    for ext in ("png", "svg"):
        path = FIGURES / f"Figure_P16_cross_trait_profile_similarity.{ext}"
        fig.savefig(path, dpi=300, bbox_inches="tight")
        outputs.append(path)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.8, 5.2))
    for trait, color in zip(TRAITS, [blue, violet, teal]):
        d = consensus.query("trait == @trait")
        ax.scatter(d["mean_donor_delta"], d["beta"], label=trait, s=48, alpha=0.8, color=color)
    ax.axvline(0, color="0.75", lw=1)
    ax.axhline(0, color="0.75", lw=1)
    ax.set_xlabel("donor-aware mean delta")
    ax.set_ylabel("MAGMA gene-property beta")
    ax.set_title("P17  Donor-aware scDRS vs MAGMA broad-cell evidence")
    ax.legend()
    fig.tight_layout(pad=1.2)
    for ext in ("png", "svg"):
        path = FIGURES / f"Figure_P17_scDRS_vs_MAGMA_concordance.{ext}"
        fig.savefig(path, dpi=300, bbox_inches="tight")
        outputs.append(path)
    plt.close(fig)
    return outputs


def fmt(v, digits=4):
    if v is None or not np.isfinite(float(v)):
        return "NA"
    return f"{float(v):.{digits}f}"


def write_reports(consensus: pd.DataFrame, state: pd.DataFrame, similarity: pd.DataFrame, group_concordance: pd.DataFrame, score_qc: pd.DataFrame, magma: pd.DataFrame, verdict_name: str, next_action: str, figure_paths: list[Path]) -> None:
    eligible = state.query("state_analysis_eligibility == 'STATE_ANALYSIS_ELIGIBLE'").sort_values("priority_rank", na_position="last").head(5)
    tier_rows = {trait: consensus.loc[(consensus["trait"] == trait) & consensus["tier"].isin(["TIER1_BROAD_CLASS", "TIER2_BROAD_ASSOCIATION"])].sort_values("global_broad_FDR") for trait in TRAITS}
    supported_magma = {trait: magma.loc[(magma["trait"] == trait) & (magma["global_broad_FDR"] < 0.05), "broad_class"].tolist() for trait in TRAITS}
    top_effects = {}
    for trait in TRAITS:
        trait_rows = consensus.loc[consensus["trait"] == trait]
        estimable = trait_rows.loc[trait_rows["n_donors"] >= 5]
        top_effects[trait] = (estimable if len(estimable) else trait_rows).sort_values(
            "mean_donor_delta", ascending=False
        ).iloc[0]
    primary_text = []
    for trait in TRAITS:
        q = score_qc.loc[score_qc["trait"] == trait].iloc[0]
        rows = consensus.loc[consensus["trait"] == trait].sort_values("global_broad_FDR")
        primary_text.append(f"- {trait}: {int(q['input_cells'])} input/output cells; finite norm_score {int(q['finite_normalized_scores'])}/{int(q['output_cells'])}; matched genes {int(q['matched_genes'])}; lowest global broad FDR {fmt(rows['global_broad_FDR'].min())}.")
    report = f"""# PHASE 3A-1 BROAD CELLULAR ARCHITECTURE

## VERDICT

{verdict_name}

## REFERENCE AND LOCK

Primary reference: `HRCA_EXACT_MAPPED_HEALTHY_CORE` (`{H5AD.relative_to(ROOT)}`), all readable exact official-ID mapped healthy/normal cells, no post-mapping cap.

- cells = 98,834
- genes = 36,601
- healthy donors = 8
- samples = 8
- broad classes = 9: {', '.join(CLASSES)}
- exact mapping fraction = 98,834 / 114,267 = 0.8649
- expression matrix = sparse raw 10x counts; scDRS normalization was applied once (size-factor 1e4 + log1p)
- donor-class eligibility = at least 30 cells; broad-class inference = at least 5 eligible donors
- Phase 3A-1 lock = [`phase3a1_analysis_lock.yaml`](../config/phase3a1_analysis_lock.yaml)

## SCORE QC

{chr(10).join(primary_text)}

All POAG and RE finite normalized-score fractions were at least 99%; no technical stop rule was triggered. Technical covariate QC is descriptive and does not alter the locked scoring model.

## PRIMARY DONOR-AWARE INFERENCE

The primary effect is the within-donor mean scDRS normalized score for one broad class minus the mean among all other eligible broad classes. Donors are equally weighted. Primary P values are one-sample donor-level t tests; exact donor sign-flip P values and 5,000-donor-bootstrap intervals are reported separately.

### Tier 1

"""
    for trait in TRAITS:
        t1 = tier_rows[trait].query("tier == 'TIER1_BROAD_CLASS'")
        report += f"- {trait}: {', '.join(t1['broad_class']) if len(t1) else 'none'}\n"
    report += "\n### Tier 2\n\n"
    for trait in TRAITS:
        t2 = tier_rows[trait].query("tier == 'TIER2_BROAD_ASSOCIATION'")
        report += f"- {trait}: {', '.join(t2['broad_class']) if len(t2) else 'none'}\n"
    report += "\n### Strongest positive donor-aware estimable effect by trait\n\n"
    for trait in TRAITS:
        row = top_effects[trait]
        report += f"- {trait}: {row['broad_class']}; mean delta = {fmt(row['mean_donor_delta'])}, 95% CI = [{fmt(row['CI95_lower'])}, {fmt(row['CI95_upper'])}], global FDR = {fmt(row['global_broad_FDR'])}, donor direction = {fmt(row['donor_direction_fraction_positive'])}.\n"
    report += f"""
## OFFICIAL scDRS GROUP ANALYSIS

Official `downstream_group_analysis` results using `broad_class` are secondary. They report cell-level Monte Carlo association/heterogeneity summaries and are not donor-level replication. Donor-versus-group concordance was saved to `results/phase3a1/donor_vs_scdrs_group_concordance.tsv`.

## ORTHOGONAL MAGMA BASELINE

MAGMA gene-property analysis used donor-balanced pseudobulk log1p-CPM profiles, an equal-class `Average` covariate, and one-sided positive testing. It is an orthogonal transcriptomic baseline, not a causal analysis.

- AMD supported classes (global MAGMA FDR < 0.05): {', '.join(supported_magma['AMD']) if supported_magma['AMD'] else 'none'}
- POAG supported classes (global MAGMA FDR < 0.05): {', '.join(supported_magma['POAG']) if supported_magma['POAG'] else 'none'}
- RE supported classes (global MAGMA FDR < 0.05): {', '.join(supported_magma['RE']) if supported_magma['RE'] else 'none'}

## STATE-ANALYSIS ELIGIBILITY

No fine-state analysis was run. The following locked-rule candidates are eligible for a possible next phase, limited to five priorities:

"""
    if len(eligible):
        for row in eligible.itertuples():
            report += f"- {row.lineage}: traits = {row.qualifying_traits}; broad evidence = {row.broad_evidence}; fine-state donor support = {row.fine_state_donor_support}; ready states = {row.fine_states_ready}.\n"
    else:
        report += "- none\n"
    report += f"""
## INTERPRETATION

- The primary unit is the donor, so the cellular result is not driven by treating tens of thousands of cells as independent biological replicates.
- scDRS/MAGMA enrichment is unsigned with respect to cross-trait allele direction.
- Similar POAG and RE cellular profiles, if present, represent cellular-context convergence despite their negative global LDSC genetic correlation; they do not establish same-direction effects or shared risk alleles.
- AMD remains `AMD_CELLULAR_POWER = LIMITED` because its LDSC h2 Z = 2.4934; a nonsignificant AMD result cannot be interpreted as strong evidence of biological absence.

## NEXT ACTION

{next_action}

Phase 3A-2 was not run automatically.

## FIGURES

"""
    for path in figure_paths:
        report += f"- {path.relative_to(ROOT)}\n"
    report += """
## FILES

- reports/PHASE3A1_BROAD_CELLULAR_ARCHITECTURE.md
- reports/PHASE3A1_DONOR_AWARE_INFERENCE.md
- reports/PHASE3A1_CROSS_TRAIT_CELLULAR_CONVERGENCE.md
- results/phase3a1/scdrs_score_qc.tsv
- results/phase3a1/scdrs_technical_covariate_qc.tsv
- results/phase3a1/donor_class_effects.tsv
- results/phase3a1/scdrs_official_group_analysis.tsv
- results/phase3a1/donor_vs_scdrs_group_concordance.tsv
- results/phase3a1/magma_broad_celltype_enrichment.tsv
- results/phase3a1/BROAD_CELLULAR_CONSENSUS.tsv
- results/phase3a1/leave_one_donor_out.tsv
- results/phase3a1/mapping_selection_sensitivity.tsv
- results/phase3a1/cell_count_balance_sensitivity.tsv
- results/phase3a1/cross_trait_broad_profile_similarity.tsv
- results/phase3a1/state_analysis_eligibility.tsv
- figures/phase3a1/Figure_P14_donor_aware_effect_heatmap.png/.svg
- figures/phase3a1/Figure_P15_donor_level_forest_dot.png/.svg
- figures/phase3a1/Figure_P16_cross_trait_profile_similarity.png/.svg
- figures/phase3a1/Figure_P17_scDRS_vs_MAGMA_concordance.png/.svg

STOP.
"""
    (ROOT / "reports" / "PHASE3A1_BROAD_CELLULAR_ARCHITECTURE.md").write_text(report, encoding="utf-8")

    donor_report = f"""# PHASE 3A-1 DONOR-AWARE INFERENCE

Primary biological unit: donor × broad class. A donor-class pair required at least 30 exact-mapped cells; broad-class results required at least five eligible donors. The estimand was the within-donor class mean norm_score minus the mean of all other eligible broad-class cells in the same donor. Donors were equally weighted.

Primary one-sample donor t tests, exact sign-flip sensitivity P values, leave-one-donor-out estimates, 5,000 donor bootstrap replicates, mapping-selection sensitivities, and cell-count balance sensitivities are in the corresponding TSV files under `results/phase3a1/`.

The official scDRS group association is secondary because its association P values operate on cells and control-score Monte Carlo distributions. It does not replace donor-level inference.

AMD is reported with effect estimates, confidence intervals, donor consistency and FDR, while retaining the locked `AMD_CELLULAR_POWER = LIMITED` interpretation.
"""
    (ROOT / "reports" / "PHASE3A1_DONOR_AWARE_INFERENCE.md").write_text(donor_report, encoding="utf-8")

    sim_lines = []
    for row in similarity.itertuples():
        sim_lines.append(f"- {row.trait_a}–{row.trait_b}: rho = {fmt(row.spearman_rho)}, bootstrap 95% CI = [{fmt(row.bootstrap_CI95_lower)}, {fmt(row.bootstrap_CI95_upper)}], fraction > 0 = {fmt(row.bootstrap_fraction_gt_0)}.")
    convergence_classes = consensus.loc[consensus["cellular_convergence_classification"].isin(["PAN_TRAIT_CONVERGENT", "POAG_RE_CONVERGENT", "AMD_POAG_CONVERGENT", "AMD_RE_CONVERGENT"]), ["broad_class", "cellular_convergence_classification"]].drop_duplicates()
    convergence_text = "; ".join(f"{r.broad_class} ({r.cellular_convergence_classification})" for r in convergence_classes.itertuples()) if len(convergence_classes) else "none"
    cross_report = f"""# PHASE 3A-1 CROSS-TRAIT CELLULAR CONVERGENCE

Cross-trait comparisons used the nine-class donor-aware effect vectors after standardization within each trait. Spearman correlation was primary; cosine similarity and 5,000 donor-bootstrap intervals were secondary.

{chr(10).join(sim_lines)}

Robust cellular convergence classifications: {convergence_text}.

The POAG–RE LDSC correlation is negative (`rg = -0.1732`, FDR `3.68e-09`). Therefore, a similar broad cellular profile is described as convergence on overlapping retinal cellular contexts, not concordant genetic effects, shared risk alleles, or same-direction allele effects. AMD profile similarity is secondary because `AMD_CELLULAR_POWER = LIMITED`.

No fine-state mining was performed. The state-entry gate is recorded in `results/phase3a1/state_analysis_eligibility.tsv`.
"""
    (ROOT / "reports" / "PHASE3A1_CROSS_TRAIT_CELLULAR_CONVERGENCE.md").write_text(cross_report, encoding="utf-8")


def update_status(verdict_name: str, next_action: str) -> None:
    block = f"""

## Phase 3A-1 Update (2026-09-05)

The authorized donor-aware broad polygenic retinal cellular architecture was completed on the full `HRCA_EXACT_MAPPED_HEALTHY_CORE` reference (98,834 exact official-ID mapped healthy cells, 8 donors, 9 broad classes). Primary inference used donor × broad class, with >=30 cells per donor-class pair, >=5 eligible donors per class, one-sample donor tests, exact sign-flip sensitivity, 5,000 donor bootstraps, leave-one-donor-out analysis, mapping-selection sensitivity and cell-count balance sensitivity. Official scDRS group analysis and donor-balanced pseudobulk MAGMA gene-property analysis were secondary/orthogonal.

Phase 3A-1 verdict: `{verdict_name}`. Next action: `{next_action}`. AMD remains `AMD_CELLULAR_POWER = LIMITED`; POAG–RE cellular similarity, if observed, is interpreted as cellular-context convergence despite negative genome-wide LDSC rg, not same-direction genetic effects or shared risk alleles. Phase 3A-2 was not run automatically. See `reports/PHASE3A1_BROAD_CELLULAR_ARCHITECTURE.md`, `reports/PHASE3A1_DONOR_AWARE_INFERENCE.md`, and `reports/PHASE3A1_CROSS_TRAIT_CELLULAR_CONVERGENCE.md`.
"""
    for path in (ROOT / "PROJECT_STATUS.md", ROOT / "reports" / "CURRENT_STATUS.md"):
        with path.open("a", encoding="utf-8") as handle:
            handle.write(block)


def main() -> None:
    require_lock()
    ensure_dirs()
    adata = load_normalized_reference()
    cov = build_covariates(adata)
    dict_gs, overlaps = load_gene_sets(adata)
    score_qc, score_paths, full_paths, donor_delta_by_trait = score_traits(adata, cov, dict_gs, overlaps)
    effects, donor_long, loo = primary_effects(adata.obs, donor_delta_by_trait)
    mapping_detail, mapping_summary = mapping_sensitivity(adata.obs, donor_delta_by_trait, effects)
    balance = cell_balance_sensitivity(adata.obs, effects)
    group = official_group_analysis(adata, full_paths)
    group_concordance, effect_group = donor_group_concordance(effects, group)
    similarity, effect_matrix, standardized = cross_trait_similarity(donor_delta_by_trait)
    del adata
    covar_path, covar = make_pseudobulk_gene_covariates(pd.DataFrame())
    magma = magma_gene_property(covar_path)
    consensus, state, class_labels = consensus_and_state(effects, group, magma, group_concordance, mapping_summary)
    verdict_name, next_action = verdict(consensus)
    figure_paths = make_figures(effects, consensus, similarity, magma)
    write_reports(consensus, state, similarity, group_concordance, score_qc, magma, verdict_name, next_action, figure_paths)
    update_status(verdict_name, next_action)
    software = pd.DataFrame([
        {"software": "scDRS", "version": "1.0.2", "path": str(SCDRS_PYTHON)},
        {"software": "MAGMA", "version": "v1.10 (custom, self-reported)", "path": str(MAGMA)},
        {"software": "Python", "version": sys.version.split()[0], "path": sys.executable},
        {"software": "reference", "version": "HRCA_EXACT_MAPPED_HEALTHY_CORE", "path": str(H5AD)},
    ])
    software.to_csv(OUT / "software_versions.tsv", sep="\t", index=False)
    print(f"PHASE3A1_VERDICT\t{verdict_name}")
    print(f"NEXT_ACTION\t{next_action}")
    print(f"REPORT\t{ROOT / 'reports' / 'PHASE3A1_BROAD_CELLULAR_ARCHITECTURE.md'}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"PHASE3A1_FAILURE\t{type(exc).__name__}: {exc}", file=sys.stderr)
        raise
