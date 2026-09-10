#!/usr/bin/env python3
"""Phase 3A-2 donor-aware Müller-glia state and program analysis.

State discovery is expression-only and is completed before any scDRS score is
loaded. All primary trait tests use donor-level summaries with equal donor
weight. Raw counts are retained in memory for the donor-blocked pseudobulk
marker analysis.
"""
from __future__ import annotations

import gzip
import hashlib
import itertools
import json
import math
import os
import subprocess
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import linear_sum_assignment
from scipy.stats import spearmanr, ttest_1samp
from sklearn.decomposition import NMF
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

if not hasattr(np, "float_"):
    np.float_ = np.float64

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "phase3a2"
LOGS = OUT / "logs"
FIGURES = ROOT / "figures" / "phase3a2"
H5AD = ROOT / "data" / "single_cell" / "HRCA" / "HRCA_EXACT_MAPPED_HEALTHY_CORE.h5ad"
GS = ROOT / "results" / "phase3ar8" / "retinaplex_traits.gs"
P1 = ROOT / "results" / "phase3a1"
LOCK = ROOT / "config" / "phase3a2_analysis_lock.yaml"
SC_PYTHON = Path(os.environ.get("SC_PYTHON", "python"))
TRAITS = ["AMD", "POAG", "RE"]
RESOLUTIONS = [0.2, 0.4, 0.6, 0.8]
SEEDS = [0, 1, 2]
BASE_SEED = 20260905
N_HVG = 2000
N_PCS = 30
N_NEIGHBORS = 15
N_BOOT = 5000
MIN_DONORS = 5
MIN_CELLS_PER_DONOR = 30
MAX_DONOR_FRACTION = 0.50
PROGRAM_K_GRID = [4, 5, 6]
PROGRAM_FEATURES = 1000

def bh(x):
    x = np.asarray(x, dtype=float)
    out = np.full(x.shape, np.nan, dtype=float)
    ok = np.isfinite(x)
    if not ok.any():
        return out
    p = x[ok]
    order = np.argsort(p)
    q = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    tmp = np.empty_like(q)
    tmp[order] = q
    out[ok] = tmp
    return out

def safe_float(x):
    try:
        y = float(x)
        return y if np.isfinite(y) else np.nan
    except Exception:
        return np.nan

def sign_flip_p(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan
    obs = abs(float(values.mean()))
    if len(values) <= 15:
        vals = []
        for signs in itertools.product((-1.0, 1.0), repeat=len(values)):
            vals.append(abs(float(np.mean(values * np.asarray(signs)))))
        return float(np.mean(np.asarray(vals) >= obs - 1e-15))
    rng = np.random.default_rng(BASE_SEED)
    signs = rng.choice([-1.0, 1.0], size=(50000, len(values)))
    null = np.abs((signs * values[None, :]).mean(axis=1))
    return float((np.sum(null >= obs) + 1) / (len(null) + 1))

def bootstrap_mean(values, n=N_BOOT):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan, np.nan, np.nan, np.nan
    rng = np.random.default_rng(BASE_SEED)
    boot = rng.choice(values, size=(n, len(values)), replace=True).mean(axis=1)
    return (float(np.median(boot)), float(np.quantile(boot, 0.025)),
            float(np.quantile(boot, 0.975)), float(np.mean(boot > 0)))

def write_tsv(df, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, na_rep="NA")

def load_muller():
    import anndata as ad
    adata = ad.read_h5ad(H5AD, backed="r")
    obs = adata.obs.copy()
    mask = obs["harmonized_broad_class"].astype(str).eq("Muller glia").to_numpy()
    m = adata[mask].to_memory()
    try:
        adata.file.close()
    except Exception:
        pass
    if m.n_obs != 9534:
        warnings.warn("Müller cell count differs from Phase 3A-1 audit: %s" % m.n_obs)
    if not sparse.issparse(m.X):
        raise RuntimeError("Müller primary matrix is not sparse")
    m.X = m.X.tocsr().astype(np.float32)
    m.obs_names = m.obs_names.astype(str)
    m.obs["donor_id"] = m.obs["donor_id"].astype(str)
    m.obs["sample"] = m.obs["sample"].astype(str)
    return m

def audit_population(m):
    obs = m.obs.copy()
    rows = []
    for donor, g in obs.groupby("donor_id", sort=True):
        rows.append({
            "record_type": "donor", "donor": donor,
            "sample": ";".join(sorted(g["sample"].unique())), "cells": len(g),
            "n_donors": 1, "n_samples": g["sample"].nunique(),
            "region": ";".join(sorted(g["region"].dropna().astype(str).unique())) if "region" in g else "NA",
            "age": ";".join(sorted(g["age"].dropna().astype(str).unique())) if "age" in g else "NA",
            "mean_total_counts": pd.to_numeric(g["total_counts"], errors="coerce").mean(),
            "mean_n_genes": pd.to_numeric(g["n_genes"], errors="coerce").mean(),
            "mean_pct_counts_mt": pd.to_numeric(g["pct_counts_mt"], errors="coerce").mean(),
            "inclusion": "healthy_exact_official_ID_Muller_glia",
        })
    rows.append({
        "record_type": "TOTAL", "donor": "ALL",
        "sample": ";".join(sorted(obs["sample"].unique())), "cells": len(obs),
        "n_donors": obs["donor_id"].nunique(), "n_samples": obs["sample"].nunique(),
        "region": ";".join(sorted(obs["region"].dropna().astype(str).unique())) if "region" in obs else "NA",
        "age": ";".join(sorted(obs["age"].dropna().astype(str).unique())) if "age" in obs else "NA",
        "mean_total_counts": pd.to_numeric(obs["total_counts"], errors="coerce").mean(),
        "mean_n_genes": pd.to_numeric(obs["n_genes"], errors="coerce").mean(),
        "mean_pct_counts_mt": pd.to_numeric(obs["pct_counts_mt"], errors="coerce").mean(),
        "inclusion": "healthy_exact_official_ID_Muller_glia",
    })
    write_tsv(pd.DataFrame(rows), OUT / "muller_population_audit.tsv")
    return {
        "cells": int(len(obs)), "donors": int(obs["donor_id"].nunique()),
        "samples": int(obs["sample"].nunique()),
        "donor_counts": obs["donor_id"].value_counts().sort_index().to_dict(),
        "region_values": sorted(obs["region"].dropna().astype(str).unique()) if "region" in obs else [],
        "age_values": sorted(obs["age"].dropna().astype(str).unique()) if "age" in obs else [],
    }

def annotation_inventory(m):
    obs = m.obs.copy()
    fields = []
    for c in ["majorclass", "cell_type", "author_cell_type", "subtype", "cluster"]:
        if c in obs.columns:
            fields.append(c)
    for c in obs.columns:
        cl = c.lower()
        if c not in fields and any(k in cl for k in ["fine", "subtype", "cluster", "annotation"]):
            if obs[c].nunique(dropna=True) <= 100:
                fields.append(c)
    rows = []
    for field in fields:
        values = obs[field].astype(object).where(obs[field].notna(), "NA").astype(str)
        for value, g in obs.groupby(values, sort=True):
            donor_counts = g["donor_id"].value_counts()
            donor_total = obs["donor_id"].value_counts()
            frac = donor_counts / donor_total.reindex(donor_counts.index)
            enough_donors = int((donor_counts >= MIN_CELLS_PER_DONOR).sum())
            n_donors = int(g["donor_id"].nunique())
            ontology_broad = field in {"majorclass", "cell_type", "author_cell_type"} and value in {"MG", "CL:0000636", "Muller glia"}
            ready = (n_donors >= MIN_DONORS and enough_donors >= MIN_DONORS and
                     float(frac.max()) <= MAX_DONOR_FRACTION and not ontology_broad)
            rows.append({
                "field": field, "value": value, "cells": int(len(g)),
                "donors": n_donors, "samples": int(g["sample"].nunique()),
                "minimum_cells_per_donor": int(donor_counts.min()),
                "donors_with_30_cells": enough_donors,
                "maximum_donor_fraction_of_label": float(frac.max()),
                "official_state_ready": "OFFICIAL_STATE_READY" if ready else "OFFICIAL_STATE_INSUFFICIENT",
                "ontology_broad_label_excluded": bool(ontology_broad),
            })
    inv = pd.DataFrame(rows)
    write_tsv(inv, OUT / "muller_official_annotation_inventory.tsv")
    usable = inv[(inv["official_state_ready"] == "OFFICIAL_STATE_READY") &
                 (~inv["ontology_broad_label_excluded"])]
    official_ready = int(usable["value"].nunique()) >= 2
    return inv, official_ready

def write_lock(m, inv):
    digest = hashlib.sha256(H5AD.read_bytes()).hexdigest()
    if LOCK.exists():
        text = LOCK.read_text(encoding="utf-8")
        required = ["phase: Phase3A-2", "status: LOCKED",
                    "trait_blind_state_discovery: true", "primary_lineage: Muller glia", digest]
        for token in required:
            if token not in text:
                raise RuntimeError("existing Phase 3A-2 lock is incompatible: %s" % token)
        return
    LOCK.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "phase: Phase3A-2", "status: LOCKED", "authorized_date: 2026-09-06",
        "question: convergent Muller-glia transcriptional state and program vulnerability across AMD POAG and RE",
        "primary_lineage: Muller glia", "trait_blind_state_discovery: true",
        "primary_reference:",
        "  h5ad: data/single_cell/HRCA/HRCA_EXACT_MAPPED_HEALTHY_CORE.h5ad",
        "  h5ad_sha256: %s" % digest,
        "  inclusion: healthy_or_normal_exact_official_ID_harmonized_broad_class_Muller_glia_all_8_donors",
        "  excluded: disease_model_transferred_labels_unmapped_cells_other_lineages",
        "  cells: %d" % m.n_obs, "  donors: %d" % m.obs["donor_id"].nunique(),
        "official_annotation_policy:",
        "  inspect_fields: [majorclass, cell_type, author_cell_type, subtype, cluster, fine_annotation]",
        "  official_state_ready: >=5 donors and >=30 cells in >=5 donors and max donor fraction <=0.50",
        "  ontology_rule: CL:0000636 is broad Muller cell label and cannot be a fine state",
        "state_discovery:",
        "  mode: blind_de_novo_if_two_replicated_official_substates_absent",
        "  normalized_counts: library_size_1e4_then_log1p", "  raw_counts_retained: true",
        "  hvg_rule: Muller_specific_cell_ranger_HVG_top_2000_excluding_union_of_AMD_POAG_RE_top_1000_and_MT_ribosomal_MALAT1",
        "  pca_dimensions: 30", "  neighbors: 15", "  leiden_resolution_grid: [0.2, 0.4, 0.6, 0.8]",
        "  seeds: [0, 1, 2]", "  subsampling: repeated_80_percent_and_donor_balanced",
        "  stability: median_ARI_NMI_donor_representation_smallest_cluster",
        "  primary_resolution_selection: stability_donor_representation_tiny_cluster_avoidance_only",
        "state_eligibility:", "  minimum_donors: 5", "  minimum_cells_in_donor: 30",
        "  max_donor_fraction: 0.50",
        "  tier1_requires: global_FDR_lt_0.05_positive_ge75_percent_direction_LODO_no_composition_failure_stability_mapping_sensitivity",
        "donor_effect:", "  definition: mean_state_norm_score_minus_mean_other_eligible_Muller_cells_same_donor",
        "  donor_weighting: equal", "  primary_test: donor_level_one_sample_t_test",
        "  exact_sign_flip: two_sided_separately_reported", "  bootstrap_replicates: 5000",
        "  bootstrap_seed: 20260905",
        "multiple_testing:", "  global_family: 3_traits_times_K_eligible_states",
        "  method: Benjamini_Hochberg", "  secondary: within_trait_FDR",
        "pathways:", "  databases: [GO_Biological_Process, Reactome, Hallmark]",
        "  marker_selection: donor_blocked_edgeR_quasi_likelihood_positive_marker_FDR_and_donor_consistency",
        "programs:", "  method: unsupervised_NMF_consensus_stability_reconstruction",
        "  K_candidates: [4, 5, 6]", "  trait_association: donor_level_program_activity_vs_donor_mean_scDRS",
        "sensitivities:", "  cell_balance_cap: 500",
        "  mapping_selection: [top75_percent_mapping, exclude_lowest_mapping, median_donor_score]",
        "  disease_gene_exclusion_full_HVG: true", "  clustering_resolution: true",
        "  age_region_secondary_if_estimable: true",
        "cross_trait:", "  PAN_TRAIT_CONVERGENT: robust_positive_AMD_POAG_RE",
        "  pairwise_convergent: robust_positive_in_any_two_traits",
        "  trait_specific: robust_positive_in_one_trait_only",
        "prohibited: [score_guided_state_discovery, amacrine_analysis, Phase3B_automatic_start, cell_level_primary_inference]",
    ]
    LOCK.write_text("\n".join(lines) + "\n", encoding="utf-8")

def load_gene_sets():
    sys.path.insert(0, str(ROOT / ".vendor" / "phase3a2"))
    import scdrs
    gs = scdrs.util.load_gs(str(GS))
    return {trait: set(map(str, gs[trait][0])) for trait in TRAITS}

def prepare_features(m, gene_sets):
    import scanpy as sc
    raw = m.X.tocsr().copy().astype(np.float32)
    sc.pp.normalize_total(m, target_sum=1e4, inplace=True)
    sc.pp.log1p(m)
    # Scanpy's cell_ranger binning fails on the repeated zero-expression bin
    # edges in this sparse subset. Reproduce its mean/dispersion ranking with
    # duplicate-safe quantile bins, retaining the locked cell_ranger-like rule.
    x = m.X.tocsr().astype(np.float64)
    means = np.asarray(x.mean(axis=0)).ravel()
    vars_ = np.asarray(x.power(2).mean(axis=0)).ravel() - means ** 2
    vars_ = np.maximum(vars_, 0)
    disp = vars_ / np.maximum(means, 1e-12)
    try:
        bins = pd.qcut(means, q=20, labels=False, duplicates="drop")
        norm_disp = np.full(len(means), np.nan)
        for b in np.unique(bins.dropna() if hasattr(bins, "dropna") else bins[np.isfinite(bins)]):
            mask = np.asarray(bins == b)
            med = np.nanmedian(disp[mask])
            norm_disp[mask] = disp[mask] / med if med > 0 else 0
    except Exception:
        norm_disp = disp.copy()
    norm_disp[~np.isfinite(norm_disp)] = -np.inf
    top = np.argsort(norm_disp)[::-1][:min(N_HVG, len(norm_disp))]
    hv = np.zeros(len(norm_disp), dtype=bool)
    hv[top] = True
    m.var["means"] = means
    m.var["variances"] = vars_
    m.var["dispersions"] = disp
    m.var["dispersions_norm"] = norm_disp
    m.var["highly_variable"] = hv
    gene_names = np.asarray(m.var_names.astype(str), dtype=str)
    upper = np.char.upper(gene_names)
    disease_union = set().union(*gene_sets.values())
    disease_flag = np.array([g in disease_union for g in gene_names], dtype=bool)
    mt_flag = np.char.startswith(upper, "MT-")
    ribo_flag = np.char.startswith(upper, "RPS") | np.char.startswith(upper, "RPL")
    tech_flag = gene_names == "MALAT1"
    hvg_flag = m.var["highly_variable"].to_numpy(dtype=bool)
    keep = hvg_flag & ~disease_flag & ~mt_flag & ~ribo_flag & ~tech_flag
    if keep.sum() < 500:
        raise RuntimeError("too few state-discovery genes after locked exclusion")
    reason = np.full(len(gene_names), "not_HVG", dtype=object)
    reason[hvg_flag] = "HVG_retained"
    reason[hvg_flag & disease_flag] = "HVG_excluded_disease_gene"
    reason[hvg_flag & mt_flag] = "HVG_excluded_mitochondrial"
    reason[hvg_flag & ribo_flag] = "HVG_excluded_ribosomal"
    reason[hvg_flag & tech_flag] = "HVG_excluded_technical"
    audit = pd.DataFrame({
        "gene": gene_names, "hvg_candidate": hvg_flag,
        "excluded_disease_union": disease_flag, "excluded_mitochondrial": mt_flag,
        "excluded_ribosomal": ribo_flag, "excluded_technical": tech_flag,
        "retained_for_state_discovery": keep, "reason": reason,
    })
    for trait in TRAITS:
        audit["in_%s_top1000" % trait] = np.array([g in gene_sets[trait] for g in gene_names], dtype=bool)
    write_tsv(audit, OUT / "state_discovery_gene_audit.tsv")
    return m, raw, gene_names[keep].tolist(), audit

def cluster_from_pca(pca, idx, resolution, seed):
    import anndata as ad
    import scanpy as sc
    x = np.asarray(pca[idx], dtype=np.float32)
    a = ad.AnnData(x)
    a.obs_names = np.asarray(idx).astype(str)
    sc.pp.neighbors(a, n_neighbors=N_NEIGHBORS, n_pcs=min(N_PCS, x.shape[1]), random_state=BASE_SEED)
    sc.tl.leiden(a, resolution=resolution, random_state=seed, flavor="igraph", directed=False)
    return a.obs["leiden"].astype(str).to_numpy()

def donor_metrics(labels, donors):
    labs = np.asarray(labels).astype(str)
    d = np.asarray(donors).astype(str)
    counts = pd.crosstab(labs, d)
    donor_rep = int((counts > 0).sum(axis=1).min()) if len(counts) else 0
    fractions = counts.div(counts.sum(axis=1), axis=0)
    max_frac = float(fractions.to_numpy().max()) if fractions.size else np.nan
    return int(len(np.unique(labs))), int(pd.Series(labs).value_counts().min()), donor_rep, max_frac

def discover_states(m, features):
    import scanpy as sc
    proc = m[:, features].copy()
    sc.pp.scale(proc, max_value=10)
    sc.tl.pca(proc, n_comps=N_PCS, svd_solver="arpack", random_state=BASE_SEED)
    sc.pp.neighbors(proc, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, random_state=BASE_SEED)
    sc.tl.umap(proc, random_state=BASE_SEED, min_dist=0.3)
    m.obsm["X_umap"] = np.asarray(proc.obsm["X_umap"], dtype=np.float32)
    m.obsm["X_pca_state"] = np.asarray(proc.obsm["X_pca"], dtype=np.float32)
    donors = m.obs["donor_id"].astype(str).to_numpy()
    rng = np.random.default_rng(BASE_SEED)
    records, primary_labels, resolution_summary = [], {}, []
    all_idx = np.arange(m.n_obs)
    for res in RESOLUTIONS:
        ref = cluster_from_pca(proc.obsm["X_pca"], all_idx, res, SEEDS[0])
        primary_labels[res] = ref
        nc, sm, dr, mf = donor_metrics(ref, donors)
        records.append({"resolution": res, "mode": "full_seed", "replicate": 0,
                        "ARI_vs_reference": 1.0, "NMI_vs_reference": 1.0,
                        "cluster_count": nc, "smallest_cluster_size": sm,
                        "min_donor_representation": dr, "max_donor_fraction": mf})
        for rep, seed in enumerate(SEEDS[1:], start=1):
            lab = cluster_from_pca(proc.obsm["X_pca"], all_idx, res, seed)
            nc, sm, dr, mf = donor_metrics(lab, donors)
            records.append({"resolution": res, "mode": "full_seed", "replicate": rep,
                            "ARI_vs_reference": adjusted_rand_score(ref, lab),
                            "NMI_vs_reference": normalized_mutual_info_score(ref, lab),
                            "cluster_count": nc, "smallest_cluster_size": sm,
                            "min_donor_representation": dr, "max_donor_fraction": mf})
        for rep in range(3):
            idx = np.sort(rng.choice(m.n_obs, size=int(m.n_obs * 0.8), replace=False))
            lab = cluster_from_pca(proc.obsm["X_pca"], idx, res, int(100 + rep))
            nc, sm, dr, mf = donor_metrics(lab, donors[idx])
            records.append({"resolution": res, "mode": "80_percent_subsample", "replicate": rep,
                            "ARI_vs_reference": adjusted_rand_score(ref[idx], lab),
                            "NMI_vs_reference": normalized_mutual_info_score(ref[idx], lab),
                            "cluster_count": nc, "smallest_cluster_size": sm,
                            "min_donor_representation": dr, "max_donor_fraction": mf})
        donor_counts = pd.Series(donors).value_counts().sort_index()
        n_each = int(donor_counts.min() * 0.8)
        for rep in range(2):
            idxs = []
            for donor in sorted(donor_counts.index):
                pool = np.where(donors == donor)[0]
                idxs.extend(rng.choice(pool, size=n_each, replace=False))
            idx = np.sort(np.asarray(idxs, dtype=int))
            lab = cluster_from_pca(proc.obsm["X_pca"], idx, res, int(200 + rep))
            nc, sm, dr, mf = donor_metrics(lab, donors[idx])
            records.append({"resolution": res, "mode": "donor_balanced_80_percent", "replicate": rep,
                            "ARI_vs_reference": adjusted_rand_score(ref[idx], lab),
                            "NMI_vs_reference": normalized_mutual_info_score(ref[idx], lab),
                            "cluster_count": nc, "smallest_cluster_size": sm,
                            "min_donor_representation": dr, "max_donor_fraction": mf})
        sub = [r for r in records if r["resolution"] == res and
               not (r["mode"] == "full_seed" and r["replicate"] == 0)]
        med_ari = float(np.nanmedian([r["ARI_vs_reference"] for r in sub]))
        med_nmi = float(np.nanmedian([r["NMI_vs_reference"] for r in sub]))
        min_cluster = int(min(r["smallest_cluster_size"] for r in sub))
        min_rep = int(min(r["min_donor_representation"] for r in sub))
        max_dom = float(max(r["max_donor_fraction"] for r in sub))
        eligible = bool(med_ari >= 0.70 and med_nmi >= 0.70 and min_cluster >= 100 and
                        min_rep >= MIN_DONORS and max_dom <= MAX_DONOR_FRACTION)
        resolution_summary.append({"resolution": res, "median_ARI": med_ari, "median_NMI": med_nmi,
                                   "minimum_cluster_size": min_cluster, "minimum_donor_representation": min_rep,
                                   "maximum_donor_fraction": max_dom, "stability_eligible": eligible})
    stab = pd.DataFrame(records)
    write_tsv(stab, OUT / "muller_cluster_stability.tsv")
    summ = pd.DataFrame(resolution_summary)
    write_tsv(summ, OUT / "muller_resolution_selection.tsv")
    candidates = summ[summ["stability_eligible"]]
    if len(candidates):
        chosen = float(candidates.sort_values(["median_ARI", "median_NMI", "resolution"],
                                              ascending=[False, False, True]).iloc[0]["resolution"])
        selection = "highest stability among predefined resolutions meeting locked criteria"
    else:
        chosen = float(summ.sort_values(["minimum_donor_representation", "minimum_cluster_size",
                                         "median_ARI", "resolution"],
                                        ascending=[False, False, False, True]).iloc[0]["resolution"])
        selection = "best donor representation and smallest-cluster profile among predefined resolutions; no resolution met all locked criteria"
    labels = primary_labels[chosen]
    state = np.array(["MG_state_%02d" % (i + 1) for i in pd.Categorical(labels).codes], dtype=object)
    m.obs["muller_state"] = state
    write_tsv(pd.DataFrame([{"primary_resolution": chosen, "selection_rule": selection,
                             "n_states": int(len(np.unique(state))),
                             "state_discovery_gene_count": len(features)}]),
              OUT / "muller_primary_resolution.tsv")
    return m, proc, stab, summ, chosen, features, labels

def state_composition(m):
    obs = m.obs.copy()
    obs["state"] = obs["muller_state"].astype(str)
    obs["donor_id"] = obs["donor_id"].astype(str)
    state_tot = obs["state"].value_counts()
    donor_tot = obs["donor_id"].value_counts()
    rows = []
    for (state, donor), g in obs.groupby(["state", "donor_id"], sort=True):
        rows.append({
            "state": state, "donor": donor, "cells": int(len(g)),
            "state_total_cells": int(state_tot[state]),
            "donor_total_muller_cells": int(donor_tot[donor]),
            "donors_in_state": int(obs.loc[obs["state"] == state, "donor_id"].nunique()),
            "donor_fraction_of_state": float(len(g) / state_tot[state]),
            "state_fraction_of_donor": float(len(g) / donor_tot[donor]),
            "age": ";".join(sorted(g["age"].dropna().astype(str).unique())) if "age" in g else "NA",
            "region": ";".join(sorted(g["region"].dropna().astype(str).unique())) if "region" in g else "NA",
            "mean_total_counts": safe_float(pd.to_numeric(g["total_counts"], errors="coerce").mean()),
            "mean_n_genes": safe_float(pd.to_numeric(g["n_genes"], errors="coerce").mean()),
            "mean_pct_counts_mt": safe_float(pd.to_numeric(g["pct_counts_mt"], errors="coerce").mean()),
            "state_donor_status": "INSUFFICIENT_STATE_CELLS" if len(g) < MIN_CELLS_PER_DONOR else "PASS",
        })
    comp = pd.DataFrame(rows)
    eligible_rows = []
    for state, g in comp.groupby("state", sort=True):
        good = g[g["cells"] >= MIN_CELLS_PER_DONOR]
        eligible = bool(len(good) >= MIN_DONORS and float(g["donor_fraction_of_state"].max()) <= MAX_DONOR_FRACTION)
        eligible_rows.append({"state": state, "cells": int(g["state_total_cells"].iloc[0]),
                              "donors": int(g["donors_in_state"].iloc[0]),
                              "donors_with_30_cells": int(len(good)),
                              "max_donor_fraction": float(g["donor_fraction_of_state"].max()),
                              "donor_dominated": bool(g["donor_fraction_of_state"].max() > MAX_DONOR_FRACTION),
                              "state_eligible": eligible})
    eligibility = pd.DataFrame(eligible_rows)
    comp = comp.merge(eligibility[["state", "state_eligible"]], on="state", how="left")
    write_tsv(comp, OUT / "muller_state_composition.tsv")
    write_tsv(eligibility, OUT / "muller_state_eligibility.tsv")
    return comp, eligibility

def load_scores(m):
    scores = {}
    for trait in TRAITS:
        df = pd.read_csv(P1 / ("%s.score.gz" % trait), sep="\t", index_col=0)
        if not m.obs_names.isin(df.index).all():
            raise RuntimeError("score index does not cover all Müller cells for %s" % trait)
        aligned = df.reindex(m.obs_names)
        scores[trait] = aligned["norm_score"].to_numpy(dtype=float)
        m.obs["%s_norm_score" % trait] = scores[trait]
    return scores

def heterogeneity():
    x = pd.read_csv(P1 / "scdrs_official_group_analysis.tsv", sep="\t")
    x = x[(x["broad_class"] == "Muller glia") & x["trait"].isin(TRAITS)].copy()
    out = x[["trait", "broad_class", "n_cell", "hetero_mcp", "hetero_mcz"]].copy()
    out["analysis_role"] = "secondary_official_within_Muller_heterogeneity"
    write_tsv(out, OUT / "muller_scdrs_heterogeneity.tsv")
    return out

def compute_donor_deltas(m, scores, eligibility):
    states = sorted(eligibility.loc[eligibility["state_eligible"], "state"].astype(str))
    donors = sorted(m.obs["donor_id"].astype(str).unique())
    state_arr = m.obs["muller_state"].astype(str).to_numpy()
    donor_arr = m.obs["donor_id"].astype(str).to_numpy()
    rows = []
    for trait, vals in scores.items():
        for state in states:
            deltas = []
            for donor in donors:
                sm = (state_arr == state) & (donor_arr == donor)
                other = (state_arr != state) & (donor_arr == donor)
                if sm.sum() < MIN_CELLS_PER_DONOR or other.sum() == 0:
                    continue
                delta = float(np.mean(vals[sm]) - np.mean(vals[other]))
                deltas.append({"trait": trait, "state": state, "donor": donor,
                                "state_cells": int(sm.sum()), "other_cells": int(other.sum()),
                                "donor_delta": delta})
            rows.extend(deltas)
    donor_df = pd.DataFrame(rows, columns=["trait", "state", "donor", "state_cells", "other_cells", "donor_delta"])
    return donor_df

def compute_effect_table(donor_df, eligibility):
    if len(donor_df) == 0:
        cols = ["trait", "state", "eligible_donors", "cells", "mean_donor_delta", "SE",
                "CI95_low_t", "CI95_high_t", "median_donor_delta", "donors_positive",
                "donors_negative", "direction_fraction_positive", "primary_P_t_test",
                "sign_flip_P", "bootstrap_median", "bootstrap_CI95_low", "bootstrap_CI95_high",
                "bootstrap_fraction_gt0", "state_eligible", "donor_composition_status",
                "global_FDR", "within_trait_FDR", "direction_ge75pct"]
        out = pd.DataFrame(columns=cols)
        write_tsv(out, OUT / "muller_state_trait_effects.tsv")
        return out
    records = []
    for (trait, state), g in donor_df.groupby(["trait", "state"], sort=True):
        x = g["donor_delta"].to_numpy(dtype=float)
        n = len(x)
        mean = float(np.mean(x)) if n else np.nan
        se = float(np.std(x, ddof=1) / np.sqrt(n)) if n > 1 else np.nan
        ci_low = mean - 1.96 * se if np.isfinite(se) else np.nan
        ci_high = mean + 1.96 * se if np.isfinite(se) else np.nan
        tp = float(ttest_1samp(x, 0).pvalue) if n >= 2 and np.std(x, ddof=1) > 0 else np.nan
        med, boot_low, boot_high, boot_pos = bootstrap_mean(x)
        elig = eligibility.loc[eligibility["state"] == state].iloc[0]
        records.append({
            "trait": trait, "state": state, "eligible_donors": n,
            "cells": int(g["state_cells"].sum()), "mean_donor_delta": mean, "SE": se,
            "CI95_low_t": ci_low, "CI95_high_t": ci_high, "median_donor_delta": float(np.median(x)),
            "donors_positive": int((x > 0).sum()), "donors_negative": int((x < 0).sum()),
            "direction_fraction_positive": float((x > 0).mean()),
            "primary_P_t_test": tp, "sign_flip_P": sign_flip_p(x),
            "bootstrap_median": med, "bootstrap_CI95_low": boot_low,
            "bootstrap_CI95_high": boot_high, "bootstrap_fraction_gt0": boot_pos,
            "state_eligible": bool(elig["state_eligible"]),
            "donor_composition_status": "DONOR_DOMINATED" if bool(elig["donor_dominated"]) else "PASS",
        })
    eff = pd.DataFrame(records)
    if len(eff):
        eff["global_FDR"] = bh(eff["primary_P_t_test"].to_numpy())
        eff["within_trait_FDR"] = eff.groupby("trait", sort=False)["primary_P_t_test"].transform(
            lambda x: bh(x.to_numpy()))
        eff["direction_ge75pct"] = eff["direction_fraction_positive"] >= 0.75
    else:
        eff["global_FDR"] = []
        eff["within_trait_FDR"] = []
        eff["direction_ge75pct"] = []
    write_tsv(eff, OUT / "muller_state_trait_effects.tsv")
    return eff

def leave_one_out(donor_df):
    if len(donor_df) == 0:
        out = pd.DataFrame(columns=["trait", "state", "removed_donor", "remaining_donors",
                                    "leave_one_out_mean_delta", "direction_after_removal", "status"])
        write_tsv(out, OUT / "muller_leave_one_donor_out.tsv")
        return out
    rows = []
    for (trait, state), g in donor_df.groupby(["trait", "state"], sort=True):
        for donor in sorted(g["donor"].unique()):
            x = g.loc[g["donor"] != donor, "donor_delta"].to_numpy(dtype=float)
            value = float(np.mean(x)) if len(x) else np.nan
            rows.append({"trait": trait, "state": state, "removed_donor": donor,
                         "remaining_donors": len(x), "leave_one_out_mean_delta": value,
                         "direction_after_removal": "POSITIVE" if value > 0 else "NEGATIVE" if value < 0 else "ZERO",
                         "status": "PASS" if np.isfinite(value) and value > 0 else "DIRECTION_REVERSED_OR_NOT_ESTIMABLE"})
    out = pd.DataFrame(rows)
    write_tsv(out, OUT / "muller_leave_one_donor_out.tsv")
    return out

def mapping_sensitivity(m, scores, effects):
    donors = sorted(m.obs["donor_id"].astype(str).unique())
    states = sorted(effects["state"].unique()) if len(effects) else []
    state_arr = m.obs["muller_state"].astype(str).to_numpy()
    donor_arr = m.obs["donor_id"].astype(str).to_numpy()
    rows = []
    for trait, vals in scores.items():
        for state in states:
            primary = effects.loc[(effects["trait"] == trait) & (effects["state"] == state), "mean_donor_delta"]
            primary = float(primary.iloc[0]) if len(primary) else np.nan
            for label, keep, rule in [
                ("A_TOP75_MAPPING", [d for d in donors if d in MAPPING_TOP75_DONORS], "top_75_percent_mapping"),
                ("B_EXCLUDE_LOWEST_MAPPING", [d for d in donors if d != "D027-13"], "exclude_D027-13"),
                ("C_MEDIAN_DONOR_SCORE", donors, "median_instead_of_mean"),
            ]:
                ds = []
                for donor in keep:
                    sm = (state_arr == state) & (donor_arr == donor)
                    other = (state_arr != state) & (donor_arr == donor)
                    if sm.sum() < MIN_CELLS_PER_DONOR or not other.any():
                        continue
                    if label == "C_MEDIAN_DONOR_SCORE":
                        ds.append(float(np.median(vals[sm]) - np.median(vals[other])))
                    else:
                        ds.append(float(np.mean(vals[sm]) - np.mean(vals[other])))
                value = float(np.mean(ds)) if len(ds) else np.nan
                rows.append({"trait": trait, "state": state, "sensitivity": label,
                             "rule": rule, "donors": len(ds), "effect": value,
                             "primary_effect": primary,
                             "direction_reversal": bool(np.isfinite(value) and np.isfinite(primary) and value * primary < 0),
                             "status": "PASS" if len(ds) >= MIN_DONORS and np.isfinite(value) else "NOT_ESTIMABLE"})
    out = pd.DataFrame(rows, columns=["trait", "state", "sensitivity", "rule", "donors", "effect",
                                      "primary_effect", "direction_reversal", "status"])
    write_tsv(out, OUT / "muller_mapping_sensitivity.tsv")
    return out

def balance_sensitivity(m, scores, eligibility):
    cap = 500
    rng = np.random.default_rng(BASE_SEED)
    state_arr = m.obs["muller_state"].astype(str).to_numpy()
    donor_arr = m.obs["donor_id"].astype(str).to_numpy()
    eligible = sorted(eligibility.loc[eligibility["state_eligible"], "state"].astype(str))
    rows = []
    for trait, vals in scores.items():
        for state in eligible:
            ds = []
            for donor in sorted(np.unique(donor_arr)):
                sm_idx = np.where((state_arr == state) & (donor_arr == donor))[0]
                ot_idx = np.where((state_arr != state) & (donor_arr == donor))[0]
                if len(sm_idx) < MIN_CELLS_PER_DONOR or len(ot_idx) == 0:
                    continue
                ns = min(cap, len(sm_idx), len(ot_idx))
                sm_take = rng.choice(sm_idx, size=ns, replace=False)
                ot_take = rng.choice(ot_idx, size=ns, replace=False)
                ds.append(float(np.mean(vals[sm_take]) - np.mean(vals[ot_take])))
            rows.append({"trait": trait, "state": state, "cap_per_donor_state_vs_other": cap,
                         "eligible_donors": len(ds), "balanced_mean_delta": np.mean(ds) if ds else np.nan,
                         "direction": "POSITIVE" if ds and np.mean(ds) > 0 else "NEGATIVE_OR_NA"})
    out = pd.DataFrame(rows, columns=["trait", "state", "cap_per_donor_state_vs_other",
                                      "eligible_donors", "balanced_mean_delta", "direction"])
    write_tsv(out, OUT / "muller_cell_count_balance_sensitivity.tsv")
    return out

def cluster_sensitivity(m, proc, features, chosen, labels, stab):
    import scanpy as sc
    pca = proc.obsm["X_pca"]
    primary = np.asarray(labels).astype(str)
    rows = []
    for res in RESOLUTIONS:
        lab = cluster_from_pca(pca, np.arange(m.n_obs), res, 0)
        rows.append({"comparison": "resolution_vs_primary", "resolution": res,
                     "primary_resolution": chosen, "ARI": adjusted_rand_score(primary, lab),
                     "NMI": normalized_mutual_info_score(primary, lab),
                     "cluster_count": len(np.unique(lab))})
    all_hvg = m.var["highly_variable"].to_numpy(dtype=bool)
    full_genes = np.asarray(m.var_names.astype(str))[all_hvg].tolist()
    full_hvg = m[:, full_genes].copy()
    sc.pp.scale(full_hvg, max_value=10)
    sc.tl.pca(full_hvg, n_comps=N_PCS, svd_solver="arpack", random_state=BASE_SEED)
    full_lab = cluster_from_pca(full_hvg.obsm["X_pca"], np.arange(m.n_obs), chosen, 0)
    rows.append({"comparison": "full_HVG_disease_gene_inclusion", "resolution": chosen,
                 "primary_resolution": chosen, "ARI": adjusted_rand_score(primary, full_lab),
                 "NMI": normalized_mutual_info_score(primary, full_lab),
                 "cluster_count": len(np.unique(full_lab))})
    out = pd.DataFrame(rows)
    write_tsv(out, OUT / "muller_cluster_sensitivity.tsv")
    return out

def write_pseudobulk(m, raw):
    states = m.obs["muller_state"].astype(str).to_numpy()
    donors = m.obs["donor_id"].astype(str).to_numpy()
    groups = [(d, s) for d in sorted(np.unique(donors)) for s in sorted(np.unique(states))]
    mats, meta = [], []
    for donor, state in groups:
        idx = np.where((donors == donor) & (states == state))[0]
        if len(idx) == 0:
            continue
        mats.append(np.asarray(raw[idx].sum(axis=0)).ravel())
        meta.append({"sample_id": "%s__%s" % (donor, state), "donor": donor,
                     "state": state, "n_cells": len(idx)})
    count_mat = np.asarray(mats, dtype=np.int64).T
    genes = np.asarray(m.var_names.astype(str))
    path = OUT / "muller_pseudobulk_counts.tsv.gz"
    pd.DataFrame(count_mat, index=genes, columns=[x["sample_id"] for x in meta]).to_csv(path, sep="\t", compression="gzip")
    meta_df = pd.DataFrame(meta)
    write_tsv(meta_df, OUT / "muller_pseudobulk_metadata.tsv")
    return path, OUT / "muller_pseudobulk_metadata.tsv"

def run_r(script, args, log_name):
    cmd = ["Rscript", str(ROOT / "scripts" / "08_scrna_state" / script)] + [str(x) for x in args]
    log = LOGS / log_name
    with log.open("w", encoding="utf-8") as h:
        h.write("COMMAND\t" + " ".join(cmd) + "\n")
        p = subprocess.run(cmd, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        h.write(p.stdout or "")
        h.write("\nEXIT_CODE\t%d\n" % p.returncode)
    if p.returncode != 0:
        raise RuntimeError("R analysis failed; see %s" % log)
    return p

def markers_and_pathways(m, raw, eligibility):
    count_path, meta_path = write_pseudobulk(m, raw)
    marker_path = OUT / "muller_state_markers.tsv"
    if not bool(eligibility["state_eligible"].any()):
        marker_cols = ["gene", "state", "logFC", "SE", "statistic", "P", "FDR",
                       "donor_consistency", "n_eligible_donors", "n_state_cells", "n_other_cells"]
        pathway_cols = ["state", "pathway_source", "pathway", "overlap", "gene_set_size",
                        "selected_gene_count", "background_gene_count", "P", "selection_rule",
                        "leading_genes", "FDR"]
        write_tsv(pd.DataFrame(columns=marker_cols), marker_path)
        write_tsv(pd.DataFrame(columns=pathway_cols), OUT / "muller_state_pathways.tsv")
        return pd.DataFrame(columns=marker_cols), pd.DataFrame(columns=pathway_cols)
    run_r("phase3a2_edgeR_markers.R", [count_path, meta_path, marker_path], "edgeR_markers.log")
    pathway_path = OUT / "muller_state_pathways.tsv"
    run_r("phase3a2_pathways.R", [marker_path, pathway_path], "pathways.log")
    markers = pd.read_csv(marker_path, sep="\t")
    pathways = pd.read_csv(pathway_path, sep="\t")
    return markers, pathways

def align_component_similarity(a, b):
    aa = a / (np.linalg.norm(a, axis=1, keepdims=True) + 1e-12)
    bb = b / (np.linalg.norm(b, axis=1, keepdims=True) + 1e-12)
    sim = aa @ bb.T
    ri, ci = linear_sum_assignment(-sim)
    return float(np.mean(sim[ri, ci]))

def program_analysis(m, features, scores):
    gene_names = np.asarray(m.var_names.astype(str))
    prog_genes = [g for g in features if g in set(gene_names)][:PROGRAM_FEATURES]
    x = m[:, prog_genes].X
    x = x.toarray().astype(np.float32) if sparse.issparse(x) else np.asarray(x, dtype=np.float32)
    donors = m.obs["donor_id"].astype(str).to_numpy()
    rng = np.random.default_rng(BASE_SEED)
    balanced = []
    for donor in sorted(np.unique(donors)):
        pool = np.where(donors == donor)[0]
        balanced.extend(rng.choice(pool, size=min(400, len(pool)), replace=False))
    balanced = np.sort(np.asarray(balanced, dtype=int))
    selection = []
    model_cache = {}
    for k in PROGRAM_K_GRID:
        errs, comps = [], []
        for seed in SEEDS:
            model = NMF(n_components=k, init="nndsvda", random_state=seed,
                        max_iter=500, tol=1e-4, solver="cd")
            model.fit(x[balanced])
            errs.append(float(model.reconstruction_err_ / np.sqrt(np.prod(x[balanced].shape))))
            comps.append(model.components_)
        sim = [align_component_similarity(comps[0], comps[i]) for i in range(1, len(comps))]
        selection.append({"K": k, "mean_normalized_reconstruction_error": float(np.mean(errs)),
                          "sd_normalized_reconstruction_error": float(np.std(errs)),
                          "mean_component_cosine_stability": float(np.mean(sim)),
                          "minimum_component_cosine_stability": float(np.min(sim))})
        model_cache[k] = (errs, comps)
    sel = pd.DataFrame(selection)
    # Expression-only rule: select the stable K within 1 SD of the best
    # reconstruction, breaking ties toward the smaller model.
    best_err = float(sel["mean_normalized_reconstruction_error"].min())
    candidates = sel[sel["mean_normalized_reconstruction_error"] <= best_err * 1.05]
    if len(candidates) == 0:
        candidates = sel
    chosen = int(candidates.sort_values(["mean_component_cosine_stability", "K"],
                                        ascending=[False, True]).iloc[0]["K"])
    sel["selected"] = sel["K"] == chosen
    sel["selection_rule"] = "expression-only: reconstruction within 5% of best, then component stability, then smaller K"
    write_tsv(sel, OUT / "muller_program_number_selection.tsv")
    model = NMF(n_components=chosen, init="nndsvda", random_state=BASE_SEED,
                max_iter=700, tol=1e-4, solver="cd")
    w = model.fit_transform(x)
    h = model.components_
    rows = []
    donor_activity = {}
    for pidx in range(chosen):
        order = np.argsort(h[pidx])[::-1][:30]
        top = gene_names[np.array([list(gene_names).index(g) for g in prog_genes])[order]]
        top_text = ";".join(top.tolist())
        for trait in TRAITS:
            dmeans, smeans = [], []
            for donor in sorted(np.unique(donors)):
                ix = donors == donor
                dmeans.append(float(np.mean(w[ix, pidx])))
                smeans.append(float(np.mean(scores[trait][ix])))
            rho, p = spearmanr(dmeans, smeans)
            rows.append({
                "program": "MG_program_%02d" % (pidx + 1), "program_index": pidx + 1,
                "trait": trait, "n_genes_used": len(prog_genes),
                "top_genes": top_text, "normalized_reconstruction_error": float(model.reconstruction_err_ / np.sqrt(np.prod(x.shape))),
                "donor_activity_mean": float(np.mean(dmeans)), "donor_activity_sd": float(np.std(dmeans, ddof=1)),
                "donor_score_rho": float(rho), "donor_score_P": float(p),
                "program_generation": "unsupervised_expression_only",
            })
            donor_activity[(pidx + 1, trait)] = (np.asarray(dmeans), np.asarray(smeans))
    out = pd.DataFrame(rows)
    out["donor_score_FDR"] = out.groupby("trait", sort=False)["donor_score_P"].transform(lambda q: bh(q.to_numpy()))
    out["robust_positive_program_trait"] = (out["donor_score_rho"] > 0) & (out["donor_score_FDR"] < 0.05)
    write_tsv(out, OUT / "muller_gene_programs.tsv")
    return out, w, h, prog_genes, chosen

def gene_correlation(m, scores, markers, program_table, program_h, program_genes):
    gene_names = np.asarray(m.var_names.astype(str))
    candidates = set(program_genes)
    if len(markers):
        candidates.update(markers.loc[markers["FDR"] < 0.10, "gene"].astype(str).tolist())
    candidates.update(g for g in gene_names if g in set().union(*[set(x) for x in load_gene_sets().values()]))
    candidates = [g for g in gene_names if g in candidates]
    # Limit to a reproducible, variance-rich panel if marker output is very broad.
    if len(candidates) > 6000:
        var = np.asarray(m[:, candidates].X.power(2).mean(axis=0) - np.square(m[:, candidates].X.mean(axis=0))).ravel()
        candidates = [g for _, g in sorted(zip(var, candidates), reverse=True)[:6000]]
    x = m[:, candidates].X
    x = x.toarray().astype(np.float32) if sparse.issparse(x) else np.asarray(x, dtype=np.float32)
    donors = m.obs["donor_id"].astype(str).to_numpy()
    marker_support = set(markers.loc[(markers["FDR"] < 0.05) & (markers["logFC"] > 0), "gene"].astype(str)) if len(markers) else set()
    program_support = {}
    for pidx in range(program_h.shape[0]):
        order = np.argsort(program_h[pidx])[::-1][:100]
        program_support.update({program_genes[i]: "MG_program_%02d" % (pidx + 1) for i in order})
    rows = []
    for trait, score in scores.items():
        donor_scores = np.array([np.mean(score[donors == d]) for d in sorted(np.unique(donors))])
        for j, gene in enumerate(candidates):
            cell_rho, cell_p = spearmanr(x[:, j], score)
            donor_expr = np.array([np.mean(x[donors == d, j]) for d in sorted(np.unique(donors))])
            donor_rho, donor_p = spearmanr(donor_expr, donor_scores)
            rows.append({"trait": trait, "gene": gene, "cell_spearman_rho": float(cell_rho),
                         "cell_P_descriptive": float(cell_p), "donor_spearman_rho": float(donor_rho),
                         "donor_P": float(donor_p), "marker_support": gene in marker_support,
                         "program_support": program_support.get(gene, "NONE")})
    out = pd.DataFrame(rows)
    out["cell_FDR_descriptive"] = out.groupby("trait", sort=False)["cell_P_descriptive"].transform(lambda q: bh(q.to_numpy()))
    out["donor_FDR"] = out.groupby("trait", sort=False)["donor_P"].transform(lambda q: bh(q.to_numpy()))
    out["candidate_support"] = np.where(out["marker_support"] | out["program_support"].ne("NONE"),
                                        "MARKER_OR_PROGRAM_SUPPORTED", "SCORE_CORRELATION_ONLY")
    out["analysis_note"] = "cell-level correlation is descriptive; donor-level support is required for primary candidates"
    out = out.sort_values(["trait", "donor_FDR", "cell_FDR_descriptive"], na_position="last")
    out = out.groupby("trait", sort=False).head(2000).reset_index(drop=True)
    write_tsv(out, OUT / "muller_trait_correlated_genes.tsv")
    return out

def state_program_consensus(m, effects, program_table, w, chosen):
    donors = m.obs["donor_id"].astype(str).to_numpy()
    states = m.obs["muller_state"].astype(str).to_numpy()
    eligible_states = sorted(effects["state"].unique()) if len(effects) else []
    rows = []
    for state in eligible_states:
        for pidx in range(chosen):
            for trait in TRAITS:
                state_d = effects.loc[(effects["trait"] == trait) & (effects["state"] == state), "mean_donor_delta"]
                p_row = program_table.loc[(program_table["program_index"] == pidx + 1) &
                                           (program_table["trait"] == trait)]
                if len(state_d) == 0 or len(p_row) == 0:
                    continue
                state_val = float(state_d.iloc[0])
                p_val = float(p_row["donor_score_rho"].iloc[0])
                per_donor_state = []
                per_donor_program = []
                for donor in sorted(np.unique(donors)):
                    ix = donors == donor
                    sm = ix & (states == state)
                    ot = ix & (states != state)
                    if sm.sum() >= MIN_CELLS_PER_DONOR and ot.sum() > 0:
                        per_donor_state.append(float(np.mean(m.obs.loc[sm, "%s_norm_score" % trait].to_numpy(dtype=float)) -
                                                      np.mean(m.obs.loc[ot, "%s_norm_score" % trait].to_numpy(dtype=float))))
                        per_donor_program.append(float(np.mean(w[ix, pidx])))
                rows.append({"state": state, "program": "MG_program_%02d" % (pidx + 1), "trait": trait,
                             "state_effect": state_val, "program_donor_score_rho": p_val,
                             "state_positive": state_val > 0, "program_positive": p_val > 0,
                             "state_and_program_positive": state_val > 0 and p_val > 0,
                             "state_donor_effects": ";".join("%.5g" % x for x in per_donor_state),
                             "program_donor_activity": ";".join("%.5g" % x for x in per_donor_program)})
    out = pd.DataFrame(rows, columns=["state", "program", "trait", "state_effect",
                                      "program_donor_score_rho", "state_positive", "program_positive",
                                      "state_and_program_positive", "state_donor_effects",
                                      "program_donor_activity"])
    out["consensus_flag"] = np.where(out["state_and_program_positive"], "CONCORDANT_DIRECTION", "NO_CONSENSUS")
    write_tsv(out, OUT / "MULLER_STATE_PROGRAM_CONSENSUS.tsv")
    return out

def tier_and_verdict(effects, lodo, mapping, cluster_sens, program_table, eligibility):
    if len(effects) == 0:
        return effects, "MULLER_STATE_ANALYSIS_WEAK", []
    out = effects.copy()
    lodo_pass = {}
    for (trait, state), g in lodo.groupby(["trait", "state"]):
        lodo_pass[(trait, state)] = bool((g["leave_one_out_mean_delta"] > 0).all())
    map_fail = {}
    for (trait, state), g in mapping.groupby(["trait", "state"]):
        map_fail[(trait, state)] = bool((g["direction_reversal"] == True).any() or
                                        (g["status"] == "NOT_ESTIMABLE").any())
    stability_ok = {}
    for state in out["state"].unique():
        overlaps = cluster_sens.loc[cluster_sens["comparison"] == "resolution_vs_primary", "ARI"]
        stability_ok[state] = bool(len(overlaps) and float(np.nanmedian(overlaps)) >= 0.60)
    out["LODO_all_positive"] = [lodo_pass.get((r.trait, r.state), False) for r in out.itertuples()]
    out["mapping_sensitivity_pass"] = [not map_fail.get((r.trait, r.state), True) for r in out.itertuples()]
    out["cluster_sensitivity_pass"] = [stability_ok.get(s, False) for s in out["state"]]
    out["tier"] = np.where(
        (out["global_FDR"] < 0.05) & (out["mean_donor_delta"] > 0) &
        (out["direction_fraction_positive"] >= 0.75) & out["LODO_all_positive"] &
        out["mapping_sensitivity_pass"] & out["cluster_sensitivity_pass"] &
        (out["donor_composition_status"] == "PASS"), "TIER1_STATE",
        np.where((out["within_trait_FDR"] < 0.05) | ((out["primary_P_t_test"] < 0.05) &
                    out["LODO_all_positive"] & out["mapping_sensitivity_pass"]), "TIER2_STATE", "NONE"))
    def cross(row):
        hit = out[(out["state"] == row["state"]) & (out["mean_donor_delta"] > 0) &
                  (out["tier"].isin(["TIER1_STATE", "TIER2_STATE"]))]
        traits = set(hit["trait"])
        if traits == set(TRAITS):
            return "PAN_TRAIT_CONVERGENT"
        if len(traits) >= 2:
            return "_".join(sorted(traits)) + "_CONVERGENT"
        if len(traits) == 1:
            return "TRAIT_SPECIFIC"
        return "NONE"
    out["cross_trait_class"] = [cross(r) for _, r in out.iterrows()]
    write_tsv(out, OUT / "muller_state_trait_effects.tsv")
    t1 = out[out["tier"] == "TIER1_STATE"]
    prog = program_table[(program_table["donor_score_rho"] > 0) & (program_table["donor_score_FDR"] < 0.05)]
    convergent_state = int(t1.groupby("state")["trait"].nunique().max()) if len(t1) else 0
    convergent_prog = int(prog.groupby("program")["trait"].nunique().max()) if len(prog) else 0
    if convergent_state >= 2 and convergent_prog >= 2:
        verdict = "MULLER_MIXED_STATE_PROGRAM_GO"
    elif convergent_state >= 2:
        verdict = "MULLER_STATE_CONVERGENCE_GO"
    elif convergent_prog >= 2:
        verdict = "MULLER_PROGRAM_CONVERGENCE_GO"
    elif len(t1):
        verdict = "MULLER_TRAIT_SPECIFIC_STATE_GO"
    elif len(effects) and float(effects["global_FDR"].min()) >= 0.05:
        verdict = "MULLER_BROAD_UNIFORM"
    else:
        verdict = "MULLER_STATE_ANALYSIS_WEAK"
    return out, verdict, t1.to_dict("records")

def plot_figures(m, effects, markers, pathways, consensus, program_table, chosen_state, verdict):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    FIGURES.mkdir(parents=True, exist_ok=True)
    colors = {"AMD": "#C44E52", "POAG": "#4C72B0", "RE": "#55A868"}
    obs = m.obs
    umap = np.asarray(m.obsm["X_umap"])
    def save(fig, stem):
        fig.savefig(FIGURES / (stem + ".png"), dpi=300, bbox_inches="tight")
        fig.savefig(FIGURES / (stem + ".svg"), bbox_inches="tight")
        plt.close(fig)
    fig, axes = plt.subplots(1, 5, figsize=(22, 4.2))
    plots = [
        (obs["muller_state"].astype(str).to_numpy(), "Müller state", True),
        (obs["donor_id"].astype(str).to_numpy(), "Donor", True),
        (obs["AMD_norm_score"].to_numpy(float), "AMD scDRS", False),
        (obs["POAG_norm_score"].to_numpy(float), "POAG scDRS", False),
        (obs["RE_norm_score"].to_numpy(float), "RE scDRS", False),
    ]
    for ax, (value, title, categorical) in zip(axes, plots):
        if categorical:
            cats = sorted(np.unique(value))
            cmap = plt.get_cmap("tab10")
            for i, cat in enumerate(cats):
                q = value == cat
                ax.scatter(umap[q, 0], umap[q, 1], s=2, alpha=0.55, color=cmap(i % 10), label=cat, rasterized=True)
            ax.legend(frameon=False, markerscale=3, fontsize=7, loc="best")
        else:
            sca = ax.scatter(umap[:, 0], umap[:, 1], c=value, s=2, alpha=0.55, cmap="viridis", rasterized=True)
            fig.colorbar(sca, ax=ax, fraction=0.046, pad=0.04)
        ax.set_title(title)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("P18. Müller-only expression state map", y=1.03)
    save(fig, "Figure_P18_Muller_state_map")

    if len(effects):
        states = sorted(effects["state"].unique())
        mat = effects.pivot(index="state", columns="trait", values="mean_donor_delta").reindex(columns=TRAITS)
        fdr = effects.pivot(index="state", columns="trait", values="global_FDR").reindex(columns=TRAITS)
        fig, ax = plt.subplots(figsize=(6.5, max(3.2, 0.55 * len(states) + 1.5)))
        sns.heatmap(mat, annot=mat.round(2), fmt="", cmap="RdBu_r", center=0, linewidths=.5,
                    cbar_kws={"label": "Mean donor state delta"}, ax=ax)
        for i, st in enumerate(mat.index):
            for j, tr in enumerate(mat.columns):
                if np.isfinite(fdr.loc[st, tr]) and fdr.loc[st, tr] < 0.05:
                    ax.text(j + 0.82, i + 0.18, "*", ha="center", va="center", color="black", fontsize=14)
        ax.set_title("P19. Donor-aware Müller state effects; * global FDR < 0.05")
        ax.set_xlabel("Trait")
        ax.set_ylabel("State")
        save(fig, "Figure_P19_Muller_state_effect_heatmap")
        # Individual donor display for the strongest multi-trait state.
        counts = effects[effects["tier"] == "TIER1_STATE"].groupby("state")["trait"].nunique()
        state = str(counts.sort_values(ascending=False).index[0]) if len(counts) else str(states[0])
        d = donor_df_global[(donor_df_global["state"] == state)].copy()
        fig, ax = plt.subplots(figsize=(8, 4.8))
        if len(d):
            sns.stripplot(data=d, x="trait", y="donor_delta", hue="donor", dodge=True, size=6, ax=ax)
            ax.axhline(0, color="black", lw=.8)
            ax.set_title("P20. Individual-donor effects for " + state)
            ax.set_ylabel("State minus other Müller cells")
            ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", title="Donor")
            save(fig, "Figure_P20_Muller_candidate_donor_effects")
        else:
            plt.close(fig)
    else:
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, "No Müller state passed donor-level eligibility", ha="center", va="center")
        ax.set_axis_off()
        ax.set_title("P19. Donor-aware Müller state effects")
        save(fig, "Figure_P19_Muller_state_effect_heatmap")
        fig, ax = plt.subplots(figsize=(8, 4.8))
        ax.text(0.5, 0.5, "No candidate state passed donor-level eligibility", ha="center", va="center")
        ax.set_axis_off()
        ax.set_title("P20. Individual-donor Müller state effects")
        save(fig, "Figure_P20_Muller_candidate_donor_effects")

    if len(pathways):
        p = pathways.copy()
        p["minus_log10_FDR"] = -np.log10(np.clip(pd.to_numeric(p["FDR"], errors="coerce"), 1e-300, 1))
        p = p.sort_values(["FDR", "P"]).head(15)
        p["label"] = p["pathway"].astype(str).str.replace("_", " ", regex=False).str.slice(0, 55)
        fig, ax = plt.subplots(figsize=(9, max(4.5, 0.34 * len(p) + 1.5)))
        sns.barplot(data=p, x="minus_log10_FDR", y="label", hue="pathway_source", dodge=False, ax=ax)
        ax.set_title("P21. Unbiased Müller marker/pathway programs")
        ax.set_xlabel("-log10 pathway FDR")
        ax.set_ylabel("")
        ax.spines[["top", "right"]].set_visible(False)
        save(fig, "Figure_P21_Muller_marker_pathway_programs")
    elif len(markers):
        p = markers.sort_values(["FDR", "P"]).head(20).copy()
        p["label"] = p["gene"].astype(str)
        p["minus_log10_FDR"] = -np.log10(np.clip(p["FDR"], 1e-300, 1))
        fig, ax = plt.subplots(figsize=(8, 5))
        sns.barplot(data=p, x="minus_log10_FDR", y="label", hue="state", dodge=False, ax=ax)
        ax.set_title("P21. Müller state marker programs")
        save(fig, "Figure_P21_Muller_marker_pathway_programs")
    else:
        fig, ax = plt.subplots(figsize=(8, 3.5))
        ax.text(0.5, 0.5, "No eligible donor-replicated states; marker/pathway analysis not estimable",
                ha="center", va="center", wrap=True)
        ax.set_axis_off()
        ax.set_title("P21. Unbiased Müller marker/pathway programs")
        save(fig, "Figure_P21_Muller_marker_pathway_programs")

    if len(consensus):
        q = consensus.pivot_table(index="state", columns="program", values="state_effect", aggfunc="mean")
        fig, ax = plt.subplots(figsize=(7, max(3.2, .55 * len(q) + 1.5)))
        sns.heatmap(q, annot=True, fmt=".2f", cmap="RdBu_r", center=0, linewidths=.5, ax=ax,
                    cbar_kws={"label": "State effect"})
        ax.set_title("P22. Müller state × continuous program convergence")
        ax.set_xlabel("Unsupervised expression program")
        ax.set_ylabel("Müller state")
        save(fig, "Figure_P22_Muller_state_program_convergence")
    else:
        fig, ax = plt.subplots(figsize=(7, 3))
        ax.text(0.5, 0.5, "No eligible discrete-state × program rows", ha="center", va="center")
        ax.set_axis_off()
        ax.set_title("P22. Müller state × continuous program convergence")
        save(fig, "Figure_P22_Muller_state_program_convergence")

def report_files(m, audit, inv, stab, summary, comp, hetero, effects, markers,
                 pathways, programs, consensus, gene_corr, verdict, chosen, official_ready,
                 next_action):
    def fmt(x, digits=4):
        if x is None or (isinstance(x, float) and not np.isfinite(x)):
            return "NA"
        try:
            return ("%." + str(digits) + "f") % float(x)
        except Exception:
            return str(x)
    eligible = comp[comp["state_eligible"]].copy()
    tier1 = effects[effects["tier"] == "TIER1_STATE"].copy() if len(effects) else pd.DataFrame()
    robust_programs = programs[(programs["donor_score_rho"] > 0) &
                               (programs["donor_score_FDR"] < 0.05)].copy()
    top_path = pathways.sort_values(["FDR", "P"]).head(10) if len(pathways) else pathways
    state_rows = []
    for st in sorted(eligible["state"].unique()):
        q = effects[effects["state"] == st] if len(effects) else pd.DataFrame()
        labs = ";".join(q.loc[q["tier"].isin(["TIER1_STATE", "TIER2_STATE"]), "trait"].astype(str)) if len(q) else ""
        state_rows.append("%s: %d cells, %d donors; robust traits=%s" % (
            st, int(eligible.loc[eligible["state"] == st, "cells"].iloc[0]),
            int(eligible.loc[eligible["state"] == st, "donors_with_30_cells"].iloc[0]), labs or "none"))
    pathway_text = "\n".join("- %s | %s | FDR=%s | genes=%s" % (
        r.pathway_source, r.pathway, fmt(r.FDR, 3), str(r.leading_genes)[:180])
        for r in top_path.itertuples()) if len(top_path) else "- No pathway rows passed the unbiased marker input."
    effect_text = "\n".join("- %s / %s: delta=%s, global FDR=%s, direction=%d/%d, tier=%s" % (
        r.trait, r.state, fmt(r.mean_donor_delta), fmt(r.global_FDR, 3),
        int(r.donors_positive), int(r.eligible_donors), r.tier)
        for r in effects.itertuples()) if len(effects) else "- No state met primary eligibility."
    program_text = "\n".join("- %s / %s: donor rho=%s, FDR=%s; top genes=%s" % (
        r.program, r.trait, fmt(r.donor_score_rho), fmt(r.donor_score_FDR, 3),
        str(r.top_genes)[:180])
        for r in robust_programs.itertuples()) if len(robust_programs) else "- No program-trait association passed FDR < 0.05."
    common = [
        "# Phase 3A-2 Müller state and program analysis",
        "",
        "This report follows the locked Phase 3A-2 analysis. State discovery used expression only and was completed before loading AMD, POAG, or RE scDRS scores.",
        "",
        "## Decision",
        "",
        "Verdict: %s" % verdict,
        "Primary Leiden resolution: %s. Official fine substates ready: %s. Primary de novo state count: %d." %
        (chosen, official_ready, int(m.obs["muller_state"].nunique())),
        "",
        "## Müller reference",
        "",
        "Cells=%d; donors=%d; samples=%d. Broad class was the exact official-ID mapped healthy/reference Müller population. The only ontology label CL:0000636 was treated as broad Müller identity, not a fine state." %
        (m.n_obs, m.obs["donor_id"].nunique(), m.obs["sample"].nunique()),
        "",
        "Eligible state audit:",
        "",
        "\n".join("- " + x for x in state_rows) if state_rows else "- No de novo state passed eligibility.",
        "",
        "## Within-Müller scDRS heterogeneity",
        "",
        hetero.to_string(index=False) if len(hetero) else "No heterogeneity rows.",
        "",
        "## State effects",
        "",
        effect_text,
        "",
        "Primary inference uses equal-weight donor deltas, exact sign-flip P separately from t-test P, 5,000 donor bootstrap replicates, global BH FDR across three traits and eligible states, leave-one-donor-out direction checks, and mapping-selection sensitivity.",
        "",
        "## Tier 1 states",
        "",
        "\n".join("- %s / %s: delta=%s, global FDR=%s, cross-trait=%s" % (
            r.state, r.trait, fmt(r.mean_donor_delta), fmt(r.global_FDR, 3), r.cross_trait_class)
            for r in tier1.itertuples()) if len(tier1) else "- None.",
        "",
        "## Limitations",
        "",
        "- Eight donors provide the donor-level replication unit; donor-level estimates remain small-n.",
        "- AMD cellular power is limited, so a null AMD state result is weak evidence of absence.",
        "- One sample per donor means donor and sample are confounded.",
        "- Primary cells are healthy/reference cells; state names are descriptive transcriptional labels and do not establish a disease state.",
        "- UMAP, cell-level correlations, and pathway enrichment are secondary descriptive evidence; regulatory or causal validation was not performed.",
        "",
        "## Next action",
        "",
        next_action,
    ]
    (ROOT / "reports" / "PHASE3A2_MULLER_STATE_VULNERABILITY.md").write_text("\n".join(common), encoding="utf-8")
    audit_text = [
        "# Phase 3A-2 Müller state discovery audit", "",
        "## Locked design", "",
        "The analysis used exact official-ID mapped healthy/reference Müller cells from the Phase 3A-1 H5AD. State discovery was blind to all three scDRS score files. The disease-gene union was excluded only from the 2,000-gene cell-ranger HVG state space; the genes remained in the expression object for marker, program, and correlation analyses.",
        "",
        "## Official annotation inventory", "",
        inv.to_string(index=False),
        "",
        "Official replicated fine states were %s. Because CL:0000636 is a broad ontology label and no two replicated official fine substates were available, de novo states were used as primary." % official_ready,
        "",
        "## Stability and selection", "",
        summary.to_string(index=False),
        "",
        stab.to_string(index=False),
        "",
        "Selection used predefined stability, donor representation, and tiny-cluster criteria only; trait association did not enter resolution selection.",
        "",
        "## Composition", "",
        comp.to_string(index=False),
        "",
        "The one-sample-per-donor design prevents separating donor from sample effects. Age and region fields were audited; no estimable variation was available in the primary H5AD.",
    ]
    (ROOT / "reports" / "PHASE3A2_MULLER_STATE_DISCOVERY_AUDIT.md").write_text("\n".join(audit_text), encoding="utf-8")
    program_text_full = [
        "# Phase 3A-2 continuous Müller programs", "",
        "Continuous programs were generated by unsupervised NMF on expression-only genes after state-discovery disease-gene exclusion. The program number was selected from K=4,5,6 using reconstruction and component stability before donor-level scDRS associations were evaluated.",
        "",
        programs.to_string(index=False),
        "",
        "## State-program consensus", "",
        consensus.to_string(index=False) if len(consensus) else "No consensus rows.",
        "",
        "## Gene-scDRS candidates", "",
        gene_corr.head(100).to_string(index=False) if len(gene_corr) else "No candidates.",
        "",
        "Program associations are donor-level summaries across eight donors; cell-level gene correlations are descriptive and require marker or program support for biological prioritization.",
    ]
    (ROOT / "reports" / "PHASE3A2_MULLER_PROGRAM_ANALYSIS.md").write_text("\n".join(program_text_full), encoding="utf-8")

def write_status(verdict, chosen, m, effects, programs):
    tier1 = effects[effects["tier"] == "TIER1_STATE"] if len(effects) else effects
    robust = programs[(programs["donor_score_rho"] > 0) & (programs["donor_score_FDR"] < 0.05)] if len(programs) else programs
    status = [
        "# Current status",
        "",
        "Updated: 2026-09-06",
        "",
        "## Phase 3A-2",
        "",
        "Completed donor-aware Müller-glia state and unsupervised program analysis.",
        "",
        "- Verdict: %s" % verdict,
        "- Reference: %d Müller cells, %d donors, %d samples" % (m.n_obs, m.obs["donor_id"].nunique(), m.obs["sample"].nunique()),
        "- Primary de novo resolution: %s; states: %d" % (chosen, m.obs["muller_state"].nunique()),
        "- Tier 1 state-trait rows: %d" % len(tier1),
        "- Robust program-trait rows: %d" % len(robust),
        "- Phase 3B and RE amacrine Phase 3A-2B were not started.",
        "",
        "Required reports and results are under retinaplex_state/reports and retinaplex_state/results/phase3a2.",
    ]
    (ROOT / "reports" / "CURRENT_STATUS.md").write_text("\n".join(status) + "\n", encoding="utf-8")
    (ROOT / "PROJECT_STATUS.md").write_text("\n".join(status) + "\n", encoding="utf-8")

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    start = time.time()
    m = load_muller()
    audit = audit_population(m)
    inv, official_ready = annotation_inventory(m)
    write_lock(m, inv)
    gene_sets = load_gene_sets()
    m, raw, features, gene_audit = prepare_features(m, gene_sets)
    m, proc, stab, summary, chosen, features, labels = discover_states(m, features)
    comp, eligibility = state_composition(m)
    # No trait score is loaded until the blind state-discovery stage is frozen.
    scores = load_scores(m)
    hetero = heterogeneity()
    donor_df = compute_donor_deltas(m, scores, eligibility)
    write_tsv(donor_df, OUT / "muller_state_trait_deltas.tsv")
    effects = compute_effect_table(donor_df, eligibility)
    lodo = leave_one_out(donor_df)
    mapping = mapping_sensitivity(m, scores, effects)
    balance = balance_sensitivity(m, scores, eligibility)
    cluster_sens = cluster_sensitivity(m, proc, features, chosen, labels, stab)
    markers, pathways = markers_and_pathways(m, raw, eligibility)
    programs, w, h, prog_genes, nprograms = program_analysis(m, features, scores)
    gene_corr = gene_correlation(m, scores, markers, programs, h, prog_genes)
    effects, verdict, tier1 = tier_and_verdict(effects, lodo, mapping, cluster_sens, programs, eligibility)
    consensus = state_program_consensus(m, effects, programs, w, nprograms)
    global donor_df_global
    donor_df_global = donor_df
    candidate_state = tier1[0]["state"] if tier1 else (effects.sort_values("global_FDR").iloc[0]["state"] if len(effects) else "NA")
    if verdict in {"MULLER_MIXED_STATE_PROGRAM_GO", "MULLER_STATE_CONVERGENCE_GO", "MULLER_PROGRAM_CONVERGENCE_GO"}:
        next_action = "A. RUN PHASE 3B ATAC/REGULATORY VALIDATION"
    elif verdict == "MULLER_BROAD_UNIFORM":
        next_action = "C. KEEP MÜLLER AT BROAD-CELL LEVEL"
    else:
        next_action = "D. EXPAND/VALIDATE MÜLLER STATE ANALYSIS"
    plot_figures(m, effects, markers, pathways, consensus, programs, candidate_state, verdict)
    report_files(m, audit, inv, stab, summary, comp, hetero, effects, markers, pathways,
                 programs, consensus, gene_corr, verdict, chosen, official_ready, next_action)
    write_status(verdict, chosen, m, effects, programs)
    software = []
    for name in ["anndata", "scanpy", "scdrs", "numpy", "pandas", "scikit-learn", "scipy"]:
        try:
            from importlib.metadata import version
            software.append({"package": name, "version": version(name)})
        except Exception:
            software.append({"package": name, "version": "NA"})
    software.append({"package": "igraph/leidenalg", "version": "vendor phase3a2"})
    write_tsv(pd.DataFrame(software), OUT / "software_versions.tsv")
    print(json.dumps({"verdict": verdict, "cells": int(m.n_obs), "donors": int(m.obs["donor_id"].nunique()),
                      "states": int(m.obs["muller_state"].nunique()), "eligible_states": int(eligibility["state_eligible"].sum()),
                      "programs": int(nprograms), "tier1_rows": int(len(tier1)),
                      "runtime_min": round((time.time() - start) / 60, 2)}, ensure_ascii=False))

if __name__ == "__main__":
    main()
