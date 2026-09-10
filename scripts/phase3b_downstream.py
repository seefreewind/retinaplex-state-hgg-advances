#!/usr/bin/env python3
"""Build Phase 3B consensus tables, figures, and audit-ready reports."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr, t as tdist


TRAITS = ["AMD", "POAG", "RE"]
RNA_CLASS = ["Muller glia", "amacrine cells", "bipolar cells", "cones", "rods", "horizontal cells", "RGC", "astrocytes", "microglia"]
ATAC_CLASS = ["MG", "AC", "BC", "Cone", "Rod", "HC", "RGC", "Astrocyte", "Microglia"]
CLASS_TO_RNA = dict(zip(ATAC_CLASS, RNA_CLASS))
RNA_TO_ATAC = dict(zip(RNA_CLASS, ATAC_CLASS))


def read_tsv(path):
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fnum(x, default=float("nan")):
    try:
        if x in (None, "", "NA", "nan", "NaN"): return default
        return float(x)
    except (TypeError, ValueError):
        return default


def bh(values):
    out = [float("nan")] * len(values)
    good = sorted((i, v) for i, v in enumerate(values) if np.isfinite(v))
    m = len(good); prev = 1.0
    for rank in range(m, 0, -1):
        i, p = good[rank - 1]; prev = min(prev, p * m / rank); out[i] = prev
    return out


def zscore(vals):
    x = np.asarray(vals, dtype=float)
    ok = np.isfinite(x)
    out = np.full(len(x), np.nan)
    if ok.sum() >= 2:
        sd = np.std(x[ok], ddof=1)
        if sd > 0: out[ok] = (x[ok] - np.mean(x[ok])) / sd
    return out


def write_tsv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(header); w.writerows(rows)


def ols(y, x):
    y = np.asarray(y, float); x = np.asarray(x, float)
    n, p = x.shape
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    resid = y - x @ beta
    dof = max(1, n - p)
    s2 = float(resid @ resid) / dof
    cov = s2 * np.linalg.pinv(x.T @ x)
    se = np.sqrt(np.maximum(0, np.diag(cov)))
    tstat = np.divide(beta, se, out=np.full_like(beta, np.nan), where=se > 0)
    pval = 2 * tdist.sf(np.abs(tstat), dof)
    return beta, se, tstat, pval, n, dof


def load_magma_gene_universe(project: Path):
    scores = {r["gene_symbol"]: r for r in read_tsv(project / "results/phase3ar8/magma_gene_scores.tsv") if r.get("gene_symbol")}
    ctx = {r["gene_symbol"]: r for r in read_tsv(project / "results/phase3b/official_extracted/magma_gene_muller_nearest_context.tsv") if r.get("gene_symbol")}
    cov = {r["GENE"]: r for r in read_tsv(project / "results/phase3a1/magma_broad_gene_covariates.tsv") if r.get("GENE")}
    nsnps = {t: {} for t in TRAITS}
    for trait in TRAITS:
        p = project / f"results/phase3ar8/magma/{trait}.genes.out"
        if not p.exists(): continue
        with p.open(encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("GENE") or line.startswith("#") or not line.strip(): continue
                f = line.split()
                if len(f) >= 9: nsnps[trait][f[0]] = {"nsnps": f[4], "length": str(max(1, int(f[3]) - int(f[2]) + 1))}
    return scores, ctx, cov, nsnps


def parse_published_cre_genes(project: Path):
    out = {}
    for row in read_tsv(project / "results/phase3b/hrca_cre_gene_link_audit.tsv"):
        if row.get("resource_layer") != "published_finemap_context": continue
        note = row.get("notes", "")
        m = re.search(r"genes=(.*)$", note)
        if m:
            # Extractor stores semicolon-delimited official gene symbols.
            out.setdefault(row.get("metric", ""), set()).update(g for g in m.group(1).split(";") if g)
    # The audit has one row per trait, so recover trait order from the source
    # row position when the metric name is repeated.
    rows = read_tsv(project / "results/phase3b/hrca_cre_gene_link_audit.tsv")
    vals = []
    for row in rows:
        if row.get("resource_layer") == "published_finemap_context":
            m = re.search(r"genes=(.*)$", row.get("notes", "")); vals.append(set(m.group(1).split(";")) if m else set())
    return dict(zip(TRAITS, vals))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--project", required=True); args = ap.parse_args()
    project = Path(args.project); results = project / "results/phase3b"; figdir = project / "figures/phase3b"; figdir.mkdir(parents=True, exist_ok=True)
    phase3a = {(r["trait"], r["broad_class"]): r for r in read_tsv(project / "results/phase3a1/BROAD_CELLULAR_CONSENSUS.tsv")}
    magma = {(r["trait"], r["broad_class"]): r for r in read_tsv(project / "results/phase3a1/magma_broad_celltype_enrichment.tsv")}
    sldsc_path = results / "sldsc_retinal_ocr_enrichment.tsv"
    sldsc = {(r["trait"], r["cell_class"]): r for r in read_tsv(sldsc_path)} if sldsc_path.exists() else {}

    # RNA--ATAC cellular consensus.
    consensus = []
    for trait in TRAITS:
        for atac in ATAC_CLASS:
            rna_name = CLASS_TO_RNA[atac]; rr = phase3a.get((trait, rna_name), {}); mm = magma.get((trait, rna_name), {}); aa = sldsc.get((trait, atac), {})
            rna_effect = fnum(rr.get("mean_donor_delta")); rna_fdr = fnum(rr.get("global_broad_FDR")); magma_p = fnum(mm.get("p")); magma_fdr = fnum(mm.get("global_broad_FDR")); atac_coef = fnum(aa.get("coefficient")); atac_fdr = fnum(aa.get("global_primary_fdr")); atac_enrich = fnum(aa.get("enrichment"))
            rna_tier1 = rr.get("tier", "") == "TIER1_BROAD_CLASS"
            atac_tier1 = np.isfinite(atac_coef) and atac_coef > 0 and np.isfinite(atac_fdr) and atac_fdr < 0.05
            if rna_tier1 and atac_tier1: category = "MULTIOMIC_HIGH_CONFIDENCE"
            elif rna_tier1 and not np.isfinite(atac_coef): category = "RNA_ONLY"
            elif rna_tier1: category = "RNA_ONLY"
            elif atac_tier1: category = "ATAC_ONLY"
            else: category = "UNSUPPORTED"
            consensus.append([trait, rna_name, atac, rna_effect, rna_fdr, rr.get("tier", ""), rr.get("confidence", ""), mm.get("beta", ""), magma_p, magma_fdr, atac_coef, atac_enrich, atac_fdr, aa.get("annotation_type", "OFFICIAL_MAJOR_CLASS_OCR_SURROGATE"), category])
    write_tsv(results / "RNA_ATAC_CELLULAR_CONSENSUS.tsv", ["trait", "rna_broad_class", "atac_broad_class", "phase3a_donor_aware_effect", "phase3a_global_fdr", "phase3a_tier", "phase3a_confidence", "magma_beta", "magma_p", "magma_fdr", "atac_sldsc_coefficient", "atac_sldsc_enrichment", "atac_global_fdr", "atac_annotation_type", "consensus_classification"], consensus)

    # Standardized profile concordance within each trait.
    concordance = []; regulatory_vectors = {}
    for trait in TRAITS:
        rna = []; atac = []; labels = []
        for cls in ATAC_CLASS:
            rn = CLASS_TO_RNA[cls]; rr = phase3a.get((trait, rn), {}); aa = sldsc.get((trait, cls), {})
            rna.append(fnum(rr.get("mean_donor_delta"))); atac.append(fnum(aa.get("coefficient"))); labels.append(cls)
        zr, za = zscore(rna), zscore(atac); ok = np.isfinite(zr) & np.isfinite(za)
        if ok.sum() >= 3:
            rho, p = spearmanr(zr[ok], za[ok])
        else: rho, p = float("nan"), float("nan")
        concordance.append([trait, int(ok.sum()), rho, p, ";".join(labels), "Spearman on within-modality standardized broad-class effects; no donor bootstrap because pooled official OCR was the only auditable ATAC annotation"])
        regulatory_vectors[trait] = np.asarray(atac)
    write_tsv(results / "rna_atac_profile_concordance.tsv", ["trait", "n_classes", "spearman_rho", "p_value", "classes", "method_note"], concordance)

    # Gene-level regulatory context regression, using the frozen genome-wide
    # MAGMA universe. This is explicitly nearest-gene context, not a formal
    # CRE--gene link analysis.
    scores, ctx, cov, nsnps, = load_magma_gene_universe(project)
    regression_rows = []; gene_context = {}
    for trait in TRAITS:
        genes = []; y = []; xrows = []
        for gene, score in scores.items():
            z = fnum(score.get(f"{trait}_ZSTAT")); c = ctx.get(gene, {}); co = cov.get(score.get("gene_id", ""), {})
            if not np.isfinite(z): continue
            mg = fnum(c.get("muller_nearest_peak_count"), 0); total = fnum(c.get("all_class_nearest_peak_count"), 0); spec = fnum(c.get("muller_specificity_ratio"), 0); avg = fnum(co.get("Average"), 0)
            raw = nsnps.get(trait, {}).get(score.get("gene_id", ""), {})
            n_snp = fnum(raw.get("nsnps"), 1); length = fnum(raw.get("length"), 1)
            genes.append(gene); y.append(z); xrows.append([1, spec, math.log1p(total), avg, math.log1p(n_snp), math.log1p(length)])
        y = np.asarray(y); x = np.asarray(xrows)
        beta, se, tstat, pval, n, dof = ols(y, x)
        terms = ["intercept", "muller_specificity_ratio", "log1p_all_class_nearest_peak_count", "average_retinal_expression", "log1p_magma_nsnp", "log1p_magma_gene_length"]
        for term, b, s, tval, pv in zip(terms, beta, se, tstat, pval):
            regression_rows.append([trait, term, b, s, tval, pv, n, dof, "MAGMA gene Z ~ MG-specific nearest-gene OCR ratio + total retinal nearest-gene OCR burden + Average expression + MAGMA SNP/gene covariates", "EXPLORATORY_NEAREST_GENE_CONTEXT_NOT_FORMAL_CRE_GENE_LINK"])
        gene_context[trait] = {g: (fnum(scores[g].get(f"{trait}_ZSTAT")), ctx.get(g, {})) for g in genes}
    write_tsv(results / "gene_regulatory_context_regression.tsv", ["trait", "term", "beta", "se", "t_stat", "p_value", "n_genes", "residual_dof", "model", "status"], regression_rows)

    # Published locus-level CRE context candidates, ranked by frozen MAGMA Z.
    published = parse_published_cre_genes(project); prioritized = []
    for trait in TRAITS:
        for gene in sorted(published.get(trait, set()), key=lambda g: fnum(scores.get(g, {}).get(f"{trait}_ZSTAT"), -999), reverse=True)[:15]:
            s = scores.get(gene, {}); c = ctx.get(gene, {})
            prioritized.append([trait, gene, fnum(s.get(f"{trait}_ZSTAT")), fnum(s.get(f"{trait}_P")), fnum(c.get("muller_nearest_peak_count"), 0), fnum(c.get("all_class_nearest_peak_count"), 0), fnum(c.get("muller_specificity_ratio"), 0), "PUBLISHED_TABLE_S20_CRE_CONTEXT", "prioritized regulatory candidate; not a causal gene and not a formal Müller CRE-gene link"])
    write_tsv(results / "muller_regulatory_trait_genes.tsv", ["trait", "gene", "magma_z", "magma_p", "muller_nearest_ocr_count", "all_class_nearest_ocr_count", "muller_specificity_ratio", "official_regulatory_evidence", "interpretation"] , prioritized)

    # Cross-trait regulatory profile similarity, calculated on coefficients and
    # therefore independent of raw P-value ranking.
    cross = []
    for i, t1 in enumerate(TRAITS):
        for t2 in TRAITS[i+1:]:
            a, b = regulatory_vectors[t1], regulatory_vectors[t2]; ok = np.isfinite(a) & np.isfinite(b)
            if ok.sum() >= 3: rho, p = spearmanr(a[ok], b[ok])
            else: rho, p = float("nan"), float("nan")
            cross.append([t1, t2, int(ok.sum()), rho, p, "Spearman on broad-class s-LDSC coefficient profiles; sign retained"])
    write_tsv(results / "cross_trait_regulatory_profile_similarity.tsv", ["trait_1", "trait_2", "n_classes", "spearman_rho", "p_value", "method_note"], cross)

    # Figures are compact and auditable; they use NA panels where the official
    # resource does not support the requested claim.
    make_p23(figdir, sldsc)
    make_p24(figdir, consensus)
    make_p25(figdir, prioritized)
    make_p26(figdir, project)
    make_p27(figdir)
    print(f"consensus_rows={len(consensus)} regression_rows={len(regression_rows)} prioritized_rows={len(prioritized)}")


def savefig(fig, path):
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def make_p23(figdir, sldsc):
    mat = np.full((len(ATAC_CLASS), len(TRAITS)), np.nan); labels = np.full(mat.shape, "NA", dtype=object)
    for j, trait in enumerate(TRAITS):
        for i, cls in enumerate(ATAC_CLASS):
            r = sldsc.get((trait, cls), {}); mat[i, j] = fnum(r.get("coefficient"));
            if r: labels[i, j] = fnum(r.get("global_primary_fdr"))
    fig, ax = plt.subplots(figsize=(6.5, 6.8)); im = ax.imshow(mat, cmap="coolwarm", aspect="auto")
    ax.set_xticks(range(3), TRAITS); ax.set_yticks(range(len(ATAC_CLASS)), ATAC_CLASS); ax.set_title("P23  HRCA broad-class OCR s-LDSC")
    for i in range(len(ATAC_CLASS)):
        for j in range(3):
            txt = "NA" if not np.isfinite(mat[i,j]) else f"{mat[i,j]:.3g}\nFDR {labels[i,j]:.2g}" if np.isfinite(labels[i,j]) else f"{mat[i,j]:.3g}\nFDR NA"
            ax.text(j, i, txt, ha="center", va="center", fontsize=6)
    fig.colorbar(im, ax=ax, label="conditional LDSC coefficient"); savefig(fig, figdir / "Figure_P23_trait_class_ATAC_sLDSC_heatmap")


def make_p24(figdir, consensus):
    vals = np.zeros((len(TRAITS)*len(ATAC_CLASS), 3)); text = []
    for k, r in enumerate(consensus):
        # consensus is kept as the row-list representation used by
        # write_tsv(), so use the fixed output-column positions here.
        vals[k, 0] = 1 if r[5] == "TIER1_BROAD_CLASS" else 0
        vals[k, 1] = 1 if np.isfinite(fnum(r[9])) and fnum(r[9]) < .05 else 0
        vals[k, 2] = 1 if np.isfinite(fnum(r[12])) and fnum(r[12]) < .05 else 0
        text.append(f"{r[0]} × {r[1]}")
    fig, ax = plt.subplots(figsize=(6.8, 8.0)); ax.imshow(vals, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(3), ["RNA donor Tier1", "MAGMA FDR<0.05", "ATAC FDR<0.05"]); ax.set_yticks(range(len(text)), text, fontsize=6); ax.set_title("P24  RNA versus ATAC evidence matrix")
    for i in range(vals.shape[0]):
        for j in range(3): ax.text(j, i, "✓" if vals[i,j] else "–", ha="center", va="center", fontsize=7)
    savefig(fig, figdir / "Figure_P24_RNA_vs_ATAC_evidence_matrix")


def make_p25(figdir, prioritized):
    fig, ax = plt.subplots(figsize=(8, 5)); ax.axis("off"); ax.set_title("P25  Müller regulatory context: formal CRE–gene network audit")
    ax.text(.5, .60, "Formal broad-class HRCA CRE–gene edge matrix\nwas not present in the downloaded official workbook", ha="center", va="center", fontsize=13, weight="bold")
    ax.text(.5, .38, "Nearest-gene OCR annotations and Table S20 locus-level CRE context\nare retained in the TSV audit but are not drawn as causal links.", ha="center", va="center", fontsize=10)
    ax.text(.5, .16, f"Published locus-context candidate rows available: {len(prioritized)}", ha="center", va="center", fontsize=9, color="#555555")
    savefig(fig, figdir / "Figure_P25_Muller_CRE_gene_network")


def make_p26(figdir, project):
    rows = read_tsv(project / "results/phase3b/official_extracted/hrca_muller_top_rss.tsv")
    names = [r["official_regulon_or_TF"] for r in rows[:12]][::-1]; vals = [fnum(r["MG_RSS"]) for r in rows[:12]][::-1]
    fig, ax = plt.subplots(figsize=(7, 5.4)); ax.barh(range(len(names)), vals, color="#7b3294"); ax.set_yticks(range(len(names)), names, fontsize=7); ax.set_xlabel("official Müller regulon specificity score (RSS)"); ax.set_title("P26  Official Müller RSS inventory")
    ax.text(.98, .02, "Trait enrichment not estimable:\n target-gene matrix absent", transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color="#555555")
    savefig(fig, figdir / "Figure_P26_Muller_regulon_enrichment")


def make_p27(figdir):
    fig, ax = plt.subplots(figsize=(10, 3.7)); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_title("P27  Phase 3B conceptual architecture", pad=14)
    boxes = [(0.03, .56, .17, .2, "Global genetic\narchitecture"), (0.25, .56, .17, .2, "Diffuse polygenic\nliability"), (0.47, .56, .18, .2, "Broad retinal\ncellular context"), (0.72, .56, .22, .2, "Orthogonal regulatory\nvalidation")]
    for x,y,w,h,t in boxes:
        ax.add_patch(plt.Rectangle((x,y), w,h, facecolor="#e8f1f8", edgecolor="#2166ac", linewidth=1.5)); ax.text(x+w/2,y+h/2,t,ha="center",va="center",fontsize=10)
    for x in [.20,.42,.67]: ax.annotate("", xy=(x+.04,.66), xytext=(x,.66), arrowprops=dict(arrowstyle="->",lw=1.5,color="#2166ac"))
    ax.plot([.56,.56],[.56,.33],color="#4d9221",lw=2); ax.text(.56,.27,"Müller shared branch\n(AMD · POAG · RE)",ha="center",va="center",color="#276419",fontsize=10)
    ax.plot([.83,.83],[.56,.16],color="#d01c8b",lw=2); ax.text(.83,.09,"RE amacrine comparator",ha="center",va="center",color="#a50f67",fontsize=9)
    savefig(fig, figdir / "Figure_P27_conceptual_architecture")


if __name__ == "__main__": main()
