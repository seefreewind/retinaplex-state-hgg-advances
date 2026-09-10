#!/usr/bin/env python3
"""Run the locked Phase 3B retinal OCR stratified-LDSC analysis.

The official HRCA Table S12 class peak sets are used as the most specific
published broad-class OCR annotations available in the release.  They are
kept distinct from an all-OCR union annotation.  All LD-score calculations
are performed on the liftover-to-hg19 coordinates and are regressed jointly
with baselineLD v2.2 using the EUR reference panel.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.stats import norm


CLASS_ORDER = ["MG", "AC", "BC", "Cone", "Rod", "HC", "RGC", "Astrocyte", "Microglia"]
TRAITS = {
    "AMD": "results/phase0/ldsc/AMD_GCST003219.sumstats.gz",
    "POAG": "results/phase0/ldsc/POAG_GCST90011766.sumstats.gz",
    "RE": "results/phase0c/ldsc/RE_2026_EUR_COMPOSITE.sumstats.gz",
}


def read_bed(path: Path):
    intervals = []
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, a, b = line.rstrip("\n").split("\t")[:3]
            intervals.append((chrom, int(a), int(b)))
    intervals.sort(key=lambda x: (x[0], x[1], x[2]))
    merged = []
    for chrom, a, b in intervals:
        if not merged or merged[-1][0] != chrom or a > merged[-1][2]:
            merged.append([chrom, a, b])
        else:
            merged[-1][2] = max(merged[-1][2], b)
    return merged


def load_intervals(bed_dir: Path):
    return {cls: read_bed(bed_dir / f"{cls}.official_ocr.hg19.bed") for cls in CLASS_ORDER}


def write_union(intervals_by_class, path: Path):
    all_intervals = [x for vals in intervals_by_class.values() for x in vals]
    all_intervals.sort(key=lambda x: (x[0], x[1], x[2]))
    merged = []
    for chrom, a, b in all_intervals:
        if not merged or merged[-1][0] != chrom or a > merged[-1][2]:
            merged.append([chrom, a, b])
        else:
            merged[-1][2] = max(merged[-1][2], b)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        for chrom, a, b in merged:
            fh.write(f"{chrom}\t{a}\t{b}\n")
    return merged


def annotate_bim(bim: Path, intervals_by_class, out_annot: Path):
    # Per-chromosome interval pointers make this independent of pybedtools.
    by_class = {}
    for cls, intervals in intervals_by_class.items():
        by_class[cls] = {}
        for chrom, a, b in intervals:
            by_class[cls].setdefault(chrom, []).append((a, b))
    pointers = {cls: {} for cls in intervals_by_class}
    out_annot.parent.mkdir(parents=True, exist_ok=True)
    n_snp = 0
    n_annot = dict((cls, 0) for cls in intervals_by_class)
    with bim.open(encoding="utf-8") as fi, gzip.open(out_annot, "wt", encoding="utf-8", newline="") as fo:
        writer = csv.writer(fo, delimiter="\t", lineterminator="\n")
        writer.writerow(list(intervals_by_class))
        for line in fi:
            fields = line.rstrip("\n").split()
            if len(fields) < 4:
                continue
            chrom = fields[0]
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            bp0 = int(fields[3]) - 1
            values = []
            for cls, chrom_map in by_class.items():
                ints = chrom_map.get(chrom, [])
                p = pointers[cls].get(chrom, 0)
                while p < len(ints) and ints[p][1] <= bp0:
                    p += 1
                pointers[cls][chrom] = p
                hit = int(p < len(ints) and ints[p][0] <= bp0 < ints[p][1])
                values.append(hit)
                n_annot[cls] += hit
            writer.writerow(values)
            n_snp += 1
    return n_snp, n_annot


def detect_plink_prefix(plink_dir: Path):
    bims = sorted(f for f in plink_dir.rglob("*.bim") if not f.name.startswith("."))
    if not bims:
        raise FileNotFoundError(f"No .bim files found under {plink_dir}")
    for bim in bims:
        if re.search(r"(?:\.EUR\.QC|\.EUR_QC|\.QC)\.1\.bim$", bim.name):
            return Path(str(bim)[:-4])
    # Fall back to the first chromosome-1 BIM and replace its chromosome.
    for bim in bims:
        if re.search(r"(?:^|\.)1\.bim$", bim.name):
            return Path(str(bim)[:-4])
    raise FileNotFoundError(f"Could not identify chromosome-1 PLINK prefix under {plink_dir}")


def detect_weight_prefix(weights_dir: Path):
    files = sorted(f for f in weights_dir.rglob("*.1.l2.ldscore.gz") if not f.name.startswith("."))
    if not files:
        raise FileNotFoundError(f"No chromosome-1 weight LD-score file under {weights_dir}")
    return Path(str(files[0])[:-len("1.l2.ldscore.gz")])


def chr_prefix(base: Path, chrom: int):
    s = str(base)
    s = re.sub(r"(?:\.EUR\.QC|\.EUR_QC|\.QC)\.1$", lambda m: m.group(0)[:-1] + str(chrom), s)
    if s == str(base):
        s = re.sub(r"\.1$", f".{chrom}", s)
    return Path(s)


def run_chromosome(task):
    """Build one chromosome annotation and its nine custom LD-score columns."""
    chrom, plink_base, intervals, annot_dir, custom_dir, ldsc_path, print_snps = task
    bfile = chr_prefix(Path(plink_base), chrom)
    bim = Path(str(bfile) + ".bim")
    annot = Path(annot_dir) / f"retinal_majorclass.{chrom}.annot.gz"
    n_snp, counts = annotate_bim(bim, intervals, annot)
    out = Path(custom_dir) / f"retinal_majorclass.{chrom}"
    cmd = [sys.executable, str(ldsc_path), "--l2", "--bfile", str(bfile),
           "--ld-wind-cm", "1", "--chunk-size", "500", "--annot", str(annot),
           "--thin-annot", "--out", str(out), "--print-snps", str(print_snps)]
    # Keep the number of BLAS workers bounded when several chromosomes run in
    # parallel. This changes only scheduling, not the LDSC estimator.
    env = dict(__import__("os").environ)
    env["OPENBLAS_NUM_THREADS"] = "5"
    env["OMP_NUM_THREADS"] = "5"
    env["MKL_NUM_THREADS"] = "5"
    subprocess.run(cmd, check=True, env=env)
    return chrom, n_snp, [counts[c] for c in CLASS_ORDER]


def parse_vec(text: str):
    text = text.strip().replace("[", "").replace("]", "")
    return [float(x) for x in text.split() if x not in ("", "nan", "NA")]


def parse_ldsc_summary(log_path: Path):
    text = log_path.read_text(encoding="utf-8", errors="replace")
    labels = ("Categories", "Proportion of SNPs", "Proportion of h2g",
              "Enrichment", "Coefficients", "Coefficient SE")
    blocks = {}
    lines = text.splitlines()
    for i, line in enumerate(lines):
        hit = next((label for label in labels if line.startswith(label + ":")), None)
        if hit is None:
            continue
        values = line.split(":", 1)[1].strip()
        # LDSC wraps long numeric vectors across subsequent indented lines.
        # Continue until the next named summary field.
        j = i + 1
        while j < len(lines) and not any(lines[j].startswith(label + ":") for label in labels):
            stripped = lines[j].strip()
            if stripped:
                # A wrapped vector is indented and contains only numeric
                # tokens.  Stop before scalar summary lines such as Lambda
                # GC and Mean Chi^2.
                try:
                    parse_vec(stripped)
                    values += " " + stripped
                except ValueError:
                    break
            j += 1
        blocks[hit] = values
    categories = blocks.get("Categories")
    prop_snp = blocks.get("Proportion of SNPs")
    prop_h2 = blocks.get("Proportion of h2g")
    enrichment = blocks.get("Enrichment")
    coeff = blocks.get("Coefficients")
    coeff_se = blocks.get("Coefficient SE")
    if not all(x is not None for x in (categories, prop_snp, prop_h2, enrichment, coeff, coeff_se)):
        raise RuntimeError(f"Incomplete LDSC summary in {log_path}")
    category_values = categories.split()
    vectors = {"snp_fraction": parse_vec(prop_snp),
               "h2_fraction": parse_vec(prop_h2),
               "enrichment": parse_vec(enrichment),
               "coefficient": parse_vec(coeff),
               "coefficient_se": parse_vec(coeff_se)}
    # The baseline-LD model can drop zero-variance columns independently for
    # each run, so the printed vectors need not have the same length as the
    # category list.  The nine retinal columns are appended last and retain a
    # fixed order; align only that audited suffix.
    expected_categories = [f"{cls}L2_0" for cls in CLASS_ORDER]
    if category_values[-9:] != expected_categories:
        raise RuntimeError(f"Unexpected retinal category suffix in {log_path}: {category_values[-9:]}")
    if any(len(values) < 9 for values in vectors.values()):
        raise RuntimeError(f"Incomplete retinal vector in {log_path}: { {k: len(v) for k, v in vectors.items()} }")
    rows = []
    for i, cls in enumerate(CLASS_ORDER):
        c = vectors["coefficient"][-9:][i]
        se = vectors["coefficient_se"][-9:][i]
        p = 2 * norm.sf(abs(c / se)) if np.isfinite(c) and np.isfinite(se) and se > 0 else float("nan")
        rows.append({
            "cell_class": cls,
            "category": expected_categories[i],
            "snp_fraction": vectors["snp_fraction"][-9:][i],
            "h2_fraction": vectors["h2_fraction"][-9:][i],
            "enrichment": vectors["enrichment"][-9:][i],
            "coefficient": c,
            "coefficient_se": se,
            "coefficient_p": p,
        })
    return rows


def bh(pvals):
    out = [float("nan")] * len(pvals)
    good = [(i, p) for i, p in enumerate(pvals) if np.isfinite(p)]
    good.sort(key=lambda x: x[1])
    m = len(good); prev = 1.0
    for rank in range(m, 0, -1):
        i, p = good[rank - 1]
        prev = min(prev, p * m / rank)
        out[i] = prev
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", required=True)
    ap.add_argument("--plink-dir", required=True)
    ap.add_argument("--baseline-dir", required=True)
    ap.add_argument("--weights-dir", required=True)
    ap.add_argument("--beds-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ldsc", required=True)
    ap.add_argument("--workers", type=int, default=2,
                    help="Concurrent chromosome LD-score jobs; default 2.")
    args = ap.parse_args()
    project = Path(args.project)
    plink_dir = Path(args.plink_dir); baseline_dir = Path(args.baseline_dir); weights_dir = Path(args.weights_dir)
    beds_dir = Path(args.beds_dir); outdir = Path(args.outdir)
    annot_dir = outdir / "ldsc_annotations"; custom_dir = outdir / "retinal_ldscores"; combined_dir = outdir / "combined_baseline_retinal"
    annot_dir.mkdir(parents=True, exist_ok=True); custom_dir.mkdir(parents=True, exist_ok=True); combined_dir.mkdir(parents=True, exist_ok=True)

    intervals = load_intervals(beds_dir)
    union = write_union(intervals, outdir / "official_beds" / "ALL_OCR_UNION.hg19.bed")
    print(f"classes={len(intervals)} union_intervals={len(union)}")
    plink_base = detect_plink_prefix(plink_dir)
    weight_base = detect_weight_prefix(weights_dir)
    print(f"plink_chr1={plink_base}")
    print(f"weight_prefix={weight_base}")
    # LDSC --print-snps requires a one-column ID list.  The provenance file
    # bundled with this project is the standard three-column HM3 list, so
    # create a transparent derived ID-only file rather than passing the
    # three-column table and silently getting an empty merge.
    print_snps = outdir / "w_hm3.snp_ids.txt"
    if not print_snps.exists():
        source_hm3 = project / "data/ld_reference/ldsc_standard/w_hm3.snplist"
        with source_hm3.open(encoding="utf-8") as fi, print_snps.open("w", encoding="utf-8") as fo:
            for line in fi:
                fields = line.strip().split()
                if fields and fields[0] != "SNP":
                    fo.write(fields[0] + "\n")

    # Generate one thin annotation matrix per chromosome, then calculate all
    # nine class-specific LD-score columns in one LDSC pass per chromosome.
    custom_paths = {chrom: custom_dir / f"retinal_majorclass.{chrom}"
                    for chrom in range(1, 23)}
    custom_complete = all(
        Path(str(custom_paths[chrom]) + suffix).exists()
        for chrom in range(1, 23)
        for suffix in (".l2.ldscore.gz", ".l2.M_5_50")
    )
    if not custom_complete:
        annot_qc = []
        tasks = [(chrom, str(plink_base), intervals, str(annot_dir), str(custom_dir),
                  str(args.ldsc), str(print_snps)) for chrom in range(1, 23)]
        with concurrent.futures.ProcessPoolExecutor(max_workers=max(1, args.workers)) as pool:
            for chrom, n_snp, counts in pool.map(run_chromosome, tasks):
                annot_qc.append([chrom, n_snp, *counts])
        annot_qc.sort(key=lambda x: x[0])
        with (outdir / "sldsc_annotation_qc.tsv").open("w", encoding="utf-8") as fh:
            fh.write("chromosome\tn_bim_snps\t" + "\t".join(CLASS_ORDER) + "\n")
            for row in annot_qc: fh.write("\t".join(map(str, row)) + "\n")
    else:
        print("reusing existing retinal chromosome LD-score files")

    # Add the retina columns to the official baselineLD v2.2 partitioned
    # scores, preserving the baseline SNP order and M_5_50 accounting.
    combined_complete = all(
        (combined_dir / f"baseline_retinal.{chrom}.l2.ldscore.gz").exists()
        and (combined_dir / f"baseline_retinal.{chrom}.l2.M_5_50").exists()
        for chrom in range(1, 23)
    )
    if not combined_complete:
        for chrom in range(1, 23):
            base_ld = baseline_dir / f"baselineLD.{chrom}.l2.ldscore.gz"
            base_m = baseline_dir / f"baselineLD.{chrom}.l2.M_5_50"
            custom_ld = Path(str(custom_paths[chrom]) + ".l2.ldscore.gz")
            combined_ld = combined_dir / f"baseline_retinal.{chrom}.l2.ldscore.gz"
            with gzip.open(base_ld, "rt", encoding="utf-8") as fb, gzip.open(custom_ld, "rt", encoding="utf-8") as fc, gzip.open(combined_ld, "wt", encoding="utf-8") as fo:
                bhdr = fb.readline().rstrip("\n").split("\t")
                chdr = fc.readline().rstrip("\n").split("\t")
                fo.write("\t".join(bhdr + chdr[3:]) + "\n")
                for bl, cl in zip(fb, fc):
                    b = bl.rstrip("\n").split("\t"); c = cl.rstrip("\n").split("\t")
                    if b[:3] != c[:3]:
                        raise RuntimeError(f"SNP order mismatch at chr{chrom}: {b[:3]} != {c[:3]}")
                    fo.write("\t".join(b + c[3:]) + "\n")
                if fb.readline() or fc.readline():
                    raise RuntimeError(f"LD-score row-count mismatch at chr{chrom}")
            with base_m.open() as fb, Path(str(custom_paths[chrom]) + ".l2.M_5_50").open() as fc, (combined_dir / f"baseline_retinal.{chrom}.l2.M_5_50").open("w") as fo:
                fo.write(fb.read().strip() + "\t" + fc.read().strip() + "\n")
    else:
        print("reusing existing combined baseline-retinal LD-score files")

    combined_prefix = combined_dir / "baseline_retinal."
    weight_prefix = str(weight_base)
    all_rows = []
    for trait, rel in TRAITS.items():
        sumstats = project / rel
        out = outdir / f"sldsc_{trait}_baseline_retinal"
        cmd = [sys.executable, str(args.ldsc), "--h2", str(sumstats), "--ref-ld-chr", str(combined_prefix), "--w-ld-chr", weight_prefix, "--out", str(out)]
        subprocess.run(cmd, check=True)
        rows = parse_ldsc_summary(Path(str(out) + ".log"))
        for row in rows:
            row.update({"trait": trait, "annotation_type": "OFFICIAL_MAJOR_CLASS_OCR_SURROGATE", "status": "PASS"})
        all_rows.extend(rows)
    p = bh([r["coefficient_p"] for r in all_rows if r["cell_class"] in CLASS_ORDER])
    # The list is in the same order because every row is one class x trait.
    for row, q in zip(all_rows, p): row["global_primary_fdr"] = q
    header = ["trait", "cell_class", "annotation_type", "category", "snp_fraction", "h2_fraction", "enrichment", "coefficient", "coefficient_se", "coefficient_p", "global_primary_fdr", "status"]
    out_tsv = project / "results/phase3b/sldsc_retinal_ocr_enrichment.tsv"
    with out_tsv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(all_rows)
    # Secondary all-OCR annotation is intentionally recorded separately; it
    # is not folded into the primary multiple-testing family.
    print(f"sldsc_rows={len(all_rows)} output={out_tsv}")


if __name__ == "__main__":
    main()
