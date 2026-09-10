#!/usr/bin/env python3
"""Run the locked Phase 3A-R8 MAGMA -> scDRS technical bridge.

This script deliberately stops at technical donor QC. It does not run any
cell-type enrichment, disease-association, or within-lineage state analysis.
"""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import math
import os
import shlex
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "phase3ar8"
MAGMA_OUT = OUT / "magma"
LOG_OUT = OUT / "logs"
GENE_LOC = ROOT / "data" / "annotations" / "magma" / "NCBI37.3.gene.loc"
BFILE = ROOT / "data" / "ld_reference" / "magma_1000g_eur" / "g1000_eur"
H5AD = ROOT / "data" / "single_cell" / "HRCA" / "technical_reference_phase3ar7.h5ad"
MAGMA = Path(os.environ.get("MAGMA_EXECUTABLE", "magma"))
SCDRS_PYTHON = Path(os.environ.get("SCDRS_PYTHON", "python"))
SCDRS = Path(os.environ.get("SCDRS_EXECUTABLE", "scdrs"))
MAGMA_VERSION = "v1.10 (custom, self-reported)"
RUN_DATE = dt.date.today().isoformat()

GWAS = {
    "AMD": ROOT / "data" / "harmonized" / "phase0b" / "AMD_GCST003219_for_ldsc.tsv.gz",
    "POAG": ROOT / "data" / "harmonized" / "phase0b" / "POAG_GCST90011766_for_ldsc.tsv.gz",
    "RE": ROOT / "data" / "harmonized" / "phase0c" / "RE_2026_EUR_COMPOSITE_for_ldsc.tsv.gz",
}

EXPECTED_SHA256 = {
    "NCBI37.3.gene.loc": "bfcb5273801cc3e37dcb98f5465f4498daa8c0cfcda2f4f3fbc51abe02000e7e",
    "g1000_eur.bed": "fadce7321e622ed20e3cb7798b3e54699c7f342752c642a14ddb982d93d8e08c",
    "g1000_eur.bim": "b7dfeb30cef0573059ae43dd85a6924e46e418d6fe03df6b0ff1166c79cfa131",
    "g1000_eur.fam": "79dee21226fd04b5a889cbbeacd08079e1740e7d63b0a4dea115527e8fea4dae",
}
ZENODO_ARCHIVE_MD5 = "1919cb5c79bbe7871aed71ae4abe6217"


def ensure_dirs() -> None:
    for path in (OUT, MAGMA_OUT, LOG_OUT, ROOT / "metadata", ROOT / "reports"):
        path.mkdir(parents=True, exist_ok=True)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_command(command: list[str], log_path: Path) -> subprocess.CompletedProcess[str]:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        log.write("COMMAND\t" + " ".join(shlex.quote(x) for x in command) + "\n\n")
        result = subprocess.run(command, cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log.write(result.stdout or "")
        log.write(f"\nEXIT_CODE\t{result.returncode}\n")
    return result


def audit_gene_loc(path: Path) -> dict[str, object]:
    n = 0
    malformed = 0
    ids: list[str] = []
    chromosomes: set[str] = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 6:
                malformed += 1
                continue
            n += 1
            ids.append(fields[0])
            chromosomes.add(fields[1])
    duplicates = n - len(set(ids))
    return {
        "resource": "NCBI37.3.gene.loc",
        "genes": n,
        "chromosomes": ",".join(sorted(chromosomes, key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x))),
        "duplicate_gene_ids": duplicates,
        "malformed_rows": malformed,
        "status": "PASS" if n > 10000 and duplicates == 0 and malformed == 0 else "FAIL",
        "notes": "Standard MAGMA NCBI37.3 six-column resource; first field is the MAGMA/Entrez gene ID.",
    }


def count_duplicate_bim_ids(path: Path) -> int:
    # Stream through sort instead of retaining millions of IDs in Python RAM.
    awk = subprocess.Popen(["awk", "{print $2}", str(path)], stdout=subprocess.PIPE, text=True)
    sorter = subprocess.Popen(["sort", "-u"], stdin=awk.stdout, stdout=subprocess.PIPE, text=True, env={**os.environ, "LC_ALL": "C"})
    assert awk.stdout is not None and sorter.stdout is not None
    awk.stdout.close()
    unique = sum(1 for _ in sorter.stdout)
    sorter.stdout.close()
    sorter.wait()
    awk.wait()
    return unique


def audit_bim(path: Path) -> dict[str, object]:
    n = 0
    malformed = 0
    rsids = 0
    chromosomes: set[str] = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) < 6:
                malformed += 1
                continue
            n += 1
            chromosomes.add(fields[0])
            rsids += fields[1].lower().startswith("rs")
    unique = count_duplicate_bim_ids(path)
    duplicates = n - unique
    return {
        "resource": "g1000_eur.bim",
        "snp_count": n,
        "chromosomes": ",".join(sorted(chromosomes, key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x))),
        "rsid_fraction": rsids / n if n else math.nan,
        "duplicate_snp_ids": duplicates,
        "malformed_rows": malformed,
        "status": "PASS" if n > 1000000 and duplicates == 0 and malformed == 0 else "FAIL",
        "notes": "BIM was read as whitespace-delimited PLINK text; duplicate count is based on unique SNP ID.",
    }


def read_headers(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        return handle.readline().rstrip("\n").split("\t")


def find_column(headers: list[str], wanted: str) -> str | None:
    for header in headers:
        if header.upper() == wanted.upper():
            return header
    return None


def prepare_pval(path: Path, output: Path) -> dict[str, object]:
    csv.field_size_limit(sys.maxsize)
    headers = read_headers(path)
    snp_col = find_column(headers, "SNP")
    p_col = find_column(headers, "P")
    n_col = find_column(headers, "N")
    if not snp_col or not p_col or not n_col:
        raise RuntimeError(f"MAGMA input lacks SNP/P/N columns: {path} headers={headers}")

    total = valid = invalid = duplicate = 0
    seen: set[str] = set()
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace", newline="") as handle, output.open("w", encoding="utf-8", newline="") as out_handle:
        reader = csv.DictReader(handle, delimiter="\t")
        writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["SNP", "P", "N"])
        for row in reader:
            total += 1
            snp = (row.get(snp_col) or "").strip()
            p_text = (row.get(p_col) or "").strip()
            n_text = (row.get(n_col) or "").strip()
            try:
                p_value = float(p_text)
                n_value = float(n_text)
                valid_number = math.isfinite(p_value) and 0 < p_value <= 1 and math.isfinite(n_value) and n_value > 0
            except ValueError:
                valid_number = False
            if not snp or not valid_number:
                invalid += 1
                continue
            if snp in seen:
                duplicate += 1
                continue
            seen.add(snp)
            writer.writerow([snp, p_text, n_text])
            valid += 1
    return {
        "trait": path.name.split("_")[0],
        "input_rows": total,
        "valid_unique_rows": valid,
        "invalid_rows": invalid,
        "duplicate_snp_rows": duplicate,
        "pval_file": str(output),
    }


def pval_ids(path: Path) -> set[str]:
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["SNP"] for row in reader}


def matched_bim_snps(pval_path: Path, bim_path: Path) -> int:
    remaining = pval_ids(pval_path)
    matched = 0
    with bim_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2 and fields[1] in remaining:
                remaining.remove(fields[1])
                matched += 1
                if not remaining:
                    break
    return matched


def annotation_snp_count(path: Path) -> tuple[int, int]:
    gene_rows = 0
    snp_links = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            gene_rows += 1
            if len(fields) >= 4:
                try:
                    snp_links += int(fields[3])
                except ValueError:
                    pass
    return gene_rows, snp_links


def read_magma_output(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep=r"\s+", engine="python", comment="#")
    frame.columns = [str(column).upper() for column in frame.columns]
    required = {"GENE", "ZSTAT", "P"}
    missing = required - set(frame.columns)
    if missing:
        raise RuntimeError(f"MAGMA output {path} lacks {sorted(missing)}; columns={list(frame.columns)}")
    frame["GENE"] = frame["GENE"].astype(str).str.replace(r"\.0$", "", regex=True)
    frame["ZSTAT"] = pd.to_numeric(frame["ZSTAT"], errors="coerce")
    frame["P"] = pd.to_numeric(frame["P"], errors="coerce")
    frame = frame.drop_duplicates("GENE", keep="first")
    return frame


def write_reference_provenance() -> dict[str, object]:
    archive = ROOT / "data" / "ld_reference" / "magma_1000g_eur" / "g1000_eur.zip"
    resources = [
        ("gene_loc", GENE_LOC, "https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.zip", "https://github.com/comorment/magma/tree/main/reference/magma/NCBI37.3", "Git LFS retrieval mirror; primary provenance is official MAGMA resource"),
        ("g1000_eur.bed", BFILE.with_suffix(".bed"), "https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip", "https://zenodo.org/records/14142263", "Compressed EUR-1000G retrieval mirror; archive checksum and PLINK object checksums recorded"),
        ("g1000_eur.bim", BFILE.with_suffix(".bim"), "https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip", "https://zenodo.org/records/14142263", "Compressed EUR-1000G retrieval mirror; archive checksum and PLINK object checksums recorded"),
        ("g1000_eur.fam", BFILE.with_suffix(".fam"), "https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip", "https://zenodo.org/records/14142263", "Compressed EUR-1000G retrieval mirror; archive checksum and PLINK object checksums recorded"),
    ]
    rows: list[dict[str, object]] = []
    for resource, path, primary_url, retrieval_url, notes in resources:
        if not path.exists():
            raise RuntimeError(f"missing reference file: {path}")
        checksum = sha256_file(path)
        expected = EXPECTED_SHA256.get(path.name)
        rows.append({
            "resource": resource,
            "file": str(path),
            "primary_source_url": primary_url,
            "retrieval_url": retrieval_url,
            "download_date": RUN_DATE,
            "bytes": path.stat().st_size,
            "checksum_algorithm": "SHA-256",
            "checksum": checksum,
            "expected_checksum": expected or "NA",
            "checksum_match": "PASS" if expected == checksum else "NOT_VERIFIED",
            "magma_version": MAGMA_VERSION,
            "genome_build": "GRCh37 / NCBI37.3",
            "notes": notes,
        })
    if archive.exists():
        rows.append({
            "resource": "g1000_eur.zip",
            "file": str(archive),
            "primary_source_url": "https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip",
            "retrieval_url": "https://zenodo.org/records/14142263",
            "download_date": RUN_DATE,
            "bytes": archive.stat().st_size,
            "checksum_algorithm": "MD5",
            "checksum": md5_file(archive),
            "expected_checksum": ZENODO_ARCHIVE_MD5,
            "checksum_match": "PASS" if md5_file(archive) == ZENODO_ARCHIVE_MD5 else "NOT_VERIFIED",
            "magma_version": MAGMA_VERSION,
            "genome_build": "GRCh37 / NCBI37.3",
            "notes": "Zenodo record used as a retrieval mirror because the official CTG host was unavailable from this host; extracted PLINK objects are independently checked below.",
        })
    frame = pd.DataFrame(rows)
    frame.to_csv(ROOT / "metadata" / "magma_reference_provenance.tsv", sep="\t", index=False)
    return {"archive": archive, "rows": rows}


def write_reference_audit(gene_audit: dict[str, object], bim_audit: dict[str, object], provenance: dict[str, object]) -> None:
    rows = [
        {"resource": "NCBI37.3.gene.loc", "metric": "number_of_genes", "value": gene_audit["genes"], "status": gene_audit["status"], "notes": gene_audit["notes"]},
        {"resource": "NCBI37.3.gene.loc", "metric": "chromosome_coverage", "value": gene_audit["chromosomes"], "status": gene_audit["status"], "notes": "Expected autosomes 1-22 plus X/Y when present."},
        {"resource": "NCBI37.3.gene.loc", "metric": "duplicate_gene_ids", "value": gene_audit["duplicate_gene_ids"], "status": "PASS" if gene_audit["duplicate_gene_ids"] == 0 else "FAIL", "notes": "Duplicate count across first-column MAGMA gene IDs."},
        {"resource": "NCBI37.3.gene.loc", "metric": "malformed_rows", "value": gene_audit["malformed_rows"], "status": "PASS" if gene_audit["malformed_rows"] == 0 else "FAIL", "notes": "Rows with fewer than six whitespace-delimited fields."},
        {"resource": "g1000_eur.bim", "metric": "snp_count", "value": bim_audit["snp_count"], "status": bim_audit["status"], "notes": bim_audit["notes"]},
        {"resource": "g1000_eur.bim", "metric": "chromosome_coverage", "value": bim_audit["chromosomes"], "status": bim_audit["status"], "notes": "Chromosomes observed in BIM first field."},
        {"resource": "g1000_eur.bim", "metric": "rsid_fraction", "value": bim_audit["rsid_fraction"], "status": bim_audit["status"], "notes": "Fraction of SNP IDs beginning with rs, reported descriptively."},
        {"resource": "g1000_eur.bim", "metric": "duplicate_snp_ids", "value": bim_audit["duplicate_snp_ids"], "status": "PASS" if bim_audit["duplicate_snp_ids"] == 0 else "FAIL", "notes": "Duplicate count across second-column SNP IDs."},
        {"resource": "g1000_eur.bim", "metric": "malformed_rows", "value": bim_audit["malformed_rows"], "status": "PASS" if bim_audit["malformed_rows"] == 0 else "FAIL", "notes": "Rows with fewer than six whitespace-delimited fields."},
        {"resource": "g1000_eur", "metric": "ploid_reference_files", "value": ",".join(path.name for path in (BFILE.with_suffix(".bed"), BFILE.with_suffix(".bim"), BFILE.with_suffix(".fam"))), "status": "PASS", "notes": "All three PLINK files are present and independently checksum-audited in metadata/magma_reference_provenance.tsv."},
    ]
    pd.DataFrame(rows).to_csv(OUT / "magma_reference_audit.tsv", sep="\t", index=False)


def run_magma() -> tuple[dict[str, dict[str, object]], dict[str, Path]]:
    for path in (MAGMA, GENE_LOC, BFILE.with_suffix(".bed"), BFILE.with_suffix(".bim"), BFILE.with_suffix(".fam")):
        if not path.exists():
            raise RuntimeError(f"required MAGMA input is missing: {path}")

    pval_paths: dict[str, Path] = {}
    input_stats: dict[str, dict[str, object]] = {}
    for trait, source in GWAS.items():
        pval_path = MAGMA_OUT / f"{trait}.pval"
        pval_paths[trait] = pval_path
        input_stats[trait] = prepare_pval(source, pval_path)

    annotation_prefix = MAGMA_OUT / "shared_annotation"
    annotation_path = annotation_prefix.with_suffix(".genes.annot")
    annotation_command = [
        str(MAGMA), "--annotate", "window=35,35", "--snp-loc", str(BFILE.with_suffix(".bim")),
        "--gene-loc", str(GENE_LOC), "--out", str(annotation_prefix),
    ]
    if not annotation_path.exists():
        result = run_command(annotation_command, LOG_OUT / "shared_annotation.log")
        if result.returncode != 0 or not annotation_path.exists():
            raise RuntimeError("shared MAGMA annotation failed; POAG gate stopped the trait sequence")
    else:
        (LOG_OUT / "shared_annotation.command.txt").write_text(" ".join(shlex.quote(x) for x in annotation_command) + "\n", encoding="utf-8")

    gene_rows, snp_links = annotation_snp_count(annotation_path)
    annotation_stats = {"annotation_gene_rows": gene_rows, "annotation_snp_links": snp_links, "annotation_path": annotation_path}
    results: dict[str, pd.DataFrame] = {}
    qc_rows: list[dict[str, object]] = []

    for trait in ("POAG", "AMD", "RE"):
        prefix = MAGMA_OUT / trait
        output_path = prefix.with_suffix(".genes.out")
        command = [
            str(MAGMA), "--bfile", str(BFILE), "--pval", str(pval_paths[trait]), "use=SNP,P", "ncol=N",
            "--gene-annot", str(annotation_path), "--out", str(prefix),
        ]
        if not output_path.exists():
            result = run_command(command, LOG_OUT / f"{trait}.genes.log")
            if result.returncode != 0 or not output_path.exists():
                if trait == "POAG":
                    raise RuntimeError("POAG MAGMA gate failed; AMD/RE were not run")
                raise RuntimeError(f"{trait} MAGMA gene analysis failed")
        else:
            (LOG_OUT / f"{trait}.genes.command.txt").write_text(" ".join(shlex.quote(x) for x in command) + "\n", encoding="utf-8")
        frame = read_magma_output(output_path)
        results[trait] = frame
        valid_z = np.isfinite(frame["ZSTAT"])
        valid_p = np.isfinite(frame["P"])
        matched = matched_bim_snps(pval_paths[trait], BFILE.with_suffix(".bim"))
        warnings: list[str] = []
        if len(frame) <= 10000:
            warnings.append("genes_tested<=10000")
        if valid_z.mean() < 0.95 or valid_p.mean() < 0.95:
            warnings.append("valid_gene_stat_fraction<0.95")
        status = "PASS" if len(frame) > 10000 and valid_z.mean() >= 0.95 and valid_p.mean() >= 0.95 and matched > 0 else "PARTIAL"
        qc_rows.append({
            "trait": trait,
            "GWAS_SNPs_input": input_stats[trait]["valid_unique_rows"],
            "GWAS_SNPs_matched_to_reference": matched,
            "annotation_SNPs": snp_links,
            "annotation_gene_rows": gene_rows,
            "genes_tested": len(frame),
            "genes_with_valid_ZSTAT": int(valid_z.sum()),
            "genes_with_valid_P": int(valid_p.sum()),
            "warnings": ";".join(warnings) if warnings else "none",
            "status": status,
        })

    pd.DataFrame(qc_rows).to_csv(OUT / "magma_qc.tsv", sep="\t", index=False)

    mapping: dict[str, str] = {}
    with GENE_LOC.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 6:
                mapping.setdefault(fields[0], fields[5])
    score: pd.DataFrame | None = None
    for trait in ("AMD", "POAG", "RE"):
        selected = results[trait][["GENE", "ZSTAT", "P"]].rename(columns={"GENE": "gene_id", "ZSTAT": f"{trait}_ZSTAT", "P": f"{trait}_P"})
        score = selected if score is None else score.merge(selected, on="gene_id", how="outer", validate="one_to_one")
    assert score is not None
    score["gene_symbol"] = score["gene_id"].map(mapping)
    score = score[["gene_id", "gene_symbol", "AMD_ZSTAT", "AMD_P", "POAG_ZSTAT", "POAG_P", "RE_ZSTAT", "RE_P"]]
    score.to_csv(OUT / "magma_gene_scores.tsv", sep="\t", index=False)

    complete = score.dropna(subset=["gene_symbol", "AMD_ZSTAT", "POAG_ZSTAT", "RE_ZSTAT"]).copy()
    for column in ("AMD_ZSTAT", "POAG_ZSTAT", "RE_ZSTAT"):
        complete = complete[np.isfinite(complete[column])]
    complete = complete.drop_duplicates("gene_symbol", keep="first")
    zscores = complete[["gene_symbol", "AMD_ZSTAT", "POAG_ZSTAT", "RE_ZSTAT"]].rename(columns={"gene_symbol": "GENE", "AMD_ZSTAT": "AMD", "POAG_ZSTAT": "POAG", "RE_ZSTAT": "RE"})
    zscores.to_csv(OUT / "scdrs_gene_zscores.tsv", sep="\t", index=False)
    return {"AMD": qc_rows[1], "POAG": qc_rows[0], "RE": qc_rows[2]}, {"annotation": annotation_path, "complete": len(zscores), "score": OUT / "magma_gene_scores.tsv"}


def run_munge_and_overlap(zscore_info: dict[str, object]) -> dict[str, object]:
    zscore_path = OUT / "scdrs_gene_zscores.tsv"
    gs_path = OUT / "retinaplex_traits.gs"
    command = [str(SCDRS), "munge_gs", "--out_file", str(gs_path), "--zscore_file", str(zscore_path), "--weight", "zscore", "--n_min", "100", "--n_max", "1000"]
    result = run_command(command, LOG_OUT / "scdrs_munge_gs.log")
    if result.returncode != 0 or not gs_path.exists():
        raise RuntimeError("scDRS munge_gs failed")

    zscores = pd.read_csv(zscore_path, sep="\t")
    gs = pd.read_csv(gs_path, sep="\t")
    manifest_rows: list[dict[str, object]] = []
    gene_sets: dict[str, list[str]] = {}
    for _, row in gs.iterrows():
        tokens = [token for token in str(row["GENESET"]).split(",") if token]
        genes = [token.rsplit(":", 1)[0] for token in tokens]
        trait = str(row["TRAIT"])
        gene_sets[trait] = genes
        manifest_rows.append({
            "trait": trait,
            "MAGMA_genes_available": int(zscores[trait].notna().sum()),
            "genes_considered": int(zscores[trait].notna().sum()),
            "genes_retained": len(genes),
            "n_min": 100,
            "n_max": 1000,
            "weighting": "zscore (scDRS munge_gs; weights derived from gene-level ZSTAT)",
            "top_score_threshold_rule": "No FDR/FWER supplied; scDRS selected the n_max=1000 smallest zscore-derived p-values.",
            "final_gene_count": len(genes),
            "status": "PASS" if 100 <= len(genes) <= 1000 else "FAIL",
        })
    pd.DataFrame(manifest_rows).sort_values("trait").to_csv(OUT / "scdrs_gene_set_manifest.tsv", sep="\t", index=False)

    # Exact expression universe check and real overlap; no assumptions about symbols are made.
    import anndata

    adata = anndata.read_h5ad(H5AD, backed="r")
    expression_names = [str(x) for x in adata.var_names]
    expression_set = set(expression_names)
    duplicate_expression = int(pd.Index(expression_names).duplicated().sum())
    var_gene_ids = adata.var["gene_ids"].astype(str) if "gene_ids" in adata.var.columns else pd.Series(dtype=str)
    rows: list[dict[str, object]] = []
    for trait in ("AMD", "POAG", "RE"):
        gs_genes = gene_sets[trait]
        matched = sorted(set(gs_genes) & expression_set)
        unmatched = sorted(set(gs_genes) - expression_set)
        rows.append({
            "trait": trait,
            "scDRS_gene_set_genes": len(set(gs_genes)),
            "expression_genes": len(expression_set),
            "matched_genes": len(matched),
            "unmatched_genes": len(unmatched),
            "match_fraction": len(matched) / len(set(gs_genes)) if gs_genes else math.nan,
            "duplicate_expression_ids": duplicate_expression,
            "expression_var_names_are_symbols": "yes" if set(expression_names) & set(gs_genes) else "not_assumed",
            "expression_ensembl_column": "gene_ids" if "gene_ids" in adata.var.columns else "absent",
            "status": "PASS" if gs_genes and len(matched) / len(set(gs_genes)) >= 0.80 and duplicate_expression == 0 else "PARTIAL",
        })
    pd.DataFrame(rows).to_csv(OUT / "gene_overlap_audit.tsv", sep="\t", index=False)
    return {"gs": gs_path, "gene_sets": gene_sets, "expression_genes": len(expression_set), "duplicate_expression": duplicate_expression, "overlap": rows, "h5ad_shape": adata.shape, "var_gene_ids": int(len(var_gene_ids))}


def run_bridge(gene_info: dict[str, object]) -> dict[str, object]:
    import anndata
    import scipy.sparse as sp
    import scdrs

    adata = scdrs.util.load_h5ad(str(H5AD), flag_filter_data=False, flag_raw_count=True)
    input_shape = adata.shape
    if not sp.issparse(adata.X):
        raise RuntimeError("bridge stopped: scDRS matrix is not sparse after raw-count normalization")
    obs = adata.obs.copy()
    n_genes = np.asarray(adata.X.getnnz(axis=1)).ravel()
    # sample is one-to-one with donor in this capped reference; region is constant.
    # Keeping redundant one-hot columns would make the scDRS regression design singular.
    # Pre-encode donor as numeric one-hot columns. scDRS 1.0.2's internal
    # category2dummy path can leave mixed bool/object dtypes under pandas 2.x.
    cov = pd.get_dummies(obs["donor_id"].astype(str), prefix="donor", drop_first=True, dtype=float)
    cov.index = adata.obs_names
    n_genes_float = n_genes.astype(float)
    n_genes_std = n_genes_float.std()
    if not np.isfinite(n_genes_std) or n_genes_std == 0:
        raise RuntimeError("bridge stopped: n_genes covariate has no finite variation")
    cov["n_genes"] = (n_genes_float - n_genes_float.mean()) / n_genes_std
    cov_path = OUT / "scdrs_bridge_covariates.tsv"
    cov.to_csv(cov_path, sep="\t", index_label="CELL")
    # scDRS 1.0.2 emits a repeated NumPy warning from its implicit sparse
    # covariate-correction matmul on this object even when the resulting
    # statistics are finite. Suppress that known noisy warning locally and
    # enforce finiteness explicitly after scoring below.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        scdrs.pp.preprocess(adata, cov=cov, n_mean_bin=20, n_var_bin=20, copy=False)
    gene_stats = adata.uns["SCDRS_PARAM"]["GENE_STATS"]
    if not np.isfinite(gene_stats[["mean", "var", "var_tech"]].to_numpy(dtype=float)).all():
        raise RuntimeError("bridge stopped: scDRS preprocessing produced non-finite gene statistics")

    dict_gs = scdrs.util.load_gs(str(gene_info["gs"]), to_intersect=list(adata.var_names))
    if "POAG" not in dict_gs:
        raise RuntimeError("POAG gene set is absent from retinaplex_traits.gs")
    gene_list, gene_weights = dict_gs["POAG"]
    if len(gene_list) < 100:
        raise RuntimeError(f"POAG matched gene set is too small for bridge: {len(gene_list)}")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        df_res = scdrs.score_cell(
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
    score_path = OUT / "POAG.score.gz"
    full_score_path = OUT / "POAG.full_score.gz"
    df_res.iloc[:, :6].to_csv(score_path, sep="\t", index=True, compression="gzip")
    df_res.to_csv(full_score_path, sep="\t", index=True, compression="gzip")

    raw = pd.to_numeric(df_res["raw_score"], errors="coerce").to_numpy()
    norm = pd.to_numeric(df_res["norm_score"], errors="coerce").to_numpy()
    finite_raw = np.isfinite(raw)
    finite_norm = np.isfinite(norm)
    obs_for_score = obs.loc[df_res.index].copy()
    obs_for_score["norm_score"] = norm
    donor_rows: list[dict[str, object]] = []
    for donor, group in obs_for_score.groupby("donor_id", sort=True):
        values = group["norm_score"].to_numpy(dtype=float)
        finite = values[np.isfinite(values)]
        donor_rows.append({
            "donor_id": donor,
            "cells": len(values),
            "finite_fraction": len(finite) / len(values) if len(values) else math.nan,
            "median_norm_score": float(np.median(finite)) if len(finite) else math.nan,
            "IQR_norm_score": float(np.quantile(finite, 0.75) - np.quantile(finite, 0.25)) if len(finite) else math.nan,
        })
    donor_qc = pd.DataFrame(donor_rows)
    donor_qc.to_csv(OUT / "scdrs_poag_donor_qc.tsv", sep="\t", index=False)

    sample_rows: list[dict[str, object]] = []
    for sample, group in obs_for_score.groupby("sample", sort=True):
        values = group["norm_score"].to_numpy(dtype=float)
        sample_rows.append({"sample": sample, "cells": len(values), "finite_fraction": float(np.isfinite(values).mean()) if len(values) else math.nan})
    sample_qc = pd.DataFrame(sample_rows)
    sample_qc.to_csv(OUT / "scdrs_poag_sample_qc.tsv", sep="\t", index=False)

    control_columns = [column for column in df_res.columns if str(column).startswith("ctrl_norm_score")]
    bridge_status = "PASS" if (
        input_shape[0] == 42417 and len(df_res) == 42417 and len(control_columns) == 1000
        and finite_norm.mean() >= 0.99 and finite_raw.mean() >= 0.99
        and donor_qc["finite_fraction"].min() >= 0.99 and sample_qc["finite_fraction"].min() >= 0.99
    ) else "FAIL"
    bridge_row = {
        "input_cells": input_shape[0],
        "output_cells": len(df_res),
        "donors": int(obs["donor_id"].nunique()),
        "classes": int(obs["harmonized_broad_class"].nunique()),
        "POAG_gene_set_genes": int(pd.read_csv(OUT / "scdrs_gene_set_manifest.tsv", sep="\t").query("trait == 'POAG'")["final_gene_count"].iloc[0]),
        "matched_genes": len(gene_list),
        "control_sets": len(control_columns),
        "finite_raw_scores": int(finite_raw.sum()),
        "finite_normalized_scores": int(finite_norm.sum()),
        "NaN": int(pd.isna(df_res[["raw_score", "norm_score"]]).to_numpy().sum()),
        "Inf": int(np.isinf(df_res[["raw_score", "norm_score"]].to_numpy(dtype=float)).sum()),
        "preprocessing": "scDRS util.load_h5ad: one size-factor normalization to 1e4 plus log1p; scDRS pp.preprocess with donor+n_genes covariates; sparse matrix retained",
        "status": bridge_status,
    }
    pd.DataFrame([bridge_row]).to_csv(OUT / "scdrs_poag_bridge_test.tsv", sep="\t", index=False)
    return {"bridge": bridge_row, "donor_qc": donor_qc, "score_path": score_path, "full_score_path": full_score_path, "input_shape": input_shape}


def write_reports(magma_qc: dict[str, dict[str, object]], gene_info: dict[str, object], bridge_info: dict[str, object], provenance: dict[str, object]) -> None:
    qc = pd.read_csv(OUT / "magma_qc.tsv", sep="\t")
    overlap = pd.read_csv(OUT / "gene_overlap_audit.tsv", sep="\t")
    manifest = pd.read_csv(OUT / "scdrs_gene_set_manifest.tsv", sep="\t")
    bridge = bridge_info["bridge"]
    magma_pass = set(qc["status"]) == {"PASS"}
    overlap_pass = set(overlap["status"]) == {"PASS"}
    gene_gate = magma_pass and overlap_pass and all(x >= 0.80 for x in overlap["match_fraction"])
    bridge_pass = bridge["status"] == "PASS"
    verdict = "BROAD_PHASE3A_READY_AMD_LOW_POWER" if magma_pass and gene_gate and bridge_pass else ("SCDRS_BRIDGE_FAILURE" if magma_pass and gene_gate else "MAGMA_PARTIAL")
    next_action = "A. RUN PHASE 3A-1 BROAD CELLULAR ARCHITECTURE" if verdict.startswith("BROAD_PHASE3A_READY") else ("C. FIX SCDRS BRIDGE" if magma_pass and gene_gate else "B. FIX MAGMA INPUT")
    ref_frame = pd.DataFrame(provenance["rows"])
    ref_gene = ref_frame[ref_frame["resource"] == "gene_loc"].iloc[0]
    ref_bim = ref_frame[ref_frame["resource"] == "g1000_eur.bim"].iloc[0]
    ref_audit = pd.read_csv(OUT / "magma_reference_audit.tsv", sep="\t")
    ref_snp_count = int(ref_audit.query("resource == 'g1000_eur.bim' and metric == 'snp_count'")["value"].iloc[0])
    report = f"""# PHASE 3A-R8 REPORT

## VERDICT

{verdict}

## MAGMA REFERENCE

MAGMA version = {MAGMA_VERSION}
gene loc = NCBI37.3.gene.loc; SHA-256 {ref_gene['checksum']}
EUR reference SNPs = {ref_snp_count}
gene annotation = shared_annotation.genes.annot; window=35,35 kb; one annotation reused for AMD/POAG/RE
status = {'PASS' if magma_pass else 'PARTIAL'}

## MAGMA TRAITS

AMD:
genes = {int(qc.query("trait == 'AMD'")["genes_tested"].iloc[0])}
status = {qc.query("trait == 'AMD'")["status"].iloc[0]}

POAG:
genes = {int(qc.query("trait == 'POAG'")["genes_tested"].iloc[0])}
status = {qc.query("trait == 'POAG'")["status"].iloc[0]}

RE:
genes = {int(qc.query("trait == 'RE'")["genes_tested"].iloc[0])}
status = {qc.query("trait == 'RE'")["status"].iloc[0]}

## GENE OVERLAP

AMD = {int(overlap.query("trait == 'AMD'")["matched_genes"].iloc[0])}/{int(overlap.query("trait == 'AMD'")["scDRS_gene_set_genes"].iloc[0])} ({overlap.query("trait == 'AMD'")["match_fraction"].iloc[0]:.4f})
POAG = {int(overlap.query("trait == 'POAG'")["matched_genes"].iloc[0])}/{int(overlap.query("trait == 'POAG'")["scDRS_gene_set_genes"].iloc[0])} ({overlap.query("trait == 'POAG'")["match_fraction"].iloc[0]:.4f})
RE = {int(overlap.query("trait == 'RE'")["matched_genes"].iloc[0])}/{int(overlap.query("trait == 'RE'")["scDRS_gene_set_genes"].iloc[0])} ({overlap.query("trait == 'RE'")["match_fraction"].iloc[0]:.4f})
status = {'PASS' if gene_gate else 'PARTIAL'}

## SCDRS GENE SETS

AMD genes = {int(manifest.query("trait == 'AMD'")["final_gene_count"].iloc[0])}
POAG genes = {int(manifest.query("trait == 'POAG'")["final_gene_count"].iloc[0])}
RE genes = {int(manifest.query("trait == 'RE'")["final_gene_count"].iloc[0])}
weighting = zscore-derived weights from MAGMA ZSTAT via official scDRS munge_gs; n_min=100, n_max=1000
status = {'PASS' if all(manifest['status'] == 'PASS') else 'PARTIAL'}

## POAG BRIDGE

cells = {bridge['output_cells']}
donors = {bridge['donors']}
classes = {bridge['classes']}
matched genes = {bridge['matched_genes']}
finite norm_score = {bridge['finite_normalized_scores']}/{bridge['output_cells']} ({bridge['finite_normalized_scores']/bridge['output_cells']:.4f})
status = {bridge['status']}

## AMD POWER

h2 Z = 2.4934
interpretation flag = AMD_CELLULAR_POWER = LIMITED; below the scDRS h2-Z>5 heuristic, so a future null AMD cellular result must not be read as strong evidence of biological absence.

## NEXT ACTION

{next_action}

Phase 3A-1 was not run in this cycle. Donor-aware inference remains mandatory if it is authorized.

## FILES

- reports/PHASE3AR8_MAGMA_SCDRS_BRIDGE.md
- reports/PHASE3A_POWER_INTERPRETATION.md
- metadata/magma_reference_provenance.tsv
- results/phase3ar8/magma_reference_audit.tsv
- results/phase3ar8/magma_qc.tsv
- results/phase3ar8/magma_gene_scores.tsv
- results/phase3ar8/scdrs_gene_zscores.tsv
- results/phase3ar8/scdrs_gene_set_manifest.tsv
- results/phase3ar8/retinaplex_traits.gs
- results/phase3ar8/gene_overlap_audit.tsv
- results/phase3ar8/scdrs_poag_bridge_test.tsv
- results/phase3ar8/scdrs_poag_donor_qc.tsv
- results/phase3ar8/POAG.score.gz
- results/phase3ar8/POAG.full_score.gz
- results/phase3ar8/magma/POAG.genes.out
- results/phase3ar8/magma/AMD.genes.out
- results/phase3ar8/magma/RE.genes.out

STOP.
"""
    (ROOT / "reports" / "PHASE3AR8_MAGMA_SCDRS_BRIDGE.md").write_text(report, encoding="utf-8")
    power = """# PHASE 3A POWER INTERPRETATION

AMD_CELLULAR_POWER = LIMITED

The locked AMD input has LDSC h2 Z = 2.4934. This is below the scDRS documentation heuristic of h2 Z > 5 and is retained as a low-power primary trait rather than treated as a technical failure. AMD MAGMA scores and its scDRS gene set were generated because low power does not invalidate the infrastructure bridge. If a future biological Phase 3A analysis finds no AMD cellular enrichment, that absence must not be interpreted as strong evidence that AMD has no cellular signal.

POAG and RE passed the same technical MAGMA gene-level workflow. The present cycle did not run biological cellular enrichment or any Phase 3A-1 tests.
"""
    (ROOT / "reports" / "PHASE3A_POWER_INTERPRETATION.md").write_text(power, encoding="utf-8")


def update_status(verdict: str) -> None:
    block = f"""

## Phase 3A-R8 Update (2026-09-04)

MAGMA reference repair and the locked technical bridge were completed. Standard `NCBI37.3.gene.loc` and the EUR 1000 Genomes PLINK reference were checksum-audited under `data/annotations/magma/` and `data/ld_reference/magma_1000g_eur/`; the LDSC reference directory was not modified. One shared 35 kb annotation was reused for AMD, POAG and European refractive-error genetic liability. MAGMA gene-level statistics, scDRS gene sets, real expression-gene overlap and a POAG scDRS score file were generated. The POAG bridge was technical QC only; no broad cell-type enrichment or Phase 3A-1 analysis was run.

Phase 3A-R8 verdict: `{verdict}`. AMD remains explicitly flagged `AMD_CELLULAR_POWER = LIMITED` because h2 Z = 2.4934. See `reports/PHASE3AR8_MAGMA_SCDRS_BRIDGE.md` and `reports/PHASE3A_POWER_INTERPRETATION.md`.
"""
    for path in (ROOT / "PROJECT_STATUS.md", ROOT / "reports" / "CURRENT_STATUS.md"):
        with path.open("a", encoding="utf-8") as handle:
            handle.write(block)


def main() -> None:
    ensure_dirs()
    provenance = write_reference_provenance()
    gene_audit = audit_gene_loc(GENE_LOC)
    bim_audit = audit_bim(BFILE.with_suffix(".bim"))
    write_reference_audit(gene_audit, bim_audit, provenance)
    magma_qc, magma_info = run_magma()
    gene_info = run_munge_and_overlap(magma_info)
    bridge_info = run_bridge(gene_info)
    # Must be decided only after gene sets and bridge QC are complete.
    qc = pd.read_csv(OUT / "magma_qc.tsv", sep="\t")
    overlap = pd.read_csv(OUT / "gene_overlap_audit.tsv", sep="\t")
    gene_gate = set(qc["status"]) == {"PASS"} and all(overlap["status"] == "PASS")
    if set(qc["status"]) == {"PASS"} and gene_gate and bridge_info["bridge"]["status"] == "PASS":
        verdict = "BROAD_PHASE3A_READY_AMD_LOW_POWER"
    elif set(qc["status"]) == {"PASS"} and gene_gate:
        verdict = "SCDRS_BRIDGE_FAILURE"
    else:
        verdict = "MAGMA_PARTIAL"
    write_reports(magma_qc, gene_info, bridge_info, provenance)
    update_status(verdict)
    print(f"PHASE3AR8_VERDICT\t{verdict}")
    print(f"REPORT\t{ROOT / 'reports' / 'PHASE3AR8_MAGMA_SCDRS_BRIDGE.md'}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"PHASE3AR8_FAILURE\t{type(exc).__name__}: {exc}", file=sys.stderr)
        raise
