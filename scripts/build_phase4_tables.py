#!/usr/bin/env python3
"""Build the three Phase 4 manuscript tables from frozen Phase 0–3B outputs.

This script is a deterministic formatting step only. It does not run a new
biological or statistical analysis.
"""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "tables"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fmt(value: str, digits: int = 4, sign: bool = False) -> str:
    if value in {"", "NA", "nan", "NaN", "NULL", "None"}:
        return "NA"
    try:
        number = float(value)
    except ValueError:
        return value
    if number != 0 and abs(number) < 0.001:
        return f"{number:.3e}"
    if sign:
        return f"{number:+.{digits}f}"
    return f"{number:.{digits}f}"


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_table1() -> None:
    fields = [
        "Trait",
        "GWAS",
        "Ancestry",
        "Effective/sample size",
        "h2",
        "h2 SE",
        "h2 Z",
        "Primary interpretation",
        "Power flag",
    ]
    rows = [
        {
            "Trait": "AMD",
            "GWAS": "GCST003219 / IAMDGC-Fritsche 2016",
            "Ancestry": "European",
            "Effective/sample size": "16,144 cases + 17,832 controls; N_eff=33,892.14",
            "h2": "0.3947",
            "h2 SE": "0.1583",
            "h2 Z": "2.4934",
            "Primary interpretation": "Low-precision conditional primary",
            "Power flag": "Lower power; signed-Z reconstruction from direction and P",
        },
        {
            "Trait": "POAG",
            "GWAS": "GCST90011766 / Gharahkhani 2021 first-stage meta-analysis",
            "Ancestry": "European",
            "Effective/sample size": "16,677 cases + 199,580 controls; N_eff=61,563.71",
            "h2": "0.2220",
            "h2 SE": "0.0168",
            "h2 Z": "13.2143",
            "Primary interpretation": "Primary GWAS",
            "Power flag": "Strong h2 precision; no-UKBB sensitivity available",
        },
        {
            "Trait": "RE",
            "GWAS": "RE_2026_EUR_COMPOSITE / Cheng et al. 2026",
            "Ancestry": "European",
            "Effective/sample size": "Per-variant N; headline N=1,495,159",
            "h2": "0.1095",
            "h2 SE": "0.0048",
            "h2 Z": "22.8125",
            "Primary interpretation": "Primary alternative composite liability GWAS",
            "Power flag": "Strong h2 precision; not pure quantitative MSE",
        },
    ]
    write_tsv(TABLES / "TABLE1_GWAS_ARCHITECTURE.tsv", fields, rows)


def build_table2() -> None:
    source = read_tsv(ROOT / "results/phase3a1/BROAD_CELLULAR_CONSENSUS.tsv")
    fields = [
        "Trait",
        "Cell class",
        "Donors",
        "Donor-aware effect",
        "95% CI",
        "Global FDR",
        "Direction consistency",
        "MAGMA support",
        "Evidence tier",
    ]
    rows: list[dict[str, str]] = []
    for row in source:
        low = fmt(row["CI95_lower"])
        high = fmt(row["CI95_upper"])
        if low == "NA" or high == "NA":
            ci = "NA"
        else:
            ci = f"[{low}, {high}]"
        positive = row["donors_delta_gt_0"]
        negative = row["donors_delta_lt_0"]
        donors = row["n_donors"]
        if positive not in {"", "NA"} and negative not in {"", "NA"}:
            direction = f"{positive}/{donors} positive; {negative}/{donors} negative"
        else:
            direction = row.get("leave_one_donor_stability", "NA") or "NA"
        magma = (
            f"beta={fmt(row['beta'], sign=True)}; P={fmt(row['magma_p'])}; "
            f"FDR={fmt(row['magma_FDR'])}"
        )
        rows.append(
            {
                "Trait": row["trait"],
                "Cell class": row["broad_class"],
                "Donors": donors,
                "Donor-aware effect": fmt(row["mean_donor_delta"], sign=True),
                "95% CI": ci,
                "Global FDR": fmt(row["global_broad_FDR"]),
                "Direction consistency": direction,
                "MAGMA support": magma,
                "Evidence tier": row["tier"],
            }
        )
    write_tsv(TABLES / "TABLE2_BROAD_CELLULAR_ASSOCIATIONS.tsv", fields, rows)


def build_table3() -> None:
    genetic = read_tsv(ROOT / "results/phase1/pairwise_ldsc_rg.tsv")
    profiles = read_tsv(ROOT / "results/phase3a1/cross_trait_broad_profile_similarity.tsv")
    profile_map = {(row["trait_a"], row["trait_b"]): row for row in profiles}
    fields = [
        "Pair",
        "Genetic rg",
        "rg FDR",
        "Cellular-profile rho",
        "Bootstrap CI",
        "RNA convergence",
        "State evidence",
        "ATAC evidence",
        "Interpretation",
    ]
    pair_labels = {
        ("AMD", "POAG"): "AMD–POAG",
        ("AMD", "RE"): "AMD–RE",
        ("POAG", "RE"): "POAG–RE",
    }
    rows: list[dict[str, str]] = []
    for row in genetic:
        a = "AMD" if row["trait_1"].startswith("AMD") else "POAG" if row["trait_1"].startswith("POAG") else "RE"
        b = "AMD" if row["trait_2"].startswith("AMD") else "POAG" if row["trait_2"].startswith("POAG") else "RE"
        profile = profile_map.get((a, b)) or profile_map.get((b, a))
        if profile is None:
            raise RuntimeError(f"Missing broad-profile row for {a}-{b}")
        if {a, b} == {"POAG", "RE"}:
            rna = "Shared Müller glia Tier 1; RE also has Tier 1 amacrine association"
            interpretation = "Cellular-context convergence despite oppositely signed genome-wide architecture"
        elif {a, b} == {"AMD", "RE"}:
            rna = "Shared Müller glia Tier 1; RE also has Tier 1 amacrine association"
            interpretation = "Shared Müller context with weaker cross-trait profile similarity and RE-specific amacrine component"
        else:
            rna = "Shared Müller glia Tier 1"
            interpretation = "Partial genetic relationship with shared broad Müller context"
        rows.append(
            {
                "Pair": pair_labels[(a, b)] if (a, b) in pair_labels else pair_labels[(b, a)],
                "Genetic rg": fmt(row["rg"], sign=True),
                "rg FDR": fmt(row["FDR"]),
                "Cellular-profile rho": fmt(profile["spearman_rho"], sign=True),
                "Bootstrap CI": f"[{fmt(profile['bootstrap_CI95_lower'])}, {fmt(profile['bootstrap_CI95_upper'])}]",
                "RNA convergence": rna,
                "State evidence": "No donor-replicated Müller state; no consensus state",
                "ATAC evidence": "No FDR-supported class-level OCR enrichment",
                "Interpretation": interpretation,
            }
        )
    write_tsv(TABLES / "TABLE3_CROSS_TRAIT_CONVERGENCE.tsv", fields, rows)


if __name__ == "__main__":
    TABLES.mkdir(parents=True, exist_ok=True)
    build_table1()
    build_table2()
    build_table3()
    print("Wrote TABLE1_GWAS_ARCHITECTURE.tsv, TABLE2_BROAD_CELLULAR_ASSOCIATIONS.tsv, and TABLE3_CROSS_TRAIT_CONVERGENCE.tsv")
