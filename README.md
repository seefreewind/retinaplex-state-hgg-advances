# RetinaPlex-State HGG Advances release package

This directory is a distributable release candidate for the manuscript “Cross-trait genetic architecture reveals shared retinal cellular contexts despite opposing genetic direction.” It contains the analysis configuration, software/resource manifest, provenance notes, derived result tables, figure source data and selected scripts needed to inspect or reassemble the reported outputs.

## Scope and data-use boundary

The package contains derived tables and metadata only. Raw or restricted GWAS summary-statistics files, individual-level genotype files, H5AD objects, PLINK files, OCR matrices and other large source datasets are intentionally excluded. Users should retrieve public inputs from the accessions and source publications listed in `metadata/` and should follow each resource’s access terms.

The reported manuscript values are read from the retained result tables. The package does not claim that every included script can be rerun without the original public inputs, local resource downloads or author-supplied paths.

## Contents

- `config/`: analysis configuration, thresholds, trait definitions and software/resource manifest.
- `metadata/`: source accessions, atlas provenance, mapping records and resource manifests.
- `provenance/`: LD-reference and MAGMA-reference provenance records.
- `results/main/`: main analysis tables used for the manuscript.
- `results/supplementary/`: supplementary result tables and machine-readable summaries.
- `figure_source_data/`: source tables for Figures 1–5 and Tables 1–3.
- `scripts/`: selected analysis, table-assembly and display-only figure-assembly scripts.
- `environment/`: available environment specifications and requirements files.
- `LICENSE_PLACEHOLDER.txt`: license decision required before public release.

## Reproduction and inspection

1. Review `config/FINAL_ANALYSIS_FREEZE_v1.yaml`, `config/SOFTWARE_AND_RESOURCES.md` and the provenance tables.
2. Create an environment from the file appropriate to the host system in `environment/`.
3. Supply local paths to the public inputs and references described in `metadata/`; do not commit those inputs unless their terms permit redistribution.
4. For portable script execution, set `SCDRS_PYTHON`, `SC_PYTHON`, `SCDRS_EXECUTABLE` and `MAGMA_EXECUTABLE` as needed. The scripts otherwise fall back to commands on `PATH`.
5. Use `scripts/build_phase5_figures.py` only for deterministic display assembly from retained outputs. It does not constitute a new biological analysis.

## Release status

This release candidate is publicly available at https://github.com/seefreewind/retinaplex-state-hgg-advances. The archived DOI/accession, final license and exact input-file manifest must still be confirmed by the authors before the manuscript is submitted. No DOI, accession or license is invented in this package.
