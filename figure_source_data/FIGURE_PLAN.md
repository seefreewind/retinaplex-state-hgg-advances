# Final figure plan

## Figure 1 — Study architecture and genetic landscape

**Panel A: study design.** Flow from three European GWAS through global genetic architecture, HDL-L local resolution, donor-aware retinal RNA mapping, Müller state testing, and HRCA OCR assessment. Mark the final interpretive boundary as RNA-defined broad cellular context.

**Panel B: pairwise genome-wide genetic correlations.** Show signed `rg` with confidence intervals for AMD–POAG, AMD–RE and POAG–RE. Annotate FDR and the POAG no-UKBB sensitivity.

**Panel C: global versus local architecture.** Summarize the absence of FDR-significant HDL-L local hotspots across all tested pairs. Use a compact count plot or schematic; keep the exact regional audit in the supplement.

**Main message.** The traits have partial and oppositely signed relationships, but no robust local hotspot dominates the architecture.

**Primary sources.** `results/phase1/pairwise_ldsc_rg.tsv`, `results/phase1/poag_noukbb_sensitivity.tsv`, `results/phase2a/hdll_local_rg_all.tsv`, `tables/TABLE1_GWAS_ARCHITECTURE.tsv`.

## Figure 2 — Broad retinal cellular architecture

**Panel A: donor-aware heatmap.** Traits on rows and the nine locked broad retinal classes on columns. Use the donor-aware mean effect as color and global FDR as an overlaid symbol or label.

**Panel B: Müller donor distributions.** Show donor-level effects for AMD, POAG and RE with the mean and 95% interval. Annotate 8/8 positive donors.

**Panel C: RE amacrine component.** Show the RE amacrine donor-level distribution beside the Müller result. Add MAGMA beta/FDR as a compact orthogonal-support annotation.

**Main message.** Müller glia are the shared broad RNA context, while RE has an additional amacrine component.

**Primary sources.** `results/phase3a1/BROAD_CELLULAR_CONSENSUS.tsv`, `results/phase3a1/magma_broad_celltype_enrichment.tsv`, `tables/TABLE2_BROAD_CELLULAR_ASSOCIATIONS.tsv`.

## Figure 3 — Genetic direction versus cellular context

**Panel A: pair-level comparison.** Plot signed genome-wide `rg` against standardized broad-cellular profile similarity for the three pairs. Label the POAG–RE point prominently.

**Panel B: POAG–RE profile.** Display the nine standardized broad-cell effects for POAG and RE with a connecting line or paired dot plot. Include `rho=0.8833` and bootstrap 95% CI `[0.5714, 0.9524]`.

**Panel C: conceptual distinction.** A two-axis schematic separates signed genetic-effect direction from similarity of cellular-context enrichment.

**Main message.** Negative genome-wide genetic correlation can coexist with positive similarity of broad cellular profiles; this does not imply shared risk alleles or same-direction effects.

**Primary sources.** `results/phase1/pairwise_ldsc_rg.tsv`, `results/phase3a1/cross_trait_broad_profile_similarity.tsv`, `results/phase3a1/cross_trait_standardized_effect_profiles.tsv`, `tables/TABLE3_CROSS_TRAIT_CONVERGENCE.tsv`.

## Figure 4 — Robustness and orthogonal evidence

**Panel A: mapping-selection sensitivity.** Display the principal Müller effects under the three locked mapping-selection rules.

**Panel B: leave-one-donor-out.** Display the range of Müller effects after removing each donor, grouped by trait.

**Panel C: scDRS versus MAGMA.** Use a compact comparison of broad-cell effect direction and FDR; annotate that MAGMA is orthogonal associative support and not a causal validation.

**Main message.** The broad RNA pattern is directionally stable across prespecified mapping and donor-sensitivity checks, with the expected lower precision for AMD.

**Primary sources.** `results/phase3a1/mapping_selection_sensitivity.tsv`, `results/phase3a1/leave_one_donor_out.tsv`, `results/phase3a1/BROAD_CELLULAR_CONSENSUS.tsv`, `results/phase3a1/magma_broad_celltype_enrichment.tsv`.

## Figure 5 — Boundary of cellular convergence

**Panel A: broad RNA tier.** Show the three Müller Tier 1 broad-cell results and the RE amacrine Tier 1 result.

**Panel B: Müller state resolution.** Show the state eligibility outcome and NMF program summary as unsupported, with donor/state criteria explicitly labelled.

**Panel C: regulatory resolution.** Show class-level OCR s-LDSC FDR for Müller and the other major classes. Mark all traits as no FDR-supported class-level enrichment. Include a note that the official OCR resource is pooled.

**Panel D: evidence ladder.** A compact ladder from broad RNA `SUPPORTED` to Müller state `UNSUPPORTED` to class-level ATAC `UNSUPPORTED`; use “RNA-defined broad cellular convergence” as the endpoint.

**Main message.** The evidence supports a broad transcriptional context, not a discrete Müller state or multiomic regulatory convergence.

**Primary sources.** `results/phase3a2/MULLER_STATE_PROGRAM_CONSENSUS.tsv`, `results/phase3a2/muller_state_trait_effects.tsv`, `results/phase3a2/muller_gene_programs.tsv`, `results/phase3b/sldsc_retinal_ocr_enrichment.tsv`, `results/phase3b/RNA_ATAC_CELLULAR_CONSENSUS.tsv`.

## Supplementary figures

Move technical detail to supplementary figures: GWAS QC and LDSC diagnostics; GenomicSEM matrices and saturation note; full HDL-L regional results; exact-ID mapping recovery and donor coverage; all broad-class donor effects; mapping and leave-one-donor-out details; Müller clustering grids and NMF diagnostics; all class-level OCR s-LDSC results; HRCA donor/sample and genome-build audits; CRE–gene resource audit; and descriptive RNA–ATAC profile concordance.

## Assembly rules

- Do not delete preliminary P23–P27 figures; retain them as audit/provenance outputs.
- Use self-contained legends with sample/donor counts, effect definition, FDR family, and abbreviation definitions.
- Use a consistent sign convention across Figures 1–4.
- Never label Figure 5 as “multiomic validation,” “open-chromatin validation,” or “Müller state discovery.”
- Final artwork has not been generated in Phase 4; this file is the authoritative construction plan.
