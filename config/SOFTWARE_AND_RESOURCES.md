# Software and resources

## Analysis software

| Component | Version / locked setting | Role |
|---|---|---|
| LDSC | v3.0.1; project Python 3.9 environment | SNP heritability and genome-wide genetic correlation |
| HDL-L | v1.4.3; R 4.4.3 | Local genetic covariance and correlation |
| MAGMA | v1.10 (custom, self-reported) | Gene-property enrichment |
| scDRS | v1.0.2 | Single-cell polygenic gene-score mapping |
| GenomicSEM | locked project environment | Descriptive three-indicator factor model |
| Scanpy | software citation recorded in reference list | Normalization, PCA and graph preparation |
| Leiden | locked primary resolution 0.4 | Graph-based state discovery |
| s-LDSC | baselineLD v2.2 | HRCA pooled OCR annotation assessment |

## Key datasets and resources

- AMD GWAS: GCST003219; source study DOI 10.1038/ng.3448.
- POAG GWAS: GCST90011766; source study DOI 10.1038/s41467-020-20851-4.
- RE: European composite refractive-error resource RE_2026_EUR_COMPOSITE; source study DOI 10.1038/s41588-026-02576-0.
- Retinal RNA: GSE265774, exact-ID-mapped healthy reference, 98,834 cells, eight donors and eight samples.
- HRCA pooled RNA/ATAC resource: source publication DOI 10.1038/s41588-025-02454-1; associated accessions include GSE265774, GSE265801, PRJNA1093457 and SRP499089.
- Reference genome builds: RNA/ATAC source GRCh38/hg38; GWAS regulatory analysis hg19 after locked coordinate conversion.
- EUR LD reference: 22,665,064 reference SNPs for MAGMA and 335,272 individuals for HDL-L as specified in the frozen analysis configuration.

## Reproducibility notes

Primary FDR families, thresholds, random seed (20260905 for donor bootstrap), cell and donor eligibility rules, gene-set overlaps and coordinate-conversion settings are preserved in the analysis-plan manifest. Derived source-data tables for the main figures are the three doc-ready TSVs in the manuscript package.
