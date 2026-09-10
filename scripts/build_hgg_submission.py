from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
SUB = ROOT / "submission"
MAN = ROOT / "manuscript"


def read_body(path: Path) -> str:
    text = path.read_text(encoding="utf-8").strip()
    lines = text.splitlines()
    if lines and lines[0].startswith("# "):
        text = "\n".join(lines[1:]).strip()
    return text


OLD_TO_NEW = {
    1: 1,   # Fritsche AMD
    2: 2,   # Gharahkhani POAG
    3: 3,   # Cheng RE
    4: 7,   # LDSC method
    5: 5,   # genetic-correlation atlas
    6: 16,  # HDL
    7: 17,  # enhanced local correlation
    8: 8,   # Lukowski retina atlas
    9: 10,  # Wang retina multiome
    10: 11, # HRCA
    11: 13, # MAGMA
    12: 12, # scDRS
    13: 19, # GenomicSEM
    14: 20, # Scanpy
    15: 21, # Leiden
    16: 22, # NMF
    17: 23, # UCSC/liftOver
    18: 14, # Müller biology
    19: 24, # s-LDSC
    20: 15, # Benjamini-Hochberg
    21: 18, # 1000 Genomes
}


def parse_citation_token(token: str):
    token = token.strip()
    if token.isdigit():
        return [int(token)]
    if re.fullmatch(r"\d+[-–]\d+", token):
        a, b = re.split(r"[-–]", token)
        return list(range(int(a), int(b) + 1))
    return []


def remap_citations(text: str) -> str:
    pattern = re.compile(r"\[([0-9][0-9,;\s–-]*)\]")

    def replace(match):
        body = match.group(1)
        old_numbers = []
        for token in re.split(r"[,;\s]+", body.strip()):
            old_numbers.extend(parse_citation_token(token))
        if not old_numbers:
            return match.group(0)
        new_numbers = sorted({OLD_TO_NEW[n] for n in old_numbers if n in OLD_TO_NEW})
        if not new_numbers:
            return match.group(0)
        return "[" + ", ".join(str(n) for n in new_numbers) + "]"

    return pattern.sub(replace, text)


def tsv_to_markdown(path: Path) -> str:
    rows = []
    lines = path.read_text(encoding="utf-8").strip().splitlines()
    for line in lines:
        rows.append(line.split("\t"))
    header = rows[0]
    body = rows[1:]
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    out.extend("| " + " | ".join(row) + " |" for row in body)
    return "\n".join(out)


def clean_editor_language(text: str) -> str:
    replacements = [
        ("analysis freeze", "prespecified analysis plan"),
        ("locked preparation rule", "prespecified preparation rule"),
        ("locked reference convention", "prespecified reference convention"),
        ("locked LDSC v3.0.1 environment", "LDSC v3.0.1"),
        ("locked high-definition likelihood (HDL-L) workflow", "high-definition likelihood (HDL-L) analysis"),
        ("under the locked rules", "under the prespecified rules"),
        ("nine locked broad retinal classes", "nine prespecified broad retinal classes"),
        ("locked covariate specification", "prespecified covariate specification"),
        ("nine locked classes", "nine prespecified classes"),
        ("locked reconstruction and stability rule", "prespecified reconstruction and stability rule"),
        ("locked hg38-to-hg19 chain", "specified hg38-to-hg19 chain"),
        ("locked 1-cM resolution", "1-cM resolution"),
        ("under the locked design", "under the prespecified design"),
        ("from the locked analysis", "from the prespecified analysis"),
        ("the locked donor-direction rule", "the prespecified donor-direction criterion"),
        ("the locked broad-class result", "the prespecified broad-class result"),
        ("locked mapping-selection rules", "prespecified mapping-selection rules"),
        ("locked replication and composition criteria", "prespecified replication and composition criteria"),
        ("passed the locked donor, FDR and direction-consistency criteria", "passed the donor, FDR and direction-consistency criteria"),
        ("the remaining estimates were retained in the audit as non-interpretable", "the remaining estimates were retained as non-interpretable"),
        ("no region met the promotion criteria for a robust local hotspot", "no region met the prespecified threshold for retaining a robust local hotspot"),
        ("no region met the prespecified criteria for promotion", "no region met the prespecified threshold for retaining a robust hotspot"),
        ("the broad-class conclusions were based on the frozen method-specific outputs", "the broad-class conclusions were based on the method-specific outputs"),
        ("frozen LDSC estimates", "LDSC estimates from the primary analysis"),
        ("frozen s.e.", "corresponding s.e."),
        ("a frozen signed genome-wide rᵍ", "a signed genome-wide rᵍ"),
        ("frozen HRCA pooled-annotation FDR values", "HRCA pooled-annotation FDR values"),
        ("The supplementary package should preserve the frozen audit trail", "The supplementary package should preserve the complete provenance record"),
        ("retained in the audit", "retained in the supplementary record"),
        ("auditable CRE–gene target framework", "traceable CRE–gene target framework"),
        ("audited and included where estimable", "checked and included where estimable"),
    ]
    for old, new in replacements:
        text = text.replace(old, new)
    text = text.replace("authoritative", "corresponding")
    text = text.replace("frozen", "prespecified")
    text = text.replace("freeze", "analysis plan")
    text = text.replace("locked", "prespecified")
    text = text.replace("promotion criteria", "retention threshold")
    text = text.replace("project-level", "study-specific")
    text = text.replace("auditable", "traceable")
    text = text.replace("audited", "checked")
    text = text.replace("audit", "quality-control record")
    text = text.replace("workflow", "analysis")
    return text


def optimize_discussion(text: str) -> str:
    text = text.replace(
        "This study was designed to keep two questions separate: how polygenic effects relate across traits and where trait-associated gene scores are concentrated in a human retinal reference. Across AMD, POAG and European composite refractive-error genetic liability, the strongest evidence was a donor-consistent broad Müller-glial RNA context. The principal conceptual result was that POAG and RE had a negative genome-wide genetic correlation but a highly similar retinal cellular-context profile. RE also showed an amacrine-cell component. Local genetic analysis, Müller-state analysis and pooled OCR analysis did not support a more specific resolution.",
        "This study separates two dimensions of polygenic architecture: the signed relationship between traits and the retinal cellular contexts in which trait-associated genes are concentrated. Across AMD, POAG and European composite refractive-error genetic liability, the strongest evidence was a donor-consistent broad Müller-glial RNA context. POAG and RE combined a negative genome-wide genetic correlation with a highly similar retinal cellular-context profile, while RE also showed an amacrine-cell component. Local, Müller-state and pooled OCR analyses defined the resolution of the supported signal."
    )
    text = text.replace(
        "This distinction changes how cross-trait pleiotropy should be described. A positive cellular-context comparison does not establish shared risk alleles, concordant allelic effects or a common signed genetic factor. Conversely, a negative genetic correlation does not require the traits to occupy different cellular contexts. The POAG–RE result would have been incompletely described by either statistic alone: r_g captures direction, whereas the cellular-context profile captures context. Reporting both makes the complementary interpretation explicit.",
        "Reporting both statistics provides a fuller description of cross-trait architecture: r_g captures genome-wide direction, whereas cellular-context similarity captures the relative distribution of trait-associated gene scores across cell classes. A positive cellular-context comparison identifies localization similarity; it does not establish shared risk alleles or concordant allelic effects. The POAG–RE result demonstrates why the two summaries should be interpreted together."
    )
    text = text.replace(
        "The same logic can be applied outside ocular genetics, but the present data do not establish that the POAG–RE pattern will generalize to all complex traits. These results illustrate how signed genetic relationships and cellular localization can provide complementary views of cross-trait architecture when the relevant reference and donor structure are available.",
        "The distinction is relevant beyond ocular genetics. Its use in other complex traits will depend on phenotype definitions, ancestry scope, the cellular reference and the available donor structure."
    )
    text = text.replace(
        "These results place Müller glia within a shared broad transcriptional context, while they do not identify a disease-specific Müller process.",
        "Together, these results place Müller glia within a shared broad transcriptional context in the tested healthy reference."
    )
    text = text.replace(
        "The present analysis does not show that Müller glia initiate disease, mediate the genetic effects or represent the only relevant retinal compartment. It shows that trait-associated gene scores are enriched in genes expressed by Müller glia in the tested healthy reference.",
        "The analysis localizes trait-associated gene scores to Müller glia in this reference; it does not identify disease initiation, mediation or exclusivity."
    )
    text = text.replace(
        "This negative result indicates that the current healthy reference and replication requirements do not support selecting a Müller substate. It does not establish that Müller heterogeneity is irrelevant; disease, age, perturbation or developmental references may contain structure that is absent from a healthy static reference.",
        "This result sets the current resolution boundary: the healthy reference and replication requirements do not support selecting a Müller substate. Disease, age, perturbation or developmental references may contain additional structure."
    )
    text = text.replace("The layers should therefore not be collapsed into a single regulatory claim.", "These layers represent distinct evidence resolutions and cannot be combined into a single regulatory claim.")
    text = text.replace(
        "The absence of state and OCR support is scientifically useful because it identifies where the evidence stops. The manuscript supports RNA-defined broad cellular convergence. It does not support a discrete disease state, a continuous Müller program, an open-chromatin mechanism or a linked causal regulatory element. Future work could test these possibilities with disease or perturbation references, donor-resolved chromatin profiles and a complete auditable CRE–gene target framework.",
        "The state and OCR results set the current resolution boundary at broad RNA-defined cellular context. Future work could test finer resolution with disease or perturbation references, donor-resolved chromatin profiles and a traceable CRE–gene target framework."
    )
    limitations_pattern = re.compile(r"Several limitations remain\..*?These constraints limit mechanistic, causal and clinical interpretation, but they do not remove the supported broad RNA association\.", re.S)
    limitations = (
        "The limitations fall into four themes. **Genetic inputs and phenotype definition.** AMD had the least precise heritability estimate and limited power for some cross-trait and gene-level comparisons. The RE phenotype is a composite European liability measure, the three-indicator GenomicSEM model was saturated and descriptive, and local estimates depend on the selected regional partition, LD reference, eigenvalue threshold and estimability rules. **RNA reference limitations.** The reference comprised healthy donors, one sample per RNA donor created donor/sample confounding, and exact-ID mapping retained 86.49% of source cells, which may have introduced selection related to identifier recovery. **Associative gene-to-cell inference.** scDRS and MAGMA localize gene scores or gene properties to cell classes but do not identify a causal cell of origin; the healthy reference and gene-level inputs also limited AMD power relative to POAG and RE. **State and regulatory resolution.** No donor-replicated Müller state or expression program was observed, and the pooled ATAC analysis had limited donor-level regulatory sensitivity. These themes define the scope of the broad RNA association and guide the next validation steps."
    )
    return limitations_pattern.sub(limitations, text)


def main():
    methods = remap_citations(read_body(MAN / "METHODS_v1.md"))
    methods = methods.replace(
        "This study is an integrative analysis of public GWAS summary statistics and public human retinal single-cell or single-nucleus resources. The analysis was organized around three estimands: signed genome-wide genetic relationships, broad retinal cellular enrichment, and the resolution at which a cellular signal remains reproducible. All primary numerical results reported here were taken from the frozen analysis package and its authoritative result tables. No new biological or statistical analysis was introduced during manuscript assembly.",
        "This study is an integrative analysis of public GWAS summary statistics and public human retinal single-cell or single-nucleus resources. The analysis was organized around three estimands: signed genome-wide genetic relationships, broad retinal cellular enrichment, and the resolution at which a cellular signal remains reproducible. The reported estimates were taken from the completed analysis outputs, and no additional biological or statistical analyses were performed for this manuscript."
    )
    methods = methods.replace(
        "The analysis freeze specified the input accessions, software versions, genome builds, gene annotations, eligibility thresholds, multiplicity families and random seeds before manuscript assembly. It also specified the claim hierarchy used to distinguish Tier 1 broad-class evidence from secondary or unsupported resolution. Authoritative result tables were treated as immutable source data for the present draft; display rounding was applied only after values had been traced to those tables.",
        "The analysis plan prespecified the input accessions, software versions, genome builds, gene annotations, eligibility thresholds, multiplicity families and random seeds."
    )
    methods = methods.replace(
        "MAGMA was run with the locked v1.10-custom executable and the NCBI37.3 gene-location file.",
        "MAGMA v1.10 was used for gene-property analysis."
    )
    methods = methods.replace(
        "Software versions, source accessions, coordinate builds, thresholds and unresolved submission metadata are recorded in the accompanying audit files.",
        "Software versions, source accessions, coordinate builds and thresholds are reported in the Supplementary Information."
    )
    methods = methods.replace(
        "The prepared inputs underwent the project-level access, allele, sample-size and summary-statistic checks recorded in the analysis freeze.",
        "The prepared inputs underwent prespecified access, allele, sample-size and summary-statistic checks."
    )
    methods = methods.replace("the frozen local eligibility rules", "the prespecified local eligibility rules")
    methods = methods.replace("retained in the audit but were not treated", "retained in the supplementary output but were not treated")
    methods = clean_editor_language(methods)
    methods = re.sub(
        r"The reported estimates were taken from the completed analysis outputs, and no additional biological or statistical analyses were performed for this manuscript\.\s*",
        "",
        methods,
    )
    methods = methods.replace("using the LDSC v3.0.1.", "using LDSC v3.0.1.")
    methods = methods.replace("with an effective sample size of 33,892.14", "with an effective sample size of approximately 33,892")
    methods = methods.replace("with an effective sample size of 61,563.71", "with an effective sample size of approximately 61,564")
    ancestry_note = (" The analysis used European-ancestry summary-statistics inputs for all three traits. "
                     "The POAG and RE source publications describe multi-ancestry studies, whereas the POAG file and RE composite used here were the European-ancestry inputs.")
    methods = methods.replace("The source study was the multi-ancestry refractive-error GWAS by Cheng and colleagues. [3]", "The source study was the multi-ancestry refractive-error GWAS by Cheng and colleagues. [3]" + ancestry_note)
    intro = clean_editor_language(read_body(SUB / 'HGG_INTRODUCTION_v1.md'))
    results = clean_editor_language(read_body(SUB / 'HGG_RESULTS_v1.md'))
    discussion = clean_editor_language(optimize_discussion(read_body(SUB / 'HGG_DISCUSSION_v1.md')))
    figures = clean_editor_language(read_body(MAN / "FIGURE_LEGENDS_v1.md"))
    figures = figures.replace(
        "Counts summarize the 107 interpretable AMD–POAG regional estimates and the absence of FDR-supported local hotspots across the tested pairwise analyses.",
        "Counts show interpretable regional estimates for AMD–POAG (107), AMD–RE (8) and POAG–RE (6); FDR-supported local hotspots = 0 for all three pairs."
    )
    figures = figures.replace(
        "The POAG–RE point is highlighted.",
        "The POAG–RE point is highlighted; bootstrap 95% CIs are shown for all three profile similarities: AMD–POAG 0.5476–0.8167, AMD–RE 0.1905–0.6833 and POAG–RE 0.5714–0.9524."
    )
    figures = figures.replace(
        "The displayed ρ = 0.8833 has a donor-bootstrap 95% CI of 0.5714–0.9524 from 5,000 replicates.",
        "The displayed rᵍ = −0.1732 and ρ = +0.8833 have a donor-bootstrap 95% CI for ρ of 0.5714–0.9524 from 5,000 replicates."
    )
    figures = figures.replace(
        "Signed allelic architecture and profile-level cellular context are separate dimensions; similarity in cellular context does not encode shared alleles or same-direction effects.",
        "Signed allelic direction ≠ cellular-context localization."
    )
    figures = figures.replace("The inferential unit is the donor, with a minimum", "The inferential unit is the donor (n = 8 donors for the focal Müller-glial panels), with a minimum")
    figures = figures.replace("The two methods have different estimands and are shown as orthogonal associative evidence, not as causal evidence.", "The two methods have different estimands and are shown as orthogonal associative evidence.")
    figures = figures.replace("No state was eligible because the prespecified replication and composition criteria were not met.", "No state was eligible under the prespecified replication and composition criteria.")
    tables = read_body(MAN / "TABLE_LEGENDS_v1.md")
    supplements = clean_editor_language(read_body(MAN / "SUPPLEMENTARY_INDEX_v1.md"))
    supplements += "\n\nSupplementary Table 8 records the MAGMA executable banner as `v1.10 (custom)`, together with the executable checksum and the available provenance limits."
    table1 = tsv_to_markdown(MAN / "Table1.doc-ready.tsv")
    table2 = tsv_to_markdown(MAN / "Table2.doc-ready.tsv")
    table3 = tsv_to_markdown(MAN / "Table3.doc-ready.tsv")
    table1 = clean_editor_language(table1).replace("N_eff=33,892.14", "N_eff≈33,892").replace("N_eff=61,563.71", "N_eff≈61,564")
    table2 = clean_editor_language(table2)
    table3 = clean_editor_language(table3).replace("Cellular-profile rho", "Cellular-context similarity ρ")
    abstract = clean_editor_language(read_body(SUB / 'HGG_ABSTRACT_v1.md'))
    abstract = abstract.replace("10^-4", "10⁻⁴").replace("10^-5", "10⁻⁵").replace("10^-6", "10⁻⁶").replace("10^-7", "10⁻⁷").replace("10^-9", "10⁻⁹")
    author_contributions = (
        "Da Lin: Conceptualization, methodology, software, formal analysis, visualization and writing – original draft. "
        "Ying Chen: Data curation, investigation, validation and writing – review and editing. "
        "Yue Liu: Data curation, investigation, validation and writing – review and editing. "
        "Yu Zhang: Conceptualization, supervision, project administration, resources and writing – review and editing."
    )

    references = """1. Fritsche LG et al. A large genome-wide association study of age-related macular degeneration highlights contributions of rare and common variants. Nature Genetics. 2016;48:134-143. doi:10.1038/ng.3448. PMID:26691988.
2. Gharahkhani P et al. Genome-wide meta-analysis identifies 127 open-angle glaucoma loci with consistent effect across ancestries. Nature Communications. 2021;12:1258. doi:10.1038/s41467-020-20851-4. PMID:33627673.
3. Cheng F-F et al. Multi-ancestry genome-wide association analyses of refractive error augment genetic discovery and polygenic prediction. Nature Genetics. 2026;58:1030-1039. doi:10.1038/s41588-026-02576-0.
4. Visscher PM, Yang J. A plethora of pleiotropy across complex traits. Nature Genetics. 2016;48:707-708. doi:10.1038/ng.3604.
5. Bulik-Sullivan B et al. An atlas of genetic correlations across human diseases and traits. Nature Genetics. 2015;47:1236-1241. doi:10.1038/ng.3406.
6. O'Connor LJ, Price AL. Distinguishing genetic correlation from causation across 52 diseases and complex traits. Nature Genetics. 2018;50:1728-1734. doi:10.1038/s41588-018-0255-0.
7. Bulik-Sullivan B et al. LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. Nature Genetics. 2015;47:291-295. doi:10.1038/ng.3211.
8. Lukowski SW et al. A single-cell transcriptome atlas of the adult human retina. The EMBO Journal. 2019;38:e100811. doi:10.15252/embj.2018100811. PMID:31436334.
9. Menon M et al. Single-cell transcriptomic atlas of the human retina identifies cell types associated with age-related macular degeneration. Nature Communications. 2019;10:4902. doi:10.1038/s41467-019-12780-8. PMID:31653841.
10. Wang SK et al. Single-cell multiome of the human retina and deep learning nominate causal variants in complex eye diseases. Cell Genomics. 2022;2:100164. doi:10.1016/j.xgen.2022.100164. PMID:36277849.
11. Li J et al. Single-cell atlas of the transcriptome and chromatin accessibility in the human retina. Nature Genetics. 2026;58:418-433. doi:10.1038/s41588-025-02454-1. PMID:41578023.
12. Zhang MJ et al. Polygenic enrichment distinguishes disease associations of individual cells in single-cell RNA-seq data. Nature Genetics. 2022;54:1572-1580. doi:10.1038/s41588-022-01167-z. PMID:36216977.
13. de Leeuw CA, Mooij JM, Heskes T, Posthuma D. MAGMA: generalized gene-set analysis of GWAS data. PLoS Computational Biology. 2015;11:e1004219. doi:10.1371/journal.pcbi.1004219. PMID:25906354.
14. Bringmann A et al. Müller cells in the healthy and diseased retina. Progress in Retinal and Eye Research. 2006;25(4):397-424. doi:10.1016/j.preteyeres.2006.05.003. PMID:16839797.
15. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B. 1995;57:289-300. doi:10.1111/j.2517-6161.1995.tb02031.x.
16. Ning Z, Pawitan Y, Shen X. High-definition likelihood inference of genetic correlations across human complex traits. Nature Genetics. 2020;52:859-864. doi:10.1038/s41588-020-0653-y.
17. Li Y, Pawitan Y, Shen X. An enhanced framework for local genetic correlation analysis. Nature Genetics. 2025;57:1053-1058. doi:10.1038/s41588-025-02123-3.
18. Auton A et al. A global reference for human genetic variation. Nature. 2015;526:68-74. doi:10.1038/nature15393. PMID:26432245.
19. Grotzinger AD et al. Genomic structural equation modelling provides insights into the multivariate genetic architecture of complex traits. Nature Human Behaviour. 2019;3:513-525. doi:10.1038/s41562-019-0566-x. PMID:31068683.
20. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biology. 2018;19(1):15. doi:10.1186/s13059-017-1382-0. PMID:29409532.
21. Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. Scientific Reports. 2019;9(1):5233. doi:10.1038/s41598-019-41695-z. PMID:30914743.
22. Lee DD, Seung HS. Learning the parts of objects by non-negative matrix factorization. Nature. 1999;401:788-791. doi:10.1038/44565. PMID:10553949.
23. Kent WJ et al. The human genome browser at UCSC. Genome Research. 2002;12:996-1006. doi:10.1101/gr.229102. PMID:12045153.
24. Finucane HK et al. Partitioning heritability by functional annotation using genome-wide association summary statistics. Nature Genetics. 2015;47:1228-1235. doi:10.1038/ng.3404. PMID:26414678."""

    manuscript = f"""# Cross-trait genetic architecture reveals shared retinal cellular contexts despite opposing genetic direction

**Article**

**Authors:** Da Lin¹, Ying Chen², Yue Liu², Yu Zhang¹

**Affiliations:**
¹ Department of Ophthalmology, The Second Affiliated Hospital of Wenzhou Medical University, No. 109 Xueyuan West Road, Lucheng District, Wenzhou, Zhejiang Province, China
² Wenzhou Medical University, Wenzhou, Zhejiang Province, China

**Correspondence:** Yu Zhang, Department of Ophthalmology, The Second Affiliated Hospital of Wenzhou Medical University, No. 109 Xueyuan West Road, Lucheng District, Wenzhou, Zhejiang Province, China; zhangyu1@wzhealth.com; ORCID: 0000-0001-8579-3692

**Short title:** Ocular genetic direction and retinal cellular context

## Abstract

{abstract}

## Keywords

Human genetics; genetic correlation; polygenic liability; single-cell genomics; retinal cellular context; Müller glia

## Introduction

{intro}

## Materials and methods

{methods}

## Results

{results}

## Discussion

{discussion}

## Data and code availability

The GWAS summary statistics and retinal transcriptomic and chromatin-accessibility resources used in this study are publicly available through the accessions and source publications listed in Table 1 and the Web resources section. Derived analysis tables and figure source data will be deposited in a versioned public repository before submission; the permanent accession and access conditions will be added to this statement when the deposit is complete.

The analysis code and figure-assembly code will be deposited in a versioned public repository with an archived release before submission; the public URL, release tag, archival DOI and license will be added when the release is complete.

## Web resources

- GWAS Catalog: GCST003219 and GCST90011766.
- GEO: GSE265774 and GSE265801; the source study also reports GSE281526.
- NCBI BioProject: PRJNA1104225; HCA Data Portal eye network.
- CELLxGENE HRCA collection; UCSC Retina Cell Browser; Single Cell Portal accessions SCP2805-SCP2808.
- RE_2026_EUR_COMPOSITE source publication: doi:10.1038/s41588-026-02576-0.
- Project-derived tables and code: versioned public repository release to be added before submission.

## Supplemental Information

Supplemental Information should contain the complete GWAS access and QC audit, pairwise LDSC and no-UK Biobank sensitivity outputs, HDL-L regional results, donor coverage and cell-balance checks, all donor-level scDRS results, Müller-state and NMF outputs, pooled HRCA OCR results, RNA-ATAC concordance, and the software/resource manifest. The detailed index is provided below.

{supplements}

## Acknowledgments

None.

## Author contributions

{author_contributions}

## Funding

This work received no dedicated funding.

## Ethics approval and consent to participate

This secondary analysis used publicly available GWAS summary statistics and public, de-identified human retinal single-cell and single-nucleus resources. No new human participants or specimens were enrolled, and no new specimens were collected for this study. Ethics approval and informed-consent procedures were conducted by the original contributing studies as described in their source publications and data records. No new institutional approval is claimed for the present computational analysis.

## Consent for publication

Not applicable to this analysis of public, de-identified resources.

## Declaration of interests

The authors declare no competing interests.

## References

{references}

## Figure titles and legends

{figures}

## Table titles and legends

{tables}

## Tables

### Table 1 | GWAS and cross-trait genetic architecture

{table1}

### Table 2 | Donor-aware broad retinal cellular associations

{table2}

### Table 3 | Cross-trait convergence and resolution boundary

{table3}
"""

    manuscript = clean_editor_language(manuscript)
    manuscript = manuscript.replace("cellular-profile similarity", "cellular-context similarity")
    manuscript = manuscript.replace("cellular-profile comparison", "cellular-context comparison")
    manuscript = manuscript.replace("cellular profiles", "cellular-context profiles")
    manuscript = manuscript.replace("cellular profile", "cellular-context profile")
    manuscript = manuscript.replace("Cellular-profile rho", "Cellular-context similarity rho")

    manuscript = manuscript.replace("10^-4", "10⁻⁴").replace("10^-5", "10⁻⁵").replace("10^-6", "10⁻⁶").replace("10^-7", "10⁻⁷").replace("10^-9", "10⁻⁹")
    manuscript = manuscript.replace("Opposing ocular genetic architectures converge on shared retinal cellular contexts", "Cross-trait genetic architecture reveals shared retinal cellular contexts despite opposing genetic direction")
    out = SUB / "HGG_ADVANCES_MANUSCRIPT_v2.md"
    out.write_text(manuscript.strip() + "\n", encoding="utf-8")
    print(out)


if __name__ == "__main__":
    main()
