# Genomic hallmarks of DMPA-associated meningiomas — analysis code

Analysis code accompanying *"Genomic hallmarks of depot medroxyprogesterone
acetate–associated meningiomas"* (*Neuro-Oncology*, 2026; doi:10.1093/neuonc/noag136; PMID: 42258618).

The repository contains the two self-contained R pipelines used in the study:
DNA methylation and RNA-seq. Each script runs top to bottom, writes its figures
and tables to a results directory you set at the top of the file, and ends with a
`sessionInfo()` capture.

## Scripts

| File | Analysis |
|------|----------|
| `01_methylation_analysis.R` | EPIC/450k methylation pipeline (Sections 1–23): cohort assembly, ChAMP filtering, ComBat batch correction, MAD probe selection, consensus clustering, t-SNE/PCA, random-forest group assignment, SeSAMe copy number, PGR and 11q22.1 differential methylation, GlioSeq lollipop/oncoprint figures, and progesterone-pathway methylation. |
| `02_rnaseq_TRAF7_analysis.R` | TRAF7-mutant vs. TRAF7-wildtype differential expression within the DMPA cohort: edgeR/limma-voom, fgsea (MSigDB Hallmark and the progesterone GO terms GO:0032570 and GO:0050847), volcano, pathway barplot, and PGR/pathway expression boxplots. |

## Cohorts and patient labels

Patients are labeled **DMPA-1 … DMPA-10** throughout, matching the manuscript.
Two external reference cohorts are used in the methylation analysis:

- **Baylor** — Bayley et al., *Science Advances* 2022 (GEO **GSE189521**), EPIC v1
- **Heidelberg** — Capper et al., *Nature* 2018 (GEO **GSE109381**), 450k

The DMPA cohort was profiled on EPIC v2.

## Data availability

DMPA-cohort data are deposited in GEO under SuperSeries **GSE335177**
(BioProject **PRJNA1477512**), comprising two SubSeries:

- DMPA RNA-seq — FASTQs and `featureCounts` matrix: GEO **GSE335175** (→ SRA)
- DMPA DNA methylation — IDATs and processed beta matrix: GEO **GSE335176**

The two reference cohorts are already public (Baylor GSE189521, Heidelberg
GSE109381). For a citable, frozen snapshot of this code, archive a tagged
release on Zenodo and cite the resulting DOI.

## Running the code

Each script has a short configuration block at the top. In most cases you only
need to set one path.

- **`01_methylation_analysis.R`** — set `base_dir` (Section 1) to the project
  root containing the cohort data folders. Every other path is derived from it.
  The script runs end-to-end, or can be restarted at Section 10 by reading the
  saved beta/pd objects written by Sections 1–9.
- **`02_rnaseq_TRAF7_analysis.R`** — set `rnaseq_dir` (output) and `counts_file`
  (the raw `featureCounts` matrix). The pipeline begins from the count matrix;
  FASTQ-to-count generation followed the published Supplementary Methods. TRAF7
  status is defined by the `TRAF7_MUT` / `TRAF7_WT` vectors in Section 1.

Expected project layout for the methylation script (relative to `base_dir`):

```
<base_dir>/
├── Beta values/DMPA Idat files/        # DMPA EPIC v2 IDATs
├── Baylor/GSE189521_RAW/               # Baylor EPIC v1 IDATs
├── Baylor/GSE189521_Clinical_data (Bayley et al).xlsx
├── GSE109381_MNG_2018_Nature_paper/    # Heidelberg 450k IDATs
└── Combined Analysis_WHO I_Final/      # outputs (created automatically)
```

## Dependencies

R ≥ 4.3 with Bioconductor. Each script loads (and, where applicable, installs)
the packages it needs on first run. Principal packages:

- **Methylation:** `minfi`, `ChAMP`, `sva` (ComBat), `sesame`, `ConsensusClusterPlus`,
  `Rtsne`, `randomForest`, `ComplexHeatmap`, `limma`, `missMethyl`, `org.Hs.eg.db`
- **RNA-seq:** `edgeR`, `limma`, `fgsea`, `msigdbr` (v10+; uses the `collection=`
  argument), `org.Hs.eg.db`

Each script ends with `sessionInfo()`. For exact reproducibility, capture package
versions with `renv` (`renv::snapshot()` to write an `renv.lock`).

## De-identification

These scripts contain no protected health information — no names, dates, record
numbers, or local file paths. Patient identifiers are the manuscript labels
`DMPA-1` … `DMPA-10` only.
