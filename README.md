# Reproducible GEO Reanalysis of Allergen-Induced Airway Remodeling in Asthma

This repository contains a reproducible reanalysis workflow for GEO dataset `GSE193816`, a human single-cell RNA-seq study of endobronchial brush samples from allergic asthmatics and allergic non-asthmatic controls before and after airway challenge.

The goal is to provide a clean, independent, GitHub-ready analysis that can be run from processed GEO `.h5ad` files. It focuses on sample metadata curation, cell composition, airway epithelial remodeling signatures, T-cell and mononuclear phagocyte programs, and pseudo-bulk differential expression inputs.

## Biological Question

How does allergen exposure change airway epithelial, T-cell, and mononuclear phagocyte programs differently in allergic asthmatics compared with allergic non-asthmatic controls?

Main hypothesis tested:

> In response to allergen, asthmatic airways show stronger epithelial remodeling and type-2 inflammatory programs, while allergic non-asthmatic controls show more repair/pro-resolution responses.

Published findings are treated as hypotheses to reproduce and evaluate, not as assumed results.

## Study Design

This is a challenge-response single-cell reanalysis, not a simple asthma-versus-healthy comparison. The expected sample labels encode:

- `AA`: allergic asthma
- `ANA`: allergic non-asthmatic control
- `Pre`: baseline/pre-challenge
- `Dil`: diluent/saline challenge
- `Ag`: allergen challenge

Primary contrasts are interpreted at the donor/sample level whenever possible:

- `AA_Ag` vs `AA_Pre`: allergen response in asthmatics
- `ANA_Ag` vs `ANA_Pre`: allergen response in allergic non-asthmatic controls
- `AA_Ag` vs `ANA_Ag`: asthma-associated post-allergen state
- `AA_Pre` vs `ANA_Pre`: baseline asthma-associated state
- `(AA_Ag - AA_Pre) - (ANA_Ag - ANA_Pre)`: asthma-specific allergen response, when the available sample structure supports it

The scripts infer `subject`, `group`, and `condition` from sample metadata and write a clean sample-level table before downstream analysis.

## What This Reanalysis Adds

The original study authors provide public resources and code. This repository is not intended to copy their full workflow. It is a compact, independent reanalysis designed for reproducibility and portfolio review:

- starts from GEO-provided processed `.h5ad` files rather than unpublished local intermediates
- creates explicit sample metadata and sample-level composition summaries
- separates all-cell overview, epithelial, T-cell, MNP, and pseudo-bulk steps into readable scripts
- uses curated biological signatures to test epithelial remodeling and immune-state hypotheses
- exports pseudo-bulk matrices so differential expression can be modeled with donor/sample-aware designs
- documents limitations and avoids treating cells from the same donor as independent biological replicates

## Repository Status

This repository is intentionally safe for public GitHub upload:

- Large GEO `.h5ad` data files are ignored by Git.
- Generated figures, tables, and rendered reports are ignored by default.
- Placeholder `.gitkeep` files preserve the expected output folders.
- Download, analysis, and report scripts are tracked.
- Lightweight GitHub Actions checks validate Python syntax and R script parsing.

## Project Layout

```text
GSE193816_asthma_reanalysis/
├── .github/workflows/
├── data/
│   ├── raw_geo/
│   ├── processed/
│   └── metadata/
├── docs/
├── scripts/
├── notebooks/
├── results/
│   ├── figures/
│   ├── tables/
│   └── reports/
├── reports/
├── config.yaml
├── environment.yml
├── requirements.txt
└── README.md
```

## Data

GEO accession: `GSE193816`

The workflow uses processed supplementary `.h5ad` files from GEO. Raw FASTQ files are not included in GEO for this series, and the processed data are not committed to this repository because they are large.

The main files expected in `data/raw_geo/` are:

```text
GSE193816_all_cells_data.h5ad
GSE193816_AEC_data.h5ad
GSE193816_t_cell_data.h5ad
GSE193816_mnp_data.h5ad
GSE193816_all_data_raw_counts.h5ad
GSE193816_all_data_log_adjusted_counts.h5ad
```

## QC and Annotation Strategy

This first version uses the curated cell annotations and embeddings distributed with the processed GEO objects. It does not recluster the entire dataset from raw FASTQs because raw sequencing files are not available from GEO for this series.

The workflow therefore emphasizes:

- metadata validation and sample-level summaries
- inspection of available `.obs`, `.var`, and `.obsm` fields
- cell count and composition summaries by sample
- lineage-specific analysis using the authors' processed AEC, T-cell, and MNP objects
- pseudo-bulk aggregation by `sample x cell_type` for donor-aware downstream testing

If raw count-level reprocessing is needed, `GSE193816_all_data_raw_counts.h5ad` is used for pseudo-bulk export where available. Any additional filtering thresholds should be reported in the final rendered report before biological interpretation.

## Installation

```bash
conda env create -f environment.yml
conda activate gse193816
```

Or with pip:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Reproduce the Analysis

From the repository root:

```bash
python scripts/00_download_geo.py
bash scripts/run_all.sh
```

Or run each step manually:

```bash
python scripts/01_download_and_metadata.py
python scripts/02_all_cells_overview.py
python scripts/03_epithelial_signature_analysis.py
python scripts/04_tcell_analysis.py
python scripts/05_mnp_analysis.py
python scripts/06_pseudobulk_export.py
```

If R/DESeq2 is available:

```bash
Rscript scripts/06_pseudobulk_de.R
```

Render the report:

```bash
quarto render reports/08_final_report.qmd --output-dir ../results/reports
```

## Primary Outputs

The analysis creates the following local outputs:

```text
data/metadata/
├── gse193816_all_cells_cell_metadata.csv
├── gse193816_all_cells_gene_metadata.csv
└── sample_metadata_clean.csv

results/tables/
├── cell_counts_by_sample.csv
├── celltype_composition_by_sample.csv
├── aec_signature_scores_by_sample.csv
├── tcell_signature_scores_by_sample.csv
├── mnp_signature_scores_by_sample.csv
├── pseudobulk_counts_sample_by_celltype.csv
└── pseudobulk_metadata_sample_by_celltype.csv

results/figures/
├── all_cells_umap_group.png
├── all_cells_umap_condition.png
├── celltype_composition_stacked_bar.png
├── aec_signature_scores_by_condition.png
└── additional lineage-specific UMAP/signature plots
```

Generated outputs are ignored by Git so the repository stays small. For a portfolio release, selected final figures can be copied into `docs/figures/` and explicitly committed.

## Analysis Modules

- `01_download_and_metadata.py`: loads all-cell object, infers sample/group/condition fields, exports metadata.
- `02_all_cells_overview.py`: UMAP overview, cell counts, and composition tables.
- `03_epithelial_signature_analysis.py`: airway epithelial remodeling signatures and marker feature plots.
- `04_tcell_analysis.py`: TH2/IL9 and activation signature summaries.
- `05_mnp_analysis.py`: DC2, CCR2 monocyte, macrophage repair, and inflammatory MNP summaries.
- `06_pseudobulk_export.py`: exports sample-by-cell-type pseudo-bulk count matrices.
- `06_pseudobulk_de.R`: DESeq2 template for patient-aware contrasts where replicate structure allows.

## Statistical Interpretation

Single-cell measurements from the same donor are not independent biological replicates. For that reason:

- composition is summarized at the sample level
- signature scores are averaged to sample-level summaries before group-level interpretation
- differential expression is prepared as pseudo-bulk `sample x cell_type` count matrices
- paired formulas such as `~ subject + condition` are preferred when repeated measures exist
- interaction-style contrasts are treated cautiously because the cohort is small and some sample-condition combinations are missing

## Limitations

The cohort is small and GEO lists 21 samples, so not every subject has every condition. Composition and signature results should be interpreted as exploratory. Pseudo-bulk contrasts should be checked for sufficient biological replication before strong biological claims are made.

## Author

**Srilaxmi Nerella**  
UCSF profile: [profiles.ucsf.edu/srilaxmi.nerella](https://profiles.ucsf.edu/srilaxmi.nerella)  
Google Scholar: [scholar.google.com/citations?user=wjN338cAAAAJ](https://scholar.google.com/citations?user=wjN338cAAAAJ&hl=en)  
LinkedIn: [linkedin.com/in/srilaxmi-nerella-90000146](https://www.linkedin.com/in/srilaxmi-nerella-90000146)

## References

- GEO accession: [GSE193816](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193816)
- Original project resources: [Villani Lab asthma browser](https://villani.mgh.harvard.edu/allergy-asthma/)
- Original code resource: [villani-lab/airway_allergic_asthma](https://github.com/villani-lab/airway_allergic_asthma)
