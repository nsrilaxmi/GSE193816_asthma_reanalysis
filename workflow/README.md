# Snakemake Workflow

This folder provides an optional Snakemake entry point for the asthma reanalysis.

The script-based workflow remains supported:

```bash
bash scripts/run_all.sh
```

## Install Snakemake

Using conda or mamba:

```bash
conda install -c conda-forge -c bioconda snakemake
```

## Run the Python Analysis Workflow

From the repository root:

```bash
snakemake --snakefile workflow/Snakefile --cores 4
```

This target downloads the processed GEO `.h5ad` files if needed, creates cleaned metadata, runs all-cell overview, epithelial, T-cell, MNP, and pseudo-bulk export steps.

## Dry Run

```bash
snakemake --snakefile workflow/Snakefile --cores 1 --dry-run
```

The Snakemake workflow documents dependencies between the existing scripts. Optional R/DESeq2 and Quarto report-rendering steps remain available through `scripts/run_all.sh` because those depend on local R and Quarto installations.
