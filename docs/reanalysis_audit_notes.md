# Reanalysis Audit Notes

This document tracks scientific and portfolio-quality checks for the `GSE193816` reanalysis.

## Scientific Framing

- Treat the dataset as a challenge-response design: allergic asthma (`AA`) and allergic non-asthmatic controls (`ANA`) across `Pre`, `Dil`, and `Ag` conditions.
- Avoid describing the analysis as a simple disease-versus-healthy comparison.
- Keep conclusions focused on allergen-response programs in airway epithelial, T-cell, and mononuclear phagocyte compartments.

## Reproducibility

- Record exactly which GEO `.h5ad` objects are used.
- Do not commit large downloaded `.h5ad` files.
- Keep scripts runnable from repository root.
- Save inferred cell metadata, gene metadata, and sample metadata.

## Donor-Aware Analysis

- Composition summaries should use samples as the unit of observation.
- Signature scores should be summarized by sample before group-level interpretation.
- Differential expression should use pseudo-bulk counts aggregated by `sample x cell_type`.
- Cell-level tests may be used for exploration, but they should not be presented as donor-level evidence.

## Independence From Original Public Code

This repository should make clear what is original to this reanalysis:

- simplified reproducible workflow from GEO inputs
- explicit metadata curation
- sample-level composition outputs
- curated signature scoring across epithelial, T-cell, and MNP compartments
- pseudo-bulk export for donor-aware modeling
- concise Quarto report structure

If any code is copied or adapted from the original authors' repository, document it in the relevant script and README.

## Current Limitations To State

- small cohort size
- incomplete subject-by-condition coverage
- reliance on processed GEO objects and curated annotations
- pathway/signature results are exploratory until pseudo-bulk contrasts are reviewed for adequate replicate structure

