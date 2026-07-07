# Validation Against Original Study

This repository is a compact independent reanalysis of `GSE193816`, not a full reproduction of the original authors' complete workflow.

## What Matches The Original Study Direction

The original study reports that allergen challenge in allergic asthma is associated with airway epithelial remodeling, including type-2 response, mucus/metaplasia, matrix remodeling, and altered metabolism. In the current reanalysis:

- sample metadata recovers the expected `AA` and `ANA` groups and `Pre`, `Dil`, and `Ag` conditions
- GEO reports 21 samples; the curated sample metadata table also contains 21 sample-condition profiles
- airway epithelial pseudo-bulk `AA_Ag` vs `AA_Pre` top genes include remodeling/type-2-associated genes such as `POSTN`, `SERPINB4`, `SLC26A4`, `TNC`, and mucin-associated genes
- epithelial signature scores show higher sample-level type-2/IL13 and glycolysis scores in `AA_Ag` than in `AA_Pre`

These checks support the view that the first-pass epithelial analysis is directionally consistent with the original study.

## What Should Not Be Overclaimed Yet

The original paper emphasizes specific cellular states, including IL9-expressing pathogenic TH2 cells and DC2/CCR2 monocyte-derived cells after allergen challenge in asthma. The current preview analysis does not yet fully reproduce those claims because:

- T-cell and MNP pseudo-bulk previews are aggregated at broad lineage/cell-type levels
- rare or subtype-specific populations can be diluted in broad pseudo-bulk contrasts
- the current signature plots are exploratory sample-level summaries, not formal subtype-resolved tests
- the cohort is small and some donor-condition combinations are missing

Therefore, the repository should currently be described as a reproducible first-pass reanalysis and validation scaffold, with strongest support in the epithelial compartment.

## Recommended Next Validation Steps

To make the project scientifically stronger:

- run subtype-resolved pseudo-bulk contrasts within epithelial subclusters, T-cell subclusters, and MNP subclusters
- explicitly quantify `IL9`, `IL13`, `GATA3`, and TH2 signatures in T-cell subsets by sample
- quantify `CD1C`, `FCER1A`, `CLEC10A`, `CCR2`, `FCN1`, and monocyte/DC2 signatures in MNP subsets by sample
- add gene-level dot plots or violin plots for the specific marker genes emphasized by the original study
- compare top pseudo-bulk genes against the authors' published gene-contrast tables where available

## Portfolio Interpretation

As posted, this is useful as a reproducibility and workflow project. To become a stronger independent analysis project, it should add either:

- a more rigorous subtype-resolved validation of the original findings, or
- a different biological angle that was not the main focus of the original paper.

