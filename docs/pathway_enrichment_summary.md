# Pathway Enrichment Summary

This first-pass enrichment layer tests whether committed pseudo-bulk top-gene lists overlap curated airway, immune, epithelial, and myeloid gene modules.

Because the committed inputs are compact top-50 preview tables rather than full ranked differential-expression results, p-values should be treated as exploratory overlap statistics rather than formal genome-wide enrichment claims.

## Top Overlaps

- **MNP AA Ag vs Pre (down in allergen)**: Macrophage lipid/complement program overlaps 7 genes (`APOC1;APOE;C1QA;C1QB;C1QC;MARCO;SCD`).
- **AEC ANA Ag vs Pre (up in allergen)**: Type 2 epithelial response overlaps 4 genes (`POSTN;SERPINB2;SERPINB4;SLC26A4`).
- **AEC AA Ag vs Pre (up in allergen)**: Type 2 epithelial response overlaps 3 genes (`POSTN;SERPINB4;SLC26A4`).
- **AEC AA Ag vs Pre (up in allergen)**: Epithelial injury/repair overlaps 2 genes (`RARRES1;TNC`).
- **AEC AA Ag vs Pre (up in allergen)**: Matrix remodeling overlaps 2 genes (`POSTN;TNC`).
- **AEC ANA Ag vs Pre (up in allergen)**: Epithelial injury/repair overlaps 2 genes (`CD55;RARRES1`).
- **AEC AA Ag vs Pre (down in allergen)**: Ciliated epithelial program overlaps 2 genes (`PIFO;SNTN`).
- **MNP AA Ag vs Pre (up in allergen)**: Type 2 epithelial response overlaps 1 gene (`ALOX15`).
- **AEC AA Ag vs Pre (up in allergen)**: Mucus/metaplasia overlaps 1 gene (`MUC2`).
- **AEC ANA Ag vs Pre (up in allergen)**: Macrophage lipid/complement program overlaps 1 gene (`SCD`).
- **AEC ANA Ag vs Pre (up in allergen)**: Matrix remodeling overlaps 1 gene (`POSTN`).
- **AEC ANA Ag vs Pre (down in allergen)**: Ciliated epithelial program overlaps 1 gene (`DNAH5`).

## Recommended Interpretation

- Use these results to prioritize biological themes for follow-up.
- Re-run enrichment on full ranked pseudo-bulk DE results before making strong pathway-level claims.
- Interpret broad-lineage T-cell and MNP results cautiously because rare subtype programs can be diluted.
