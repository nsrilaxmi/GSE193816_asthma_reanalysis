suppressPackageStartupMessages({
  library(DESeq2)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
  project_root <- normalizePath(file.path(dirname(script_path), ".."))
} else {
  project_root <- normalizePath(".")
}
tables_dir <- file.path(project_root, "results", "tables")

counts_path <- file.path(tables_dir, "pseudobulk_counts_sample_by_celltype.csv")
metadata_path <- file.path(tables_dir, "pseudobulk_metadata_sample_by_celltype.csv")

counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
meta <- read.csv(metadata_path, row.names = 1, check.names = FALSE)

run_simple_contrast <- function(celltype_pattern, group_keep, condition_a, condition_b, output_name) {
  keep <- grepl(celltype_pattern, meta$cell_type_clean, ignore.case = TRUE) &
    meta$group == group_keep &
    meta$condition %in% c(condition_a, condition_b) &
    meta$n_cells >= 20

  if (sum(keep) < 4) {
    message("Skipping ", output_name, ": not enough pseudo-bulk profiles.")
    return(invisible(NULL))
  }

  sub_meta <- droplevels(meta[keep, , drop = FALSE])
  sub_counts <- t(counts[rownames(sub_meta), , drop = FALSE])
  sub_counts <- sub_counts[rowSums(sub_counts) >= 10, , drop = FALSE]
  sub_meta$subject <- factor(sub_meta$subject)
  sub_meta$group <- factor(sub_meta$group)
  sub_meta$cell_type_clean <- factor(sub_meta$cell_type_clean)
  sub_meta$condition <- relevel(factor(sub_meta$condition), ref = condition_b)

  design_formula <- if (length(unique(sub_meta$subject)) < nrow(sub_meta)) {
    ~ subject + condition
  } else {
    ~ condition
  }
  design_matrix <- model.matrix(design_formula, data = sub_meta)
  if (qr(design_matrix)$rank < ncol(design_matrix)) {
    message("Using unpaired design for ", output_name, ": paired design is not full rank.")
    design_formula <- ~ condition
  }

  dds <- DESeqDataSetFromMatrix(
    countData = round(sub_counts),
    colData = sub_meta,
    design = design_formula
  )
  dds <- DESeq(dds, quiet = TRUE)
  res <- as.data.frame(results(dds, contrast = c("condition", condition_a, condition_b)))
  res$gene <- rownames(res)
  res <- res[, c("gene", setdiff(colnames(res), "gene"))]
  res <- res[order(res$padj, na.last = TRUE), ]

  write.csv(res, file.path(tables_dir, output_name), row.names = FALSE)
}

run_simple_contrast("epithelial|basal|secretory|ciliated|goblet", "AA", "Ag", "Pre", "pseudobulk_AEC_AA_Ag_vs_Pre.csv")
run_simple_contrast("epithelial|basal|secretory|ciliated|goblet", "ANA", "Ag", "Pre", "pseudobulk_AEC_ANA_Ag_vs_Pre.csv")
run_simple_contrast("T|CD4|CD8|TH", "AA", "Ag", "Pre", "pseudobulk_Tcell_AA_Ag_vs_Pre.csv")
run_simple_contrast("monocyte|macrophage|dendritic|DC|MNP|MPS", "AA", "Ag", "Pre", "pseudobulk_MNP_AA_Ag_vs_Pre.csv")
