
# =============================================================================
#Script 2 — 02_scPBMC_reference_quantification
#panel_windows and panel_promoters are required by
# 02_scPBMC_reference_quantification.R if run in the same session.
# =============================================================================
## PURPOSE
#   Quantifies single-cell PBMC reference chromatin accessibility at the
#   same gene-centered +/- 10 kb windows used for bulk samples (defined in
#   01_quantify_accessibility.R), stratified by broad immune cell group.
#   Produces the per-gene scPBMC-max value used as the normalization
#   denominator in 03_pbmc_normalization_and_ranking.R.
#
# PRECONDITION (not included in this deposit)
#   This script expects a preprocessed Signac/Seurat object `pbmc` with:
#     - DefaultAssay(pbmc) == "ATAC", a ChromatinAssay built from the public
#       10x Genomics PBMC multiome dataset (filtered feature matrix +
#       fragments file)
#     - QC filtering applied using the following thresholds:
#         nCount_ATAC < 100000
#         nCount_RNA < 25000
#         nCount_ATAC > 1800
#         nCount_RNA > 1000
#         nucleosome_signal < 2
#         TSS.enrichment > 1
#     - Cell-type annotation via Azimuth-based label transfer against a
#       multimodal PBMC reference (TransferData on celltype.l2)
#     - A metadata column `broad_plot_group` assigning each cell to one of:
#       "B", "CD4_T", "CD8_T", "NK", "Mono", "DC", "Other"
#
#   This preprocessing follows the standard published Signac and Seurat
#   workflows (Stuart et al., 2021; Hao et al., 2021) and is not reproduced
#   here, as it does not constitute original code.
#
# INPUTS
#   - pbmc            : preprocessed Signac/Seurat object (see precondition)
#   - panel_windows   : GRanges of gene-centered +/- 10 kb windows with a
#                       `gene` metadata column (from 01_quantify_accessibility.R)
#
# OUTPUTS
#   - panel_scPBMC_summary : per-gene mean/max accessibility across the six
#                            retained immune groups, excluding "Other"
#                            (Supplementary Table S8)
#
# DEPENDENCIES
#   Signac, Seurat, ggplot2, GenomicRanges
# =============================================================================

library(Signac)
library(Seurat)
library(ggplot2)
library(GenomicRanges)

DefaultAssay(pbmc) <- "ATAC"

# =============================================================================
# 1. Per-cell-group signal extraction
# =============================================================================
#
# get_scPBMC_group_signal():
#   For a given genomic region and a given broad immune group, generates a
#   Signac CoveragePlot restricted to that group, then extracts the
#   underlying numeric signal values (y / ymax) from the built ggplot
#   object. Returns the mean and maximum of those values as a summary of
#   that group's accessibility across the region.
#
#   This is a workaround for obtaining per-group coverage summary statistics
#   directly from CoveragePlot's rendered output, rather than from a
#   dedicated accessor.

get_scPBMC_group_signal <- function(seu, region_string, group_name) {
  p <- CoveragePlot(
    object   = seu,
    region   = region_string,
    group.by = "broad_plot_group",
    idents   = group_name,
    annotation = FALSE,
    peaks      = FALSE,
    links      = FALSE
  )
  
  gb <- ggplot_build(p)
  
  y_vals <- c()
  for (i in seq_along(gb$data)) {
    d <- gb$data[[i]]
    if ("y" %in% colnames(d))    y_vals <- c(y_vals, d$y)
    if ("ymax" %in% colnames(d)) y_vals <- c(y_vals, d$ymax)
  }
  
  y_vals <- y_vals[is.finite(y_vals)]
  
  if (length(y_vals) == 0) {
    return(c(mean_signal = NA, max_signal = NA))
  }
  
  c(
    mean_signal = mean(y_vals, na.rm = TRUE),
    max_signal  = max(y_vals, na.rm = TRUE)
  )
}

# =============================================================================
# 2. Quantify scPBMC signal across all panel genes x all immune groups
# =============================================================================

scPBMC_groups <- sort(unique(pbmc$broad_plot_group))

panel_scPBMC_detail_list <- list()

for (i in seq_len(length(panel_windows))) {
  gene_name <- mcols(panel_windows)$gene[i]
  region_string <- paste0(
    as.character(seqnames(panel_windows)[i]), "-",
    start(panel_windows)[i], "-",
    end(panel_windows)[i]
  )
  
  for (grp in scPBMC_groups) {
    sig <- get_scPBMC_group_signal(pbmc, region_string, grp)
    
    panel_scPBMC_detail_list[[length(panel_scPBMC_detail_list) + 1]] <- data.frame(
      gene               = gene_name,
      group              = grp,
      scPBMC_mean_signal = sig["mean_signal"],
      scPBMC_max_signal  = sig["max_signal"]
    )
  }
}

panel_scPBMC_detail <- do.call(rbind, panel_scPBMC_detail_list)

# =============================================================================
# 3. Summarize scPBMC reference excluding the "Other" group
# =============================================================================

panel_scPBMC_no_other <- subset(panel_scPBMC_detail, group != "Other")

panel_scPBMC_summary <- aggregate(
  cbind(scPBMC_mean_signal, scPBMC_max_signal) ~ gene,
  data = panel_scPBMC_no_other,
  FUN  = function(x) c(mean = mean(x, na.rm = TRUE), max = max(x, na.rm = TRUE))
)

panel_scPBMC_summary <- data.frame(
  gene                  = panel_scPBMC_summary$gene,
  scPBMC_mean_no_other  = panel_scPBMC_summary$scPBMC_mean_signal[, "mean"],
  scPBMC_max_no_other   = panel_scPBMC_summary$scPBMC_max_signal[, "max"]
)

write.csv(
  panel_scPBMC_summary,
  file = "Supplementary_Table_S8_scPBMC_summary.csv",
  row.names = FALSE
)
# =============================================================================
# panel_scPBMC_summary is required by
# 03_pbmc_normalization_and_ranking.R.
# =============================================================================
