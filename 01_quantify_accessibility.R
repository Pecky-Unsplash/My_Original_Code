Script 1 — 01_quantify_accessibility
# =============================================================================
# PURPOSE
#   Defines the 52-gene placenta-associated panel, maps each gene to its
#   GRCh38/hg38 genomic coordinates (EnsDb.Hsapiens.v86), constructs
#   gene-centered +/- 10 kb windows and TSS-centered +/- 3 kb promoter
#   windows, and quantifies bulk maternal PBMC ATAC-seq accessibility
#   across both window types using a width-weighted continuous-signal
#   approach applied to normalized bigWig coverage tracks.
#
# INPUTS
#   - Normalized bulk ATAC-seq bigWig files (one per sample), referenced
#     via `bigwig_files` below.
#
# OUTPUTS
#   - panel_bulkPBMC_signal_df    : per-sample, per-gene +/- 10 kb signal
#                                   (Supplementary Table S4)
#   - promoter_bulkPBMC_signal_df : per-sample, per-gene +/- 3 kb promoter
#                                   signal (feeds Supplementary Table S5)
#
# DEPENDENCIES
#   EnsDb.Hsapiens.v86, GenomicRanges, GenomeInfoDb, rtracklayer
# =============================================================================

library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)

# ---- User-configurable paths -----------------------------------------------
data_dir <- "."  # directory containing the bigWig files

bigwig_files <- c(
  file.path(data_dir, "ATAC_Batch2_023B_ocs.bigWig"),
  file.path(data_dir, "ATAC_Batch2_039B_ocs.bigWig"),
  file.path(data_dir, "ATAC_Batch2_040B_ocs.bigWig"),
  file.path(data_dir, "ATAC_Batch2_049B_ocs.bigWig")
)

# =============================================================================
# 1. Define the 52-gene placenta-associated panel
# =============================================================================

panel_genes <- c(
  "KRT7","GATA3","TFAP2C","ELF5","TP63","TEAD4","ITGA6","EGFR","PEG10",
  "HLA-G","PAPPA2","ITGA1","ITGA5","MMP2","FLT1","NOTUM","PLAC8","HTRA4",
  "DIO2","TAC3","FN1","ASCL2",
  "CGA","CGB3","CGB5","CGB7","CGB8","ERVFRD-1","ERVW-1","GCM1","INSL4","PLAC1",
  "LGALS13","PSG1","PSG2","PSG3","PSG5","PSG6","PSG9",
  "CYP19A1","CSH1","CSH2","CSHL1","KISS1","SDC1","INHA",
  "NOTCH1","NOTCH2","VEGFA","EPAS1","ENG","SNAI1"
)

stopifnot(length(panel_genes) == 52)

# =============================================================================
# 2. Map panel genes to genomic coordinates
# =============================================================================

g <- genes(EnsDb.Hsapiens.v86)
seqlevels(g) <- paste0("chr", seqlevels(g))
seqlevels(g)[seqlevels(g) == "chrMT"] <- "chrM"
g <- keepStandardChromosomes(g, pruning.mode = "coarse")

panel_gene_gr <- g[mcols(g)$gene_name %in% panel_genes]
panel_gene_gr <- panel_gene_gr[match(panel_genes, mcols(panel_gene_gr)$gene_name)]

stopifnot(all(!is.na(panel_gene_gr)))

# =============================================================================
# 3. Build gene-centered +/- 10 kb windows
# =============================================================================

panel_windows <- GRanges(
  seqnames = seqnames(panel_gene_gr),
  ranges = IRanges(
    start = start(panel_gene_gr) - 10000,
    end   = end(panel_gene_gr) + 10000
  ),
  strand = strand(panel_gene_gr)
)
mcols(panel_windows)$gene <- mcols(panel_gene_gr)$gene_name

# =============================================================================
# 4. Build TSS-centered +/- 3 kb promoter windows
# =============================================================================

panel_promoters <- promoters(
  panel_gene_gr,
  upstream   = 3000,
  downstream = 3000
)
mcols(panel_promoters)$gene <- mcols(panel_gene_gr)$gene_name

# =============================================================================
# 5. Width-weighted continuous-signal quantification functions
# =============================================================================
#
# mean_score_over_window():
#   Computes a per-base width-weighted mean signal across a window from
#   bigWig intervals overlapping that window. Each interval's contribution
#   is weighted by the width of its overlap with the window, then the sum
#   is divided by the total window width. Windows with no overlapping
#   signal return 0.
#
# max_score_over_window():
#   Returns the maximum bigWig signal value within the imported region.

mean_score_over_window <- function(gr, window_start, window_end) {
  if (length(gr) == 0) return(0)
  df <- as.data.frame(gr)
  overlap_start <- pmax(df$start, window_start)
  overlap_end   <- pmin(df$end, window_end)
  widths <- overlap_end - overlap_start + 1
  widths[widths < 0] <- 0
  if (sum(widths) == 0) return(0)
  sum(df$score * widths) / (window_end - window_start + 1)
}

max_score_over_window <- function(gr) {
  if (length(gr) == 0) return(0)
  max(as.data.frame(gr)$score)
}

# =============================================================================
# 6. Quantify bulk PBMC signal across +/- 10 kb gene-centered windows
# =============================================================================

panel_bulkPBMC_signal_list <- list()

for (i in seq_len(length(panel_windows))) {
  gene_name <- mcols(panel_windows)$gene[i]
  region_gr <- panel_windows[i]
  
  for (bw in bigwig_files) {
    bw_import   <- import(bw, which = region_gr, format = "BigWig")
    sample_name <- sub("_ocs.bigWig$", "", basename(bw))
    
    panel_bulkPBMC_signal_list[[length(panel_bulkPBMC_signal_list) + 1]] <- data.frame(
      gene        = gene_name,
      sample      = sample_name,
      mean_signal = mean_score_over_window(bw_import, start(region_gr), end(region_gr)),
      max_signal  = max_score_over_window(bw_import)
    )
  }
}

panel_bulkPBMC_signal_df <- do.call(rbind, panel_bulkPBMC_signal_list)

write.csv(
  panel_bulkPBMC_signal_df,
  file = "Supplementary_Table_S4_panel_bulkPBMC_10kb_signal.csv",
  row.names = FALSE
)

# =============================================================================
# 7. Quantify bulk PBMC signal across +/- 3 kb promoter windows
# =============================================================================

promoter_bulkPBMC_signal_list <- list()

for (i in seq_len(length(panel_promoters))) {
  gene_name <- mcols(panel_promoters)$gene[i]
  region_gr <- panel_promoters[i]
  
  for (bw in bigwig_files) {
    bw_import   <- import(bw, which = region_gr, format = "BigWig")
    sample_name <- sub("_ocs.bigWig$", "", basename(bw))
    
    promoter_bulkPBMC_signal_list[[length(promoter_bulkPBMC_signal_list) + 1]] <- data.frame(
      gene        = gene_name,
      sample      = sample_name,
      mean_signal = mean_score_over_window(bw_import, start(region_gr), end(region_gr)),
      max_signal  = max_score_over_window(bw_import)
    )
  }
}

promoter_bulkPBMC_signal_df <- do.call(rbind, promoter_bulkPBMC_signal_list)

write.csv(
  promoter_bulkPBMC_signal_df,
  file = "promoter_bulkPBMC_signal_per_sample.csv",
  row.names = FALSE
)

# =============================================================================
# panel_windows and panel_promoters are required by
# 02_scPBMC_reference_quantification.R
# =============================================================================
