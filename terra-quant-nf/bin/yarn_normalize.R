#!/usr/bin/env Rscript

# yarn_normalize.R
# YARN normalization of gene + TERRA counts with tissue-aware qsmooth.

suppressPackageStartupMessages({
  library(optparse)
  library(Biobase)
  library(yarn)
})

option_list <- list(
  make_option(c("-c", "--counts"), type = "character", default = NULL,
              help = "Gene counts CSV (without TERRA features)"),
  make_option(c("-t", "--terra_counts"), type = "character", default = NULL,
              help = "TERRA repeat counts CSV"),
  make_option(c("-a", "--annotation"), type = "character", default = NULL,
              help = "Sample annotation CSV/XLSX with 'Sample' column"),
  make_option(c("-k", "--target"), type = "character", default = NULL,
              help = "Column name in annotation for normalization grouping"),
  make_option(c("-o", "--outdir"), type = "character", default = ".",
              help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$counts) || is.null(opt$terra_counts) ||
    is.null(opt$annotation) || is.null(opt$target)) {
  stop("All options --counts, --terra_counts, --annotation, --target are required.")
}

# ---- Load data ----
cat("Loading gene counts...\n")
counts_df <- read.csv(opt$counts, row.names = 1, check.names = FALSE)
counts_mat <- as.matrix(counts_df)

cat("Loading TERRA counts...\n")
terra_df <- read.csv(opt$terra_counts, row.names = 1, check.names = FALSE)
terra_mat <- as.matrix(terra_df)

cat("Loading annotation...\n")
if (grepl("\\.xlsx$", opt$annotation)) {
  annotation <- readxl::read_xlsx(opt$annotation)
  annotation <- as.data.frame(annotation)
} else {
  annotation <- read.csv(opt$annotation, stringsAsFactors = FALSE)
}
rownames(annotation) <- annotation$Sample

# ---- Align sample order ----
common_samples <- intersect(colnames(counts_mat), rownames(annotation))
if (length(common_samples) == 0) {
  stop("No matching samples between counts and annotation. Check 'Sample' column.")
}
cat("Found", length(common_samples), "samples in common\n")

counts_mat <- counts_mat[, common_samples, drop = FALSE]
terra_mat <- terra_mat[, common_samples, drop = FALSE]
annotation <- annotation[common_samples, , drop = FALSE]

# ---- Build ExpressionSet ----
metadata <- data.frame(
  labelDescription = colnames(annotation),
  row.names = colnames(annotation)
)
anno_data <- new("AnnotatedDataFrame", data = annotation, varMetadata = metadata)
eset <- ExpressionSet(assayData = counts_mat, phenoData = anno_data)

# ---- Filter low-expression genes ----
cat("Filtering low-expression genes...\n")
eset_filtered <- filterLowGenes(eset, opt$target)
cat("Genes before filter:", nrow(eset), "\n")
cat("Genes after filter:", nrow(eset_filtered), "\n")

# Density plots
png(file.path(opt$outdir, "density_unfiltered.png"), width = 800, height = 600)
plotDensity(eset, opt$target, main = paste(opt$target, "- Unfiltered"))
dev.off()

png(file.path(opt$outdir, "density_filtered.png"), width = 800, height = 600)
plotDensity(eset_filtered, opt$target, main = paste(opt$target, "- Filtered"))
dev.off()

# ---- Combine filtered genes + TERRA ----
filtered_mat <- exprs(eset_filtered)
combined_mat <- rbind(terra_mat, filtered_mat)

# Write filtered counts with TERRA
combined_df <- data.frame(Gene = rownames(combined_mat), combined_mat, check.names = FALSE)
write.csv(combined_df, file.path(opt$outdir, "filtered_counts_with_TERRA.csv"), row.names = FALSE)

# ---- YARN normalization ----
cat("Running YARN qsmooth normalization...\n")
combined_eset <- ExpressionSet(assayData = combined_mat, phenoData = anno_data)
norm_eset <- normalizeTissueAware(combined_eset, opt$target, normalizationMethod = "qsmooth")

png(file.path(opt$outdir, "density_normalized.png"), width = 800, height = 600)
plotDensity(norm_eset, opt$target, normalized = TRUE, main = paste(opt$target, "- Normalized"))
dev.off()

# ---- Export results ----
norm_mat <- norm_eset@assayData$normalizedMatrix
norm_df <- data.frame(Gene = rownames(norm_mat), norm_mat, check.names = FALSE)
write.csv(norm_df, file.path(opt$outdir, "YARN_normalized_count_withTERRA.csv"), row.names = FALSE)

# Extract just TERRA
terra_norm <- norm_df[grep("^TERRA", norm_df$Gene), ]
write.csv(terra_norm, file.path(opt$outdir, "YARN_normalized_TERRA_count.csv"), row.names = FALSE)

cat("Done! Outputs written to:", opt$outdir, "\n")
