#!/usr/bin/env Rscript

# summarize_counts.R
# Merges individual HTSeq count files into a single matrix,
# then splits out TERRA repeat and subtelomeric counts.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
})

option_list <- list(
  make_option(c("-d", "--counts_dir"), type = "character", default = ".",
              help = "Directory containing *.count.txt files"),
  make_option(c("-o", "--output_prefix"), type = "character", default = ".",
              help = "Output directory/prefix for result files")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Find count files
count_files <- list.files(opt$counts_dir, pattern = "\\.count\\.txt$", full.names = TRUE)
if (length(count_files) == 0) {
  stop("No .count.txt files found in: ", opt$counts_dir)
}

cat("Found", length(count_files), "count files\n")

# Read and merge
read_counts <- function(f) {
  df <- read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                   col.names = c("Gene", "Count"))
  # Remove HTSeq summary lines
  df <- df[!grepl("^__", df$Gene), ]
  sample_name <- sub("\\.count\\.txt$", "", basename(f))
  colnames(df)[2] <- sample_name
  return(df)
}

counts_list <- lapply(count_files, read_counts)

# Merge all into one data frame
merged <- counts_list[[1]]
if (length(counts_list) > 1) {
  for (i in 2:length(counts_list)) {
    merged <- merge(merged, counts_list[[i]], by = "Gene", all = TRUE)
  }
}

# Replace NA with 0
merged[is.na(merged)] <- 0

# Write full raw counts matrix
outdir <- opt$output_prefix
write.csv(merged, file.path(outdir, "raw_counts_matrix.csv"), row.names = FALSE)
cat("Wrote raw_counts_matrix.csv\n")

# Identify TERRA features
terra_repeat <- merged %>% filter(grepl("TERRA.*repeat", Gene, ignore.case = TRUE))
terra_subtelo <- merged %>% filter(grepl("TERRA.*subtelo", Gene, ignore.case = TRUE))

# Also grab ITS features
its_features <- merged %>% filter(grepl("_ITS", Gene))

# Combine TERRA subtelo + ITS into one table
terra_subtelo_all <- bind_rows(terra_subtelo, its_features)

# Gene counts without any TERRA/ITS features
terra_its_genes <- c(terra_repeat$Gene, terra_subtelo_all$Gene)
gene_counts <- merged %>% filter(!Gene %in% terra_its_genes)

# Write outputs
write.csv(terra_repeat, file.path(outdir, "TERRA_repeat_counts.csv"), row.names = FALSE)
write.csv(terra_subtelo_all, file.path(outdir, "TERRA_subtelo_counts.csv"), row.names = FALSE)
write.csv(gene_counts, file.path(outdir, "gene_counts_no_TERRA.csv"), row.names = FALSE)

cat("Wrote TERRA_repeat_counts.csv\n")
cat("Wrote TERRA_subtelo_counts.csv\n")
cat("Wrote gene_counts_no_TERRA.csv\n")
cat("Done!\n")
