#!/usr/bin/env Rscript

# summarize_counts.R
# Merges individual HTSeq count files into a single matrix,
# then splits out TERRA repeat and subtelomeric counts.
#
# Usage: summarize_counts.R <counts_dir> [output_dir]

args <- commandArgs(trailingOnly = TRUE)
counts_dir <- if (length(args) >= 1) args[1] else "."
outdir     <- if (length(args) >= 2) args[2] else "."

# Find count files
count_files <- list.files(counts_dir, pattern = "\\.count\\.txt$", full.names = TRUE)
if (length(count_files) == 0) {
  stop("No .count.txt files found in: ", counts_dir)
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
write.csv(merged, file.path(outdir, "raw_counts_matrix.csv"), row.names = FALSE)
cat("Wrote raw_counts_matrix.csv\n")

# Identify TERRA features
terra_repeat <- merged[grepl("TERRA.*repeat", merged$Gene, ignore.case = TRUE), ]
terra_subtelo <- merged[grepl("TERRA.*subtelo", merged$Gene, ignore.case = TRUE), ]

# Also grab ITS features
its_features <- merged[grepl("_ITS", merged$Gene), ]

# Combine TERRA subtelo + ITS into one table
terra_subtelo_all <- rbind(terra_subtelo, its_features)

# Gene counts without any TERRA/ITS features
terra_its_genes <- c(terra_repeat$Gene, terra_subtelo_all$Gene)
gene_counts <- merged[!merged$Gene %in% terra_its_genes, ]

# Write outputs
write.csv(terra_repeat, file.path(outdir, "TERRA_repeat_counts.csv"), row.names = FALSE)
write.csv(terra_subtelo_all, file.path(outdir, "TERRA_subtelo_counts.csv"), row.names = FALSE)
write.csv(gene_counts, file.path(outdir, "gene_counts_no_TERRA.csv"), row.names = FALSE)

cat("Wrote TERRA_repeat_counts.csv\n")
cat("Wrote TERRA_subtelo_counts.csv\n")
cat("Wrote gene_counts_no_TERRA.csv\n")
cat("Done!\n")
