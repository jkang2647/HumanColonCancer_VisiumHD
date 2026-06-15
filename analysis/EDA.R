# ── EDA.R ─────────────────────────────────────────────────────────────────────
# Exploratory Data Analysis — P1CRC Visium HD Dataset
# Run this to understand the structure of your data before any analysis
# ──────────────────────────────────────────────────────────────────────────────

library(arrow)   # read parquet files
library(dplyr)   # data manipulation
library(rhdf5)   # read H5 count matrix

# ── MASTER PATHS (update only here if files move) ─────────────────────────────
meta_path    <- "/Users/jacobkang/Desktop/Outputs/MetaData/P1CRC_Metadata.parquet"
h5_path      <- "/Users/jacobkang/Desktop/VisiumHD/PatientCRC1/outs/binned_outputs/square_008um/filtered_feature_bc_matrix.h5"
pos_path     <- "/Users/jacobkang/Desktop/VisiumHD/PatientCRC1/outs/binned_outputs/square_008um/spatial/tissue_positions.parquet"
barcodes_path <- "/Users/jacobkang/Desktop/VisiumHD/PatientCRC1/outs/binned_outputs/square_008um/filtered_feature_bc_matrix/barcodes.tsv.gz"
features_path <- "/Users/jacobkang/Desktop/VisiumHD/PatientCRC1/outs/binned_outputs/square_008um/filtered_feature_bc_matrix/features.tsv.gz"

# ═══════════════════════════════════════════════════════════════════════════════
# PART 1: METADATA PARQUET
# This is the paper's pre-computed table — cell type calls, cluster labels, etc.
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n══════════════════════════════════════\n")
cat("PART 1: P1CRC_Metadata.parquet\n")
cat("══════════════════════════════════════\n")

meta <- read_parquet(meta_path)

cat(sprintf("\nShape: %d rows × %d columns\n", nrow(meta), ncol(meta)))
cat("\nColumn names and types:\n")
for (col in colnames(meta)) {
  cat(sprintf("  %-35s %s\n", col, class(meta[[col]])))
}

cat("\nFirst 3 rows:\n")
print(head(meta, 3))

cat("\nUnique values per categorical column:\n")
for (col in colnames(meta)) {
  vals <- unique(meta[[col]])
  if (length(vals) <= 25 && !is.numeric(meta[[col]])) {
    cat(sprintf("  %s:\n    %s\n", col, paste(sort(vals), collapse=", ")))
  }
}

# ═══════════════════════════════════════════════════════════════════════════════
# PART 2: SPATIAL POSITIONS
# Maps each bin barcode to its x/y location on the tissue
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n══════════════════════════════════════\n")
cat("PART 2: tissue_positions.parquet (square_008um)\n")
cat("══════════════════════════════════════\n")

pos <- read_parquet(pos_path)
cat(sprintf("\nShape: %d rows × %d columns\n", nrow(pos), ncol(pos)))
cat("\nColumn names:", paste(colnames(pos), collapse=", "), "\n")
cat("\nFirst 3 rows:\n")
print(head(pos, 3))
cat(sprintf("\nSpatial range — pxl_col_in_fullres (X): %.0f to %.0f\n",
            min(pos$pxl_col_in_fullres), max(pos$pxl_col_in_fullres)))
cat(sprintf("Spatial range — pxl_row_in_fullres (Y): %.0f to %.0f\n",
            min(pos$pxl_row_in_fullres), max(pos$pxl_row_in_fullres)))
cat(sprintf("Bins on tissue: %d\n", sum(pos$in_tissue == 1)))

# ═══════════════════════════════════════════════════════════════════════════════
# PART 3: GENE/BARCODE LISTS
# Features = genes. Barcodes = bin IDs.
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n══════════════════════════════════════\n")
cat("PART 3: Features (genes) and Barcodes (bins)\n")
cat("══════════════════════════════════════\n")

barcodes <- read.table(barcodes_path, header=FALSE)$V1
features <- read.table(features_path, header=FALSE, sep="\t")

cat(sprintf("\nNumber of bins (barcodes): %d\n", length(barcodes)))
cat(sprintf("Number of genes (features): %d\n", nrow(features)))
cat("\nFirst 5 barcodes:\n")
print(head(barcodes, 5))
cat("\nFirst 5 genes (ID, name, type):\n")
print(head(features, 5))

# ═══════════════════════════════════════════════════════════════════════════════
# PART 4: H5 COUNT MATRIX STRUCTURE
# The actual gene expression counts — genes × bins
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n══════════════════════════════════════\n")
cat("PART 4: filtered_feature_bc_matrix.h5 structure\n")
cat("══════════════════════════════════════\n")

h5_contents <- h5ls(h5_path)
print(h5_contents)

cat("\n══════════════════════════════════════\n")
cat("DONE — paste this output back to Claude\n")
cat("══════════════════════════════════════\n")