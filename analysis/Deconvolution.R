library(arrow)
library(Seurat)
library(spacexr)

# Load the pre-made metadata from GitHub (skips needing to run FlexSingleCell.R)
MetaData <- read_parquet('~/Outputs/MetaData/P1CRC_Metadata.parquet')

# Load the single-cell count matrix
FlexOutPath <- "~/VisiumHD/PatientCRC1/outs"
FlexRef <- Read10X_h5(paste0(FlexOutPath, "/binned_outputs/square_008um/filtered_feature_bc_matrix.h5"))

# Filter to cell types with >25 cells
KpIdents <- names(which(table(MetaData$UnsupervisedL2) > 25))
MetaData <- MetaData[MetaData$UnsupervisedL2 %in% KpIdents, ]

# Fix labels
CTRef <- MetaData$UnsupervisedL2
CTRef <- gsub("/", "_", CTRef)
CTRef <- as.factor(CTRef)
names(CTRef) <- MetaData$barcode

# Build reference
reference <- Reference(FlexRef[, names(CTRef)], CTRef, colSums(FlexRef))

# Load spatial data
counts <- Read10X_h5("~/VisiumHD/PatientCRC1/outs/binned_outputs/square_008um/filtered_feature_bc_matrix.h5")
coords <- read_parquet("~/VisiumHD/PatientCRC1/outs/binned_outputs/square_008um/spatial/tissue_positions.parquet", as_data_frame = TRUE)
rownames(coords) <- coords$barcode
coords <- coords[colnames(counts), ]
coords <- coords[, 3:4]
nUMI <- colSums(counts)

# Run deconvolution
puck <- SpatialRNA(coords, counts, nUMI)
myRCTD <- create.RCTD(puck, reference, max_cores = 4)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')

# Save output
dir.create("~/Outputs/Deconvolution", recursive = TRUE, showWarnings = FALSE)
saveRDS(myRCTD, file = "~/Outputs/Deconvolution/PatientCRC1_Deconvolution_HD.rds")