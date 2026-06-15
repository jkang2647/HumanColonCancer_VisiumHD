library(arrow)
library(Seurat)
library(tidyverse)
library(scattermore)
library(patchwork)
library(RColorBrewer)
library(paletteer)
library(rjson)
library(readbitmap)
library(grid)


# Source AuxFunctions exactly as the paper does
source("/Users/jacobkang/VisiumHD/analysis/AuxFunctions.R")

PATH <- "/Users/jacobkang/VisiumHD/PatientCRC1/outs"

# Load spatial data using paper's own GenerateSampleData function
cat("Loading spatial data...\n")
SampleData    <- GenerateSampleData(PATH, size = "008um")
bcs           <- SampleData$bcs
images_tibble <- SampleData$images_tibble

# Load pre-computed metadata
cat("Loading metadata...\n")
MetaData <- read_parquet('/Users/jacobkang/Outputs/MetaData/P1CRC_Metadata.parquet')

# Merge metadata into bcs by barcode
bcs <- bcs %>%
  left_join(MetaData %>% select(barcode, DeconvolutionClass, DeconvolutionLabel1,
                                 DeconvolutionLabel2, Periphery,
                                 UnsupervisedL1, UnsupervisedL2,
                                 MacrophageSubtype, GobletSubcluster),
            by = "barcode")

# Paper's exact color palette
Colors <- ColorPalette()

# Fix label mismatches to match actual data values
periphery_colors <- c(
  "Tumor"      = "#E41A1C",
  "50 micron"  = "#377EB8",
  "Tissue"     = "#AAAAAA"
)

macro_colors <- c(
  "Macrophage-SELENOP+" = "#E41A1C",
  "Macrophage-SPP1+"    = "#377EB8"
)

# Dynamic palette for unsupervised clusters
pal          <- as.vector(paletteer::paletteer_d("ggsci::default_igv"))
bcs_tissue   <- bcs %>% filter(tissue == 1)
unsup_types  <- sort(unique(na.omit(bcs_tissue$UnsupervisedL2)))
unsup_colors <- setNames(pal[seq_along(unsup_types)], unsup_types)

# Deconvolution colors using paper's palette
deconv_types  <- sort(unique(na.omit(bcs_tissue$DeconvolutionLabel1)))
deconv_colors <- Colors[names(Colors) %in% deconv_types]
missing       <- deconv_types[!deconv_types %in% names(deconv_colors)]
if (length(missing) > 0) {
  extra         <- setNames(pal[seq_len(length(missing))], missing)
  deconv_colors <- c(deconv_colors, extra)
}

# Helper: base plot with H&E background using paper's coordinate system
img_h <- images_tibble$height[[1]]
img_w <- images_tibble$width[[1]]

PlotBase <- function() {
  ggplot() +
    annotation_custom(
      grob  = images_tibble$grob[[1]],
      xmin  = 0,     xmax  = img_w,
      ymin  = -img_h, ymax = 0
    ) +
    coord_cartesian(
      xlim   = c(0, img_w),
      ylim   = c(-img_h, 0),
      expand = FALSE
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.text     = element_text(size = 6),
      legend.key.size = unit(0.3, "cm"),
      plot.title      = element_text(size = 10, face = "bold")
    )
}

# Plot 1: Unsupervised clustering — Figure 3a
cat("Making plot 1: unsupervised cell types...\n")
p1 <- PlotBase() +
  geom_scattermore(
    data    = bcs_tissue %>% filter(!is.na(UnsupervisedL2)),
    mapping = aes(x = imagecol_scaled, y = -imagerow_scaled, color = UnsupervisedL2),
    pointsize = 1.2, pixels = rep(2000, 2), alpha = 0.85
  ) +
  scale_color_manual(values = unsup_colors, na.value = "transparent") +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  ggtitle("Unsupervised clustering (level 2) — P1CRC")

# Plot 2: Deconvolution labels — Figure 3b
cat("Making plot 2: deconvolution labels...\n")
p2 <- PlotBase() +
  geom_scattermore(
    data    = bcs_tissue %>% filter(!is.na(DeconvolutionLabel1)),
    mapping = aes(x = imagecol_scaled, y = -imagerow_scaled, color = DeconvolutionLabel1),
    pointsize = 1.2, pixels = rep(2000, 2), alpha = 0.85
  ) +
  scale_color_manual(values = deconv_colors, na.value = "transparent") +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  ggtitle("Deconvolution cell type labels — P1CRC")

# Plot 3: Tumor periphery — Figure 4a
cat("Making plot 3: tumor periphery...\n")
p3 <- PlotBase() +
  geom_scattermore(
    data    = bcs_tissue %>% filter(!is.na(Periphery)),
    mapping = aes(x = imagecol_scaled, y = -imagerow_scaled, color = Periphery),
    pointsize = 1.5, pixels = rep(2000, 2), alpha = 0.9
  ) +
  scale_color_manual(values = periphery_colors, na.value = "transparent") +
  ggtitle("Tumor periphery (50 µm boundary) — P1CRC")

# Plot 4: Macrophage subtypes — Figure 4d
cat("Making plot 4: macrophage subtypes...\n")
p4 <- PlotBase() +
  geom_scattermore(
    data    = bcs_tissue %>% filter(!is.na(MacrophageSubtype)),
    mapping = aes(x = imagecol_scaled, y = -imagerow_scaled, color = MacrophageSubtype),
    pointsize = 2.5, pixels = rep(2000, 2), alpha = 0.95
  ) +
  scale_color_manual(values = macro_colors, na.value = "transparent") +
  ggtitle("Macrophage subtypes in tumor periphery — P1CRC")

# Save
cat("Saving figures...\n")
dir.create("/Users/jacobkang/Outputs/Figures", recursive = TRUE, showWarnings = FALSE)
ggsave("/Users/jacobkang/Outputs/Figures/P1CRC_Fig3a_celltypes.png",
       plot = p1, width = 12, height = 10, dpi = 200)
ggsave("/Users/jacobkang/Outputs/Figures/P1CRC_Fig3b_deconvolution.png",
       plot = p2, width = 12, height = 10, dpi = 200)
ggsave("/Users/jacobkang/Outputs/Figures/P1CRC_Fig4a_periphery.png",
       plot = p3, width = 10, height = 8, dpi = 200)
ggsave("/Users/jacobkang/Outputs/Figures/P1CRC_Fig4d_macrophages.png",
       plot = p4, width = 10, height = 8, dpi = 200)



# ── Plot 5: Gene expression spatial plots — Figure 3c ─────────────────────────
cat("Making gene expression plots...\n")

library(rhdf5)
library(Matrix)

# Read count matrix using rhdf5 directly (avoids hdf5r dependency)
h5_path <- "/Users/jacobkang/VisiumHD/PatientCRC1/outs/binned_outputs/square_008um/filtered_feature_bc_matrix.h5"

# Read barcodes and gene names from H5 file
barcodes_h5 <- h5read(h5_path, "matrix/barcodes")
genes_h5    <- h5read(h5_path, "matrix/features/name")
data_h5     <- h5read(h5_path, "matrix/data")
indices_h5  <- h5read(h5_path, "matrix/indices")
indptr_h5   <- h5read(h5_path, "matrix/indptr")
shape_h5    <- h5read(h5_path, "matrix/shape")

cat("Matrix shape:", shape_h5, "\n")
cat("Genes:", length(genes_h5), "\n")
cat("Barcodes:", length(barcodes_h5), "\n")

# Build sparse matrix
counts <- sparseMatrix(
  i = indices_h5 + 1,
  p = indptr_h5,
  x = as.numeric(data_h5),
  dims = shape_h5,
  dimnames = list(genes_h5, barcodes_h5)
)

cat("Sparse matrix built\n")

# Genes to plot — key markers from the paper Figure 3c
genes_to_plot <- c("SELENOP", "SPP1", "CD3E", "CEACAM6", "COL1A1", "PIGR")

# Match bcs barcodes to matrix columns
bc_idx <- match(bcs_tissue$barcode, colnames(counts))

for (gene in genes_to_plot) {
  if (!gene %in% rownames(counts)) {
    cat("Warning: gene", gene, "not found\n")
    next
  }
  
  cat("Plotting", gene, "...\n")
  
  # Get expression for tissue bins
  expr_vals <- as.numeric(counts[gene, bc_idx])
  
  plot_df <- bcs_tissue %>%
    mutate(expression = expr_vals)
  
  # Use paper's PlotExpression approach — only show bins with expression > 0
  p_gene <- PlotBase() +
    geom_scattermore(
      data    = plot_df %>% filter(expression > 0),
      mapping = aes(x     = imagecol_scaled,
                    y     = -imagerow_scaled,
                    color = log1p(expression)),
      pointsize = 1.5,
      pixels    = rep(2000, 2)
    ) +
    scale_color_viridis_c(
      option = "B",
      name   = "log(UMI+1)",
      limits = c(0, 3)
    ) +
    ggtitle(paste0(gene, " expression — P1CRC"))
  
  ggsave(
    paste0("/Users/jacobkang/Outputs/Figures/P1CRC_", gene, "_expression.png"),
    plot   = p_gene,
    width  = 10,
    height = 8,
    dpi    = 200
  )
}

cat("Gene expression plots done\n")

cat("Done!\n")