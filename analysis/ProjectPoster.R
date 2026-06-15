# ── ProjectPoster.R ───────────────────────────────────────────────────────────
# Generates a visual HTML poster recap of the P1CRC Visium HD project
# Output: /Users/jacobkang/Desktop/Outputs/ProjectPoster.html
# ──────────────────────────────────────────────────────────────────────────────

library(arrow)
library(dplyr)

# ── Paths ─────────────────────────────────────────────────────────────────────
meta_path   <- "/Users/jacobkang/Desktop/Outputs/MetaData/P1CRC_Metadata.parquet"
figures_dir <- "/Users/jacobkang/Desktop/Outputs/Figures"
out_path    <- "/Users/jacobkang/Desktop/Outputs/ProjectPoster.html"

# ── Pull live stats from the parquet ──────────────────────────────────────────
cat("Reading metadata...\n")
meta <- read_parquet(meta_path)

n_total_bins   <- nrow(meta)
n_tissue_bins  <- sum(meta$tissue == 1, na.rm = TRUE)
n_genes        <- 18085  # from EDA
n_cols         <- ncol(meta)

deconv_counts <- meta %>%
  filter(!is.na(DeconvolutionLabel1)) %>%
  count(DeconvolutionLabel1, sort = TRUE)

unsup_types <- meta %>%
  filter(!is.na(UnsupervisedL1)) %>%
  pull(UnsupervisedL1) %>%
  unique() %>% sort()

periphery_counts <- meta %>%
  filter(!is.na(Periphery)) %>%
  count(Periphery)

get_pct <- function(label) {
  n <- periphery_counts %>% filter(Periphery == label) %>% pull(n)
  if (length(n) == 0) return("—")
  sprintf("%.1f%%", 100 * n / n_tissue_bins)
}

tumor_pct  <- get_pct("Tumor")
border_pct <- get_pct("50 micron")
tissue_pct <- get_pct("Tissue")

deconv_rows_html <- paste(
  apply(deconv_counts, 1, function(r) {
    pct <- round(100 * as.numeric(r["n"]) / n_tissue_bins, 1)
    sprintf(
      '<tr><td class="label-cell">%s</td><td class="num-cell">%s</td><td class="pct-cell">
        <div class="bar-bg"><div class="bar-fill" style="width:%.1f%%"></div></div>
      </td></tr>',
      r["DeconvolutionLabel1"],
      formatC(as.numeric(r["n"]), format = "d", big.mark = ","),
      min(pct, 100)
    )
  }),
  collapse = "\n"
)

figures <- list(
  list(name = "Fig 3b — Deconvolution map",    file = "P1CRC_Fig3b_deconvolution.png",  col = "#2dd4bf", desc = "RCTD cell type predictions across all tissue bins"),
  list(name = "Fig 4a — Tumor periphery",       file = "P1CRC_Fig4a_periphery.png",       col = "#2dd4bf", desc = "Tumor / 50µm border / stromal tissue zones"),
  list(name = "Fig 4d — Macrophage subtypes",   file = "P1CRC_Fig4d_macrophages.png",     col = "#2dd4bf", desc = "SPP1+ vs SELENOP+ macrophage spatial distribution"),
  list(name = "COL1A1 expression",              file = "P1CRC_COL1A1_expression.png",     col = "#2dd4bf", desc = "Stromal collagen marker"),
  list(name = "CD3E expression",                file = "P1CRC_CD3E_expression.png",        col = "#2dd4bf", desc = "Pan-T cell marker"),
  list(name = "CEACAM6 expression",             file = "P1CRC_CEACAM6_expression.png",    col = "#2dd4bf", desc = "Tumor epithelial marker"),
  list(name = "SELENOP expression",             file = "P1CRC_SELENOP_expression.png",    col = "#2dd4bf", desc = "Macrophage-SELENOP+ marker"),
  list(name = "SPP1 expression",                file = "P1CRC_SPP1_expression.png",        col = "#2dd4bf", desc = "Macrophage-SPP1+ marker"),
  list(name = "PIGR expression",                file = "P1CRC_PIGR_expression.png",        col = "#2dd4bf", desc = "Goblet / epithelial marker"),
  list(name = "Fig 3a — Unsupervised clusters", file = "P1CRC_Fig3a_celltypes.png",       col = "#f59e0b", desc = "⚠ Color palette mismatch — deprioritized")
)

fig_cards_html <- paste(
  lapply(figures, function(f) {
    img_path <- file.path(figures_dir, f$file)
    has_file <- file.exists(img_path)
    if (has_file) {
      img_b64  <- base64enc::base64encode(img_path)
      img_mime <- if (grepl("\\.png$", f$file)) "image/png" else "image/jpeg"
      img_tag  <- sprintf('<img src="data:%s;base64,%s" alt="%s">', img_mime, img_b64, f$name)
    } else {
      img_tag <- sprintf('<div class="img-placeholder">%s</div>', f$file)
    }
    status_dot <- if (f$col == "#f59e0b") "dot-warn" else "dot-done"
    sprintf('
      <div class="fig-card">
        <div class="fig-img-wrap">%s</div>
        <div class="fig-meta">
          <span class="fig-dot %s"></span>
          <span class="fig-name">%s</span>
          <span class="fig-desc">%s</span>
        </div>
      </div>', img_tag, status_dot, f$name, f$desc)
  }),
  collapse = "\n"
)

unsup_pills <- paste(
  sapply(unsup_types, function(t)
    sprintf('<span class="pill">%s</span>', t)),
  collapse = ""
)

dots_html <- paste(sapply(1:200, function(i) {
  colors <- c("#2dd4bf","#2dd4bf","#2dd4bf","#f87171","#f59e0b","#a78bfa","#34d399","#60a5fa")
  sprintf('<span style="background:%s"></span>', sample(colors, 1))
}), collapse="")

html <- '<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>P1CRC Project Poster</title>
<style>
  @import url("https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&family=IBM+Plex+Sans:wght@300;400;500;600&display=swap");

  *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }

  :root {
    --bg:       #080d14;
    --surface:  #0f1923;
    --border:   #1e2d3d;
    --teal:     #2dd4bf;
    --teal-dim: #0f766e;
    --amber:    #f59e0b;
    --coral:    #f87171;
    --text:     #e2e8f0;
    --muted:    #64748b;
    --mono:     "IBM Plex Mono", monospace;
    --sans:     "IBM Plex Sans", sans-serif;
  }

  body {
    background: var(--bg);
    color: var(--text);
    font-family: var(--sans);
    font-size: 14px;
    line-height: 1.6;
    padding: 48px 56px;
    max-width: 1200px;
    margin: 0 auto;
  }

  /* ── HEADER ── */
  .header { margin-bottom: 40px; }
  .eyebrow {
    font-family: var(--mono);
    font-size: 11px;
    letter-spacing: 0.12em;
    color: var(--teal);
    text-transform: uppercase;
    margin-bottom: 10px;
  }
  h1 {
    font-family: var(--mono);
    font-size: 32px;
    font-weight: 600;
    color: #fff;
    line-height: 1.15;
    margin-bottom: 6px;
  }
  .subtitle {
    font-size: 14px;
    color: var(--muted);
    font-weight: 300;
  }
  .header-rule {
    margin-top: 24px;
    height: 1px;
    background: linear-gradient(90deg, var(--teal) 0%, var(--teal-dim) 30%, var(--border) 100%);
  }

  /* ── DOT BANNER ── */
  .dot-banner {
    display: flex;
    gap: 3px;
    flex-wrap: wrap;
    margin: 28px 0;
    padding: 16px;
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 6px;
  }
  .dot-banner span {
    width: 8px; height: 8px;
    border-radius: 50%;
    display: inline-block;
    opacity: 0.85;
  }
  .dot-banner .label {
    width: auto; height: auto;
    border-radius: 0;
    font-family: var(--mono);
    font-size: 10px;
    color: var(--muted);
    align-self: center;
    margin: 0 8px;
    opacity: 1;
  }

  /* ── STAT ROW ── */
  .stat-row {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 16px;
    margin-bottom: 40px;
  }
  .stat-card {
    background: var(--surface);
    border: 1px solid var(--border);
    border-top: 2px solid var(--teal);
    padding: 20px 20px 16px;
    border-radius: 4px;
  }
  .stat-num {
    font-family: var(--mono);
    font-size: 28px;
    font-weight: 600;
    color: var(--teal);
    line-height: 1;
    margin-bottom: 6px;
  }
  .stat-label {
    font-size: 11px;
    color: var(--muted);
    text-transform: uppercase;
    letter-spacing: 0.08em;
    font-family: var(--mono);
  }

  /* ── SECTION HEADER ── */
  .section-head {
    font-family: var(--mono);
    font-size: 11px;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    color: var(--teal);
    margin-bottom: 16px;
    padding-bottom: 8px;
    border-bottom: 1px solid var(--border);
  }

  /* ── TWO-COL LAYOUT ── */
  .two-col {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 32px;
    margin-bottom: 40px;
  }

  /* ── METADATA TABLE ── */
  .meta-block {
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 4px;
    overflow: hidden;
  }
  .meta-table { width: 100%; border-collapse: collapse; }
  .meta-table th {
    background: #111e2b;
    font-family: var(--mono);
    font-size: 10px;
    color: var(--muted);
    text-transform: uppercase;
    letter-spacing: 0.08em;
    padding: 8px 14px;
    text-align: left;
    font-weight: 500;
  }
  .meta-table td {
    padding: 9px 14px;
    font-family: var(--mono);
    font-size: 12px;
    border-bottom: 1px solid var(--border);
    vertical-align: top;
  }
  .meta-table tr:last-child td { border-bottom: none; }
  .col-name { color: var(--teal); }
  .col-type { color: var(--muted); font-size: 11px; }
  .col-desc { color: var(--text); font-size: 11px; }

  /* ── DECONV TABLE ── */
  .deconv-block {
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 4px;
    overflow: hidden;
  }
  .deconv-table { width: 100%; border-collapse: collapse; }
  .deconv-table th {
    background: #111e2b;
    font-family: var(--mono);
    font-size: 10px;
    color: var(--muted);
    text-transform: uppercase;
    letter-spacing: 0.08em;
    padding: 8px 14px;
    text-align: left;
    font-weight: 500;
  }
  .label-cell { font-family: var(--mono); font-size: 12px; color: var(--text); padding: 8px 14px; }
  .num-cell   { font-family: var(--mono); font-size: 12px; color: var(--muted); padding: 8px 14px; white-space: nowrap; }
  .pct-cell   { padding: 8px 14px 8px 0; width: 40%; }
  .bar-bg     { background: var(--border); border-radius: 2px; height: 6px; }
  .bar-fill   { background: var(--teal); height: 6px; border-radius: 2px; }
  .deconv-table tr { border-bottom: 1px solid var(--border); }
  .deconv-table tr:last-child { border-bottom: none; }

  /* ── PERIPHERY ── */
  .zone-row {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 12px;
    margin-bottom: 40px;
  }
  .zone-card {
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 4px;
    padding: 18px 20px;
  }
  .zone-pct {
    font-family: var(--mono);
    font-size: 26px;
    font-weight: 600;
    line-height: 1;
    margin-bottom: 4px;
  }
  .zone-name { font-size: 11px; color: var(--muted); font-family: var(--mono); text-transform: uppercase; letter-spacing: 0.08em; }

  /* ── UNSUPERVISED PILLS ── */
  .pills-block {
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 4px;
    padding: 20px;
    margin-bottom: 40px;
    display: flex;
    flex-wrap: wrap;
    gap: 8px;
    align-items: center;
  }
  .pill {
    font-family: var(--mono);
    font-size: 11px;
    color: var(--teal);
    background: rgba(45, 212, 191, 0.08);
    border: 1px solid var(--teal-dim);
    padding: 4px 10px;
    border-radius: 3px;
    letter-spacing: 0.04em;
  }

  /* ── FIGURE GALLERY ── */
  .fig-grid {
    display: grid;
    grid-template-columns: repeat(5, 1fr);
    gap: 12px;
    margin-bottom: 40px;
  }
  .fig-card {
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 4px;
    overflow: hidden;
  }
  .fig-img-wrap { width: 100%; aspect-ratio: 1; overflow: hidden; background: #060b11; }
  .fig-img-wrap img { width: 100%; height: 100%; object-fit: cover; display: block; }
  .img-placeholder {
    width: 100%; height: 100%;
    display: flex; align-items: center; justify-content: center;
    font-family: var(--mono); font-size: 9px; color: var(--muted);
    text-align: center; padding: 8px;
  }
  .fig-meta {
    padding: 8px 10px;
    display: flex;
    flex-direction: column;
    gap: 2px;
  }
  .fig-dot {
    width: 6px; height: 6px;
    border-radius: 50%;
    display: inline-block;
    margin-bottom: 2px;
  }
  .dot-done { background: var(--teal); }
  .dot-warn { background: var(--amber); }
  .fig-name { font-family: var(--mono); font-size: 9px; color: var(--text); line-height: 1.3; }
  .fig-desc { font-size: 9px; color: var(--muted); line-height: 1.3; }

  /* ── NEXT STEPS ── */
  .next-grid {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 16px;
    margin-bottom: 40px;
  }
  .next-card {
    background: var(--surface);
    border: 1px solid var(--border);
    border-left: 3px solid var(--amber);
    border-radius: 4px;
    padding: 18px 20px;
  }
  .next-tag {
    font-family: var(--mono);
    font-size: 10px;
    color: var(--amber);
    text-transform: uppercase;
    letter-spacing: 0.1em;
    margin-bottom: 8px;
  }
  .next-title {
    font-size: 13px;
    font-weight: 500;
    color: var(--text);
    margin-bottom: 6px;
  }
  .next-desc { font-size: 12px; color: var(--muted); line-height: 1.5; }

  /* ── FOOTER ── */
  .footer {
    border-top: 1px solid var(--border);
    padding-top: 20px;
    display: flex;
    justify-content: space-between;
    font-family: var(--mono);
    font-size: 10px;
    color: var(--muted);
  }
</style>
</head>
<body>

<div class="header">
  <div class="eyebrow">UCSD Bioinformatics · Dr. Binfeng Lu Lab · Summer 2026</div>
  <h1>P1CRC Visium HD<br>Spatial Transcriptomics</h1>
  <div class="subtitle">Replication of Oliveira et al., Nature Genetics 2025 · Patient CRC1 · 10x Genomics Visium HD</div>
  <div class="header-rule"></div>
</div>

<!-- dot banner representing bins on tissue -->
<div class="dot-banner">
  <span class="label">507,684 bins on tissue →</span>
  @@DOTS@@
</div>

<!-- stat row -->
<div class="stat-row">
  <div class="stat-card">
    <div class="stat-num">@@TOTAL_BINS@@</div>
    <div class="stat-label">Total bins</div>
  </div>
  <div class="stat-card">
    <div class="stat-num">@@TISSUE_BINS@@</div>
    <div class="stat-label">On-tissue bins</div>
  </div>
  <div class="stat-card">
    <div class="stat-num">18,085</div>
    <div class="stat-label">Genes measured</div>
  </div>
  <div class="stat-card">
    <div class="stat-num">12</div>
    <div class="stat-label">Metadata columns</div>
  </div>
</div>

<!-- metadata columns + deconvolution side by side -->
<div class="section-head">Metadata table — P1CRC_Metadata.parquet</div>
<div class="two-col">

  <div class="meta-block">
    <table class="meta-table">
      <tr><th>Column</th><th>Type</th><th>Description</th></tr>
      <tr><td class="col-name">barcode</td><td class="col-type">chr</td><td class="col-desc">Unique bin ID (s_008um_RRRRR_CCCCC)</td></tr>
      <tr><td class="col-name">tissue</td><td class="col-type">factor</td><td class="col-desc">0 = off tissue, 1 = on tissue</td></tr>
      <tr><td class="col-name">X, Y</td><td class="col-type">num</td><td class="col-desc">Spatial coordinates in microns</td></tr>
      <tr><td class="col-name">DeconvolutionClass</td><td class="col-type">chr</td><td class="col-desc">singlet / doublet / reject</td></tr>
      <tr><td class="col-name">DeconvolutionLabel1</td><td class="col-type">chr</td><td class="col-desc">Broad RCTD cell type call</td></tr>
      <tr><td class="col-name">DeconvolutionLabel2</td><td class="col-type">chr</td><td class="col-desc">Fine RCTD cell type call</td></tr>
      <tr><td class="col-name">Periphery</td><td class="col-type">chr</td><td class="col-desc">Tumor / 50 micron / Tissue zone</td></tr>
      <tr><td class="col-name">UnsupervisedL1</td><td class="col-type">chr</td><td class="col-desc">Unsupervised cluster (broad)</td></tr>
      <tr><td class="col-name">UnsupervisedL2</td><td class="col-type">chr</td><td class="col-desc">Unsupervised cluster (fine)</td></tr>
      <tr><td class="col-name">MacrophageSubtype</td><td class="col-type">chr</td><td class="col-desc">SELENOP+ or SPP1+ macrophages</td></tr>
      <tr><td class="col-name">GobletSubcluster</td><td class="col-type">chr</td><td class="col-desc">Goblet cell subtypes (0–6)</td></tr>
    </table>
  </div>

  <div class="deconv-block">
    <table class="deconv-table">
      <tr><th>Cell type (DeconvolutionLabel1)</th><th>Bins</th><th>Proportion</th></tr>
      @@DECONV_ROWS@@
    </table>
  </div>

</div>

<!-- periphery zones -->
<div class="section-head">Tissue zone breakdown — Periphery column</div>
<div class="zone-row">
  <div class="zone-card">
    <div class="zone-pct" style="color:#f87171">@@TUMOR_PCT@@</div>
    <div class="zone-name">Tumor core</div>
  </div>
  <div class="zone-card">
    <div class="zone-pct" style="color:#f59e0b">@@BORDER_PCT@@</div>
    <div class="zone-name">50 µm border zone</div>
  </div>
  <div class="zone-card">
    <div class="zone-pct" style="color:#2dd4bf">@@TISSUE_PCT@@</div>
    <div class="zone-name">Stromal / normal tissue</div>
  </div>
</div>

<!-- unsupervised types -->
<div class="section-head">Unsupervised cell types identified — UnsupervisedL1</div>
<div class="pills-block">@@UNSUP_PILLS@@</div>

<!-- figure gallery -->
<div class="section-head">Figures generated — Outputs/Figures/</div>
<div class="fig-grid">@@FIG_CARDS@@</div>

<!-- next steps -->
<div class="section-head">Next directions — PI meeting notes</div>
<div class="next-grid">
  <div class="next-card">
    <div class="next-tag">Priority 1</div>
    <div class="next-title">BPCells — disk-backed matrix</div>
    <div class="next-desc">Load the 18,085 × 507,684 count matrix via BPCells for memory-efficient access. Required for neighborhood analysis compute.</div>
  </div>
  <div class="next-card">
    <div class="next-tag">Priority 2</div>
    <div class="next-title">Neighborhood analysis</div>
    <div class="next-desc">Build a spatial neighbor graph from X/Y coordinates. Identify recurring cell-type co-occurrence patterns (neighborhoods) across the tissue.</div>
  </div>
  <div class="next-card">
    <div class="next-tag">Priority 3</div>
    <div class="next-title">Sub-neighborhood analysis</div>
    <div class="next-desc">Find finer spatial structure within each neighborhood. Identify zones like immune-infiltrated tumor vs excluded T cell regions.</div>
  </div>
  <div class="next-card">
    <div class="next-tag">Priority 4</div>
    <div class="next-title">T cell + stromal spatial analysis</div>
    <div class="next-desc">Map CD3E+ T cells relative to Fibroblast, Smooth Muscle, and Endothelial bins. Are T cells infiltrating or excluded from stroma?</div>
  </div>
  <div class="next-card">
    <div class="next-tag">Priority 5</div>
    <div class="next-title">ACCESS server setup</div>
    <div class="next-desc">Connect with Atharva to set up compute server access. Needed for heavy workloads: neighborhood graph construction, BPCells at scale.</div>
  </div>
  <div class="next-card">
    <div class="next-tag">Open item</div>
    <div class="next-title">Unsupervised color palette fix</div>
    <div class="next-desc">Fig 3a failed due to mismatch between "Tumor-0" style labels and ColorPalette() in AuxFunctions.R. Deprioritized — deconvolution plot is higher quality.</div>
  </div>
</div>

<div class="footer">
  <span>P1CRC · Visium HD · square_008um · 10x Genomics</span>
  <span>Oliveira et al., Nature Genetics 2025 · Replication study</span>
  <span>Generated by ProjectPoster.R</span>
</div>

</body>
</html>'

html <- sub("@@DOTS@@",        dots_html,                                         html, fixed = TRUE)
html <- sub("@@TOTAL_BINS@@",  formatC(n_total_bins,  format = "d", big.mark = ","), html, fixed = TRUE)
html <- sub("@@TISSUE_BINS@@", formatC(n_tissue_bins, format = "d", big.mark = ","), html, fixed = TRUE)
html <- sub("@@DECONV_ROWS@@", deconv_rows_html,                                   html, fixed = TRUE)
html <- sub("@@TUMOR_PCT@@",   tumor_pct,                                          html, fixed = TRUE)
html <- sub("@@BORDER_PCT@@",  border_pct,                                         html, fixed = TRUE)
html <- sub("@@TISSUE_PCT@@",  tissue_pct,                                         html, fixed = TRUE)
html <- sub("@@UNSUP_PILLS@@", unsup_pills,                                        html, fixed = TRUE)
html <- sub("@@FIG_CARDS@@",   fig_cards_html,                                     html, fixed = TRUE)

writeLines(html, out_path)
cat(sprintf("\nPoster saved to: %s\n", out_path))
