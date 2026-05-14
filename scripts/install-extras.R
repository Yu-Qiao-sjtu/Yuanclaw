# =====================================================
# YuanSeq - Install extra packages not in conda/system
# Run via: pixi run install-extras
# =====================================================

cat("Installing extra packages...\n\n")

# --- GseaVis (GitHub, optional for GSEA visualization) ---
cat("GseaVis (optional, for enhanced GSEA plots)...\n")
tryCatch({
  if (!requireNamespace("GseaVis", quietly = TRUE)) {
    remotes::install_github("junjunlab/GseaVis", upgrade = "never")
    cat("  Installed GseaVis from GitHub\n")
  } else {
    cat("  GseaVis already available\n")
  }
}, error = function(e) {
  message("  WARNING: GseaVis install failed: ", e$message)
  message("  GSEA module will use enrichplot as fallback\n")
})

cat("\nDone.\n")
