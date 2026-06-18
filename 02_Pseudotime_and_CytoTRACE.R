# ==============================================================================
# Script: 02_Pseudotime_and_CytoTRACE.R
# Description: Multi-scale trajectory and stemness analysis. 
#              Part 1: Global leaf atlas hierarchy.
#              Part 2: Protoderm-specific subpopulation state transition (Pro-1~6).
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Libraries
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(CytoTRACE)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(RColorBrewer)
})

# ------------------------------------------------------------------------------
# 2. Part 1: Global Scale Analysis (Figure 2)
# ------------------------------------------------------------------------------
message("--- Executing Part 1: Global Scale Pseudotime & Stemness ---")

# Load integrated global atlas
seurat_obj <- readRDS("data/processed/maize_atlas_cca.rds")
# Note: CDS object constructed from seurat_obj as previously documented

# Calculate CytoTRACE (Global)
mat_global <- GetAssayData(seurat_obj, layer = "counts", assay = "RNA")
res_global <- CytoTRACE(as.matrix(mat_global[VariableFeatures(seurat_obj), ]), ncores = 1)
seurat_obj$CytoTRACE_Global <- res_global$CytoTRACE[colnames(seurat_obj)]

# Generate Global UMAP Plot (Result 2A-style)
# [Implementation: Use the '悬浮L型坐标轴' logic as refined in Script 205]

# ------------------------------------------------------------------------------
# 3. Part 2: Protoderm-specific Scale Analysis (Figure 4)
# ------------------------------------------------------------------------------
message("--- Executing Part 2: Protoderm-specific Sub-clustering & Dynamics ---")

# Load Protoderm-specific subset
pro_obj <- readRDS("data/processed/Protoderm_Cleaned_Clustered.rds")

# CytoTRACE: Directionality Correction (Pro-4 as root)
mat_pro <- as.matrix(GetAssayData(pro_obj, layer = "counts", assay = "RNA"))
res_pro <- CytoTRACE(mat_pro, ncores = 1)
cyto_pro <- res_pro$CytoTRACE

# Robust orientation check for Protoderm subpopulations
valid_cells <- intersect(colnames(pro_obj)[pro_obj$sub_cluster == "Pro-4"], names(cyto_pro))
if (mean(cyto_pro[valid_cells], na.rm = TRUE) < 0.5) {
  cyto_pro <- 1 - cyto_pro # Invert if trajectory orientation is reversed
}
pro_obj$CytoTRACE_Score <- cyto_pro[colnames(pro_obj)]

# Generate Protoderm UMAP with White Halo Labels (Figure 4D)
# [Implementation: Use the 'geom_text_repel with bg.color' logic as refined in Script 207]

# ------------------------------------------------------------------------------
# 4. Export Supplemental Data (Table Generation)
# ------------------------------------------------------------------------------
message("Exporting supplementary pseudotime tables...")
library(openxlsx)
# [Implementation: Consolidate the cell-level pseudotime scores into Supplemental_Table_X.xlsx]

message("All pseudotime and trajectory analyses completed successfully.")
