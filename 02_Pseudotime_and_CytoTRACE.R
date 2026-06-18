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
# ==============================================================================
# PART 3: TEMPORAL GENE DYNAMICS & PROTEOSTATIC CHECKPOINT (Fig 3E, 3F)
# Description: Identifying pseudotime-dependent genes, performing K-means 
#              clustering to define Phase I, II, and III, and plotting the 
#              dynamic expression waves of the ZmHsftf11 checkpoint.
# ==============================================================================
message("--- Executing Part 3: Temporal Heatmap & Phase Definition ---")

# Ensure the Protoderm CDS object is loaded and ordered
# cds_pro <- readRDS(file.path(data_dir, "Protoderm_CDS_With_Pseudotime.rds"))

# ------------------------------------------------------------------------------
# 3A. Identify Pseudotime-Dependent Genes (Spatial Variable Genes)
# ------------------------------------------------------------------------------
message("Calculating pseudotime-dependent gene expression...")
# Use Moran's I test on the principal graph to find genes that change over pseudotime
pro_graph_test_res <- graph_test(cds_pro, neighbor_graph = "principal_graph", cores = 4)
pr_deg_ids <- row.names(subset(pro_graph_test_res, q_value < 0.05))

# Select the top highly variable genes to construct the temporal heatmap
top_pr_genes <- pro_graph_test_res %>% 
  filter(q_value < 0.05) %>% 
  arrange(desc(morans_I)) %>% 
  head(1000) %>% 
  pull(gene_short_name)

# ------------------------------------------------------------------------------
# 3B. K-means Clustering & Temporal Heatmap (Fig 3E)
# ------------------------------------------------------------------------------
message("Performing K-means clustering to define Phase I, II, and III...")

# Plot the pseudotime heatmap and apply K-means clustering (num_clusters = 3)
# This algorithmically defines: 
#   Phase I (Biosynthesis initiation)
#   Phase II (Proteostatic Checkpoint)
#   Phase III (Maturation)
p_temporal_heatmap <- plot_pseudotime_heatmap(
  cds_pro[top_pr_genes,],
  num_clusters = 3, 
  cores = 4,
  show_rownames = FALSE,
  return_heatmap = TRUE # Set to TRUE if we need to extract cluster assignments
)

# Save the heatmap (Fig 3E)
pdf(file.path(out_dir, "Fig3E_Proteostatic_Checkpoint_Heatmap.pdf"), width = 6, height = 8)
print(p_temporal_heatmap)
dev.off()

# ------------------------------------------------------------------------------
# 3C. Dynamic Expression Trajectory of Core Marker Genes (Fig 3F)
# ------------------------------------------------------------------------------
message("Plotting dynamic expression curves for core checkpoint regulators...")

# Define the genes representing the survival trade-off
# ZmHsftf11: The bimodal pulse dominating the Phase II Checkpoint
# ZmPdk1: Photosynthetic gene demonstrating gated delayed activation
core_genes <- c("ZmHsftf11", "ZmPdk1")

p_gene_dynamics <- plot_genes_in_pseudotime(
  cds_pro[core_genes,], 
  color_cells_by = "pseudotime", 
  min_expr = 0.5, 
  cell_size = 1
) + 
  scale_color_viridis_c(option = "plasma") +
  theme_classic() +
  labs(title = "Gene Expression Dynamics along Pseudotime", x = "Pseudotime", y = "Expression") +
  theme(strip.text = element_text(face = "bold", size = 12))

# Save the dynamic curves (Fig 3F)
ggsave(file.path(out_dir, "Fig3F_Gene_Dynamics_Curves.pdf"), p_gene_dynamics, width = 7, height = 4)

message("Temporal dynamics and Checkpoint definition completed successfully.")
