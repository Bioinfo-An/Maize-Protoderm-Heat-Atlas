# ==============================================================================
# Script: 04_Network_and_Spatial.R
# Description: High-dimensional Weighted Gene Co-expression Network Analysis 
#              (hdWGCNA) across global and protoderm-specific scales, including 
#              hub gene identification and topological visualizations.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Libraries & Environment Setup
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(corrplot)
  library(openxlsx)
})

# Allow hdWGCNA to use large memory blocks
options(future.globals.maxSize = 10 * 1024^3)

data_dir <- "data/processed/"
out_dir <- "results/figures/hdWGCNA/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ==============================================================================
# PART 1: GLOBAL ATLAS hdWGCNA (Strict Network)
# ==============================================================================
message("--- Executing Part 1: Global hdWGCNA Network ---")

# Load integrated global atlas
sc_obj <- readRDS(file.path(data_dir, "maize_atlas_cca.rds"))
DefaultAssay(sc_obj) <- "RNA"

# 1A. Stringent Gene Filtering (Top 5000 highly expressed + ZmHsftf11 VIP inclusion)
expr_matrix <- GetAssayData(sc_obj, assay = "RNA", layer = "data")
cells_expressing <- rowSums(expr_matrix > 0)
top5k_genes <- head(names(sort(cells_expressing, decreasing = TRUE)), 5000)

target_gene <- "Zm00001eb060680" # ZmHsftf11
if (!(target_gene %in% top5k_genes)) top5k_genes <- c(top5k_genes, target_gene)
top5k_genes <- intersect(top5k_genes, rownames(sc_obj))
rm(expr_matrix, cells_expressing); gc()

# 1B. Metacell Construction & Soft Power Testing
sc_obj <- SetupForWGCNA(sc_obj, gene_select = "custom", gene_list = top5k_genes, wgcna_name = "Global_Strict")
sc_obj <- MetacellsByGroups(sc_obj, group.by = "cell_type", k = 25, max_shared = 10, ident.group = "cell_type")
sc_obj <- NormalizeMetacells(sc_obj)

sc_obj <- SetDatExpr(sc_obj, group_name = as.character(unique(sc_obj$cell_type)), group.by = "cell_type", assay = 'RNA', slot = 'data')
# sc_obj <- TestSoftPowers(sc_obj, networkType = 'signed') # Skipped in final run

# 1C. Network Construction (Using optimal soft power = 10)
sc_obj <- ConstructNetwork(sc_obj, soft_power = 10, setDatExpr = FALSE, tom_name = "Global_Strict_TOM", overwrite_tom = TRUE)

# 1D. Compute Module Eigengenes & Visualizations
sc_obj <- ModuleEigengenes(sc_obj)
sc_obj <- ModuleConnectivity(sc_obj)

# Bubble Plot (Module-Trait Relationship)
mods <- rownames(sc_obj[["MEs"]])[rownames(sc_obj[["MEs"]]) != "MEgrey"]
DefaultAssay(sc_obj) <- "MEs"
p_bubble_global <- DotPlot(sc_obj, features = mods, group.by = "cell_type", cluster.idents = TRUE, dot.scale = 12) + 
  scale_color_gradient(low = "#A8E6CF", high = "#FF8B94") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))
ggsave(file.path(out_dir, "FigS9_Global_Bubble.pdf"), p_bubble_global, width = 12, height = 7)
DefaultAssay(sc_obj) <- "RNA"

# Global Topological Network Plot (No Text for clean PPT layout)
pdf(file.path(out_dir, "FigS9_Global_Network_Clean.pdf"), width = 12, height = 12)
ModuleNetworkPlot(sc_obj, n_inner = 15, n_outer = 35, n_conns = 400, edge.alpha = 0.15, vertex.label.cex = 0)
dev.off()


# ==============================================================================
# PART 2: PROTODERM-SPECIFIC SUB-NETWORK
# ==============================================================================
message("--- Executing Part 2: Protoderm Subpopulation hdWGCNA ---")

sc_sub <- readRDS(file.path(data_dir, "Protoderm_Cleaned_Clustered.rds"))
sc_sub <- SetupForWGCNA(sc_sub, gene_select = "fraction", fraction = 0.05, wgcna_name = "Proto_WGCNA")

# 2A. Metacell & Network Construction (Using optimal soft power = 9)
sc_sub <- MetacellsByGroups(sc_sub, group.by = "sub_cluster", k = 20, max_shared = 10, ident.group = "sub_cluster")
sc_sub <- NormalizeMetacells(sc_sub)
sc_sub <- SetDatExpr(sc_sub, group_name = as.character(unique(sc_sub$sub_cluster)), group.by = "sub_cluster")

sc_sub <- ConstructNetwork(sc_sub, soft_power = 9, networkType = 'signed', corFnc = 'bicor', overwrite_tom = TRUE, tom_name = 'Proto_Subtype_TOM')
sc_sub <- ModuleEigengenes(sc_sub)
sc_sub <- ModuleConnectivity(sc_sub)

# 2B. Subpopulation Bubble Plot (Pale color scheme)
mods_sub <- colnames(GetMEs(sc_sub, harmonized = FALSE))
mods_sub <- mods_sub[mods_sub != "MEgrey"]

sc_sub <- AddMetaData(sc_sub, metadata = GetMEs(sc_sub, harmonized = FALSE))
p_bubble_sub <- DotPlot(sc_sub, features = mods_sub, group.by = "sub_cluster") + coord_flip() + 
  scale_color_gradientn(colors = c("#A1D99B", "#FFFFE0", "#FA8072")) + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "FigS11_Proto_Bubble.pdf"), p_bubble_sub, width = 8, height = 6)

# 2C. Module Eigengene Correlation (Ellipse matrix)
MEs <- GetMEs(sc_sub, harmonized = FALSE)
if("MEgrey" %in% colnames(MEs)) MEs <- MEs[, colnames(MEs) != "MEgrey"]
colnames(MEs) <- gsub("^ME", "", colnames(MEs))
cor_mat <- cor(MEs, method = "pearson")

pdf(file.path(out_dir, "FigS11_Proto_ME_Correlation.pdf"), width = 7, height = 6)
pale_palette <- colorRampPalette(c("#A1D99B", "white", "#FA8072"))(200)
corrplot(cor_mat, method = "ellipse", type = "full", col = pale_palette, tl.col = "black", tl.srt = 45)
dev.off()

# 2D. Protoderm Module UMAP 
sc_sub <- RunModuleUMAP(sc_sub, n_neighbors = 15, min_dist = 0.1)
p_nebula <- ModuleUMAPPlot(sc_sub, edge.alpha = 0.15, vertex.size = 1) + theme_void()
ggsave(file.path(out_dir, "FigS11_Proto_ModuleUMAP.pdf"), p_nebula, width = 8, height = 8)


# ==============================================================================
# PART 3: HUB GENE EXTRACTION & EXPORT
# ==============================================================================
message("--- Executing Part 3: Top Hub Gene Extraction ---")

modules_df <- GetModules(sc_sub)
valid_colors <- as.character(unique(modules_df$module[modules_df$module != "grey"]))
all_top30_list <- list()

for (mod in valid_colors) {
  genes_in_mod <- modules_df[modules_df$module == mod, ]
  kme_col <- paste0("kME_", mod)
  if(kme_col %in% colnames(genes_in_mod)) {
    sorted_genes <- genes_in_mod[order(genes_in_mod[[kme_col]], decreasing = TRUE), ]
    top30_genes <- sorted_genes[1:min(30, nrow(sorted_genes)), ]
    all_top30_list[[mod]] <- data.frame(
      Module = mod,
      Gene_Name = if("gene_name" %in% names(top30_genes)) top30_genes$gene_name else rownames(top30_genes),
      kME_Score = top30_genes[[kme_col]],
      stringsAsFactors = FALSE
    )
  }
}

final_hub_df <- bind_rows(all_top30_list)
# Note: Further annotation using reference databases was performed locally to merge 
# gene symbols and functional descriptions before final supplementary table generation.

message("hdWGCNA pipeline and spatial network construction completed successfully.")
