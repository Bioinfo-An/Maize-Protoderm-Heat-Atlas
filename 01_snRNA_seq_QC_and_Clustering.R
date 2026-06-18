# ==============================================================================
# Script: 01_snRNA_seq_QC_and_Clustering.R
# Description: Quality control, cell type annotation, cross-omics validation, 
#              and global UMAP visualization for the maize protoderm heat atlas.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Libraries
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggrepel)
  library(readxl)
  library(pheatmap)
})

# ------------------------------------------------------------------------------
# 2. Environment Setup & Data Loading
# ------------------------------------------------------------------------------
# Define standard relative paths for an open-source repository
data_dir <- "data/processed/"
out_dir <- "results/figures/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message("Loading single-nucleus RNA-seq data...")
# Note: The raw matrices are deposited in GSA (CRA043742).
# Here we load the pre-filtered and batch-corrected Seurat object.
seurat_obj <- readRDS(file.path(data_dir, "maize_atlas_cca.rds"))

# ------------------------------------------------------------------------------
# 3. Quality Control Metrics (Supplementary Fig. S1 A-D)
# ------------------------------------------------------------------------------
message("Generating QC metrics plots...")

# S1A: Sequencing Depth & Gene Detection (Violin + Boxplot)
qc_feats <- c("nFeature_RNA", "nCount_RNA")
p_s1a <- VlnPlot(seurat_obj, features = qc_feats, group.by = "orig.ident", 
                 ncol = length(qc_feats), pt.size = 0, 
                 cols = rep("#AEC7E8", length(unique(seurat_obj$orig.ident)))) + 
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) + 
  theme(axis.title.x = element_blank())
ggsave(file.path(out_dir, "FigS1A_Violin_QC.pdf"), p_s1a, width = 8, height = 6)

# S1B: Doublet Check (Scatter)
p_s1b <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
                        group.by = "orig.ident", pt.size = 0.5) +
  geom_smooth(method = "lm", color = "red", se = FALSE) + 
  ggtitle("Doublet Check: nCount vs nFeature") +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "FigS1B_Scatter_Doublet.pdf"), p_s1b, width = 6, height = 6)

# S1C: Library Complexity (Histogram)
median_val <- median(seurat_obj$nFeature_RNA)
p_s1c <- ggplot(seurat_obj@meta.data, aes(x = nFeature_RNA)) + 
  geom_histogram(bins = 60, fill = "grey70", color = "black", alpha = 0.8) +
  geom_vline(xintercept = median_val, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = median_val * 1.1, y = Inf, label = paste("Median:", median_val), 
           color = "red", hjust = 0, vjust = 2, fontface = "bold", size = 5) +
  theme_classic() +
  labs(title = "Library Complexity Distribution", x = "Detected Genes per Nucleus", y = "Cell Count")
ggsave(file.path(out_dir, "FigS1C_Complexity_Hist.pdf"), p_s1c, width = 6, height = 5)

# S1D: Batch Effect Mitigation Check (UMAP)
p_s1d <- DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident", 
                 pt.size = 0.1, shuffle = TRUE, raster = FALSE) +
  ggtitle("Data Integration (Batch Effect Check)") +
  theme(legend.position = "right")
ggsave(file.path(out_dir, "FigS1D_Batch_UMAP.pdf"), p_s1d, width = 7, height = 6)

# ------------------------------------------------------------------------------
# 4. Cross-omics Validation: snRNA-seq vs Bulk RNA-seq (Supplementary Fig. S1E)
# ------------------------------------------------------------------------------
message("Executing cross-omics validation...")

# Load and align Bulk RNA-seq matrix
bulk_raw <- read_excel("data/raw/bulk_gene_count_matrix.xlsx")
bulk_matrix <- as.matrix(bulk_raw[, -1])
rownames(bulk_matrix) <- toupper(gsub("[^a-zA-Z0-9]", "", bulk_raw[[1]]))

# Calculate CPM and log2 transform
bulk_cpm <- apply(bulk_matrix, 2, function(x) (x / sum(x)) * 1e6)
bulk_log <- log2(bulk_cpm + 1)

# Average biological replicates
bulk_A0_mean <- rowMeans(bulk_log[, c("A0h_1", "A0h_2")])
bulk_A2_mean <- rowMeans(bulk_log[, c("A2h_1", "A2h_2")])
bulk_A24_mean <- rowMeans(bulk_log[, c("A24h_1", "A24h_2")])

# Retrieve snRNA-seq pseudo-bulk log2(CPM+1) matrix (assuming pseudo_log is pre-calculated)
# Align overlapping genes
common_genes <- intersect(rownames(pseudo_log), rownames(bulk_log))
df_T0 <- data.frame(snRNA = pseudo_log[common_genes, "T0"], Bulk = bulk_A0_mean[common_genes], Timepoint = "T0")
df_T2 <- data.frame(snRNA = pseudo_log[common_genes, "T2"], Bulk = bulk_A2_mean[common_genes], Timepoint = "T2")
df_T24 <- data.frame(snRNA = pseudo_log[common_genes, "T24"], Bulk = bulk_A24_mean[common_genes], Timepoint = "T24")

df_all <- rbind(df_T0, df_T2, df_T24)
df_all$Timepoint <- factor(df_all$Timepoint, levels = c("T0", "T2", "T24"))

# Calculate Pearson correlations
label_df <- data.frame(
  Timepoint = factor(c("T0", "T2", "T24"), levels = c("T0", "T2", "T24")),
  label = c(sprintf("R = %.2f", cor(df_T0$snRNA, df_T0$Bulk)), 
            sprintf("R = %.2f", cor(df_T2$snRNA, df_T2$Bulk)), 
            sprintf("R = %.2f", cor(df_T24$snRNA, df_T24$Bulk)))
)

p_cross_omics <- ggplot(df_all, aes(x = snRNA, y = Bulk)) +
  geom_point(alpha = 0.15, size = 0.8, color = "grey30") +
  geom_smooth(method = "lm", color = "firebrick3", fill = "lightcoral", linewidth = 1) +
  facet_wrap(~Timepoint) +
  geom_text(data = label_df, aes(x = 1, y = max(df_all$Bulk) - 1, label = label), 
            hjust = 0, size = 5.5, fontface = "bold", color = "black") +
  labs(x = "snRNA-seq Expression [log2(CPM+1)]", y = "Bulk RNA-seq Expression [log2(CPM+1)]",
       title = "Transcriptome-wide Concordance: snRNA-seq vs Bulk RNA-seq") +
  theme_bw(base_size = 15) +
  theme(strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(face = "bold", size = 14),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(out_dir, "FigS1E_Cross_omics.png"), p_cross_omics, width = 12, height = 5, dpi = 600)

# ------------------------------------------------------------------------------
# 5. Global Atlas & Cell Type Annotation (Figure 1A)
# ------------------------------------------------------------------------------
message("Generating high-resolution transparent UMAP atlas...")

umap_df <- as.data.frame(seurat_obj@reductions$umap@cell.embeddings)
umap_df$cell_type <- Idents(seurat_obj)
umap_df$cluster_id <- if("seurat_clusters" %in% colnames(seurat_obj@meta.data)) seurat_obj$seurat_clusters else ""

cols_celltype <- c("Protoderm" = "#d62728", "Epidermis" = "#9467bd", 
                   "Mesophyll" = "#2ca02c", "Bundle Sheath" = "#1f77b4", 
                   "Companion" = "#ff7f0e")

label_data <- umap_df %>% group_by(cell_type) %>% summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

# UMAP Coordinate axes geometry
x_min <- min(umap_df$UMAP_1); x_range <- max(umap_df$UMAP_1) - x_min
y_min <- min(umap_df$UMAP_2); y_range <- max(umap_df$UMAP_2) - y_min
axis_origin_x <- x_min - (x_range * 0.05); axis_origin_y <- y_min - (y_range * 0.05)

p_main <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point(size = 0.2, alpha = 0.6) +
  geom_density_2d(aes(group = cell_type), bins = 4, linewidth = 0.3, linetype = "dashed", alpha = 0.6, show.legend = FALSE) +
  scale_color_manual(values = cols_celltype) +
  geom_text_repel(data = label_data, aes(label = cell_type), color = "black", bg.color = "white", bg.r = 0.15, size = 6, fontface = "bold", seed = 123) +
  geom_segment(aes(x = axis_origin_x, xend = axis_origin_x + (x_range*0.15), y = axis_origin_y, yend = axis_origin_y), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), color = "black", size = 0.8) +
  geom_segment(aes(x = axis_origin_x, xend = axis_origin_x, y = axis_origin_y, yend = axis_origin_y + (y_range*0.15)), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), color = "black", size = 0.8) +
  annotate("text", x = axis_origin_x + (x_range*0.075), y = axis_origin_y - (y_range*0.03), label = "UMAP_1", size = 4, fontface = "bold") +
  annotate("text", x = axis_origin_x - (x_range*0.03), y = axis_origin_y + (y_range*0.075), label = "UMAP_2", size = 4, fontface = "bold", angle = 90) +
  theme_void() + theme(legend.position = "none", aspect.ratio = 1, plot.background = element_rect(fill = "transparent", color = NA))

# Customized Legend
legend_data <- umap_df %>% group_by(cell_type) %>% summarise(clusters = paste(sort(unique(as.numeric(as.character(cluster_id)))), collapse = ", ")) %>% arrange(cell_type)
legend_data$y_pos <- seq(nrow(legend_data), 1)

p_legend <- ggplot(legend_data) +
  geom_point(aes(x = 0, y = y_pos, color = cell_type), size = 10) +
  geom_text(aes(x = 0.5, y = y_pos + 0.15, label = cell_type), hjust = 0, size = 6.5, fontface = "bold", color = "black") +
  geom_text(aes(x = 0.5, y = y_pos - 0.25, label = paste0("Clusters: ", clusters)), hjust = 0, size = 5, fontface = "italic", color = "grey40") + 
  scale_color_manual(values = cols_celltype) +
  scale_x_continuous(limits = c(-0.2, 10)) + scale_y_continuous(limits = c(0.5, max(legend_data$y_pos) + 0.5)) +
  theme_void() + theme(legend.position = "none", plot.background = element_rect(fill = "transparent", color = NA))

final_umap <- p_main + p_legend + plot_layout(widths = c(2, 1.0)) & theme(plot.background = element_rect(fill = "transparent", color = NA))
ggsave(file.path(out_dir, "Fig1_Maize_Atlas_Final.pdf"), final_umap, width = 13, height = 8, device = cairo_pdf, bg = "transparent")

# ------------------------------------------------------------------------------
# 6. Global Top Marker DotPlot (Supplementary Fig. S2)
# ------------------------------------------------------------------------------
message("Calculating specific markers for all 21 algorithmically defined clusters...")
DefaultAssay(seurat_obj) <- "RNA"
Idents(seurat_obj) <- "seurat_clusters"

cluster_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
top_markers <- cluster_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top_genes_list <- unique(top_markers$gene)

p_cluster_dot <- DotPlot(seurat_obj, features = top_genes_list, dot.scale = 5) + 
  scale_color_gradientn(colors = c("lightgrey", "#E64B35FF")) + 
  labs(x = "Top 3 Marker Genes per Seurat Cluster (Algorithmically Defined)", y = "Seurat Clusters (0 - 20)") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 9, angle = 90, hjust = 1, vjust = 0.5, face = "italic", color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", color = "black"),
        panel.grid.major = element_line(colour = "grey95"), panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "FigS2_All_Clusters_Markers.png"), p_cluster_dot, width = 16, height = 6, dpi = 600)

message("Pipeline Step 01 Completed Successfully.")
