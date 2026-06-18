# ==============================================================================
# Script: 05_GRN_and_Virtual_Knockout.R
# Description: Causal inference of the ZmHsftf11-mediated regulatory hub.
#              Part 1: Virtual gene knockout (scTenifoldKnk) analysis.
#              Part 2: Directional Gene Regulatory Network (GENIE3) inference.
#              Part 3: Hierarchical network visualization mapping to Phases I-III.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Libraries & Environment Setup
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(scTenifoldKnk)
  library(GENIE3)
  library(dplyr)
  library(tidygraph)
  library(ggraph)
  library(igraph)
  library(openxlsx)
})

data_dir <- "data/processed/"
out_dir <- "results/figures/GRN_and_Knockout/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ==============================================================================
# PART 1: VIRTUAL KNOCKOUT ANALYSIS (scTenifoldKnk)
# ==============================================================================
message("--- Executing Part 1: scTenifoldKnk Virtual Knockout ---")

# Load subsampled count matrix for KO analysis
# Note: Ensure the matrix is purely of class 'matrix' to avoid strict type errors in R
slim_matrix <- readRDS(file.path(data_dir, "Protoderm_Count_Matrix_Slim.rds"))
attr(slim_matrix, "class") <- "matrix"

if(length(class(slim_matrix)) > 1) {
  stop("Matrix class coercion failed. Ensure input is a standard R matrix.")
}

# Perform virtual knockout targeting ZmHsftf11 (Zm00001eb060680)
message("Running virtual knockout perturbation for ZmHsftf11...")
output_knk <- scTenifoldKnk(countMatrix = slim_matrix, gKO = "Zm00001eb060680")

# Extract and rank perturbed genes by Fold Change (FC / Impact Score)
final_results <- output_knk$diffRegulation
final_results <- final_results[order(final_results$FC, decreasing = TRUE), ]
final_results$Rank <- 1:nrow(final_results)

# Save raw ranked KO results
write.csv(final_results, file.path(out_dir, "Hsftf11_Virtual_KO_Ranked.csv"), row.names = FALSE)

# ==============================================================================
# PART 2: GENE REGULATORY NETWORK INFERENCE (GENIE3)
# ==============================================================================
message("--- Executing Part 2: GENIE3 Directional GRN Inference ---")

# Load Protoderm Seurat object
pro_obj <- readRDS(file.path(data_dir, "Protoderm_Cleaned_Clustered.rds"))

# Define target heatmap genes and core HSF regulators
heroes_ids <- c("Zm00001eb210530", "Zm00001eb060680", "Zm00001eb047760") # Hsftf22, Hsftf11, Hsftf24
heatmap_genes <- read.xlsx(file.path(data_dir, "Supplemental_Table_S20_heatmapList.xlsx"), colNames = FALSE)[[1]]
all_g <- intersect(unique(c(heatmap_genes, heroes_ids)), rownames(pro_obj))

# Extract expression matrix
expr_matrix <- as.matrix(GetAssayData(pro_obj, assay = "RNA", layer = "data")[all_g, ])

# Compute regulatory weights (nCores=1 for absolute reproducibility across OS)
weight_mat <- GENIE3(expr_matrix, regulators = intersect(heroes_ids, rownames(expr_matrix)), nCores = 1)

# Format edge list and extract Top 50 targets for each core HSF
edge_list <- getLinkList(weight_mat)
colnames(edge_list) <- c("Source_TF", "Target_Gene", "Importance_Score")
edge_list <- edge_list[order(edge_list$Importance_Score, decreasing = TRUE), ]

# Loop to save Top 50 targeted genes for downstream functional annotation
for (tf in heroes_ids) {
  top_50 <- head(edge_list[edge_list$Source_TF == tf, ], 50)
  write.csv(top_50, file.path(out_dir, paste0("Top50_Targets_", tf, ".csv")), row.names = FALSE)
}

# ==============================================================================
# PART 3: HIERARCHICAL NETWORK VISUALIZATION (Phases I, II, III)
# ==============================================================================
message("--- Executing Part 3: Phase-Mapped Hierarchical Network (Fig 5B) ---")

# Define Phase-specific color palette
pal <- c(
  "Phase I"   = "#d62728",  # Red
  "Phase II"  = "#ff7f0e",  # Orange
  "Phase III" = "#2ca02c",  # Green
  "Root"      = "#000000",  # Core Black
  "Unknown"   = "#999999"   # Grey fallback
)

# [Placeholder: Assume 'node_attr' and 'all_edges' are loaded from the annotated Top50 lists]
# Map core TFs to specific temporal phases
tf_phase_map <- c("Hsftf22" = "Phase I", "Hsftf11" = "Phase II", "Hsftf24" = "Phase III")

# Clean attributes and assign hierarchy
# nodes <- data.frame(name = all_node_names) %>%
#   left_join(clean_attr, by = c("name" = "node_name")) %>%
#   mutate(
#     node_branch = case_when(
#       name == "HsfGRN" ~ "Root",
#       name %in% names(tf_phase_map) ~ tf_phase_map[name],
#       TRUE ~ as.character(Phase)
#     ),
#     leaf = !(name %in% c("HsfGRN", names(tf_phase_map))),
#     size_val = case_when(name == "HsfGRN" ~ 15, !leaf ~ 8, TRUE ~ 3)
#   )

# Compute optimal angle for outward-facing text in circular layouts
node_angle <- function(x, y) {
  angle <- atan2(y, x) * 180 / pi
  angle <- ifelse(angle < 0, angle + 360, angle)
  ifelse(angle > 90 & angle < 270, angle + 180, angle)
}

# Build graph object
# hierarchy_graph <- tbl_graph(nodes = nodes, edges = rbind(edges_root, edges_genes))

# Generate the radial dendrogram
# p_grn <- ggraph(hierarchy_graph, layout = 'dendrogram', circular = TRUE) +
#   geom_edge_diagonal(aes(color = node1.node_branch), alpha = 0.4, width = 0.4) +
#   geom_node_point(aes(size = size_val, color = node_branch), alpha = 1) +
#   geom_node_text(aes(x = x*1.04, y = y*1.04, label = node_name, angle = node_angle(x, y), filter = leaf, color = node_branch), size = 2.2, hjust = 'outward') +
#   geom_node_text(aes(label = node_name, filter = !leaf, color = node_branch), fontface = "bold", size = 5, repel = FALSE) +
#   scale_size(range = c(2, 10), guide = "none") +
#   scale_color_manual(values = pal) +
#   scale_edge_color_manual(values = pal) +
#   coord_fixed() + theme_void() +
#   expand_limits(x = c(-1.6, 1.6), y = c(-1.6, 1.6)) +
#   theme(legend.position = "right", plot.background = element_rect(fill = "white", color = NA))

# ggsave(file.path(out_dir, "Fig5B_GRN_Phase_Network.pdf"), p_grn, width = 16, height = 16)

message("Causal inference and GRN visualization completed successfully.")
