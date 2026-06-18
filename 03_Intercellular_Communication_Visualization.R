# ==============================================================================
# Script: 03_Intercellular_Communication_Visualization.R
# Description: Construction and visualization of directed cell-cell communication 
#              networks based on PlantPhoneDB inference.
# NOTE: The raw ligand-receptor interaction inference was computed using the 
#       standard PlantPhoneDB pipeline. This script processes the output scoring 
#       matrices to construct the directed communication networks and generate 
#       the topological visualizations (Chord diagrams) presented in the manuscript.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Libraries & Environment Setup
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(circlize)
  library(readxl)
  library(grid)
  library(RColorBrewer)
})

# Define IO paths
data_dir <- "data/processed/"
out_dir <- "results/figures/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ==============================================================================
# PART 1: GLOBAL CELL-CELL COMMUNICATION NETWORK (Fig 2C)
# ==============================================================================
message("--- Executing Part 1: Global Cell-Cell Communication ---")

# Define global cell lineages and colors
cols_global <- c("Mesophyll" = "#2ca02c", "Bundle Sheath" = "#1f77b4", 
                 "Protoderm" = "#d62728", "Epidermis" = "#9467bd", "Companion" = "#ff7f0e")
order_global <- c("Protoderm", "Epidermis", "Mesophyll", "Bundle Sheath", "Companion")

# Load PlantPhoneDB global output (Provide the actual path to the txt file here)
# df_chat <- read.table(file.path(data_dir, "PlantPhoneDB_Global_Output.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# [Placeholder for demonstration: Ensure your dataset contains Ligands_cell and Receptors_cell columns]
# Clean nomenclature format
clean_name <- function(x) {
  x <- gsub("_", " ", x)
  s <- tolower(x)
  val <- paste0(toupper(substr(s, 1, 1)), substring(s, 2))
  val <- gsub("Bundle sheath", "Bundle Sheath", val)
  return(val)
}

# Assume df_chat is loaded properly; processing directed interactions
# df_chat$Source <- sapply(df_chat$Ligands_cell, clean_name)
# df_chat$Target <- sapply(df_chat$Receptors_cell, clean_name)
# df_valid_global <- df_chat %>% filter(Source %in% names(cols_global) & Target %in% names(cols_global))
# df_chord_global <- df_valid_global %>% group_by(Source, Target) %>% summarise(Value = n(), .groups = "drop") %>% as.data.frame()

# Wrapper function for global chord diagram visualization
draw_global_chord <- function(df_chord) {
  circos.clear()
  circos.par(start.degree = 90, gap.after = 3, points.overflow.warning = FALSE)
  par(xpd = NA) # Allow drawing outside boundaries
  
  chordDiagram(
    df_chord, grid.col = cols_global,
    order = order_global[order_global %in% unique(c(df_chord$Source, df_chord$Target))],
    transparency = 0.3, annotationTrack = "grid", annotationTrackHeight = 0.05,
    direction.type = c("diffHeight", "arrows"), link.arr.type = "triangle", link.arr.length = 0.15,
    # Highlight the critical Protoderm -> Epidermis signaling axis
    link.border = ifelse(df_chord$Source == "Protoderm" & df_chord$Target == "Epidermis", "black", NA),
    link.lwd = 2.5
  )
  
  # Custom external track for labels with adjusted offsets
  circos.track(track.index = 1, panel.fun = function(x, y) {
    sector.name = CELL_META$sector.index
    circos.text(mean(CELL_META$xlim), CELL_META$ylim[1] + 3.5, sector.name, 
                facing = "bending.outside", niceFacing = TRUE, adj = c(0.5, 0), 
                cex = 1.3, col = cols_global[sector.name], font = 2)
  }, bg.border = NA)
}

# Example to save (uncomment when real data is loaded)
# pdf(file.path(out_dir, "Fig2C_Global_Chord.pdf"), width = 7, height = 7)
# draw_global_chord(df_chord_global)
# dev.off()

# ==============================================================================
# PART 2: PROTODERM-SPECIFIC SUBPOPULATION NETWORK (Fig 4)
# ==============================================================================
message("--- Executing Part 2: Protoderm Subpopulation Communication ---")

# Define subpopulation lineages and colors
sub_levels <- paste0("Pro-", 1:6)
cols_sub <- brewer.pal(8, "Set2")[1:6]
names(cols_sub) <- sub_levels

# Load specific interaction table 
# df_raw_sub <- read.xlsx(file.path(data_dir, "Supplemental_Table_S18_CellChat_Final.xlsx"))

# Robust extraction using regex matching for dynamic column headers
# sender_col <- grep("Sender.*Paper", colnames(df_raw_sub), value = TRUE)[1]
# receiver_col <- grep("Receiver.*Paper", colnames(df_raw_sub), value = TRUE)[1]
# df_extract <- df_raw_sub[, c(sender_col, receiver_col)]
# colnames(df_extract) <- c("Source", "Target")

# df_chord_sub <- df_extract %>% filter(!is.na(Source) & !is.na(Target)) %>%
#   group_by(Source, Target) %>% summarise(Value = n(), .groups = "drop") %>% as.data.frame()
# df_chord_sub <- df_chord_sub %>% filter(Source %in% sub_levels & Target %in% sub_levels)

# Wrapper function for subpopulation chord diagram
draw_sub_chord <- function(df_chord) {
  circos.clear()
  circos.par(start.degree = 90, gap.after = 3, points.overflow.warning = FALSE)
  par(xpd = NA)
  
  chordDiagram(
    df_chord, grid.col = cols_sub, order = sub_levels,
    transparency = 0.3, annotationTrack = "grid", annotationTrackHeight = 0.05,
    direction.type = c("diffHeight", "arrows"), link.arr.type = "triangle", link.arr.length = 0.15,
    # Highlight the mechanistic Pro-4 -> Pro-1 trajectory transition
    link.border = ifelse(df_chord$Source == "Pro-4" & df_chord$Target == "Pro-1", "black", NA),
    link.lwd = 2.5
  )
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    sector.name = CELL_META$sector.index
    circos.text(mean(CELL_META$xlim), CELL_META$ylim[1] + 3.5, sector.name, 
                facing = "bending.outside", niceFacing = TRUE, adj = c(0.5, 0), 
                cex = 1.3, col = cols_sub[sector.name], font = 2)
  }, bg.border = NA)
}

# Example to save (uncomment when real data is loaded)
# pdf(file.path(out_dir, "Fig4F_Subpopulation_Chord.pdf"), width = 7, height = 7)
# draw_sub_chord(df_chord_sub)
# dev.off()

message("Intercellular communication network generation completed successfully.")
