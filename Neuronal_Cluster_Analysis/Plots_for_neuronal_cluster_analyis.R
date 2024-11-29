# Neuronal Subtype scRNA-Seq Analysis based on Seurat and Bioconductor 


# Load necessary libraries ----
library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(patchwork)
library(EnhancedVolcano)
library(purrr)
library(tidyr)
library(celldex)
library(SingleR)
library(tibble)
library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db)
library(openxlsx)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)

# I need to change the plots so they show the same as they showed for Astrocytes, also maybe we first look on at 0.1 and 35 pcs 
# Load data ----
# Update these paths with the exact names of your RDS files
# Aqp4_neurons_0_1_30 <- readRDS("your/path/to/Aqp4_neurons_0.1_30.rds")
neuronsData <- readRDS('your/path/to/sub_neuronal.rds')

save_neurons_path <- "your/path/to/saves"

# Add day and a batch column: 
# based on the orig.idents column in the meta.data we can know which day and which batch each cell is coming from 
# c2307 : days = 5d, batch = 2nd, condition  = Ctrl  
# c2592 : days = 5d, batch = 2nd, condition = Aqp4  
# c2594 : days = 5d, batch = 2nd, condition = Aqp4 
# c3434 : days = 3d, batch = 1st, condition = Ctrl  
# c4769 : days = 3d, batch = 1st, condition = Aqp4 
# c3349 : days = 5d, batch = 2nd, condtion = Ctrl
# c4774 : days = 3d, batch = 1st, condition = Ctrl  
# c3436 : days = 3d, batch = 1st, condition = Aqp4

# Map day information based on orig.ident column in meta.data
neuronsData@meta.data <- neuronsData@meta.data %>%
  mutate(day = case_when(
    orig.ident == "c2307" ~ "5d",
    orig.ident == "c2592" ~ "5d",
    orig.ident == "c2594" ~ "5d",
    orig.ident == "c3434" ~ "3d",
    orig.ident == "c4769" ~ "3d",
    orig.ident == "c3349" ~ "5d",
    orig.ident == "c4774" ~ "3d",
    orig.ident == "c3436" ~ "3d",
    TRUE ~ NA_character_  # Optional: fill with NA if no match
  ))

# Verify the addition of 'day' column
head(neuronsData@meta.data)

# Create a new identity combining condition and day
neuronsData@meta.data$condition_day <- paste(neuronsData@meta.data$condition, neuronsData@meta.data$day, sep = "_")

# Define neuronal subtype marker sets ----
Excitatory.markers <- c("Slc17a7", "Nrgn", "Camk2a", "Grin1", "Grin2b")
Inhibitory.markers <- c("Gad1", "Gad2", "Slc6a1", "Vgat", "Calb2")
Dopaminergic.markers <- c("Th", "Ddc", "Slc6a3", "Drd2", "Pitx3")
Serotonergic.markers <- c("Tph2", "Slc6a4", "Htr1a", "Htr2c", "Fev")
Interneuron.markers <- c("Sst", "Pvalb", "Vip", "Calb2")
Parvalbumin.markers <- c("Pvalb", "Gad1")
SST.markers <- c("Sst", "Nxph1")
VIP.markers <- c("Vip", "Npy")
Camk2.markers <- c("Camk2a", "Camk2b")
PyramidalNeurons.markers <- c("Slc17a7", "Camk2a", "Bcl11b", "Tbr1")

# Define the object name for saving plots
object_name <- "neuronsData"

# Plots -----
# Violin plot for QC metrics
VlnPlot(neuronsData, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(paste0(save_neurons_path, "VlnPlot_QC_", object_name, ".png"), width = 10, height = 8, dpi = 300)

# UMAP plot by clusters
DimPlot(neuronsData, reduction = "umap", label = TRUE) + ggtitle(paste("Neuronal Clusters -", object_name))
ggsave(paste0(save_neurons_path, "DimPlot_Clusters_", object_name, ".png"), width = 10, height = 8, dpi = 300)

# UMAP plot by condition
DimPlot(neuronsData, reduction = "umap", group.by = "condition") + ggtitle(paste("UMAP by Condition -", object_name))
ggsave(paste0(save_neurons_path, "UMAP_Condition_", object_name, ".png"), width = 10, height = 8, dpi = 300)

# UMAP plot by condition
DimPlot(neuronsData, reduction = "umap", group.by = "batch") + ggtitle(paste("UMAP by Batch -", object_name))
ggsave(paste0(save_neurons_path, "UMAP_Batch_", object_name, ".png"), width = 10, height = 8, dpi = 300)


# UMAP plot by day
DimPlot(neuronsData, reduction = "umap", group.by = "day") + ggtitle(paste("UMAP by Day -", object_name))
ggsave(paste0(save_neurons_path, "UMAP_Day_", object_name, ".png"), width = 10, height = 8, dpi = 300)

# FeaturePlots for Excitatory Markers
available_excitatory_markers <- intersect(Excitatory.markers, rownames(neuronsData))
if (length(available_excitatory_markers) > 0) {
  feature_plots <- lapply(available_excitatory_markers, function(marker) {
    FeaturePlot(neuronsData, features = marker) + ggtitle(paste(marker, "-", object_name))
  })
  
  combined_feature_plot <- wrap_plots(feature_plots, ncol = 3) + plot_annotation(title = paste("FeaturePlots for Excitatory Markers -", object_name))
  ggsave(paste0(save_neurons_path, "FeaturePlots_ExcitatoryMarkers_", object_name, ".png"), plot = combined_feature_plot, width = 10, height = 8, dpi = 300)
} else {
  message("No available excitatory markers found in the dataset for ", object_name)
}

# Step 1: Load a Reference Dataset
# For this example, we'll use the Mouse Cell Atlas (MCA) reference dataset
reference <- MouseRNAseqData()  # Load the reference dataset from celldex

# Step 2: Get the Expression Data from Seurat Object
neurons_counts <- GetAssayData(neuronsData, slot = "data")

# Step 3: Run SingleR for Neuronal Subtype Prediction
predictions <- SingleR(
  test = neurons_counts,
  ref = reference,
  labels = reference$label.main  # or `reference$label.fine` for more specific subtypes
)

# Step 4: Add Predictions to Seurat Metadata
neuronsData <- AddMetaData(neuronsData, predictions$labels, col.name = "predicted_subtype")

# Step 5: Visualization

# 1. UMAP Plot by Predicted Subtype
umap_plot <- DimPlot(neuronsData, group.by = "predicted_subtype", label = TRUE) + 
  ggtitle("Neuronal Subtype Predictions by SingleR")
print(umap_plot)
ggsave(filename = paste0(save_neurons_path, "UMAP_PredictedSubtype.png"), plot = umap_plot, width = 10, height = 8, dpi = 300)


# VolcanoPlots ---------

# Set the combined condition_day column as the active identity
Idents(neuronsData) <- neuronsData$condition_day

# Define marker sets if applicable (update with your specific marker lists if needed)
marker_sets <- list(
  "Excitatory Markers" = c("Slc17a7", "Nrgn", "Camk2a", "Grin1", "Grin2b"),
  "Inhibitory Markers" = c("Gad1", "Gad2", "Slc6a1", "Vgat", "Calb2"),
  "Dopaminergic Markers" = c("Th", "Ddc", "Slc6a3", "Drd2", "Pitx3"),
  "Serotonergic Markers" = c("Tph2", "Slc6a4", "Htr1a", "Htr2c", "Fev")
)

# Function to perform and filter comparisons
perform_and_filter_comparisons <- function(seurat_object, ident1, ident2, title) {
  results <- FindMarkers(seurat_object, ident.1 = ident1, ident.2 = ident2, logfc.threshold = 0, test.use = "wilcox")
  filtered_results <- results[abs(results$avg_log2FC) > 0.1, ]
  
  if (nrow(filtered_results) > 0) {
    create_general_volcano_plot(filtered_results, title)
    for (marker_set_name in names(marker_sets)) {
      marker_genes <- marker_sets[[marker_set_name]]
      create_volcano_plot_with_markers(filtered_results, title, marker_set_name, marker_genes)
    }
  } else {
    message("No significant results for comparison: ", title)
  }
}

# General volcano plot function
create_general_volcano_plot <- function(results, title) {
  labeled_genes <- rownames(results)[results$p_val_adj < 0.05 & abs(results$avg_log2FC) > 0.25]
  
  p <- EnhancedVolcano(
    results,
    lab = rownames(results),
    selectLab = labeled_genes,
    x = 'avg_log2FC',
    y = 'p_val_adj',
    xlim = c(-2.5, 2.5),   
    ylim = c(0, 30),  
    title = title,
    subtitle = "All Significant Markers Labeled",
    pCutoff = 0.05,
    FCcutoff = 0.25,  
    pointSize = 2.0,
    labSize = 4.0,
    colAlpha = 0.5, 
    col = c('grey30', 'forestgreen', 'royalblue', 'red2'),  
    legendPosition = 'top',
    legendLabSize = 12,
    legendIconSize = 4.0,
    gridlines.major = TRUE,  
    gridlines.minor = TRUE,  
    border = 'full'  
  )
  
  ggsave(paste0(save_neurons_path, title, "_General.png"), plot = p, width = 8, height = 6)
}

# Volcano plot function for specific marker sets
create_volcano_plot_with_markers <- function(results, title, marker_set_name, marker_genes) {
  labeled_genes <- intersect(rownames(results), marker_genes)
  
  if (length(labeled_genes) > 0) {
    p <- EnhancedVolcano(
      results,
      lab = rownames(results),
      selectLab = labeled_genes,
      x = 'avg_log2FC',
      y = 'p_val_adj',
      xlim = c(-2.5, 2.5),   
      ylim = c(0, 30),  
      title = paste(title, "-", marker_set_name),
      subtitle = paste("Labeled Markers from", marker_set_name),
      pCutoff = 0.05,
      FCcutoff = 0.25,  
      pointSize = 2.0,
      labSize = 4.0,
      colAlpha = 0.5, 
      col = c('grey30', 'forestgreen', 'royalblue', 'red2'),  
      legendPosition = 'top',
      legendLabSize = 12,
      legendIconSize = 4.0,
      gridlines.major = TRUE,  
      gridlines.minor = TRUE,  
      border = 'full'  
    )
    
    ggsave(paste0(save_neurons_path, title, "_", marker_set_name, ".png"), plot = p, width = 8, height = 6)
  } else {
    message(paste("No markers from", marker_set_name, "found in the results for", title))
  }
}

# Run comparisons
perform_and_filter_comparisons(neuronsData, "Aqp4_5d", "Aqp4_3d", "Aqp4 5D vs Aqp4 3D")
perform_and_filter_comparisons(neuronsData, "Aqp4_3d", "Ctrl_3d", "Aqp4 3D vs Ctrl 3D")
perform_and_filter_comparisons(neuronsData, "Aqp4_5d", "Ctrl_5d", "Aqp4 5D vs Ctrl 5D")

# Perform and filter comparisons for upregulated and downregulated genes
perform_and_filter_comparisons <- function(seurat_object, ident1, ident2, title) {
  # Perform differential expression analysis
  results <- FindMarkers(seurat_object, ident.1 = ident1, ident.2 = ident2, logfc.threshold = 0, test.use = "wilcox")
  
  # Filter for significantly upregulated and downregulated genes
  significant_results <- results[results$p_val_adj < 0.05, ]
  upregulated <- significant_results[significant_results$avg_log2FC > 0.1, ]
  downregulated <- significant_results[significant_results$avg_log2FC < -0.1, ]
  
  # Save results for overlap analysis
  list(upregulated = rownames(upregulated), downregulated = rownames(downregulated))
}

# Find overlaps in upregulated and downregulated genes between two comparisons
find_gene_overlaps <- function(up1, up2, down1, down2) {
  # Find overlaps
  overlapping_up <- intersect(up1, up2)
  overlapping_down <- intersect(down1, down2)
  
  # Report counts and overlaps
  message("Number of overlapping upregulated genes: ", length(overlapping_up))
  message("Overlapping upregulated genes: ", paste(overlapping_up, collapse = ", "))
  
  message("Number of overlapping downregulated genes: ", length(overlapping_down))
  message("Overlapping downregulated genes: ", paste(overlapping_down, collapse = ", "))
  
  # Return overlaps
  list(overlapping_up = overlapping_up, overlapping_down = overlapping_down)
}

# Perform comparisons for AQP4 3 days vs Control 3 days and AQP4 5 days vs Control 5 days
comparison_3d <- perform_and_filter_comparisons(neuronsData, "Aqp4_3d", "Ctrl_3d", "Aqp4 3D vs Ctrl 3D")
comparison_5d <- perform_and_filter_comparisons(neuronsData, "Aqp4_5d", "Ctrl_5d", "Aqp4 5D vs Ctrl 5D")

# Find and report overlaps
overlap_results <- find_gene_overlaps(
  comparison_3d$upregulated,
  comparison_5d$upregulated,
  comparison_3d$downregulated,
  comparison_5d$downregulated
)

# Save overlaps to a file if needed
saveRDS(overlap_results, file = paste0(save_neurons_path, "Aqp4_Overlap_3D_5D.rds"))



# Top markers -----

# Ensure clustering information is set as active identity
Idents(neuronsData) <- "seurat_clusters"

# Initialize an empty data frame to store top 10 markers for each cluster
top_markers_df <- data.frame()

# Loop over clusters 0 to 5
for (i in 0:11) {
  # Set ident.2 to all other clusters except the current cluster (i)
  other_clusters <- setdiff(0:11, i)
  
  # Find markers distinguishing the current cluster from other clusters
  cluster_markers <- FindMarkers(neuronsData, ident.1 = i, ident.2 = other_clusters, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE)
  
  # Convert row names to a column named 'gene'
  cluster_markers <- rownames_to_column(cluster_markers, var = "gene")
  
  # Select the top 10 markers
  top10_markers <- head(cluster_markers, n = 10)
  
  # Add a column to indicate the cluster number
  top10_markers$cluster <- i
  
  # Append the top 10 markers to the data frame
  top_markers_df <- bind_rows(top_markers_df, top10_markers)
}

# Print the combined data frame to verify the results
print(top_markers_df)
View(top_markers_df)

# Extract the unique gene names (features) from the combined data frame
top_markers_genes <- unique(top_markers_df$gene)
View(top_markers_genes)

# Generate and save a heatmap for the top markers
topHeatmap <- DoHeatmap(neuronsData, features = top_markers_genes, size = 3)
ggsave(paste0(save_neurons_path, "DoHeatmap_of_top_markers.png"), plot = topHeatmap, width = 10, height = 10)

# Dot plots ----------
# Define neuronal marker sets
Excitatory.markers <- c("Slc17a7", "Nrgn", "Camk2a", "Grin1", "Grin2b")
Inhibitory.markers <- c("Gad1", "Gad2", "Slc6a1", "Vgat", "Calb2")
Dopaminergic.markers <- c("Th", "Ddc", "Slc6a3", "Drd2", "Pitx3")
Serotonergic.markers <- c("Tph2", "Slc6a4", "Htr1a", "Htr2c", "Fev")
Interneuron.markers <- c("Sst", "Pvalb", "Vip", "Calb2")
Parvalbumin.markers <- c("Pvalb", "Gad1")
SST.markers <- c("Sst", "Nxph1")
VIP.markers <- c("Vip", "Npy")
Camk2.markers <- c("Camk2a", "Camk2b")
PyramidalNeurons.markers <- c("Slc17a7", "Camk2a", "Bcl11b", "Tbr1")

# Create a named list of marker sets
marker_sets <- list(
  "Excitatory Markers" = Excitatory.markers,
  "Inhibitory Markers" = Inhibitory.markers,
  "Dopaminergic Markers" = Dopaminergic.markers,
  "Serotonergic Markers" = Serotonergic.markers,
  "Interneuron Markers" = Interneuron.markers,
  "Parvalbumin Markers" = Parvalbumin.markers,
  "SST Markers" = SST.markers,
  "VIP Markers" = VIP.markers,
  "Camk2 Markers" = Camk2.markers,
  "Pyramidal Neurons Markers" = PyramidalNeurons.markers
)

# Define a function to create and save DotPlots with a white background
save_dotplot <- function(data, features, title, filename, group_by = "seurat_clusters") {
  DefaultAssay(data) <- "RNA"  # Set the assay to RNA explicitly
  p <- DotPlot(data, features = features, group.by = group_by) + 
    RotatedAxis() + 
    ggtitle(title) + 
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white")
    )
  ggsave(filename,save_neurons_path, plot = p, width = 20, height = 10)
}

# Set the correct grouping column here (replace "seurat_clusters" with your desired grouping)
grouping_column <- "seurat_clusters"

# Loop over marker sets and save DotPlots
for (marker_set_name in names(marker_sets)) {
  marker_genes <- marker_sets[[marker_set_name]]
  filename <- paste0("DotPlot_", gsub(" ", "_", marker_set_name), ".png")
  save_dotplot(neuronsData, marker_genes, marker_set_name, filename, group_by = grouping_column)
}

# FeaturePlots and Violinplots ----------
# Define neuronal marker sets
Excitatory.markers <- c("Slc17a7", "Nrgn", "Camk2a", "Grin1", "Grin2b")
Inhibitory.markers <- c("Gad1", "Gad2", "Slc6a1", "Vgat", "Calb2")
Dopaminergic.markers <- c("Th", "Ddc", "Slc6a3", "Drd2", "Pitx3")
Serotonergic.markers <- c("Tph2", "Slc6a4", "Htr1a", "Htr2c", "Fev")
Interneuron.markers <- c("Sst", "Pvalb", "Vip", "Calb2")
Parvalbumin.markers <- c("Pvalb", "Gad1")
SST.markers <- c("Sst", "Nxph1")
VIP.markers <- c("Vip", "Npy")
Camk2.markers <- c("Camk2a", "Camk2b")
PyramidalNeurons.markers <- c("Slc17a7", "Camk2a", "Bcl11b", "Tbr1")

# Create a named list of marker sets
marker_sets <- list(
  "Excitatory Markers" = Excitatory.markers,
  "Inhibitory Markers" = Inhibitory.markers,
  "Dopaminergic Markers" = Dopaminergic.markers,
  "Serotonergic Markers" = Serotonergic.markers,
  "Interneuron Markers" = Interneuron.markers,
  "Parvalbumin Markers" = Parvalbumin.markers,
  "SST Markers" = SST.markers,
  "VIP Markers" = VIP.markers,
  "Camk2 Markers" = Camk2.markers,
  "Pyramidal Neurons Markers" = PyramidalNeurons.markers
)

# Set the default assay to RNA
DefaultAssay(neuronsData) <- "RNA"

# Define a function to create and save DotPlots with a white background
save_dotplot <- function(data, features, title, filename, group_by = "seurat_clusters") {
  p <- DotPlot(data, features = features, group.by = group_by) + 
    RotatedAxis() + 
    ggtitle(title) + 
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white")
    )
  ggsave(paste0(save_neurons_path, filename), plot = p, width = 20, height = 10)
}

# Define a function to create and save DotPlots with condition_day as x-axis
save_dotplot_condition_day <- function(data, features, title, filename) {
  p <- DotPlot(data, features = features, group.by = "condition_day") + 
    RotatedAxis() + 
    ggtitle(title) + 
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white")
    )
  ggsave(paste0(save_neurons_path, filename), plot = p, width = 20, height = 10)
}

# Loop over marker sets and save DotPlots with condition_day on x-axis
for (marker_set_name in names(marker_sets)) {
  marker_genes <- marker_sets[[marker_set_name]]
  dotplot_filename <- paste0("DotPlot_ConditionDay_", gsub(" ", "_", marker_set_name), ".png")
  save_dotplot_condition_day(neuronsData, marker_genes, marker_set_name, dotplot_filename)
}

# Loop over marker sets and save DotPlots
for (marker_set_name in names(marker_sets)) {
  marker_genes <- marker_sets[[marker_set_name]]
  dotplot_filename <- paste0("DotPlot_", gsub(" ", "_", marker_set_name), ".png")
  save_dotplot(neuronsData, marker_genes, marker_set_name, dotplot_filename, group_by = grouping_column)
}

# Define a function to create and save FeaturePlots with missing gene check
save_featureplot <- function(data, features, title, filename) {
  # Filter features to only those present in the data
  available_features <- intersect(features, rownames(data))
  
  # Check if there are any features to plot
  if (length(available_features) == 0) {
    message("No available features for ", title)
    return(NULL)
  }
  
  # Create FeaturePlots for available features only
  feature_plots <- lapply(available_features, function(feature) {
    FeaturePlot(data, features = feature, slot = "data") + ggtitle(paste(title, "-", feature))
  })
  
  # Combine plots
  combined_plot <- patchwork::wrap_plots(feature_plots, ncol = 3)
  
  # Save the combined plot
  ggsave(paste0(save_neurons_path, filename), plot = combined_plot, width = 20, height = 15)
}

# Loop over marker sets and save FeaturePlots
for (marker_set_name in names(marker_sets)) {
  marker_genes <- marker_sets[[marker_set_name]]
  featureplot_filename <- paste0("FeaturePlot_", gsub(" ", "_", marker_set_name), ".png")
  save_featureplot(neuronsData, marker_genes, marker_set_name, featureplot_filename)
}

# Define a function to create and save ViolinPlots, checking for available features
save_violinplot <- function(data, features, title, filename) {
  # Check if features exist in the RNA assay
  available_features <- intersect(features, rownames(data@assays$RNA))
  
  # If no features are available, exit the function with a message
  if (length(available_features) == 0) {
    message("No available features for ", title)
    return(NULL)
  }
  
  # Create ViolinPlot using the specified layer (e.g., "data")
  violin_plot <- VlnPlot(data, features = available_features, group.by = "seurat_clusters", ncol = 3, slot = "data") +
    ggtitle(title)
  
  # Save the ViolinPlot
  ggsave(paste0(save_neurons_path, filename), plot = violin_plot, width = 20, height = 10)
}

# Loop over marker sets and save ViolinPlots
for (marker_set_name in names(marker_sets)) {
  marker_genes <- marker_sets[[marker_set_name]]
  violinplot_filename <- paste0("ViolinPlot_", gsub(" ", "_", marker_set_name), ".png")
  save_violinplot(neuronsData, marker_genes, marker_set_name, violinplot_filename)
}

# Pheatmap to look how similar each cluster to each other is ----

# Define the save path
save_path <- "/your/path/to/save/"

# Check if the directory exists; if not, create it
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

# Define the list of marker sets
marker_sets <- list(
  Excitatory = c("Slc17a7", "Nrgn", "Camk2a", "Grin1", "Grin2b"), #glutamatergic
  Inhibitory = c("Gad1", "Gad2", "Slc6a1", "Calb2"), # gabaergic 
  Serotonergic = c("Tph2", "Slc6a4", "Htr1a", "Htr2c", "Fev"),
  Interneuron = c("Sst", "Pvalb", "Vip", "Calb2"),
  Parvalbumin = c("Pvalb", "Gad1"),
  SST = c("Sst", "Nxph1"),
  VIP = c("Vip", "Npy"),
  Camk2 = c("Camk2a", "Camk2b"),
  Pyramidal = c("Slc17a7", "Camk2a", "Bcl11b", "Tbr1")
)

# Loop through each marker set to create heatmaps
for (marker_set_name in names(marker_sets)) {
  # Extract the markers for the current set
  markers <- marker_sets[[marker_set_name]]
  
  # Calculate average expression for the marker set in the specified assay and slot
  avg_expression_candidates <- AverageExpression(
    neuronsData,  # Your Seurat object
    features = markers,
    assay = "RNA",
    slot = "data.1"  # Use the appropriate layer
  )
  
  # Extract the average expression matrix
  avg_expression_matrix <- avg_expression_candidates$RNA  # Access the RNA layer
  
  # Define the color gradient
  color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  
  # Save the heatmap to a file
  png(filename = paste0(save_path, "Heatmap_Average_", marker_set_name, ".png"), width = 10, height = 8, units = "in", res = 300)
  
  # Create the heatmap
  pheatmap(
    avg_expression_matrix,              # Use the average expression matrix
    cluster_rows = TRUE,                # Cluster the markers (rows)
    cluster_cols = TRUE,                # Cluster the clusters (columns)
    show_rownames = TRUE,               # Show marker names on the y-axis
    show_colnames = TRUE,               # Show cluster names on the x-axis
    scale = "row",                      # Scale each marker's expression
    color = color_palette,              # Apply the color palette
    main = paste("Average Expression of", marker_set_name, "Markers by Cluster")  # Add title
  )
  
  dev.off()  # Close the graphical device to save the plot
}


# KEGG per cluster ----
# Load required libraries

# Define the save path
save_path <- "/your/path/to/save"

# Check if the directory exists; if not, create it
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

# Create condition_day column in metadata and set as identities
neuronsData@meta.data$condition_day <- paste(neuronsData@meta.data$condition, neuronsData@meta.data$day, sep = "_")
Idents(neuronsData) <- neuronsData@meta.data$condition_day

# Print unique identities to confirm correct setup
print("Unique values in condition_day identities:")
print(unique(Idents(neuronsData)))

# Initialize results table
results_table <- data.frame(
  Cluster = character(), 
  Comparison = character(), 
  PathwayType = character(), 
  Pathway = character(), 
  GeneCount = integer(), 
  P.adjust = numeric(), 
  GeneNames = character(),  # Add a column for gene names
  stringsAsFactors = FALSE
)

# Function to perform GO enrichment analysis and plot with color per ontology type
perform_go_analysis <- function(go_result, cluster_number, comparison_name, type, ontology) {
  if (!is.null(go_result) && nrow(go_result) > 0) {
    # Select top 20 GO terms by p.adjust
    top_go <- go_result %>%
      as.data.frame() %>%
      arrange(p.adjust) %>%
      head(min(20, nrow(go_result)))
    
    # Choose color scheme
    color_palette <- if (ontology == "BP") "purple" else if (ontology == "MF") "blue" else "green"
    
    # Create plot
    go_plot <- ggplot(top_go, aes(x = reorder(Description, Count), y = Count, fill = -log10(p.adjust))) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_gradient(low = color_palette, high = "red") +
      labs(
        x = ontology,
        y = "Gene Count",
        title = paste("GO", ontology, "Analysis (", type, " DEGs) for Cluster", cluster_number, "(", comparison_name, ")")
      ) +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 12)
      )
    
    # Save plot
    ggsave(filename = paste0(save_path, "GO_", ontology, "_", type, "_Cluster_", cluster_number, "_", comparison_name, ".png"), plot = go_plot, width = 10, height = 8)
    
    # Append results to the results table
    results_table <<- rbind(results_table, data.frame(
      Cluster = cluster_number, 
      Comparison = comparison_name, 
      PathwayType = paste("GO", type, ontology), 
      Pathway = top_go$Description, 
      GeneCount = top_go$Count, 
      P.adjust = top_go$p.adjust,
      GeneNames = paste(rownames(top_go), collapse = ", ")  # Add gene names
    ))
  } else {
    cat("No significant GO pathways found for", type, ontology, "genes in Cluster", cluster_number, "(", comparison_name, ")\n")
  }
}

# Function to perform KEGG enrichment analysis and plot
perform_kegg_analysis <- function(kegg_result, cluster_number, comparison_name, type) {
  if (!is.null(kegg_result) && nrow(kegg_result) > 0) {
    # Clean up the pathway descriptions
    kegg_result@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", kegg_result@result$Description)
    
    # Create plot
    kegg_plot <- dotplot(kegg_result, showCategory = min(20, nrow(kegg_result))) +
      ggtitle(paste("KEGG Pathway Analysis (", type, " DEGs) for Cluster", cluster_number, "(", comparison_name, ")")) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 10))
    
    # Save plot
    ggsave(filename = paste0(save_path, "KEGG_", type, "_Cluster_", cluster_number, "_", comparison_name, ".png"), plot = kegg_plot, width = 10, height = 8)
    
    # Append results to the results table
    results_table <<- rbind(results_table, data.frame(
      Cluster = cluster_number, 
      Comparison = comparison_name, 
      PathwayType = paste("KEGG", type), 
      Pathway = kegg_result@result$Description, 
      GeneCount = kegg_result@result$Count, 
      P.adjust = kegg_result@result$p.adjust,
      GeneNames = paste(kegg_result@gene, collapse = ", ")  # Add gene names
    ))
  } else {
    cat("No significant KEGG pathways found for", type, "genes in Cluster", cluster_number, "(", comparison_name, ")\n")
  }
}

# Function to perform pathway enrichment for a specific cluster and condition comparison
perform_deg_pathway_analysis_for_cluster_condition <- function(data, cluster_number, comparison_name, condition1, condition2) {
  # Subset data for the specific cluster
  cluster_cells <- subset(data, seurat_clusters == cluster_number)
  
  # Set identities to condition_day within the cluster subset
  Idents(cluster_cells) <- cluster_cells$condition_day
  
  # Print identities in the subset for debugging
  print(paste("Processing Cluster:", cluster_number, "Comparison:", comparison_name))
  print("Identities in cluster_cells after setting to 'condition_day':")
  print(table(Idents(cluster_cells)))
  
  # Subset data to only the two conditions being compared
  cluster_cells <- subset(cluster_cells, idents = c(condition1, condition2))
  
  # Check if both conditions are present and have sufficient cells
  if (sum(Idents(cluster_cells) == condition1) < 10 | sum(Idents(cluster_cells) == condition2) < 10) {
    cat("Skipping cluster", cluster_number, "(", comparison_name, "): One or both conditions are not present or do not have enough cells\n")
    results_table <<- rbind(results_table, data.frame(
      Cluster = cluster_number, 
      Comparison = comparison_name, 
      PathwayType = "None", 
      Pathway = "Comparison skipped due to insufficient cells", 
      GeneCount = NA, 
      P.adjust = NA,
      GeneNames = NA  # No genes to report
    ))
    return(NULL)
  }
  
  # Perform differential expression analysis
  de_results <- FindMarkers(cluster_cells, ident.1 = condition1, ident.2 = condition2, logfc.threshold = 0.25, test.use = "wilcox")
  
  # Separate upregulated and downregulated genes
  upregulated_genes <- rownames(de_results[de_results$avg_log2FC > 0 & de_results$p_val_adj < 0.05, ])
  downregulated_genes <- rownames(de_results[de_results$avg_log2FC < 0 & de_results$p_val_adj < 0.05, ])
  
  # GO analysis for upregulated genes
  if (length(upregulated_genes) > 0) {
    upregulated_entrez_ids <- mapIds(org.Mm.eg.db, keys = upregulated_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    upregulated_entrez_ids <- upregulated_entrez_ids[!is.na(upregulated_entrez_ids)]
    if (length(upregulated_entrez_ids) > 0) {
      # Run GO analysis for BP, MF, CC ontologies
      go_result_bp <- enrichGO(gene = upregulated_entrez_ids, OrgDb = org.Mm.eg.db, ont = "BP", pvalueCutoff = 0.05)
      go_result_mf <- enrichGO(gene = upregulated_entrez_ids, OrgDb = org.Mm.eg.db, ont = "MF", pvalueCutoff = 0.05)
      go_result_cc <- enrichGO(gene = upregulated_entrez_ids, OrgDb = org.Mm.eg.db, ont = "CC", pvalueCutoff = 0.05)
      perform_go_analysis(go_result_bp, cluster_number, comparison_name, "Upregulated", "BP")
      perform_go_analysis(go_result_mf, cluster_number, comparison_name, "Upregulated", "MF")
      perform_go_analysis(go_result_cc, cluster_number, comparison_name, "Upregulated", "CC")
    }
  }
  
  # GO analysis for downregulated genes
  if (length(downregulated_genes) > 0) {
    downregulated_entrez_ids <- mapIds(org.Mm.eg.db, keys = downregulated_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    downregulated_entrez_ids <- downregulated_entrez_ids[!is.na(downregulated_entrez_ids)]
    if (length(downregulated_entrez_ids) > 0) {
      # Run GO analysis for BP, MF, CC ontologies
      go_result_bp <- enrichGO(gene = downregulated_entrez_ids, OrgDb = org.Mm.eg.db, ont = "BP", pvalueCutoff = 0.05)
      go_result_mf <- enrichGO(gene = downregulated_entrez_ids, OrgDb = org.Mm.eg.db, ont = "MF", pvalueCutoff = 0.05)
      go_result_cc <- enrichGO(gene = downregulated_entrez_ids, OrgDb = org.Mm.eg.db, ont = "CC", pvalueCutoff = 0.05);
      perform_go_analysis(go_result_bp, cluster_number, comparison_name, "Downregulated", "BP");
      perform_go_analysis(go_result_mf, cluster_number, comparison_name, "Downregulated", "MF");
      perform_go_analysis(go_result_cc, cluster_number, comparison_name, "Downregulated", "CC");
    }
  }
  
  # KEGG analysis for upregulated genes
  if (length(upregulated_genes) > 0) {
    upregulated_entrez_ids <- mapIds(org.Mm.eg.db, keys = upregulated_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    upregulated_entrez_ids <- upregulated_entrez_ids[!is.na(upregulated_entrez_ids)]
    if (length(upregulated_entrez_ids) > 0) {
      # Run KEGG analysis
      kegg_result_up <- enrichKEGG(gene = upregulated_entrez_ids, organism = 'mmu', pvalueCutoff = 0.05)
      perform_kegg_analysis(kegg_result_up, cluster_number, comparison_name, "Upregulated")
    }
  }
  
  # KEGG analysis for downregulated genes
  if (length(downregulated_genes) > 0) {
    downregulated_entrez_ids <- mapIds(org.Mm.eg.db, keys = downregulated_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    downregulated_entrez_ids <- downregulated_entrez_ids[!is.na(downregulated_entrez_ids)]
    if (length(downregulated_entrez_ids) > 0) {
      # Run KEGG analysis
      kegg_result_down <- enrichKEGG(gene = downregulated_entrez_ids, organism = 'mmu', pvalueCutoff = 0.05)
      perform_kegg_analysis(kegg_result_down, cluster_number, comparison_name, "Downregulated")
    }
  }
  
  # Add results to table if no DEGs were found
  if (length(upregulated_genes) == 0 && length(downregulated_genes) == 0) {
    results_table <<- rbind(results_table, data.frame(Cluster = cluster_number, Comparison = comparison_name, PathwayType = "None", Pathway = "No DEGs found", GeneCount = NA, P.adjust = NA, GeneNames = NA))
  }
}

# Define the comparisons for each condition pair with correct case
comparisons <- list(
  "Aqp4_5d_vs_Ctrl_5d" = c("Aqp4_5d", "Ctrl_5d"),
  "Aqp4_3d_vs_Ctrl_3d" = c("Aqp4_3d", "Ctrl_3d")
)

# Adjust cluster range based on actual clusters
unique_clusters <- unique(neuronsData$seurat_clusters)

# Run analysis across clusters and comparisons
for (cluster in unique_clusters) {
  for (comparison_name in names(comparisons)) {
    conditions <- comparisons[[comparison_name]]
    perform_deg_pathway_analysis_for_cluster_condition(neuronsData, cluster, comparison_name, conditions[1], conditions[2])
  }
}

# Save results table to Excel
write.xlsx(results_table, paste0(save_path, "Pathway_Analysis_Results_Neuronal_Comparisons.xlsx"), rowNames = FALSE)

cat("Cluster-based pathway analysis for condition comparisons is complete. The results are saved in", save_path)

# Annotation: ----
# based on expression of markers sets 

# Define the cluster groups
glutamatergic_clusters <- c(0, 1, 2, 4, 6, 10, 11)
gabaergic_clusters <- c(5, 8)
other_clusters <- c(3, 7, 9)

# Create a new metadata column for cluster groups
neuronsData$ClusterGroup <- "Other"  # Default group
neuronsData$ClusterGroup[neuronsData$seurat_clusters %in% glutamatergic_clusters] <- "Glutamatergic"
neuronsData$ClusterGroup[neuronsData$seurat_clusters %in% gabaergic_clusters] <- "Gabaergic"

# Plot UMAP with improved aesthetics
umap_plot <- DimPlot(neuronsData, reduction = "umap", group.by = "ClusterGroup", 
                     cols = c("Glutamatergic" = "#E69F00",  # Publication-friendly orange
                              "Gabaergic" = "#56B4E9",      # Publication-friendly blue
                              "Other" = "#999999"),         # Neutral gray
                     label = TRUE, label.size = 4) +
  labs(title = "UMAP of Neuronal Clusters",
       subtitle = "Colored by Cluster Group",
       x = "UMAP 1",  # Add axis labels
       y = "UMAP 2") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center and bold title
    plot.subtitle = element_text(hjust = 0.5, size = 14),  # Center subtitle
    axis.text = element_text(size = 12, color = "black"),  # Adjust axis text size and color
    axis.title = element_text(size = 14, color = "black", face = "bold"),  # Bold axis titles
    axis.line = element_line(color = "black"),  # Add black axis lines
    legend.title = element_text(size = 20, face = "bold"),  # Increase legend title size
    legend.text = element_text(size = 12),  # Increase legend text size
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  )

# Save the UMAP plot
ggsave(filename = paste0(save_neurons_path, "UMAP_Neuronal_Clusters_Publication.png"), 
       plot = umap_plot, width = 10, height = 8, dpi = 300)

# Print the UMAP plot
print(umap_plot)


# GO and KEGG analysis of neuronal subtypes ------
# Load required libraries
# Define the cluster groups
glutamatergic_clusters <- c(0, 1, 2, 4, 6, 10, 11)
gabaergic_clusters <- c(5, 8)
other_clusters <- c(3, 7, 9)

# Add the neuronal_subtype column to metadata
neuronsData$neuronal_subtype <- ifelse(
  neuronsData$seurat_clusters %in% glutamatergic_clusters, "Glutamatergic",
  ifelse(neuronsData$seurat_clusters %in% gabaergic_clusters, "Gabaergic", "Other")
)

# Set identities to the new neuronal_subtype column
Idents(neuronsData) <- neuronsData$neuronal_subtype

# Verify the new column and identities
table(neuronsData$neuronal_subtype)
print("Identities set to neuronal_subtype:")
print(unique(Idents(neuronsData)))

# Load required libraries
library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)
library(openxlsx)

# Define the save path
save_path <- save_neurons_path

# Check if the directory exists; if not, create it
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

# Create condition_day column in metadata and set as identities
neuronsData@meta.data$condition_day <- paste(neuronsData@meta.data$condition, neuronsData@meta.data$day, sep = "_")
Idents(neuronsData) <- neuronsData@meta.data$condition_day

# Print unique identities to confirm correct setup
print("Unique values in condition_day identities:")
print(unique(Idents(neuronsData)))

# Initialize results table
results_table <- data.frame(
  NeuronalSubtype = character(), 
  Comparison = character(), 
  PathwayType = character(), 
  Pathway = character(), 
  GeneCount = integer(), 
  P.adjust = numeric(), 
  stringsAsFactors = FALSE
)

# Function to perform GO enrichment analysis and plot with color per ontology type
perform_go_analysis <- function(go_result, neuronal_subtype, comparison_name, type, ontology) {
  if (!is.null(go_result) && nrow(go_result) > 0) {
    # Select top 20 GO terms by p.adjust
    top_go <- go_result %>%
      as.data.frame() %>%
      arrange(p.adjust) %>%
      head(min(20, nrow(go_result)))
    
    # Choose color scheme
    color_palette <- if (ontology == "BP") "purple" else if (ontology == "MF") "blue" else "green"
    
    # Create plot
    go_plot <- ggplot(top_go, aes(x = reorder(Description, Count), y = Count, fill = -log10(p.adjust))) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_gradient(low = color_palette, high = "red") +
      labs(
        x = ontology,
        y = "Gene Count",
        title = paste("GO", ontology, "Analysis (", type, " DEGs) for", neuronal_subtype, "(", comparison_name, ")")
      ) +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 12)
      )
    
    # Save plot
    ggsave(filename = paste0(save_path, "GO_", ontology, "_", type, "_", neuronal_subtype, "_", comparison_name, ".png"), plot = go_plot, width = 10, height = 8)
    
    # Append results to the results table
    results_table <<- rbind(results_table, data.frame(
      NeuronalSubtype = neuronal_subtype, 
      Comparison = comparison_name, 
      PathwayType = paste("GO", type, ontology), 
      Pathway = top_go$Description, 
      GeneCount = top_go$Count, 
      P.adjust = top_go$p.adjust
    ))
  } else {
    cat("No significant GO pathways found for", type, ontology, "genes in", neuronal_subtype, "(", comparison_name, ")\n")
  }
}

# Function to perform KEGG enrichment analysis and plot
perform_kegg_analysis <- function(kegg_result, neuronal_subtype, comparison_name, type) {
  if (!is.null(kegg_result) && nrow(kegg_result) > 0) {
    # Clean up the pathway descriptions
    kegg_result@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", kegg_result@result$Description)
    
    # Create plot
    kegg_plot <- dotplot(kegg_result, showCategory = min(20, nrow(kegg_result))) +
      ggtitle(paste("KEGG Pathway Analysis (", type, " DEGs) for", neuronal_subtype, "(", comparison_name, ")")) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 10))
    
    # Save plot
    ggsave(filename = paste0(save_path, "KEGG_", type, "_", neuronal_subtype, "_", comparison_name, ".png"), plot = kegg_plot, width = 10, height = 8)
    
    # Append results to the results table
    results_table <<- rbind(results_table, data.frame(
      NeuronalSubtype = neuronal_subtype,
      Comparison = comparison_name,
      PathwayType = paste("KEGG", type),
      Pathway = kegg_result@result$Description,
      GeneCount = kegg_result@result$Count,
      P.adjust = kegg_result@result$p.adjust
    ))
  } else {
    cat("No significant KEGG pathways found for", type, "genes in", neuronal_subtype, "(", comparison_name, ")\n")
  }
}

# Function to perform pathway enrichment for a specific neuronal subtype and condition comparison
perform_deg_pathway_analysis_for_subtype_condition <- function(data, neuronal_subtype, comparison_name, condition1, condition2) {
  # Subset data for the specific neuronal subtype
  subtype_cells <- subset(data, neuronal_subtype == neuronal_subtype)
  
  # Set identities to condition_day within the subtype subset
  Idents(subtype_cells) <- subtype_cells$condition_day
  
  # Subset data to only the two conditions being compared
  subtype_cells <- subset(subtype_cells, idents = c(condition1, condition2))
  
  # Check if both conditions are present and have sufficient cells
  if (sum(Idents(subtype_cells) == condition1) < 10 | sum(Idents(subtype_cells) == condition2) < 10) {
    cat("Skipping", neuronal_subtype, "(", comparison_name, "): One or both conditions are not present or do not have enough cells\n")
    results_table <<- rbind(results_table, data.frame(
      NeuronalSubtype = neuronal_subtype, 
      Comparison = comparison_name, 
      PathwayType = "None", 
      Pathway = "Comparison skipped due to insufficient cells", 
      GeneCount = NA, 
      P.adjust = NA
    ))
    return(NULL)
  }
  
  # Perform differential expression analysis
  de_results <- FindMarkers(subtype_cells, ident.1 = condition1, ident.2 = condition2, logfc.threshold = 0.25, test.use = "wilcox")
  
  # Separate upregulated and downregulated genes
  upregulated_genes <- rownames(de_results[de_results$avg_log2FC > 0 & de_results$p_val_adj < 0.05, ])
  downregulated_genes <- rownames(de_results[de_results$avg_log2FC < 0 & de_results$p_val_adj < 0.05, ])
  
  # GO analysis for upregulated genes
  if (length(upregulated_genes) > 0) {
    upregulated_entrez_ids <- mapIds(org.Mm.eg.db, keys = upregulated_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    upregulated_entrez_ids <- upregulated_entrez_ids[!is.na(upregulated_entrez_ids)]
    if (length(upregulated_entrez_ids) > 0) {
      # Run GO analysis for BP, MF, CC ontologies
      go_result_bp <- enrichGO(gene = upregulated_entrez_ids, OrgDb = org.Mm.eg.db, ont = "BP", pvalueCutoff = 0.05)
      go_result_mf <- enrichGO(gene = upregulated_entrez_ids, OrgDb = org.Mm.eg.db, ont = "MF", pvalueCutoff = 0.05)
      go_result_cc <- enrichGO(gene = upregulated_entrez_ids, OrgDb = org.Mm.eg.db, ont = "CC", pvalueCutoff = 0.05)
      perform_go_analysis(go_result_bp, neuronal_subtype, comparison_name, "Upregulated", "BP")
      perform_go_analysis(go_result_mf, neuronal_subtype, comparison_name, "Upregulated", "MF")
      perform_go_analysis(go_result_cc, neuronal_subtype, comparison_name, "Upregulated", "CC")
    }
  }
  
  # GO analysis for downregulated genes
  if (length(downregulated_genes) > 0) {
    downregulated_entrez_ids <- mapIds(org.Mm.eg.db, keys = downregulated_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    downregulated_entrez_ids <- downregulated_entrez_ids[!is.na(downregulated_entrez_ids)]
    if (length(downregulated_entrez_ids) > 0) {
      # Run GO analysis for BP, MF, CC ontologies
      go_result_bp <- enrichGO(gene = downregulated_entrez_ids, OrgDb = org.Mm.eg.db, ont = "BP", pvalueCutoff = 0.05)
      go_result_mf <- enrichGO(gene = downregulated_entrez_ids, OrgDb = org.Mm.eg.db, ont = "MF", pvalueCutoff = 0.05)
      go_result_cc <- enrichGO(gene = downregulated_entrez_ids, OrgDb = org.Mm.eg.db, ont = "CC", pvalueCutoff = 0.05);
      perform_go_analysis(go_result_bp, neuronal_subtype, comparison_name, "Downregulated", "BP");
      perform_go_analysis(go_result_mf, neuronal_subtype, comparison_name, "Downregulated", "MF");
      perform_go_analysis(go_result_cc, neuronal_subtype, comparison_name, "Downregulated", "CC");
    }
  }
  
  # KEGG analysis for upregulated genes
  if (length(upregulated_genes) > 0) {
    upregulated_entrez_ids <- mapIds(org.Mm.eg.db, keys = upregulated_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    upregulated_entrez_ids <- upregulated_entrez_ids[!is.na(upregulated_entrez_ids)]
    if (length(upregulated_entrez_ids) > 0) {
      # Run KEGG analysis
      kegg_result_up <- enrichKEGG(gene = upregulated_entrez_ids, organism = 'mmu', pvalueCutoff = 0.05)
      perform_kegg_analysis(kegg_result_up, neuronal_subtype, comparison_name, "Upregulated")
    }
  }
  
  # KEGG analysis for downregulated genes
  if (length(downregulated_genes) > 0) {
    downregulated_entrez_ids <- mapIds(org.Mm.eg.db, keys = downregulated_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    downregulated_entrez_ids <- downregulated_entrez_ids[!is.na(downregulated_entrez_ids)]
    if (length(downregulated_entrez_ids) > 0) {
      # Run KEGG analysis
      kegg_result_down <- enrichKEGG(gene = downregulated_entrez_ids, organism = 'mmu', pvalueCutoff = 0.05)
      perform_kegg_analysis(kegg_result_down, neuronal_subtype, comparison_name, "Downregulated")
    }
  }
  
  # Add results to table if no DEGs were found
  if (length(upregulated_genes) == 0 && length(downregulated_genes) == 0) {
    results_table <<- rbind(results_table, data.frame(NeuronalSubtype = neuronal_subtype, Comparison = comparison_name, PathwayType = "None", Pathway = "No DEGs found", GeneCount = NA, P.adjust = NA))
  }
}

# Define the comparisons for each condition pair with correct case
comparisons <- list(
  "Aqp4_5d_vs_Ctrl_5d" = c("Aqp4_5d", "Ctrl_5d"),
  "Aqp4_3d_vs_Ctrl_3d" = c("Aqp4_3d", "Ctrl_3d")
)

# Adjust group range based on actual groups
neuron_groups <- unique(neuronsData$neuronal_subtype)

# Run analysis across groups and comparisons
for (neuronal_subtype in neuron_groups) {
  for (comparison_name in names(comparisons)) {
    conditions <- comparisons[[comparison_name]]
    perform_deg_pathway_analysis_for_subtype_condition(neuronsData, neuronal_subtype, comparison_name, conditions[1], conditions[2])
  }
}

# Save results table to Excel
write.xlsx(results_table, paste0(save_path, "Pathway_Analysis_Results_Neuronal_Subtypes.xlsx"), rowNames = FALSE)

cat("Neuronal subtype-based pathway analysis for condition comparisons is complete. The results are saved in", save_path)

# GO for thesis plot -----
# Define neuronal subtypes
glutamatergic_clusters <- c(0, 1, 2, 4, 6, 10, 11)
gabaergic_clusters <- c(5, 8)
other_clusters <- c(3, 7, 9)

save_path <- save_neurons_path

# Add neuronal_subtype column to metadata
neuronsData$neuronal_subtype <- ifelse(
  neuronsData$seurat_clusters %in% glutamatergic_clusters, "Glutamatergic",
  ifelse(neuronsData$seurat_clusters %in% gabaergic_clusters, "Gabaergic", "Other")
)

# Create condition_day column in metadata
neuronsData@meta.data$condition_day <- paste(neuronsData@meta.data$condition, neuronsData@meta.data$day, sep = "_")

# Function to perform combined GO enrichment analysis
perform_go_analysis_combined <- function(go_results_list, neuronal_subtype, comparison_name, type) {
  # Check if any GO results exist
  if (!any(sapply(go_results_list, function(x) !is.null(x) && nrow(x) > 0))) {
    cat("No significant GO pathways found for", type, "genes in", neuronal_subtype, "(", comparison_name, ")\n")
    return(NULL)
  }
  
  # Combine results into a single data frame
  combined_go <- do.call(rbind, lapply(names(go_results_list), function(ontology) {
    if (!is.null(go_results_list[[ontology]]) && nrow(go_results_list[[ontology]]) > 0) {
      df <- as.data.frame(go_results_list[[ontology]])
      df$Ontology <- ontology
      df
    } else {
      NULL
    }
  }))
  
  # Select top 20 pathways for each ontology
  combined_go <- combined_go %>%
    group_by(Ontology) %>%
    arrange(p.adjust) %>%
    slice_head(n = 20) %>%
    ungroup()
  
  # Assign colors for ontologies
  ontology_colors <- c("BP" = "green", "MF" = "blue", "CC" = "red")
  
  # Create the combined plot
  combined_plot <- ggplot(combined_go, aes(
    x = Count,
    y = reorder(Description, Count),
    fill = Ontology
  )) +
    geom_bar(stat = "identity") +
    geom_text(
      aes(label = sprintf("%.2f", -log10(p.adjust))),
      hjust = -0.1, size = 3
    ) +
    facet_grid(rows = vars(Ontology), scales = "free_y", space = "free_y") +
    scale_fill_manual(values = ontology_colors, name = "Ontology") +
    labs(
      x = "Gene Count",
      y = "Pathway Description",
      title = paste("Combined GO Analysis (", type, " DEGs) for", neuronal_subtype, "(", comparison_name, ")"),
      subtitle = "Adjusted p-values are shown as -log10(p.adjust)"
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 10),
      strip.text = element_text(size = 10)
    )
  
  # Save the plot
  ggsave(
    filename = paste0(save_path, "GO_Combined_", type, "_", neuronal_subtype, "_", comparison_name, ".png"),
    plot = combined_plot,
    width = 12, height = 10
  )
  
  cat("Combined GO plot saved for", neuronal_subtype, "(", comparison_name, ")\n")
}

# Function to perform pathway enrichment analysis
perform_deg_pathway_analysis_for_subtype_condition <- function(data, neuronal_subtype, comparison_name, condition1, condition2) {
  # Subset data for the specific neuronal subtype
  subtype_cells <- subset(data, neuronal_subtype == neuronal_subtype)
  
  # Set identities to condition_day within the subtype subset
  Idents(subtype_cells) <- subtype_cells$condition_day
  
  # Subset data to only the two conditions being compared
  subtype_cells <- subset(subtype_cells, idents = c(condition1, condition2))
  
  # Check if both conditions are present and have sufficient cells
  if (sum(Idents(subtype_cells) == condition1) < 10 | sum(Idents(subtype_cells) == condition2) < 10) {
    cat("Skipping", neuronal_subtype, "(", comparison_name, "): Insufficient cells\n")
    return(NULL)
  }
  
  # Perform differential expression analysis
  de_results <- FindMarkers(subtype_cells, ident.1 = condition1, ident.2 = condition2, logfc.threshold = 0.25, test.use = "wilcox")
  
  # Separate upregulated and downregulated genes
  upregulated_genes <- rownames(de_results[de_results$avg_log2FC > 0 & de_results$p_val_adj < 0.05, ])
  downregulated_genes <- rownames(de_results[de_results$avg_log2FC < 0 & de_results$p_val_adj < 0.05, ])
  
  # GO analysis for upregulated genes
  if (length(upregulated_genes) > 0) {
    upregulated_entrez_ids <- mapIds(org.Mm.eg.db, keys = upregulated_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    upregulated_entrez_ids <- upregulated_entrez_ids[!is.na(upregulated_entrez_ids)]
    if (length(upregulated_entrez_ids) > 0) {
      go_results_list <- list(
        "BP" = enrichGO(gene = upregulated_entrez_ids, OrgDb = org.Mm.eg.db, ont = "BP", pvalueCutoff = 0.05),
        "MF" = enrichGO(gene = upregulated_entrez_ids, OrgDb = org.Mm.eg.db, ont = "MF", pvalueCutoff = 0.05),
        "CC" = enrichGO(gene = upregulated_entrez_ids, OrgDb = org.Mm.eg.db, ont = "CC", pvalueCutoff = 0.05)
      )
      perform_go_analysis_combined(go_results_list, neuronal_subtype, comparison_name, "Upregulated")
    }
  }
  
  # GO analysis for downregulated genes
  if (length(downregulated_genes) > 0) {
    downregulated_entrez_ids <- mapIds(org.Mm.eg.db, keys = downregulated_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    downregulated_entrez_ids <- downregulated_entrez_ids[!is.na(downregulated_entrez_ids)]
    if (length(downregulated_entrez_ids) > 0) {
      go_results_list <- list(
        "BP" = enrichGO(gene = downregulated_entrez_ids, OrgDb = org.Mm.eg.db, ont = "BP", pvalueCutoff = 0.05),
        "MF" = enrichGO(gene = downregulated_entrez_ids, OrgDb = org.Mm.eg.db, ont = "MF", pvalueCutoff = 0.05),
        "CC" = enrichGO(gene = downregulated_entrez_ids, OrgDb = org.Mm.eg.db, ont = "CC", pvalueCutoff = 0.05)
      )
      perform_go_analysis_combined(go_results_list, neuronal_subtype, comparison_name, "Downregulated")
    }
  }
}

# Define condition comparisons
comparisons <- list(
  "Aqp4_5d_vs_Ctrl_5d" = c("Aqp4_5d", "Ctrl_5d"),
  "Aqp4_3d_vs_Ctrl_3d" = c("Aqp4_3d", "Ctrl_3d")
)

# Run analysis for all neuronal subtypes and comparisons
for (neuronal_subtype in unique(neuronsData$neuronal_subtype)) {
  for (comparison_name in names(comparisons)) {
    conditions <- comparisons[[comparison_name]]
    perform_deg_pathway_analysis_for_subtype_condition(neuronsData, neuronal_subtype, comparison_name, conditions[1], conditions[2])
  }
}

cat("GO analysis complete. Combined plots saved in:", save_path)