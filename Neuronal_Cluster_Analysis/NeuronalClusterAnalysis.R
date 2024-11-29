###################################################################
# Script Name: Final Neuronal Cluster Analysis
# Author: Yona Perstat
# Date: 2024-10
# Description: 
#   This script performs clustering analysis on neuronal clusters 
#   derived from single-cell RNA-seq data. It includes preprocessing,
#   batch correction, clustering, and automated labeling, followed by 
#   sub-clustering of neuronal cells and pathway analysis. 
#
# Input:
#   - Seurat object with single-cell RNA-seq data
#
# Output:
#   - UMAP plots (PNG)
#   - plots of different calculations ... 
#   - Clustered Seurat object (RDS)
###################################################################

# follows the OSCA and Seurat Basics (Bioconductor and Seurat)

# 1. Load necessary libraries for analysis
# 2. Load the Seurat object containing the data

# 3. Preprocess the data: 
# data is already clustered (0-22) and normalized
# - Filter out low-quality cells based on mitochondrial gene content
# - Apply thresholds for the number of detected features (genes) via 
#   identification of highly variable features (FindVariableFeatures)
#   or check if this is already done! 
# Is the data already scaled? 

# 4. Batch integration
# - Check if batches are already integrated
# - If not, perform batch correction using a suitable method (e.g., Harmony, CCA)

# 5. Check for doublets (can I find out if that is necessary or not?)
# - Use DoubletFinder or a similar method to identify and remove potential doublets

# 6. Visualize clusters 
# - Create UMAP or t-SNE plots to visually inspect clusters

# 7. Automated cluster labeling
# - Use automated labeling tools like SingleR or CellTypist to annotate clusters
# - Try alternative labeling tools as well to compare results

# 8. Save results of the clustering and labeling

# Neuronal Clusters: 
# 9. Focus on neuronal clusters 
# filter out non-neuronal cells from the neuronal clusters 
# - Identify neuronal clusters based on neuronal markers
# - Perform sub-clustering to find neuronal subtypes (best with separate code!)

# Separate code 
# 10. Perform further analysis on neuronal subtypes
# - Find differentially expressed features (use the same plots as for Astrocytes)
# - KEGG and GO BP enrichment analysis for neuronal subtypes
# - Visualize marker gene expression for each subtype

# 11. Save all final results and visualizations



# Load Libraries ----------------------

# install BiocManager and devtools! 

# Core libraries for single-cell analysis and general data handling
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("SeuratData", quietly = TRUE)) install.packages("SeuratData")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")

# Doublet detection
if (!requireNamespace("DoubletFinder", quietly = TRUE)) install.packages("DoubletFinder")

# Batch correction and Python integration
if (!requireNamespace("reticulate", quietly = TRUE)) install.packages("reticulate")

# Clustering analysis and dimensionality reduction
if (!requireNamespace("clustree", quietly = TRUE)) install.packages("clustree")
if (!requireNamespace("dendextend", quietly = TRUE)) install.packages("dendextend")
if (!requireNamespace("mclust", quietly = TRUE)) install.packages("mclust")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")

# Install scClustViz for clustering visualization
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("BaderLab/scClustViz")

# Install tooManyCells for cell clustering validation
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("GregorySchwartz/tooManyCellsR")

# Install packages for automated cell labeling
install.packages("devtools")
install.packages("pheatmap")
install.packages("shinythemes")
remotes::install_github("ggjlab/scMCA")
remotes::install_github("RausellLab/CellID")
BiocManager::install("SingleR")
BiocManager::install("celldex")

# Load core libraries for single-cell analysis and general data handling
library(SeuratData)
library(dplyr)
library(ggplot2)
library(patchwork)

# Load libraries for doublet detection
library(DoubletFinder)

# Load libraries for batch correction and Python integration
library(reticulate)

# Load clustering analysis and dimensionality reduction tools
library(clustree)
library(dendextend)
library(mclust)
library(reshape2)
library(scater)
library(scran)

# Load clustering and visualization tools
library(scClustViz)
library(TooManyCellsR)  # Correcting the name here
library(pheatmap)

# Load cluster validation and visualization tools
library(NbClust)
library(factoextra)
library(cluster)
library(PCAtools)

# Load memory-efficient data handling tools
library(ff)

# Load statistical analysis and parallel analysis tools
library(psych)

# Load libraries for automated cell labeling
library(scMCA)
library(SingleCellExperiment)
library(CelliD)
library(SingleR)
library(celldex)
library(Azimuth)

# Load data ------------
# Load the Seurat object (adjust the path accordingly)
Aqp4 <- readRDS("C:/Users/yperstat/Documents/Studium/IDB_Master/Master_thesis/Coding_Analysis/Data/aqp4_all.rds")

# Save path to folder where you want to save your plots
save_path <- "C:/Users/yperstat/Documents/Studium/IDB_Master/Master_thesis/Coding_Analysis/NeuronClusterAnalysis/plots/"

# Filter out low quality cells ----
# Step 1: Preprocess the Data - Filter Mitochondrial Genes and Low-Quality Cells
VlnPlot(Aqp4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Aqp4[["percent.mt"]] <- PercentageFeatureSet(Aqp4, pattern = "^mt-")
Aqp4 <- subset(Aqp4, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Plot again to visualize post-filtering
VlnPlot(Aqp4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Check for performed preprocessing steps --------------- 

# Step 2: Visualize Clustering (To check if the data is already clustered)
# If the Seurat object already has clusters assigned, this plot will show them
DimPlot(Aqp4, reduction = "umap", label = TRUE)
ggsave(filename = paste0(save_path, "umap_clustering_plot.png"), 
       plot = last_plot(), 
       width = 10, height = 8, units = "in", 
       dpi = 300)

# Check the Seurat object to confirm whether clusters exist
print(Aqp4@meta.data)  # This will display if there are pre-existing clusters (e.g., "seurat_clusters" column)
Aqp4 
# An object of class Seurat 
# 18302 features across 62672 samples within 1 assay 
# Active assay: RNA (18302 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap
# so all important preprocessing steps have been  performed: normalization, variable feature selection, scaling, pca and UMAP 

# ------------- CLUSTERING:  ---------
# Recluster because Fufu originally clustered searching for Astrocytes and for Neurons the optimal clustering could be different 
# Identify the "optimal" number of PCs  -------
# Here we try to identify the "optimal" number of PCs (35?)
# the choice of PCs is always a trade-off between noise and biological signal 

# 1. Elbow Plot to visually identify the optimal number of PCs
ElbowPlot(Aqp4, ndims = 50)  # Helps locate the 'elbow', where additional PCs contribute less variance.

# 2. Variance Explained by Each PC
# This shows how much variance each PC explains, helping to quantify their significance.
pca_var <- Stdev(Aqp4, reduction = "pca")^2
percent_var <- pca_var / sum(pca_var) * 100
plot(percent_var, type = "b", log = "y", xlab = "PC", ylab = "Variance explained (%)", main = "Variance Explained by Each PC")

# 3. Save DimHeatMaps for each PC from 1 to 50
# DimHeatMaps visualize the structure captured by each PC, giving insights into biological signals.
heatmap_save_path <- paste0(save_path, "heatmaps/")
if (!dir.exists(heatmap_save_path)) dir.create(heatmap_save_path)
for (i in 1:50) {
  png(filename = paste0(heatmap_save_path, "DimHeatmap_PC_", i, ".png"), width = 1200, height = 800, res = 150)
  DimHeatmap(Aqp4, dims = i, cells = 500, balanced = TRUE)
  dev.off()
}

# 4. JackStraw Analysis to assess the statistical significance of PCs (computationally intensive).
Aqp4 <- JackStraw(Aqp4, num.replicate = 100, dims = 50)
Aqp4 <- ScoreJackStraw(Aqp4, dims = 1:50)
JackStrawPlot(Aqp4, dims = 1:50)

# 5. Extract JackStraw p-values for further analysis.
jackstraw_scores <- Aqp4@reductions$pca@jackstraw$empirical.p.values
jackstraw_df <- as.data.frame(jackstraw_scores)

# 6. Perform Horn’s Parallel Analysis to compare observed eigenvalues to those from random data.
pca_embeddings <- Embeddings(Aqp4, "pca")
fa_parallel_results <- fa.parallel(pca_embeddings, fa = "pc", n.iter = 100, show.legend = FALSE, main = "Horn's Parallel Analysis")

# Doublet Detection: -----------

# After batch correction and clustering are complete
# Step 1: Estimate the homotypic proportion (adjust for similarity between doublets and real cells)
homotypic_prop <- modelHomotypic(Aqp4@meta.data$seurat_clusters) 

# Step 2: Perform paramSweep to optimize the number of PCs for doublet detection
sweep.res <- paramSweep(Aqp4, PCs = 1:30, sct = FALSE)

# Step 3: Summarize results from paramSweep
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)

# Step 4: Find the optimal pK value
best_pK <- find.pK(sweep.stats)
print(best_pK)

# Step 5: Set the optimal pK value (based on your earlier output or find.pK results)
optimal_pK <- 0.19  # Adjust based on the actual optimal pK found

# Step 6: Calculate the number of expected doublets
nExp <- round(ncol(Aqp4) * 0.075)  # Assuming a 7.5% doublet rate

# Step 7: Run DoubletFinder to classify cells as singlets or doublets
Aqp4 <- doubletFinder(
  Aqp4, 
  PCs = 1:30,        # Number of PCs to use
  pN = 0.25,         # Default pN
  pK = optimal_pK,   # Optimal pK found earlier
  nExp = nExp,       # Number of expected doublets
  reuse.pANN = FALSE # Set to FALSE as this is the first run
)

# Step 8: Visualize doublets on UMAP
DimPlot(Aqp4, reduction = "umap", group.by = "DF.classifications_0.25")
ggsave(filename = paste0(save_path, "umap_doubletFinder_plot.png"), 
       plot = last_plot(), width = 10, height = 8, units = "in", dpi = 300)

# Step 9: Optionally subset to retain only singlets
Aqp4_singlets <- subset(Aqp4, subset = DF.classifications_0.25 == "Singlet")

# Step 10: Re-run UMAP for singlets only
DimPlot(Aqp4_singlets, reduction = "umap", label = TRUE)
ggsave(filename = paste0(save_path, "umap_singlets_only_plot.png"), 
       plot = last_plot(), width = 10, height = 8, units = "in", dpi = 300)



# Cell labeling: 
# Find "optimal" Resolution for Clustering -----
# use the top 3 different resolutions and use them in the downstream anylsis 

# Why test for multiple resolutions?
# Testing multiple resolutions allows us to find an optimal clustering of the data.
# Different resolutions in Seurat change the granularity of clusters, meaning low resolutions 
# may combine distinct cell types into one cluster, while high resolutions may over-segment clusters.
# By evaluating cluster merging and splitting behavior, biological relevance, and computational metrics (such as silhouette scores and ARI),
# we can select the resolution that provides meaningful and stable clusters.

# This code block tests resolutions from 0.01 to 1 and provides tools for evaluating the optimal resolution 
# for downstream analyses.

# We will use multiple methods for assessing clustering quality, including:
# 1. UMAP plots: Visually inspect cluster separation.
# 2. Clustree: Visualize how clusters merge and split across resolutions.
# 3. ARI: Compare cluster similarity between resolutions.
# 4. Number of clusters vs resolution: Check how cluster count changes.
# 5. Heatmaps of marker genes: See the top marker genes for clusters at each resolution.
# 6. tooManyCells: Explore the tree structure and clumpiness of clusters across resolutions.


# DEFINE RESOLUTIONS:  

# Define the range of resolutions to test (0.01 to 1.0)
resolutions <- c(0.01, 0.05, 0.1, 0.2, 0.25, 0.3, 0.5, 0.8, 1.0)

# Run UMAP using 35 PCs (do this once)
Aqp4 <- FindNeighbors(Aqp4, dims = 1:35)
Aqp4 <- RunUMAP(Aqp4, dims = 1:35)


#  STORE CLUSTER RESULTS 

# Store cluster results for different resolutions in the Seurat metadata
for (res in resolutions) {
  Aqp4 <- FindClusters(Aqp4, resolution = res)  # Find clusters with the given resolution
  Aqp4[[paste0("RNA_snn_res.", res)]] <- Idents(Aqp4)  # Save cluster IDs for each resolution
}


#  UMAP PLOTS FOR EACH RESOLUTION 

# Generate UMAP for each resolution and save the plots
for (res in resolutions) {
  umap_plot <- DimPlot(Aqp4, reduction = "umap", group.by = paste0("RNA_snn_res.", res), label = TRUE) +
    ggtitle(paste("UMAP Plot - Resolution", res))
  
  # Save UMAP plots
  ggsave(filename = paste0(save_path, "umap_res_", res, ".png"), plot = umap_plot,
         width = 10, height = 8, units = "in", dpi = 300)
}


# CLUSTREE PLOTS 

# Install clustree if not already installed
if (!requireNamespace("clustree", quietly = TRUE)) install.packages("clustree")
library(clustree)

# Create a clustree plot to visualize how clusters merge and split
clustree(Aqp4, prefix = "RNA_snn_res.") +
  ggtitle("Clustree: Clustering Tree Across Resolutions")

# Save the clustree plot
ggsave(filename = paste0(save_path, "clustree_plot.png"), width = 10, height = 8, units = "in", dpi = 300)


# DENDROGRAM PLOTS 

# Install dendextend for dendrogram visualization if not already installed
if (!requireNamespace("dendextend", quietly = TRUE)) install.packages("dendextend")
library(dendextend)

# Perform hierarchical clustering on PCA embeddings
pca_embeddings <- Embeddings(Aqp4, "pca")[, 1:30]  # Use first 30 PCs
hc <- hclust(dist(pca_embeddings), method = "ward.D2")

# Plot dendrogram
dend <- as.dendrogram(hc)
plot(dend, main = "Hierarchical Clustering Dendrogram")

# Save dendrogram plot
ggsave(filename = paste0(save_path, "dendrogram.png"), width = 10, height = 8, units = "in", dpi = 300)



#  ADJUSTED RAND INDEX (ARI) 

# Step 1: Initialize an empty matrix for ARI values
ari_matrix <- matrix(NA, nrow = length(resolutions), ncol = length(resolutions))
rownames(ari_matrix) <- resolutions
colnames(ari_matrix) <- resolutions

# Step 2: Fill the ARI matrix with calculated values
for (res1 in resolutions) {
  for (res2 in resolutions) {
    if (res1 != res2) {
      # Extract cluster identities from the metadata
      clusters_res1 <- as.numeric(Aqp4@meta.data[[paste0("RNA_snn_res.", res1)]])
      clusters_res2 <- as.numeric(Aqp4@meta.data[[paste0("RNA_snn_res.", res2)]])
      
      # Calculate Adjusted Rand Index (ARI)
      ari_value <- adjustedRandIndex(clusters_res1, clusters_res2)
      ari_matrix[as.character(res1), as.character(res2)] <- ari_value
    } else {
      # Perfect similarity for identical resolutions
      ari_matrix[as.character(res1), as.character(res2)] <- 1  
    }
  }
}

# Step 3: Convert the matrix into a long format for ggplot
ari_long <- melt(ari_matrix, na.rm = TRUE)

# Step 4: Create a heatmap using ggplot
ggplot(ari_long, aes(x = factor(Var1), y = factor(Var2), fill = value)) +
  geom_tile(color = "white") + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5, 
                       name = "ARI") +
  labs(title = "Adjusted Rand Index (ARI) Heatmap between Resolutions",
       x = "Resolution", y = "Resolution") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave(filename = paste0(save_path, "ARI_heatmap.png"), width = 10, height = 8, dpi = 300)


#  CLUSTERS VS RESOLUTION 

# Initialize a data frame to store the results
cluster_counts <- data.frame(Resolution = numeric(), Number_of_Clusters = numeric())

# Loop through each resolution, perform clustering, and count the number of clusters
for (res in resolutions) {
  # Perform clustering at the given resolution
  Aqp4 <- FindClusters(Aqp4, resolution = res)
  
  # Count the number of unique clusters
  num_clusters <- length(unique(Idents(Aqp4)))
  
  # Append the result to the data frame
  cluster_counts <- rbind(cluster_counts, data.frame(Resolution = res, Number_of_Clusters = num_clusters))
}

# Plot Number of Clusters (k) vs Resolution with finer axis grading
ggplot(cluster_counts, aes(x = Resolution, y = Number_of_Clusters)) +
  geom_line() + 
  geom_point() +
  labs(title = "Number of Clusters vs Resolution", x = "Resolution", y = "Number of Clusters") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.05)) +  # Finer grading on x-axis (Resolution)
  scale_y_continuous(breaks = seq(min(cluster_counts$Number_of_Clusters), max(cluster_counts$Number_of_Clusters), by = 1)) +  # Finer grading on y-axis (Clusters)
  theme_minimal()

# Optionally, save the plot
ggsave(filename = paste0(save_path, "clusters_vs_resolution_finer_grading.png"), width = 10, height = 8, dpi = 300)


# HEATMAPS OF TOP MARKER GENES 

# Generate heatmaps of the top marker genes for each cluster at different resolutions.

for (res in resolutions) {
  # Perform clustering at the given resolution
  Aqp4 <- FindClusters(Aqp4, resolution = res)
  
  # Identify the markers for each cluster
  markers <- FindAllMarkers(Aqp4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  # Get the top 10 markers per cluster
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  
  # Create heatmap for top 10 marker genes
  heatmap_plot <- DoHeatmap(Aqp4, features = top10$gene, label = TRUE) + 
    ggtitle(paste("Heatmap of Top 10 Markers for Resolution", res))
  
  # Save heatmap to file
  ggsave(filename = paste0(save_path, "heatmap_markers_res_", res, ".png"), plot = heatmap_plot,
         width = 10, height = 8, dpi = 300)
}

# TOO MANY CELLS ANALYSIS 

# Step 3: Loop through each resolution and run TooManyCells
for (res in resolutions) {
  
  # Step 4: Run clustering in Seurat with the current resolution
  Aqp4 <- FindClusters(Aqp4, resolution = res)
  
  # Extract cluster assignments for this resolution
  clusters <- as.numeric(Idents(Aqp4))
  
  # Save cluster results for tooManyCells (this is optional but useful for validation later)
  cluster_labels <- data.frame(
    item = colnames(mat),  # Cell barcodes
    label = clusters       # Cluster assignments from Seurat
  )
  
  # Save the cluster labels to a CSV file for tooManyCells
  labels_file <- paste0("labels_res_", res, ".csv")
  write.csv(cluster_labels, labels_file, row.names = FALSE)
  
  # Step 5: Run tooManyCells with the current resolution
  args_list <- c(
    "make-tree",
    "--smart-cutoff", "4",
    "--min-size", "1",
    "--labels-file", labels_file  # Include the cluster labels
  )
  
  # Run tooManyCells and store the result
  res_tmc <- tooManyCells(mat = mat, args = args_list)
  
  # Step 6: Plot and save the tree structure for this resolution
  png(filename = paste0("tree_plot_res_", res, ".png"), width = 1200, height = 800)
  plot(res_tmc$treePlot, axes = FALSE)
  dev.off()
  
  # Optionally, save clumpiness plot and other results
  png(filename = paste0("clumpiness_plot_res_", res, ".png"), width = 1200, height = 800)
  plot(res_tmc$clumpinessPlot, axes = FALSE)
  dev.off()
  
  # Output node information for further analysis (optional)
  write.csv(res_tmc$nodeInfo, file = paste0("node_info_res_", res, ".csv"))
  
  # Print the result for each resolution
  print(paste("Completed TooManyCells analysis for resolution:", res))
}

# GAP Statistics 
gap_stat <- clusGap(Embeddings(Aqp4, "pca"), FUN = kmeans, nstart = 25, K.max = 10, B = 50)
print(gap_stat)
plot(gap_stat)


# Evaluate the Clustering -------

# Count the number of cells per cluster
cluster_counts <- table(Idents(Aqp4))  # 'Idents' contains the cluster identities of cells
print(cluster_counts)

# Step 1: Extract PCA Embeddings and Clusters from Seurat object
pca_embeddings <- Embeddings(Aqp4, "pca")[, 1:30]  # Use the first 30 PCs

# Step 2: Calculate Silhouette Scores for the Existing Clusters
dist_matrix <- dist(pca_embeddings)  # Create the distance matrix using PCA embeddings
clusters <- Aqp4@meta.data$seurat_clusters  # Clustering results from metadata
sil_scores <- silhouette(as.numeric(clusters), dist_matrix)

# Step 3: Visualize Silhouette Scores
sil_plot <- fviz_silhouette(sil_scores) + ggtitle("Silhouette Scores for All Cells")
print(sil_plot)

# Step 4: Save the silhouette plot
ggsave(filename = paste0(save_path, "silhouette_scores_all_cells.pdf"), plot = sil_plot)

#  Optimal Cluster Determination with NbClust 
library(NbClust)
nbclust_results <- NbClust(data = pca_embeddings[, 1:30],  # Use first 30 PCs
                           distance = "euclidean", 
                           min.nc = 2, max.nc = 15, 
                           method = "kmeans")

# Step 5: Visualize NbClust Results
barplot(table(nbclust_results$Best.n[1, ]), 
        xlab = "Number of Clusters", 
        ylab = "Number of Criteria", 
        main = "NbClust Optimal Number of Clusters")

# Silhouette Scores for All Cells 
sil_scores_all <- silhouette(as.numeric(clusters), dist(pca_embeddings))

# Step 6: Visualize Silhouette Scores for All Cells
sil_plot_all <- fviz_silhouette(sil_scores_all) + ggtitle("Silhouette Scores for All Cells")
print(sil_plot_all)

# Step 7: Save the silhouette plot
ggsave(filename = paste0(save_path, "silhouette_scores_all_cells.pdf"), plot = sil_plot_all)

# Check variable features ------------

# Step 1: Check if the variable features are already selected in Seurat
variable_genes_seurat <- VariableFeatures(Aqp4)

# Step 2: If no variable features were selected, find them using Seurat
if (length(variable_genes_seurat) == 0) {
  message("No variable features found. Selecting highly variable features using Seurat.")
  Aqp4 <- FindVariableFeatures(Aqp4, selection.method = "vst", nfeatures = 2000)
  variable_genes_seurat <- VariableFeatures(Aqp4)
} else {
  message("Variable features are already selected.")
}

# Step 3: Convert Seurat object to SingleCellExperiment for scran analysis
sce <- as.SingleCellExperiment(Aqp4)

# Step 4: Use scran to select highly variable genes
dec_scran <- modelGeneVar(sce)

# Step 5: Visualize the variance against the mean
fit_scran <- metadata(dec_scran)
plot(fit_scran$mean, fit_scran$var, xlab = "Mean of log-expression",
     ylab = "Variance of log-expression")
curve(fit_scran$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)

# Step 6: Choose the highly variable genes from scran
hvg_scran <- getTopHVGs(dec_scran, prop = 0.1)  # Adjust `prop` as needed

# Step 7: Compare Seurat and scran HVGs
hvg_comparison <- data.frame(
  Gene = rownames(sce),
  Scran_HVG = rownames(sce) %in% hvg_scran,
  Seurat_HVG = rownames(sce) %in% variable_genes_seurat
)

# Step 8: Summarize results of comparison
hvg_summary <- table(hvg_comparison$Scran_HVG, hvg_comparison$Seurat_HVG)
print(hvg_summary)

# Step 9: Optionally visualize the comparison
ggplot(hvg_comparison, aes(x = Scran_HVG, fill = Seurat_HVG)) +
  geom_bar(position = "fill") +
  ylab("Proportion") +
  xlab("Scran Highly Variable Genes") +
  ggtitle("Comparison of Highly Variable Genes: scran vs Seurat") +
  scale_fill_manual(values = c("red", "blue"), labels = c("No", "Yes")) +
  theme_minimal()

# Save the plot
ggsave(filename = paste0(save_path, "hvg_comparison_plot.png"), width = 10, height = 8, units = "in", dpi = 300)

# ITS NECESSARY TO HAVE A LOOK ON THE BIOLOGY OF THE FOUND GENES BY Scran

# Step 10: Find genes unique to scran
unique_hvg_scran <- setdiff(hvg_scran, variable_genes_seurat)

# Print the unique genes found only by scran
cat("Unique HVGs identified by scran:\n")
print(unique_hvg_scran)

# Optionally, count the number of unique genes
num_unique_hvg_scran <- length(unique_hvg_scran)
cat("Number of genes identified only by scran:", num_unique_hvg_scran, "\n")

# Batch Correction and Integration using CCA and Harmony ------------

# Plot UMAP grouped by batch before correction
DimPlot(Aqp4, reduction = "umap", group.by = "batch", label = TRUE)
ggsave(filename = paste0(save_path, "umap_batch_before_integration.png"), 
       plot = last_plot(), width = 10, height = 8, units = "in", dpi = 300)

# Plot UMAP grouped by clusters before correction
DimPlot(Aqp4, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
ggsave(filename = paste0(save_path, "umap_clusters_before_integration.png"), 
       plot = last_plot(), width = 10, height = 8, units = "in", dpi = 300)

# Integration using CCA 

run_integration_cca <- function(seurat_obj, save_path, n_pcs = 20, res = 0.1) {
  message("Running integration with CCA for ", n_pcs, " PCs and resolution ", res)
  
  seurat_obj_list <- SplitObject(seurat_obj, split.by = "batch")
  
  # Find integration anchors using CCA
  anchors <- FindIntegrationAnchors(object.list = seurat_obj_list, 
                                    anchor.features = 2000, reduction = "cca", dims = 1:n_pcs)
  
  # Integrate data
  integrated <- IntegrateData(anchorset = anchors, dims = 1:n_pcs)
  integrated <- ScaleData(integrated)
  integrated <- RunPCA(integrated, npcs = n_pcs)
  
  # Find neighbors, clusters, and run UMAP
  integrated <- FindNeighbors(integrated, dims = 1:n_pcs)
  integrated <- FindClusters(integrated, resolution = res)
  integrated <- RunUMAP(integrated, dims = 1:n_pcs, reduction = "pca")
  
  # Save UMAP plots (grouped by batch and clusters)
  DimPlot(integrated, reduction = "umap", group.by = "batch", label = TRUE)
  ggsave(filename = paste0(save_path, "umap_batch_after_integration_CCA_", n_pcs, "_PCs_", res, "_res.png"), 
         plot = last_plot(), width = 10, height = 8, units = "in", dpi = 300)
  
  DimPlot(integrated, reduction = "umap", label = TRUE)
  ggsave(filename = paste0(save_path, "umap_clusters_after_integration_CCA_", n_pcs, "_PCs_", res, "_res.png"), 
         plot = last_plot(), width = 10, height = 8, units = "in", dpi = 300)
  
  return(integrated)
}

# Run multiple combinations of PCs and resolutions using CCA
integrated_Aqp4_20PCs_0_1 <- run_integration_cca(Aqp4, save_path, n_pcs = 20, res = 0.1)
integrated_Aqp4_20PCs_0_5 <- run_integration_cca(Aqp4, save_path, n_pcs = 20, res = 0.5)

# Integration using Harmony 
run_integration_harmony <- function(seurat_obj, save_path, n_pcs = 20, res = 0.1) {
  message("Running integration with Harmony for ", n_pcs, " PCs and resolution ", res)
  
  # Run Harmony batch correction
  seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "batch", reduction = "pca", dims.use = 1:n_pcs)
  
  # Scale data, run PCA on Harmony embeddings, and cluster cells
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = n_pcs, reduction = "harmony")
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_pcs)
  seurat_obj <- FindClusters(seurat_obj, resolution = res)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_pcs, reduction = "harmony")
  
  # Save UMAP plots (grouped by batch and clusters)
  DimPlot(seurat_obj, reduction = "umap", group.by = "batch", label = TRUE)
  ggsave(filename = paste0(save_path, "umap_batch_after_integration_Harmony_", n_pcs, "_PCs_", res, "_res.png"), 
         plot = last_plot(), width = 10, height = 8, units = "in", dpi = 300)
  
  DimPlot(seurat_obj, reduction = "umap", label = TRUE)
  ggsave(filename = paste0(save_path, "umap_clusters_after_integration_Harmony_", n_pcs, "_PCs_", res, "_res.png"), 
         plot = last_plot(), width = 10, height = 8, units = "in", dpi = 300)
  
  return(seurat_obj)
}

# Run multiple combinations of PCs and resolutions using Harmony
integrated_Aqp4_harmony_20PCs_0_1 <- run_integration_harmony(Aqp4, save_path, n_pcs = 20, res = 0.1)
integrated_Aqp4_harmony_20PCs_0_5 <- run_integration_harmony(Aqp4, save_path, n_pcs = 20, res = 0.5)

# Install and load the package
if (!requireNamespace("kBET", quietly = TRUE)) install.packages("kBET")
library(kBET)

# Compute kBET score
kbet_results <- kBET(pca_embeddings, batch = Aqp4$batch)

# Visualize results
plot(kbet_results$stats$kBET.observed, main = "kBET Batch Effect Assessment")

# AT THE END I NEED TO CHOOSE THE ONE THAT IS INTEGRATED BETTER FOR THE FURTHER ANALYSIS!!!! 

# ------- CELL LABELING: -----------------
# Automated Cell Labeling: ----------------

# Automated cell labeling section
# In this section, we perform automated cell type labeling using several tools 
# (scMCA, CellID, SingleR, Azimuth, and CellTypist). The results are then combined 
# using majority voting, and confidence scores are computed to assess the 
# consistency across methods.

# Step 1: Using scMCA for cell type labeling
scdata <- as.matrix(GetAssayData(Aqp4, slot = "counts"))
results <- scMCA(scdata, numbers_plot = 3)
Aqp4 <- AddMetaData(Aqp4, metadata = results$labels, col.name = "scMCA_labels")
DimPlot(Aqp4, group.by = "scMCA_labels")
ggsave(filename = paste0(save_path, "UMAP_scMCA_labels.png"))

# Step 2: Using CellID for cell type labeling
data("cellid_ref")
labels <- cellid_map(Aqp4, ref = cellid_ref)
Aqp4 <- AddMetaData(Aqp4, metadata = labels, col.name = "CellID_labels")
DimPlot(Aqp4, group.by = "CellID_labels")
ggsave(filename = paste0(save_path, "UMAP_CellID_labels.png"))

# Step 3: Using SingleR for cell type labeling
ref <- celldex::MouseRNAseqData()
expression_data <- as.matrix(GetAssayData(Aqp4, slot = "data", assay = "RNA"))
pred <- SingleR(test = expression_data, ref = ref, labels = ref$label.main)
Aqp4 <- AddMetaData(Aqp4, metadata = pred$labels, col.name = "SingleR_labels")
DimPlot(Aqp4, group.by = "SingleR_labels")
ggsave(filename = paste0(save_path, "UMAP_SingleR_labels.png"))

# Step 4: Using Azimuth for cell type labeling
Azimuth::InstallData("pbmc3k")
Aqp4 <- Azimuth::RunAzimuth(Aqp4, reference = "pbmc3k")
Aqp4 <- AddMetaData(Aqp4, metadata = Aqp4$predicted.celltype.l1, col.name = "Azimuth_labels")
DimPlot(Aqp4, group.by = "Azimuth_labels")
ggsave(filename = paste0(save_path, "UMAP_Azimuth_labels.png"))

# Step 5: Using CellTypist for cell type labeling
celltypist <- import("celltypist")
scdata <- as.matrix(GetAssayData(Aqp4, slot = "counts"))
model <- celltypist$models$Immune_All_Low()
prediction <- celltypist$annotate(scdata, model = model)
Aqp4 <- AddMetaData(Aqp4, metadata = prediction$predicted_labels, col.name = "CellTypist_labels")
DimPlot(Aqp4, group.by = "CellTypist_labels")
ggsave(filename = paste0(save_path, "UMAP_CellTypist_labels.png"))

# Step 6: Combining Labels Using Majority Voting
predictions <- data.frame(
  scMCA = Aqp4@meta.data$scMCA_labels,
  CellID = Aqp4@meta.data$CellID_labels,
  SingleR = Aqp4@meta.data$SingleR_labels,
  Azimuth = Aqp4@meta.data$Azimuth_labels,
  CellTypist = Aqp4@meta.data$CellTypist_labels
)

majority_vote <- function(labels) {
  return(names(sort(table(labels), decreasing = TRUE))[1])
}

combined_labels <- apply(predictions, 1, majority_vote)
Aqp4 <- AddMetaData(Aqp4, metadata = combined_labels, col.name = "combined_labels")

DimPlot(Aqp4, group.by = "combined_labels") + ggtitle("Majority Voting: Combined Cell Type Labels")
ggsave(filename = paste0(save_path, "UMAP_combined_labels.png"))

# Step 7: Calculating and Visualizing Confidence Scores
confidence_score <- function(labels) {
  label_count <- table(labels)
  max_count <- max(label_count)
  return(max_count / length(labels))
}

confidence_scores <- apply(predictions, 1, confidence_score)
Aqp4 <- AddMetaData(Aqp4, metadata = confidence_scores, col.name = "label_confidence")

FeaturePlot(Aqp4, features = "label_confidence", cols = c("lightblue", "red")) + 
  ggtitle("Confidence Scores")
ggsave(filename = paste0(save_path, "UMAP_confidence_scores.png"))

# Visualizing Confidence Score Distribution
ggplot(data = data.frame(Confidence = Aqp4@meta.data$label_confidence), aes(x = Confidence)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "darkblue") +
  labs(title = "Confidence Score Distribution", x = "Confidence Score", y = "Frequency") +
  theme_minimal()
ggsave(filename = paste0(save_path, "confidence_score_distribution.png"))


# Wilcoxon gene set cell label assignment ----

# Wilcoxon Gene Set Cell Label Assignment:
# In this section, we perform cell type labeling based on gene set activity using a Wilcoxon test.
# Marker sets are defined based on differential gene expression, and cell labels are assigned 
# using an Area Under the Curve (AUC) approach. The AUC scores are used to determine which gene 
# set is most active in each cell, and this is used for labeling. We also perform diagnostic 
# checks, including label comparison and threshold exploration, to validate the results.

# Step 1: Define marker sets using pairwiseWilcox
wilcox.z <- pairwiseWilcox(sce.zeisel, sce.zeisel$level1class, lfc = 1, direction = "up")
markers.z <- getTopMarkers(wilcox.z$statistics, wilcox.z$pairs, pairwise = FALSE, n = 50)

# Step 2: Create GeneSetCollection from markers
all.sets <- lapply(names(markers.z), function(x) GeneSet(markers.z[[x]], setName = x))
all.sets <- GeneSetCollection(all.sets)

# Step 3: Rank genes by expression within each cell and calculate AUC for marker sets
rankings <- AUCell_buildRankings(counts(sce.tasic), plotStats = FALSE, verbose = FALSE)
cell.aucs <- AUCell_calcAUC(all.sets, rankings)

# Step 4: Assign the top AUC marker set as the cell label
results <- t(assay(cell.aucs))
new.labels <- colnames(results)[max.col(results)]
sce.tasic$AUCell_labels <- new.labels

# Step 5: Diagnostic check by comparing new labels to known labels
tab <- table(new.labels, sce.tasic$broad_type)
print(tab)

# Step 6: Plot heatmap of label comparison
pheatmap(log2(tab + 10), color = colorRampPalette(c("white", "blue"))(101))

# Step 7: Explore thresholds for the AUC distributions
AUCell_exploreThresholds(cell.aucs, plotHist = TRUE, assign = TRUE)

# Gene Set Cell Label Assignment via GO ----
# IF I INDEED USE THIS PART I NEED TO CHANGE IT PROBABLY TO INCLUDE MORE DIFFERENT GO TERMS!!!!!!!!!
# This section uses gene set activity to assign cell type labels based on differential gene expression.
# We perform a Wilcoxon rank-sum test to identify marker genes for each cluster and use GO enrichment analysis
# to assign biological meaning to each cluster. The goal is to create biologically relevant annotations
# for each cluster using top marker genes and GO terms, which will help in better understanding the
# cell types or states present in the dataset.

# Convert Seurat object to SingleCellExperiment format for compatibility with scran and other Bioconductor tools
sce <- as.SingleCellExperiment(Aqp4)

# Step 1: Identify marker genes for all clusters using scran
markers <- scoreMarkers(sce, lfc = 1)

# Initialize a vector to store new cluster labels
new_cluster_labels <- vector("character", length = length(unique(Idents(Aqp4))))
names(new_cluster_labels) <- levels(Idents(Aqp4))

# Store GO enrichment results for all clusters
go_results_all_clusters <- list()

# Step 2: Iterate over all clusters and perform GO enrichment analysis
for (cluster_id in names(markers)) {
  
  # Get marker genes for the current cluster
  cur.markers <- markers[[cluster_id]]
  
  # Select top 100 genes based on the largest median log-fold change (Cohen's d)
  is.de <- order(cur.markers$median.logFC.cohen, decreasing = TRUE)[1:100]
  
  # Map gene symbols to Entrez IDs for GO analysis
  entrez.ids <- mapIds(org.Mm.eg.db, keys = rownames(cur.markers), 
                       column = "ENTREZID", keytype = "SYMBOL")
  
  # Step 3: Perform GO enrichment analysis using the top marker genes
  go.out <- goana(unique(entrez.ids[is.de]), species = "Mm", universe = unique(entrez.ids))
  
  # Filter for relevant biological process (BP) terms and sort by significance
  go.useful <- go.out[go.out$Ont == "BP" & go.out$N <= 200, ]
  
  # Save GO results for the current cluster
  go_results_all_clusters[[cluster_id]] <- go.useful
  
  # Step 4: Assign a new label based on the most significant GO term
  if (nrow(go.useful) > 0) {
    top_go_term <- go.useful$Term[1]
    short_label <- substr(top_go_term, 1, 30)  # Shorten the GO term for readability
    new_cluster_labels[cluster_id] <- short_label
  } else {
    new_cluster_labels[cluster_id] <- paste("Cluster", cluster_id)
  }
}

# Step 5: Add new cluster labels to the Seurat object metadata
Aqp4$GO_annotation <- factor(Idents(Aqp4), labels = new_cluster_labels[as.character(Idents(Aqp4))])

# Step 6: Plot the UMAP with new cluster labels based on GO terms
DimPlot(Aqp4, group.by = "GO_annotation", label = TRUE, repel = TRUE) + NoLegend()

# Optionally, save the UMAP plot
ggsave(filename = "UMAP_with_GO_annotations.png", width = 8, height = 6)

# Step 7: Calculate gene set scores for a specific GO term (e.g., "GO:0019432")
# Get genes for the specific GO term
go_term_id <- "GO:0019432"
tab <- select(org.Mm.eg.db, keytype = "SYMBOL", keys = rownames(sce), columns = "GOALL")
by.go <- split(tab[,1], tab[,2])
go_genes <- unique(by.go[[go_term_id]])
go_genes <- go_genes[go_genes %in% rownames(Aqp4)]  # Ensure genes are present in your dataset

# Step 8: Calculate the module score for the specific GO term
Aqp4 <- AddModuleScore(Aqp4, features = list(go_genes), name = "GO_0019432_Score")

# Step 9: Plot the gene set activity for GO:0019432 on the UMAP
FeaturePlot(Aqp4, features = "GO_0019432_Score1") + ggtitle("Activity of GO:0019432")

# Optionally, save the feature plot
ggsave(filename = "GO_0019432_Activity_UMAP.png", width = 8, height = 6)

# End of gene set cell label assignment based on Wilcoxon test

# Visual inspection of clusters based on markers --------
# Observations:
# - Clusters 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, and 20 express neuronal markers.
# - Clusters 1 and 7 should be excluded due to batch effects, cluster 10 has few cells, and cluster 4 is likely astrocytes.
# - Clusters 12 and 18 have low neuronal scores (12 = oligodendrocytes, 18 = endothelial cells).
# - Clusters 0, 2, 3, 5, 6, 14, 15, 16, and 19 have high neuronal marker scores and are candidates for neuronal clusters! 

# Define markers for neuronal and non-neuronal cell types
neuronal_markers <- c("Rbfox3", "Map2", "Tubb3", "Syt1", "Eno2")  # Rbfox3 is NeuN
oligodendrocyte_markers <- c("Mbp", "Mog", "Plp1", "Enpp6", "Sox10")
astrocyte_markers <- c("Aqp4", "Gfap", "Aldh1l1", "Slc1a3", "S100b")
microglia_markers <- c("Cx3cr1", "C1qa", "P2ry12", "Aif1", "Tmem119")
endothelial_markers <- c("Pecam1", "Cldn5", "Flt1", "Kdr")
pericyte_markers <- c("Acta2", "Tagln", "Cspg4")
ependymal_markers <- c("Foxj1", "Vim")
fibroblast_markers <- c("Col1a1", "Col1a2", "Pdgfra")
t_cell_markers <- c("Cd3", "Cd4", "Cd8", "Foxp3", "Tcf7", "Ikzf1", "IL7R", "CCR7", "S100A4")
b_cell_markers <- c("Cd19", "Ms4a1", "Cd79a", "Cd79b", "Pax5", "Mzb1", "MS4A1")
neutrophil_markers <- c("Ly6", "S100a8", "S100a9", "Mpo", "Elane", "Csf3r")

# Neuronal Markers
FeaturePlot(Aqp4, features = neuronal_markers, ncol = 2)
ggsave(filename = paste0(save_path, "FeaturePlot_NeuronalMarkers.png"), width = 10, height = 8, dpi = 300)
VlnPlot(Aqp4, features = neuronal_markers, group.by = "seurat_clusters")
ggsave(filename = paste0(save_path, "ViolinPlot_NeuronalMarkers.png"), width = 12, height = 8, dpi = 300)

# Oligodendrocyte Markers
FeaturePlot(Aqp4, features = oligodendrocyte_markers, ncol = 2)
ggsave(filename = paste0(save_path, "FeaturePlot_OligodendrocyteMarkers.png"), width = 8, height = 6, dpi = 300)
VlnPlot(Aqp4, features = oligodendrocyte_markers, group.by = "seurat_clusters")
ggsave(filename = paste0(save_path, "ViolinPlot_OligodendrocyteMarkers.png"), width = 10, height = 6, dpi = 300)

# Microglia Markers
FeaturePlot(Aqp4, features = microglia_markers, ncol = 2)
ggsave(filename = paste0(save_path, "FeaturePlot_MicrogliaMarkers.png"), width = 8, height = 6, dpi = 300)
VlnPlot(Aqp4, features = microglia_markers, group.by = "seurat_clusters")
ggsave(filename = paste0(save_path, "ViolinPlot_MicrogliaMarkers.png"), width = 10, height = 6, dpi = 300)

# Endothelial Markers
FeaturePlot(Aqp4, features = endothelial_markers, ncol = 2)
ggsave(filename = paste0(save_path, "FeaturePlot_EndothelialMarkers.png"), width = 8, height = 6, dpi = 300)
VlnPlot(Aqp4, features = endothelial_markers, group.by = "seurat_clusters")
ggsave(filename = paste0(save_path, "ViolinPlot_EndothelialMarkers.png"), width = 10, height = 6, dpi = 300)

# Pericyte Markers
VlnPlot(Aqp4, features = pericyte_markers, group.by = "seurat_clusters", ncol = 2)
ggsave(filename = paste0(save_path, "ViolinPlot_PericyteMarkers.png"), width = 10, height = 6, dpi = 300)

# Ependymal Markers
FeaturePlot(Aqp4, features = ependymal_markers)
ggsave(filename = paste0(save_path, "FeaturePlot_EpendymalMarkers.png"), width = 10, height = 6, dpi = 300)
VlnPlot(Aqp4, features = ependymal_markers, group.by = "seurat_clusters")
ggsave(filename = paste0(save_path, "ViolinPlot_EpendymalMarkers.png"), width = 10, height = 6, dpi = 300)

# Fibroblast Markers
FeaturePlot(Aqp4, features = fibroblast_markers, ncol = 2)
ggsave(filename = paste0(save_path, "FeaturePlot_FibroblastMarkers.png"), width = 8, height = 6, dpi = 300)
VlnPlot(Aqp4, features = fibroblast_markers, group.by = "seurat_clusters")
ggsave(filename = paste0(save_path, "ViolinPlot_FibroblastMarkers.png"), width = 10, height = 6, dpi = 300)

# T-cell Markers (Violin and FeaturePlot)
# Visualize the expression of T-cell markers across clusters using violin plots and feature plots.
VlnPlot(Aqp4, features = t_cell_markers, group.by = "seurat_clusters", ncol = 2)
ggsave(filename = paste0(save_path, "ViolinPlot_TcellMarkers.png"), width = 12, height = 8, dpi = 300)

# Feature plots for T-cell markers across the UMAP
FeaturePlot(Aqp4, features = t_cell_markers, cols = c("lightblue", "red"), reduction = "umap", ncol = 2)
ggsave(filename = paste0(save_path, "FeaturePlot_TcellMarkers.png"), width = 12, height = 8, dpi = 300)

# B-cell Markers (Violin and FeaturePlot)
# Visualize the expression of B-cell markers across clusters using violin plots and feature plots.
VlnPlot(Aqp4, features = b_cell_markers, group.by = "seurat_clusters", ncol = 2)
ggsave(filename = paste0(save_path, "ViolinPlot_BcellMarkers.png"), width = 14, height = 8, dpi = 300)

# Feature plots for B-cell markers across the UMAP
FeaturePlot(Aqp4, features = b_cell_markers, cols = c("lightblue", "red"), reduction = "umap", ncol = 2)
ggsave(filename = paste0(save_path, "FeaturePlot_BcellMarkers.png"), width = 14, height = 8, dpi = 300)

# Neutrophil Markers (Violin and FeaturePlot)
# Visualize the expression of Neutrophil markers across clusters using violin plots and feature plots.
VlnPlot(Aqp4, features = neutrophil_markers, group.by = "seurat_clusters", ncol = 2)
ggsave(filename = paste0(save_path, "ViolinPlot_NeutrophilMarkers.png"), width = 14, height = 8, dpi = 300)

# Feature plots for Neutrophil markers across the UMAP
FeaturePlot(Aqp4, features = neutrophil_markers, cols = c("lightblue", "red"), reduction = "umap", ncol = 2)
ggsave(filename = paste0(save_path, "FeaturePlot_NeutrophilMarkers.png"), width = 14, height = 8, dpi = 300)


# average expression of neuronal markers across clusters ------------
# Analyze Neuronal Marker Expression Across Clusters ---------

#   In this section, we calculate the average expression of neuronal markers across clusters 
#   and assign a "neuronal score" to each cell based on these markers. This helps in identifying 
#   clusters with higher neuronal activity and supports the decision on which clusters might be neuronal.

# Step 1: Calculate average expression of neuronal markers per cluster
# Get the average expression of neuronal markers across all clusters
avg_expression <- AverageExpression(Aqp4, features = neuronal_markers, group.by = "seurat_clusters")

# Print the average expression of the neuronal markers for review
print(avg_expression)

# Step 2: Add a neuronal score to each cell
# Assign a neuronal score to each cell based on the expression of the neuronal markers
Aqp4 <- AddModuleScore(Aqp4, features = list(neuronal_markers), name = "Neuronal_Score")

# Step 3: Visualize neuronal scores across clusters
# Use a violin plot to visualize the distribution of neuronal scores across clusters
cluster_neuronal_scores <- VlnPlot(Aqp4, features = "Neuronal_Score1", group.by = "seurat_clusters")

# Print the violin plot to check score distributions
print(cluster_neuronal_scores)

# Step 4: Identify clusters with neuronal markers
# Find all marker genes for each cluster
all_markers <- FindAllMarkers(Aqp4)

# Filter markers to focus on those that match neuronal markers
neuronal_markers_in_clusters <- all_markers %>%
  filter(gene %in% neuronal_markers)

# Step 5: List clusters expressing neuronal markers
# Extract the unique clusters that express neuronal markers
neuronal_clusters <- neuronal_markers_in_clusters %>%
  pull(cluster) %>%
  unique()

# Print the list of neuronal clusters for review
print(neuronal_clusters)

# Step 6: Calculate the average neuronal score per cluster
# Create a data frame of neuronal scores for each cell, grouped by cluster
neuronal_scores <- FetchData(Aqp4, vars = "Neuronal_Score1")
neuronal_scores$cluster <- Aqp4$seurat_clusters

# Calculate the average neuronal score per cluster
avg_neuronal_score_per_cluster <- aggregate(Neuronal_Score1 ~ cluster, data = neuronal_scores, FUN = mean)

# Print the average neuronal scores for each cluster
print(avg_neuronal_score_per_cluster)

# Identify Neuronal Clusters Based on Neuronal Marker Expression ---------

# In this part we determine which clusters are likely to be neuronal by examining the expression of known neuronal markers.
# This will help focus on neuronal clusters for further analyses.

# Step 1: Visualize neuronal marker expression across clusters
FeaturePlot(Aqp4, features = neuronal_markers, ncol = 2)
ggsave(filename = paste0(save_path, "FeaturePlot_NeuronalMarkers.png"), width = 10, height = 8, dpi = 300)

VlnPlot(Aqp4, features = neuronal_markers, group.by = "seurat_clusters")
ggsave(filename = paste0(save_path, "ViolinPlot_NeuronalMarkers.png"), width = 12, height = 8, dpi = 300)

# Step 2: Identify clusters expressing neuronal markers
all_markers <- FindAllMarkers(Aqp4, only.pos = TRUE)

# Filter for clusters that express neuronal markers
neuronal_clusters_info <- all_markers %>%
  filter(gene %in% neuronal_markers) %>%
  group_by(cluster) %>%
  summarize(avg_log2FC = mean(avg_log2FC)) %>%
  arrange(desc(avg_log2FC))

# Review clusters with high average expression of neuronal markers
print(neuronal_clusters_info)

# At this point, based on the visualizations and marker expression data,
# you can decide which clusters are likely to be neuronal for further analysis.


# Selection of neuronal clusters and saving  ---------

# After deciding which clusters are likely neuronal, filter those clusters and save them.
# Further filtering can be done afterwards to remove any remaining non-neuronal cells.

# Step 1: Filter the Seurat object to include only the selected neuronal clusters
# Replace neuronal_clusters with the list of clusters you've identified as neuronal.
neuronal_clusters <- c(0, 2, 3, 5, 6, 14, 15, 16, 19)  # Adjust based on your analysis.
Aqp4_neurons <- subset(Aqp4, idents = neuronal_clusters)

# Step 2: Reprocess and recluster the neuronal cells with the same parameters as before
Aqp4_neurons <- NormalizeData(Aqp4_neurons)
Aqp4_neurons <- FindVariableFeatures(Aqp4_neurons)
Aqp4_neurons <- ScaleData(Aqp4_neurons)
Aqp4_neurons <- RunPCA(Aqp4_neurons, npcs = 30)
Aqp4_neurons <- FindNeighbors(Aqp4_neurons, dims = 1:30)
Aqp4_neurons <- FindClusters(Aqp4_neurons, resolution = 0.1)
Aqp4_neurons <- RunUMAP(Aqp4_neurons, dims = 1:30)

# Step 3: Save the filtered and reclustered neuronal cells in a new RDS file
saveRDS(Aqp4_neurons, file = "neurons.rds")

# Step 4: Visualize the reclustered neuronal cells
DimPlot(Aqp4_neurons, reduction = "umap", label = TRUE)
ggsave(filename = paste0(save_path, "neuronal_clusters_umap.png"), 
       plot = last_plot(), width = 10, height = 8, units = "in", dpi = 300)

# Step 5: Optionally, perform additional filtering to remove any remaining non-neuronal cells
# You can use specific marker genes to confirm and remove unwanted cell types, e.g., based on expression of certain non-neuronal markers.
# Example: Filter based on low or absent expression of non-neuronal markers
Aqp4_neurons_filtered <- subset(Aqp4_neurons, subset = expression_level_of_non_neuronal_markers < threshold_value)

# Save the final filtered neuronal cells
saveRDS(Aqp4_neurons_filtered, file = "neurons_filtered.rds")

# Optional: Plot the final filtered neuronal cells
DimPlot(Aqp4_neurons_filtered, reduction = "umap", label = TRUE)
ggsave(filename = paste0(save_path, "neuronal_clusters_filtered_umap.png"), 
       plot = last_plot(), width = 10, height = 8, units = "in", dpi = 300)

# ------- Neuronal Subclustering: ----------
# Subclustering -----
# Subcluster with the same stats as before but now the neuronal clusters to identify subtypes 
# plot for different batches, days and conditions as done for Astrocytes 
