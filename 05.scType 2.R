#Sc-type
# load libraries and functions
library(dplyr)
library(Seurat)
library(HGNChelper)
library(openxlsx)
library(ggplot2)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# Set work directory 
setwd("/gpfs/gibbs/pi/fan/sb2723/GBM_story/RNA/06.gen5_scType")

# Define the path to the parent directory containing the sample subfolders
dir_in <- "/gpfs/gibbs/pi/fan/sb2723/GBM_story/RNA/05.rds_gen5_NICHES_20241217/"

# List all files that match the pattern '_gen4_02_noMT_on_tissue.rds' recursively in subfolders
rds_files <- list.files(dir_in, pattern = "_gen5_05_noMT_on_tissue_NICHES.rds$", recursive = TRUE, full.names = TRUE)
rds_list <- lapply(rds_files, readRDS)

names(rds_list) <- basename(dirname(rds_files))

stats_table <- sapply(names(rds_list), function(sample_name) {
  seurat_obj <- rds_list[[sample_name]]  # Get the Seurat object for this sample
  
  # Calculate the mean for both nCount_Spatial and nFeature_Spatial
  mean_nCount <- mean(seurat_obj@meta.data$nCount_Spatial, na.rm = TRUE)  # Mean of nCount_Spatial
  mean_nFeature <- mean(seurat_obj@meta.data$nFeature_Spatial, na.rm = TRUE)  # Mean of nFeature_Spatial
  
  # Return both means
  c(mean_nCount, mean_nFeature)
})

# Convert the results into a data frame and label the columns
stats_table <- as.data.frame(t(stats_table))
colnames(stats_table) <- c("Mean_nCount_Spatial", "Mean_nFeature_Spatial")


# Prepare gene sets
current_gs_list <- read.xlsx("mega_data_set.xlsx")

gs_positive <- current_gs_list
gs_negative <- lapply(gs_positive, function(x) {
  setNames(character(0), names(x))
})
gs_list_DB <- gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

names(gs_negative)[names(gs_negative) == "Endothel"] <- "Endothelial"
names(gs_positive)[names(gs_positive) == "Endothel"] <- "Endothelial"
gs_negative$Macrophage <- 'NCAM1'
gs_list <- list(gs_positive = gs_positive, gs_negative = gs_negative)

# Loop through all tumors
#rds_list[[1]]@active.assay should be SCT as active assay 
for (tumor in names(rds_list)) {
  # Assume 'sample' is the Seurat object corresponding to the current tumor
  sample <- rds_list[[tumor]]
  scRNAseqData_scaled <-  as.matrix(sample[["SCT"]]@scale.data)
  
  # Run ScType
  es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  
  # Merge by cluster
  cl_results <- do.call("rbind", lapply(unique(sample@meta.data$seurat_clusters), function(cl) {
    es.max.cl = sort(rowSums(es.max[, rownames(sample@meta.data[sample@meta.data$seurat_clusters == cl, ])]), decreasing = TRUE)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sample@meta.data$seurat_clusters == cl)), 10)
  }))
  
  sctype_scores <- cl_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
  
  # Adjust the cutoff for "unknown" classification
  cutoff_proportion <- 0.1  # Change this value to experiment with different cutoffs
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells * cutoff_proportion] <- "Unknown"
  print(sctype_scores[, 1:3])
  
  sample@meta.data$scType_SCT <- ""
  for (j in unique(sctype_scores$cluster)) {
    cl_type <- sctype_scores[sctype_scores$cluster == j, ]
    sample@meta.data$scType_SCT[sample@meta.data$seurat_clusters == j] <- as.character(cl_type$type[1])
  }
  
  # Check for unique cell types in the sample
  unique_cell_types <- unique(sample@meta.data$scType_SCT)
  print(paste("Unique cell types in", tumor, ":", paste(unique_cell_types, collapse = ", ")))
  rds_list[[tumor]] <- sample
}

saveRDS(rds_list, "gen5_rdslist_06_scType.rds")
