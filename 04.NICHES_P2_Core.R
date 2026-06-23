library(dplyr)
library(scales)
library(NICHES)
library(stringr)
library(Seurat)
library(SeuratObject)

samp <- readRDS("/gpfs/gibbs/project/fan/gil5/Transfer/P2_Core_20250409.rds")

NICHES_output <- NICHES::RunNICHES(samp,
                                   assay = 'Spatial', 
                                   LR.database = 'fantom5',
                                   species = 'human',
                                   position.x = 'x',
                                   position.y = 'y',
                                   k = 9, # I am choosing to study the edges between every cell and its 10 nearest neighbors in histologic (real) space
                                   rad.set = NULL,
                                   cell_types = "CellType_GLSB_manual",
                                   CellToCellSpatial = T,
                                   CellToNeighborhood = T,
                                   NeighborhoodToCell = T,
                                   CellToCell = F,
                                   SystemToCell = F,
                                   CellToSystem = F)


me.data <- GetAssayData(object = NICHES_output$NeighborhoodToCell, layer  = 'data')
colnames(me.data) <- NICHES_output$NeighborhoodToCell$ReceivingCell
me.assay <- CreateAssayObject(data = me.data)
samp[['NeighborhoodToCell_SCT_9']] <- me.assay

DefaultAssay(samp) <- "NeighborhoodToCell_SCT_9"

samp <- ScaleData(samp)
samp <- FindVariableFeatures(samp,selection.method = "disp")
samp <- RunPCA(samp)
ElbowPlot(samp,ndims = 50)
samp <- RunUMAP(samp,dims = 1:20)
samp <- FindNeighbors(samp,dims = 1:15)


#want metadata that keeps highlighted cells as interface, and other cells as original cell type from CellType_GLSB_manua
samp@meta.data$highlight <- as.character(samp@meta.data$highlight)
samp@meta.data$CellType_GLSB_manual <- as.character(samp@meta.data$CellType_GLSB_manual)

# Create the new column 'Interface'
samp@meta.data$Interface <- ifelse(
  samp@meta.data$highlight == "Neuron-GBM", 
  "interface", 
  samp@meta.data$CellType_GLSB_manual
)

SpatialDimPlot(samp, group.by = "Interface") 
DefaultAssay(samp) <- "NeighborhoodToCell_SCT_9"

samp <- FindClusters(samp, resolution = 0.3)
DimPlot(samp)
SpatialDimPlot(samp)
SpatialDimPlot(samp, group.by = "Interface") + SpatialDimPlot(samp)

marker_all <- FindAllMarkers(samp,logfc.threshold = 0.25, min.pct = 0.1, group.by = "Interface")
marker_all$ratio <- marker_all$pct.1/marker_all$pct.2
marker_all$power <- marker_all$ratio * marker_all$avg_log2FC



saveRDS(samp, file = "P2_Core_20250414_NICHES.rds")



DefaultAssay(samp) <- "Spatial"
p1 <- SpatialFeaturePlot(samp, "MAP2")
DefaultAssay(samp) <- "NeighborhoodToCell_SCT_9"
#p2 <- SpatialFeaturePlot(samp,"VIP—ADCYAP1R1")
p2 <- SpatialFeaturePlot(samp,"CALM2—GRM7")
p1 + p2

ph_MAP2 <- GetAssayData(samp, assay = "Spatial", layer = "data")["MAP2", ]
NtC_CALM2_GRM7 <- GetAssayData(samp, assay = "NeighborhoodToCell_SCT_9", slot = "data")["CALM2—GRM7", ]



nonzero_idx <- NtC_CALM2_GRM7 != 0

# Subset both vectors using that index
ph_MAP2_filtered <- ph_MAP2[nonzero_idx]
NtC_CALM2_GRM7_filtered <- NtC_CALM2_GRM7[nonzero_idx]
plot(ph_MAP2_filtered~NtC_CALM2_GRM7_filtered)
common_nonzero_indices <- (ph_MAP2 != 0) & (NtC_CALM2_GRM7 != 0)





#Volcano plot 
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
data <- read.table("03.25.findmarkers_neuron_v_other.txt", header = TRUE,  sep = ",")
EnhancedVolcano(
  data,
  lab = data$X,
  x = 'ratio',
  y = 'p_val',
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  xlim = c(0, 3),           # Set x-axis limits
  ylim = c(0, 7),   
  pCutoff = 0.01,           # New p-value threshold (e.g., 0.01 for stricter significance)
  FCcutoff = 1.5, 
  selectLab = c("GeneA", "GeneB", "GeneC")# Set y-axis limits
)


markers <- FindAllMarkers(samp)
markers$ratio <-markers$pct.1/markers$pct.2

SpatialFeaturePlot(samp, features = c("PIP—NPTN", "ICAM1—EGFR", "VEGFB—NRP1", "WNT5A—ROR2", "VEGFA—ITGA9", "SEMA5A—MET"))
SpatialFeaturePlot(samp, features = c("VIP—SCTR", "CXCL13—HTR2A", "HGF—MET", "SEMA5A—MET", "SEMA4D—MET", "IL16—KCNA3"))
SpatialFeaturePlot(samp, features = c("THBS2—CD47", "PON2—HTR2A", "NLGN3—NRXN1", "NLGN3—NRXN3", "NLGN3—NRXN2", "BMP5—BMPR2"))


#use iinterface to see finallmarkers for interface region not subset



p1 <- DimPlot(samp, reduction = "umap",group.by = 'seurat_clusters', label = TRUE)
p2 <- SpatialDimPlot(samp, label = TRUE,group.by = 'active.assay', label.size = 3)
p1 + p2


p1 <- DimPlot(samp, reduction = "umap",group.by = 'NeighborhoodToCell_SCT_9', label = TRUE)






# Load necessary libraries
library(ggplot2)
library(dplyr)

# Example data
colnames(mark)

mark$gene <- rownames(mark)
highlight_genes <- c("NOTCH4-5", "IFIT1", "KDM5D", "MET", "KCND2")
mark_positive <- mark[mark$avg_log2FC > 0, ]

mark$highlight <- ifelse(mark$gene %in% highlight_genes, "highlight", "other")
mark$rank <- rank(-mark$ratio, ties.method = "first")
mark_positive <- mark[mark$avg_log2FC > 0, ]


ggplot(mark, aes(x = rank, y = avg_log2FC)) +
  geom_point(aes(color = highlight), size = 1) +  # Color points based on 'highlight' column
  scale_color_manual(values = c("highlight" = "red", "other" = "grey")) +  # Highlight genes in red
  geom_text(data = subset(mark, highlight == "highlight"), aes(label = gene), vjust = -1, size = 3) +  # Label highlighted genes
  labs( 
    x = "Rank (Based on Highest Ratio)", y = "Log2 Fold Change") +
  theme_minimal() +
  theme(legend.position = "none")

