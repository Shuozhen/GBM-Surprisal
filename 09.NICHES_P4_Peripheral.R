library(dplyr)
library(scales)
library(NICHES)
library(stringr)
library(Seurat)
library(SeuratObject)

samp <- readRDS("/gpfs/gibbs/project/fan/gil5/Transfer/P4_Peripheral_20250414.rds")

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

samp <- FindClusters(samp, resolution = 0.1)
DimPlot(samp)
SpatialDimPlot(samp)
SpatialDimPlot(samp, group.by = "CellType_GLSB_manual") + SpatialDimPlot(samp)

marker_all <- FindAllMarkers(samp,logfc.threshold = 0.25, min.pct = 0.1)
marker_all$ratio <- marker_all$pct.1/marker_all$pct.2
marker_all$power <- marker_all$ratio * marker_all$avg_log2FC



saveRDS(samp, file = "P4_Peripheral_20250414_NICHES.rds")
