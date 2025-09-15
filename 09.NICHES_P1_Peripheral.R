library(patchwork)
library(colorspace)
library(ggpubr)
library(parallel)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(shiny)
library(tidyverse)
library(purrr)
library(janitor)
library(magrittr)
library(patchwork)
library(stringr)
library(R.utils)
library(layer)
library(dplyr)

#LoadP1_Peripheral_20250414_NICHES.rds
#Load the saved metadata (CSV method)
samp <- readRDS("~/project/Transfer/P1_Peripheral_20250414_NICHES.rds")

#Load the selected regions 
neuron_gbm <- read.csv("neuron_gbm_metadata.csv", row.names = 1)
samp@meta.data$neuron_gbm <- neuron_gbm$neuron_gbm

DefaultAssay(samp) <- "Spatial"
p1 <- SpatialDimPlot(samp, group.by = "neuron_gbm", cols = c( "Other" = "grey", 
                                                              "Neuron-Healthy" = "lightgreen", 
                                                              "Neuron-1A" = "plum3", 
                                                              "Neuron-1B" = "plum1", 
                                                              "GBM-1A" = "red2", 
                                                              "GBM-1B" = "yellow", 
                                                              "GBM-1C" = "blue"), image.alpha = 0, shape = 22)

p2 <- SpatialDimPlot(samp, group.by= "CellType_GLSB_manual")
p3 <- SpatialDimPlot(samp, group.by= "highlight")

Idents(samp) <- "neuron_gbm"
mark_1 <- FindMarkers(samp, ident.1 = c("GBM-1A", "GBM-1B", "GBM-1C"), ident.2 = "Neuron-Healthy")
mark_1$ratio <- mark_1$pct.1/mark_1$pct.2
gene_of_interest_1 <- rownames(mark_1[mark_1$pct.1 >= 0.3 & mark_1$pct.2 >= 0.3 & mark_1$p_val_adj < 0.05 & (mark_1$avg_log2FC > 1.5 | mark_1$avg_log2FC < -1.5), ])

mark_2 <- FindMarkers(samp, ident.1 = c("GBM-1A", "GBM-1B", "GBM-1C"), ident.2 = c("Neuron-1A", "Neuron-1B"))
mark_2$ratio <- mark_2$pct.1/mark_2$pct.2
gene_of_interest_2 <- rownames(mark_2[mark_2$pct.1 >= 0.3 & mark_2$pct.2 >= 0.3 & mark_2$p_val_adj < 0.05 & (mark_2$avg_log2FC > 1.5 | mark_2$avg_log2FC < -1.5), ])


mark_3 <- FindMarkers(samp, ident.1 = c("GBM-1A"), ident.2 = c("GBM-1B", "GBM-1C"))
mark_3$ratio <- mark_3$pct.1/mark_3$pct.2
gene_of_interest_3 <- rownames(mark_3[mark_3$pct.1 >= 0.3 & mark_3$pct.2 >= 0.3 & mark_3$p_val_adj < 0.05 & (mark_3$avg_log2FC > 1.5 | mark_3$avg_log2FC < -1.5), ])


DefaultAssay(samp) <- "NeighborhoodToCell_SCT_9"
markers_1 <- FindMarkers(samp, ident.1 = c("GBM-1A", "GBM-1B", "GBM-1C"), ident.2 = "Neuron-Healthy")
markers_1$ratio <- markers_1$pct.1/markers_1$pct.2
LRM_of_interest_1 <- rownames(markers_1[markers_1$pct.1 >= 0.3 & markers_1$pct.2 >= 0.3 & markers_1$p_val_adj < 0.05 & (markers_1$avg_log2FC > 1.5 | markers_1$avg_log2FC < -1.5), ])


DefaultAssay(samp) <- "NeighborhoodToCell_SCT_9"
markers_2 <- FindMarkers(samp, ident.1 = c("GBM-1A", "GBM-1B", "GBM-1C"), ident.2 = c("Neuron-1A", "Neuron-1B"))
markers_2$ratio <- markers_2$pct.1/markers_2$pct.2
LRM_of_interest_2 <- rownames(markers_2[markers_2$pct.1 >= 0.3 & markers_2$pct.2 >= 0.3 & markers_2$p_val_adj < 0.05 & (markers_2$avg_log2FC > 1.5 | markers_2$avg_log2FC < -1.5), ])



DefaultAssay(samp) <- "NeighborhoodToCell_SCT_9"
markers_3 <- FindMarkers(samp, ident.1 = c("GBM-1A"), ident.2 = c("GBM-1B", "GBM-1C"))
markers_3$ratio <- markers_3$pct.1/markers_3$pct.2
LRM_of_interest_3 <- rownames(markers_3[markers_3$pct.1 >= 0.3 & markers_3$pct.2 >= 0.3 & markers_3$p_val_adj < 0.05 & (markers_3$avg_log2FC > 1.5 | markers_3$avg_log2FC < -1.5), ])

setwd("./project/Transfer/P1_Peripheral")
dir_out <- "plots/"
dir.create(dir_out)
setwd(dir_out)
for (i in c(LRM_of_interest_1, LRM_of_interest_2, LRM_of_interest_3)){
  p <- SpatialFeaturePlot(samp, features = i, image.alpha = 0, shape = 22, pt.size.factor = 2)
  pdf(paste0("spatialfeatureplot_", i, ".pdf"))
  print(p)
  dev.off()
}

DefaultAssay(samp) <- "Spatial"
for (i in c(gene_of_interest_1, gene_of_interest_2, gene_of_interest_3)){
  p <- SpatialFeaturePlot(samp, features = i, image.alpha = 0, shape = 22, pt.size.factor = 2)
  pdf(paste0("spatialfeatureplot_", i, ".pdf"))
  print(p)
  dev.off()
}

write.csv(gene_of_interest_1, "gene_GBMvsNeuronHealthy.csv")
write.csv(gene_of_interest_2, "gene_GBMvsNeuronUnhealthy.csv")


#Create volcano plot
genes <- markers
#genes$genes <- genes$Gene
genes$genes <- rownames(genes)
genes$logFC <- genes$avg_log2FC 
genes$PValue <- genes$p_val_adj

AdjustedCutoff <- 0.05
LabellingCutoff <- 0.05
FCCutoff <- 1.5

genes$Significance <- "Not significant"
genes$Significance[(abs(genes$logFC) > FCCutoff)] <- "Significant fold change"
genes$Significance[(genes$PValue<AdjustedCutoff)] <- "Sigificant p-value"
genes$Significance[(genes$PValue<AdjustedCutoff) & (abs(genes$logFC)>FCCutoff)] <- "Signficant fold change & p-value"
table(genes$Significance)
genes$Significance <- factor(genes$Significance, levels=c("Not significant", 
                                                          "Significant fold change", 
                                                          "Sigificant p-value", 
                                                          "Signficant fold change & p-value"))

#genes$genes <- factor(genes$genes, levels = c(as.character(genes$genes[1:10]))) #label top 10 genes
genes$labels<- factor(genes$genes, levels = c("NLGN1—NRXN1", "SLIT3—ROBO2", "SEMA5A—MET", "HGF—MET"))

pdf(file = "volcano_plot.pdf", width = 12, height = 7) 
ggplot(genes, aes(x = logFC, y = -log10(PValue))) +
  geom_point(aes(color = Significance)) +
  scale_color_manual(values = c("black", "orange", "blue", "red")) +
  theme_bw(base_size = 10) + 
  theme(legend.title = element_blank(),
        text=element_text(size=18)) +
  xlab("log(fold change)") +
  #xlim(0, 2) + 
  ylab("-log10(p-value)") + 
  geom_vline(xintercept=0.5, linetype='dotted') + 
  geom_vline(xintercept=-.5, linetype='dotted') + 
  geom_hline(yintercept=1.3, linetype='dotted') + 
  geom_text_repel(data = genes,
                  na.rm = TRUE,
                  aes(label = labels),
                  size = 3,
                  box.padding = unit(0.7, "lines"))
dev.off()




