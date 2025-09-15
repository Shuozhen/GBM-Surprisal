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
library(plotly)

samp <- readRDS("~/project/Transfer/P2_Core_20250414_NICHES.rds")

# Load the CSV file and preserve row names
neuron_new_data <- read.csv("neuron_new_metadata_P2_core.csv", row.names = 1)
samp@meta.data$neuron_new <- neuron_new_data$neuron_new

DefaultAssay(samp) <- "Spatial"
SpatialDimPlot(
  samp, 
  group.by = "neuron_new", 
  cols = c(
    "Other" = "lightgrey", 
    "Neuron-Healthy" = "green4", 
    "Neuron-1A" = "purple3", 
    "Neuron-1B" = "deeppink", 
    "Neuron-1C" = "darkorange",
    "GBM" = "red"
  )
  , image.alpha = 0, pt.size.factor = 2, shape =22)

p2 <- SpatialDimPlot(samp, group.by = "CellType_GLSB_manual")
p1 +p2

Idents(samp) <- "neuron_new"
mark <- FindMarkers(samp, ident.1 = c("Neuron-1A", "Neuron-1B", "Neuron-1C"), ident.2 = "Neuron-Healthy")
mark$ratio <- mark$pct.1/mark$pct.2

DefaultAssay(samp) <- "NeighborhoodToCell_SCT_9"
markers <- FindMarkers(samp, ident.1 = c("Neuron-1A", "Neuron-1B", "Neuron-1C"), ident.2 = "Neuron-Healthy")
markers$ratio <- markers$pct.1/markers$pct.2

#Make volcano plot 
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
genes$labels<- factor(genes$genes, levels = c("NLGN1—NRXN1", "MATN1—PTPRD", "NXPH2—NRXN1", "EFNA5—EPHA6"))

png(filename='volcano.pdf', width=3500, height=2000, res=300)
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


SpatialFeaturePlot(samp, features = c("NLGN1—NRXN1"), image.alpha = 0, shape = 22, pt.size.factor =2)
SpatialFeaturePlot(samp, features = c("EFNA5—EPHA6"), image.alpha = 0, shape = 22, pt.size.factor =2)
SpatialFeaturePlot(samp, features = c("SEMA5A—MET"), image.alpha = 0, shape = 22, pt.size.factor =2)

SpatialFeaturePlot(samp, features = c("MAP2"), image.alpha = 0, shape = 22, pt.size.factor =2, min.cutoff = 1.5)
SpatialFeaturePlot(samp, features = c("CALM1"), image.alpha = 0, shape = 22, pt.size.factor =2, min.cutoff = 2)
SpatialFeaturePlot(samp, features = c("RBFOX3"), image.alpha = 0, shape = 22, pt.size.factor =2)


saveRDS(samp, file = "P2_Core_20250414_NICHES.rds")
