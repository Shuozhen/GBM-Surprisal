## This is scripts for zUMIs output to rds with tissue information etc, written by Shuozhen Bao
## This is for GBM NICHES project ALL TEN samples

library(Seurat)
library(magrittr)
library(cowplot)
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)

# rm(list=ls())

# Create the data frame for the whole GBM samples
mapping_df <- data.frame(
  SampleID = c("P1_Core", "P2_Core", "P3_Core", "P4_Core", "P5_Core",
               "P1_Peripheral", "P2_Peripheral", "P3_Peripheral", "P4_Peripheral", "P5_Peripheral"),
  ExperimentID = c("GBM1027C", "GBM1031C", "R1102C", "R1110C", "R0103N",
                   "GBM1027N", "GBM1031N", "R1102N", "R1110N", "R0103C"),
  TumorType = c("Recurrent", "Recurrent", "Primary", "Recurrent", "Primary",
                "Recurrent", "Recurrent", "Primary", "Recurrent", "Primary"),
  SampleType = c("Core", "Core", "Core", "Core", "Core",
                 "Peripheral", "Peripheral", "Peripheral", "Peripheral", "Peripheral")
)
# Print the data frame
print(mapping_df)

# Load the barcode file and create x.y.cord
barcodes50_file <- "/vast/palmer/pi/fan/sb2723/00.database/barcodes-AB.xls"
barcode50 <- read.table(barcodes50_file, header = FALSE, stringsAsFactors = FALSE)
colnames(barcode50) <- c("barcode", "y.cord", "x.cord")
barcode50$x.y.cord <- paste0(barcode50$y.cord, "x", barcode50$x.cord)
print(head(barcode50))

#### --------------- Variables start here

# ### Sample ID
# # Set the sample name
# samples <- list.files("/gpfs/gibbs/pi/fan/sb2723/GBM_story/RNA/00.rawdata/")
# sample <- "GBM220418A"

# Directories for input and output
dir_zUMIs_exp <- paste0("/gpfs/gibbs/pi/fan/sb2723/GBM_story/RNA/01.cutadapt_zUMIs_outs_20241204/", sample, "/zUMIs_output/expression/")
dir_out <- paste0("/gpfs/gibbs/pi/fan/sb2723/GBM_story/RNA/03.gen5_rds_20241206/", sample, "/")
dir_img <- paste0("/gpfs/gibbs/pi/fan/sb2723/GBM_story/RNA/02.10x_fmt_img/", sample, "/")
# Create output directory if it doesn't exist
if (!dir.exists(dir_out)) {dir.create(dir_out, recursive = TRUE)}

## Load Image information
img <- Read10X_Image(image.dir = paste0(dir_img, "spatial"), 
                     filter.matrix = TRUE)
image_barcodes <- img@boundaries$centroids@cells

# Check if the input experiment ID exists in the data frame
if (sample %in% mapping_df$ExperimentID) {
  # Get the corresponding SampleID
  sampleID <- mapping_df$SampleID[mapping_df$ExperimentID == sample]
  TumorType <- mapping_df$TumorType[mapping_df$ExperimentID == sample]
  SampleType <- mapping_df$SampleType[mapping_df$ExperimentID == sample]
  message <- paste("SampleID for", sample, "is", sampleID, ";")
  message <- paste("Tumor Type for", sample, "is", TumorType, ";")
  message <- paste("Sample Type for", sample, "is", SampleType, ";")
} else {
  message <- "The sample does not match any ExperimentID."
}
# Print the message
print(message)

# Load the zUMIs output data
all_files <- list.files(dir_zUMIs_exp)
rds_file <- all_files[grep("\\.rds$", all_files)]
exp_zUMIs_rds <- readRDS(paste0(dir_zUMIs_exp, rds_file))
exp_zUMIs <- as.matrix(exp_zUMIs_rds$umicount$inex$all)

### Generate Seurat Object with original matrix, no filtering, no spot selection
samp_zUMIs <-
  CreateSeuratObject(counts = exp_zUMIs,
                     project = sample,
                     assay = 'Spatial',
                     min.features = 1)

spatial_barcodes <- Cells(samp_zUMIs)
samp_zUMIs@meta.data$nucleotide_name <- spatial_barcodes
## Keep the original rownames to expriment_id
samp_zUMIs@meta.data$expriment_id <- sample
samp_zUMIs@meta.data$sample_id <- sampleID
samp_zUMIs@meta.data$TumorType <- TumorType
samp_zUMIs@meta.data$SampleType <- SampleType
samp_zUMIs@meta.data$MT_del <- "No"
common_barcodes <- intersect(spatial_barcodes,image_barcodes)
img_origin <- img[common_barcodes, ]
samp_zUMIs[[sample]] <- img_origin
names(samp_zUMIs@meta.data)

## Have to do the subset step, otherwise it will not delete the spots in the expression matrix
samp_zUMIs_ontissue <- subset(samp_zUMIs, cells = common_barcodes)
samp_zUMIs_ontissue[[sample]] <- img_origin
names(samp_zUMIs_ontissue@meta.data)

p1 <- VlnPlot(samp_zUMIs_ontissue, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
p2 <- SpatialFeaturePlot(samp_zUMIs_ontissue, features = "nCount_Spatial",pt.size.factor = 2,stroke = 0) + theme(legend.position = "right")
p3 <- VlnPlot(samp_zUMIs_ontissue, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
p4 <- SpatialFeaturePlot(samp_zUMIs_ontissue, features = "nFeature_Spatial",pt.size.factor = 2,stroke = 0) + theme(legend.position = "right")

png(filename=paste0(dir_out, sample, '_original_QC_on_tissue.png'), width = 8,height = 8,units = 'in',res=300)
plotA <- (p1 | p2)/ (p3 | p4)
plotA + plot_annotation(
  title = paste0(sample, '_original_QC_on_tissue'),
)
dev.off()

# saveRDS(samp_zUMIs, file = paste0(dir_out, sample, '_gen5_01_original.rds'))
saveRDS(samp_zUMIs_ontissue, file = paste0(dir_out, sample, '_gen5_01_original_on_tissue.rds'))

## Manual MT filtering becasue we are adding polyA tails using Patho-DBiT
df_gname_zUMIs_1 <- rownames(exp_zUMIs)
df_gname_zUMIs_filter_1 <- df_gname_zUMIs_1[!grepl('^(MT|AL|H3|EEF|RPS|RPL|FP|AC|ENSM|ENSG|RNA|RN5|RN7|RNU|RNV)', df_gname_zUMIs_1)] 
filtered_exp_zUMIs <- exp_zUMIs[df_gname_zUMIs_filter_1, ]
allZero_columns <- colSums(filtered_exp_zUMIs == 0) == nrow(filtered_exp_zUMIs)
filtered_exp_zUMIs <- filtered_exp_zUMIs[, !allZero_columns]

total_counts_original <- sum(exp_zUMIs)
total_counts_filtered <- sum(filtered_exp_zUMIs)
percentage_filtered <- (1-total_counts_filtered/total_counts_original)*100
print(paste0("MT gene percentage is ", round(percentage_filtered, 2), "%"))

samp_zUMIs_filtered <-
  CreateSeuratObject(counts = filtered_exp_zUMIs,
                     project = sample,
                     min.features = 1,
                     assay = 'Spatial')

### Same process for the filtered object
spatial_barcodes_filtered <- Cells(samp_zUMIs_filtered)
samp_zUMIs_filtered@meta.data$nucleotide_name <- spatial_barcodes_filtered
## Keep the original rownames to expriment_id
samp_zUMIs_filtered@meta.data$expriment_id <- sample
samp_zUMIs_filtered@meta.data$sample_id <- sampleID
samp_zUMIs_filtered@meta.data$TumorType <- TumorType
samp_zUMIs_filtered@meta.data$SampleType <- SampleType
samp_zUMIs_filtered@meta.data$MT_del <- "Yes"
common_barcodes_filtered <- intersect(spatial_barcodes_filtered,image_barcodes)
img_noMT <- img[common_barcodes_filtered, ]
samp_zUMIs_filtered[[sample]] <- img_noMT
names(samp_zUMIs_filtered@meta.data)

## Have to do the subset step, otherwise it will not delete the spots in the expression matrix
samp_zUMIs_filtered_ontissue <- subset(samp_zUMIs_filtered, cells = common_barcodes)
samp_zUMIs_filtered_ontissue[[sample]] <- img_origin
names(samp_zUMIs_filtered_ontissue@meta.data)

p5 <- VlnPlot(samp_zUMIs_filtered_ontissue, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
p6 <- SpatialFeaturePlot(samp_zUMIs_filtered_ontissue, features = "nCount_Spatial",pt.size.factor = 2,stroke = 0) + theme(legend.position = "right")
p7 <- VlnPlot(samp_zUMIs_filtered_ontissue, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
p8 <- SpatialFeaturePlot(samp_zUMIs_filtered_ontissue, features = "nFeature_Spatial",pt.size.factor = 2,stroke = 0) + theme(legend.position = "right")

png(filename=paste0(dir_out, sample, '_noMT_QC_on_tissue.png'), width = 8,height = 8,units = 'in',res=300)
plotB <- (p5 | p6)/(p7 | p8)
plotB + plot_annotation(
  title = paste0(sample, '_noMT_QC_on_tissue'),
)
dev.off()

# saveRDS(samp_zUMIs_filtered, file = paste0(dir_out, sample, '_gen5_02_noMT.rds'))
saveRDS(samp_zUMIs_filtered_ontissue, file = paste0(dir_out, sample, '_gen5_02_noMT_on_tissue.rds'))
