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

samples <- readRDS("~/project/manual_cell_type_20250213/07.maunally_labeled_GLSB_v2_20250218.rds")

######## P5_Core
samp <- samples[[9]]
###################
# 1) Convert coords to data frame
orig_coords <- samp@images$R0103N@boundaries$centroids@coords
coords_df <- data.frame(x = orig_coords[,1], y = orig_coords[,2])
rownames(coords_df) <- colnames(samp) 

# 2) Bin into 50 rows, 50 cols
num_bins <- 50
x_min <- min(coords_df$x); x_max <- max(coords_df$x)
y_min <- min(coords_df$y); y_max <- max(coords_df$y)

x_breaks <- seq(x_min, x_max, length.out = num_bins + 1)
y_breaks <- seq(y_min, y_max, length.out = num_bins + 1)

bin_width  <- (x_max - x_min)/num_bins
bin_height <- (y_max - y_min)/num_bins

coords_df$col_bin <- cut(coords_df$x, breaks = x_breaks, labels = FALSE, include.lowest = TRUE)
coords_df$row_bin <- cut(coords_df$y, breaks = y_breaks, labels = FALSE, include.lowest = TRUE)

# 3) Identify which row/col bins are empty
SpatialDimPlot(samp)

bins_to_remove_row <- NULL
bins_to_remove_col <- c(4, 47)
all_row_bins <- 1:50
all_col_bins <- 1:50

kept_row_bins <- setdiff(all_row_bins, bins_to_remove_row)
kept_col_bins <- setdiff(all_col_bins, bins_to_remove_col)

# 4) Filter out empty bins & re-label

coords_df <- coords_df %>%
  filter(row_bin %in% kept_row_bins,
         col_bin %in% kept_col_bins) %>%
  mutate(
    row_bin_new = match(row_bin, sort(unique(row_bin))),
    col_bin_new = match(col_bin, sort(unique(col_bin)))
  )

# 5) Shift the original coords to remove gaps
coords_df <- coords_df %>%
  mutate(
    bins_removed_below = sapply(row_bin, function(rb) sum(bins_to_remove_row < rb)),
    bins_removed_left  = sapply(col_bin, function(cb) sum(bins_to_remove_col < cb)),
    
    # Shift continuous coords properly
    x_shifted = x - bins_removed_left * bin_width,
    y_shifted = y - bins_removed_below * bin_height
  )

# 6) Save your new gap-free coords
coords_new <- as.matrix(coords_df[, c("x_shifted", "y_shifted")])
kept_spots <- rownames(coords_new)  # The barcodes that remain after filtering
samp <- subset(samp, cells = kept_spots)
samp@images$R0103N@boundaries$centroids@coords <- coords_new

SpatialDimPlot(samp)



saveRDS(samp, file = "P5_Core_20250414.rds")
