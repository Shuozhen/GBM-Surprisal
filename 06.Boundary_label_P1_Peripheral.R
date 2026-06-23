library(shiny)
library(ggplot2)
library(plotly)
library(dplyr)
library(Seurat)
library(patchwork)
library(colorspace)
library(ggpubr)
library(parallel)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(purrr)
library(janitor)
library(magrittr)
library(stringr)
library(R.utils)
library(layer)

samples <- readRDS("~/project/manual_cell_type_20250213/07.maunally_labeled_GLSB_v2_20250218.rds")


#P1 Peripheral
samp <- samples[[2]]
metadata <- samp@meta.data

#Create shiny app with lasso to define the 3 neuron regions
ui <- fluidPage(
  titlePanel("Free Draw to Select a Region"),
  sidebarLayout(
    sidebarPanel(
      h4("Selected Points"),
      verbatimTextOutput("selected_points") 
    ),
    mainPanel(
      plotlyOutput("spatial_plot")
    )
  )
)

# Define Server
server <- function(input, output, session) {
  output$spatial_plot <- renderPlotly({
    coords <- samp@meta.data[, c("x", "y")]
    clusters <- samp@meta.data$CellType_GLSB_manual
    
    gg <- ggplot(data = samp@meta.data, aes(x = x, y = y, color = clusters)) +
      geom_point(size = 1) +
      theme_minimal() +
      labs(title = "Spatial Sequencing", x = "X Coordinate", y = "Y Coordinate")

    ggplotly(gg) %>% layout(dragmode = "lasso")
  })
  
  output$selected_points <- renderPrint({
    event_data <- event_data("plotly_selected")
    
    if (!is.null(event_data)) {
      selected_data <- samp@meta.data %>%
        filter(x %in% event_data$x & y %in% event_data$y)
      
      if (nrow(selected_data) > 0) {
        return(selected_data[, c("x", "y", "seurat_clusters", "CellType_GLSB_manual")])
      } else {
        return("No points found in the selected region.")
      }
    } else {
      return("No region selected.")
    }
  })
}

# Run the App, then copy/paste the coordinates of interest
shinyApp(ui = ui, server = server)

selected_x_coords_neuron_healthy <- c(4059, 4176, 4294, 4412, 4529, 4647, 3943, 4178, 4296, 4413, 4531, 3826, 3943, 4061, 4178, 4296, 4413, 3942, 4060, 4178, 4295, 4413, 4530, 4648, 3942, 4060, 4177, 4295, 4412, 4530, 4647, 3942, 4059, 4177, 4294, 4412, 4529, 4647, 4765)
selected_y_coords_neuron_healthy <- c(6033, 6033, 6034, 6035, 6036, 6037, 6503, 6505, 6506, 6507, 6508, 6621, 6621, 6622, 6623, 6624, 6625, 6386, 6386, 6387, 6388, 6389, 6390, 6391, 6268, 6268, 6269, 6270, 6271, 6272, 6273, 6150, 6151, 6151, 6152, 6153, 6154, 6155, 6156)

selected_x_coords_neuron_1a <- c(2874, 2992, 3109, 2521, 2757, 2639, 2874, 2991, 3109, 3226, 3344, 2521, 3462, 2756, 2639, 2873, 2991, 3108, 3226, 3344, 2521, 3461, 2756, 2403, 2638, 2873, 2990, 3108, 3226, 3343, 2520, 3461, 2755, 2285, 2403, 2638, 2873, 2990, 3108, 3225, 3343, 2520, 2755, 2402, 2637, 2519)
selected_y_coords_neuron_1a <- c(3310, 3311, 3312, 3308, 3309, 3308, 3192, 3193, 3194, 3195, 3196, 3190, 3197, 3191, 3191, 3074, 3075, 3076, 3077, 3078, 3072, 3079, 3073, 3071, 3073, 2956, 2957, 2958, 2959, 2960, 2954, 2961, 2956, 2952, 2953, 2955, 2838, 2839, 2840, 2841, 2842, 2836, 2838, 2835, 2837, 2718)

selected_x_coords_neuron_1b <- c(4519, 4637, 4754, 4872, 4989, 4519, 4636, 4754, 4872, 5107, 5224, 4989, 4518, 4636, 4754, 4871, 5106, 5224, 4989, 4518, 4636, 4753, 4871, 5106, 5223, 5341, 4988, 4518, 4635, 4753, 4870, 5105, 4988, 4635, 4752, 4870)
selected_y_coords_neuron_1b <- c(3087, 3088, 3089, 3090, 3091, 2969, 2970, 2971, 2972, 2974, 2975, 2973, 2851, 2852, 2853, 2854, 2856, 2857, 2855, 2733, 2734, 2735, 2736, 2738, 2739, 2740, 2737, 2615, 2616, 2617, 2618, 2620, 2619, 2498, 2499, 2500)

#Create a dataframe with the coordinates of interest, for each region
selected_coords_neuron_healthy <- data.frame(
  selected_x_coords_neuron_healthy,  # Your x coordinates
  selected_y_coords_neuron_healthy   # Your y coordinates
)

selected_coords_neuron_1a <- data.frame(
  selected_x_coords_neuron_1a,  # Your x coordinates
  selected_y_coords_neuron_1a   # Your y coordinates
)

selected_coords_neuron_1b <- data.frame(
  selected_x_coords_neuron_1b,  # Your x coordinates
  selected_y_coords_neuron_1b   # Your y coordinates
)

selected_x_coords_invasive_1a <- c(2516, 2046, 1928, 2163, 2281, 2398, 2046, 2164, 2282, 2046, 2164, 2281, 2399, 2515, 2045, 1928, 2163, 2280, 2398, 2045, 1927, 2162, 2280, 2398, 2044, 1927, 2162, 2280)
selected_y_coords_invasive_1a <- c(1656, 1652, 1652, 1653, 1654, 1655, 1888, 1889, 1890, 1770, 1771, 1772, 1773, 1538, 1534, 1534, 1535, 1536, 1537, 1417, 1416, 1417, 1418, 1419, 1299, 1298, 1299, 1300)

selected_x_coords_invasive_1b <- c(6879, 6174, 6291, 6409, 6644, 6762, 6056, 6173, 6291, 6409, 6761, 6996, 6878, 6173, 6291, 6408, 6526, 6761, 6996, 6878, 6173, 6290, 6408, 6525, 6643, 6760, 6996, 6878, 6172, 6290, 6407, 6525, 6642, 6760, 6172, 6289, 6407, 6524, 6642, 6760, 6289, 6407, 6524, 6642)
selected_y_coords_invasive_1b<-  c(5701, 5695, 5696, 5697, 5699, 5700, 5576, 5577, 5578, 5579, 5582, 5584, 5465, 5459, 5460, 5461, 5462, 5464, 5466, 5347, 5341, 5342, 5343, 5344, 5345, 5346, 5348, 5229, 5223, 5224, 5225, 5226, 5227, 5228, 5105, 5106, 5107, 5108, 5109, 5110, 4988, 4989, 4990, 4991)

selected_x_coords_invasive_1c <- c(1589, 1471, 2059, 1942, 1706, 2177, 1588, 2059, 1941, 1706, 2176, 1588, 2058, 1941, 1706, 1588, 2058, 1940, 1705, 1590, 1472, 2060, 1943, 1708, 1589, 1472, 2059, 1942, 1707, 2177, 1590, 1472, 2060, 1942, 1707)
selected_y_coords_invasive_1c <- c(5659, 5658, 5663, 5662, 5660, 5664, 5541, 5545, 5544, 5542, 5546, 5423, 5427, 5426, 5424, 5305, 5309, 5308, 5306, 6013, 6012, 6017, 6016, 6014, 5777, 5776, 5781, 5780, 5778, 5782, 5895, 5894, 5899, 5898, 5896)

#Add metadata and then visualize to ensure correct regions were selected
samp@meta.data$neuron_gbm <- "Other"  # Set all to "Other" initially

# Update for the first group: "Neuron-Healthy"
for (i in 1:nrow(selected_coords_neuron_healthy)) {
  selected_x <- selected_coords_neuron_healthy[i, "selected_x_coords_neuron_healthy"]
  selected_y <- selected_coords_neuron_healthy[i, "selected_y_coords_neuron_healthy"]
  
  # Find matching rows in metadata
  match_idx <- which(samp@meta.data$x == selected_x & samp@meta.data$y == selected_y)
  
  # Update the neuron_new column
  samp@meta.data$neuron_gbm[match_idx] <- "Neuron-Healthy"
}

# Update for the second group: "Neuron-1A"
for (i in 1:nrow(selected_coords_neuron_1a)) {
  selected_x <- selected_coords_neuron_1a[i, "selected_x_coords_neuron_1a"]
  selected_y <- selected_coords_neuron_1a[i, "selected_y_coords_neuron_1a"]
  
  # Find matching rows in metadata
  match_idx <- which(samp@meta.data$x == selected_x & samp@meta.data$y == selected_y)
  
  # Update the neuron_new column
  samp@meta.data$neuron_gbm[match_idx] <- "Neuron-1A"
}

# Update for the third group: "Neuron-1B"
for (i in 1:nrow(selected_coords_neuron_1b)) {
  selected_x <- selected_coords_neuron_1b[i, "selected_x_coords_neuron_1b"]
  selected_y <- selected_coords_neuron_1b[i, "selected_y_coords_neuron_1b"]
  
  # Find matching rows in metadata
  match_idx <- which(samp@meta.data$x == selected_x & samp@meta.data$y == selected_y)
  
  # Update the neuron_new column
  samp@meta.data$neuron_gbm[match_idx] <- "Neuron-1B"
}

for (i in 1:nrow(selected_coords_gbm_1a)) {
  selected_x <- selected_coords_gbm_1a[i, "selected_x_coords_invasive_1a"]
  selected_y <- selected_coords_gbm_1a[i, "selected_y_coords_invasive_1a"]
  
  # Find matching rows in metadata
  match_idx <- which(samp@meta.data$x == selected_x & samp@meta.data$y == selected_y)
  
  # Update the neuron_new column
  samp@meta.data$neuron_gbm[match_idx] <- "GBM-1A"
}

for (i in 1:nrow(selected_coords_gbm_1b)) {
  selected_x <- selected_coords_gbm_1b[i, "selected_x_coords_invasive_1b"]
  selected_y <- selected_coords_gbm_1b[i, "selected_y_coords_invasive_1b"]
  
  # Find matching rows in metadata
  match_idx <- which(samp@meta.data$x == selected_x & samp@meta.data$y == selected_y)
  
  # Update the neuron_new column
  samp@meta.data$neuron_gbm[match_idx] <- "GBM-1B"
}

for (i in 1:nrow(selected_coords_gbm_1c)) {
  selected_x <- selected_coords_gbm_1c[i, "selected_x_coords_invasive_1c"]
  selected_y <- selected_coords_gbm_1c[i, "selected_y_coords_invasive_1c"]
  
  # Find matching rows in metadata
  match_idx <- which(samp@meta.data$x == selected_x & samp@meta.data$y == selected_y)
  
  # Update the neuron_new column
  samp@meta.data$neuron_gbm[match_idx] <- "GBM-1C"
}



# Visualize the spatial data
p1 <- SpatialDimPlot(
  samp, 
  group.by = "neuron_gbm", 
  cols = c(
    "Other" = "grey", 
    "Neuron-Healthy" = "green", 
    "Neuron-1A" = "blue", 
    "Neuron-1B" = "orange", 
    "GBM-1A" = "red", 
    "GBM-1B" = "red",
    "GBM-1C" = "red"
  )
)
p2 <- SpatialDimPlot(samp)
p1 + p2


#P1 Peripheral, removing empty channel
# 1) Convert coords to data frame
orig_coords <- samp@images$GBM1027N@boundaries$centroids@coords
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
bins_to_remove_col <- c(4)
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
samp@images$GBM1027N@boundaries$centroids@coords <- coords_new

SpatialDimPlot(samp)


#Save metadata labels of selected regions 
write.csv(
  data.frame(neuron_gbm = samp@meta.data$neuron_gbm), 
  file = "neuron_gbm_metadata.csv", 
  row.names = TRUE
)

