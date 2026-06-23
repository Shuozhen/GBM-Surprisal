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

samples <- readRDS("~/project/manual_cell_type_20250213/07.maunally_labeled_GLSB_v2_20250218.rds")


#P2 Core
samp <- samples[[3]]
metadata <- samp@meta.data


# Define UI for Shiny app with lasso tool 
ui <- fluidPage(
  titlePanel("Free Draw to Select a Region"),
  sidebarLayout(
    sidebarPanel(
      h4("Selected Points"),
      verbatimTextOutput("selected_points") # Display selected metadata rows
    ),
    mainPanel(
      plotlyOutput("spatial_plot") # Enable interactive selection
    )
  )
)

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
        filter(
          x %in% event_data$x & 
            y %in% event_data$y & 
            CellType_GLSB_manual == "Neuron"
        )
            if (nrow(selected_data) > 0) {
        return(selected_data[, c("x", "y", "seurat_clusters", "CellType_GLSB_manual")])
      } else {
        return("No Neuron points found in the selected region.")
      }
    } else {
      return("No region selected.")
    }
  })
}

# Run the app
shinyApp(ui, server)

#Document selected coordinates 
selected_x_coords_neuron_healthy <- c(4733, 4894, 5055, 5216, 5377, 5538, 5698, 5859, 6020, 6181, 6342, 6503, 6664, 
                                      5860, 6021, 6182, 6343, 6504, 5377, 5538, 5699, 5860, 6021, 6182, 6343, 6504, 
                                      4733, 4894, 5055, 5216, 5377, 5538, 5699, 5860, 6021, 6182, 6343, 6503, 
                                      4733, 4894, 5055, 5216, 5377, 5538, 5699, 5860, 6021, 6181, 6342, 6503, 6664, 
                                      4733, 4894, 5055, 5215, 5376, 5537, 5698, 5859, 6020, 6181, 6342, 6503, 6664, 
                                      4732, 4893, 5054, 5215, 5376, 5537, 5698, 5859, 6020, 6181, 6342, 6503, 
                                      4893, 5054, 5215, 5376, 5537, 5698, 5859, 6020, 6181, 6342)
selected_y_coords_neuron_healthy <- c(3923, 3922, 3922, 3922, 3921, 3921, 3920, 3920, 3920, 3919, 3919, 3919, 3918, 
                                      4556, 4556, 4556, 4555, 4555, 4398, 4398, 4398, 4397, 4397, 4396, 4396, 4396, 
                                      4241, 4240, 4240, 4240, 4239, 4239, 4239, 4238, 4238, 4237, 4237, 4237, 
                                      4082, 4081, 4081, 4081, 4080, 4080, 4080, 4079, 4079, 4078, 4078, 4078, 4077, 
                                      3764, 3763, 3763, 3763, 3762, 3762, 3761, 3761, 3761, 3760, 3760, 3760, 3759, 
                                      3605, 3604, 3604, 3604, 3603, 3603, 3602, 3602, 3602, 3601, 3601, 3601, 
                                      3445, 3445, 3444, 3444, 3444, 3443, 3443, 3443, 3442, 3442)


selected_x_coords_neuron_1a <- c(6669, 5864, 6186, 6990, 6025, 6186, 6347, 6668, 6346, 6507, 6829, 6990, 
                                 6346, 6507, 6829, 7151, 7312, 6668, 6185, 6507, 6829, 6990, 7151, 6668, 
                                 6507, 6829, 6668)
selected_y_coords_neuron_1a <- c(7894, 7737, 7736, 7735, 7578, 7577, 7577, 7576, 7418, 7418, 7417, 7416, 
                                 7259, 7259, 7258, 7257, 7257, 7258, 7100, 7099, 7099, 7098, 7098, 7099, 
                                 6940, 6940, 6940)

selected_x_coords_neuron_1b <- c(4253, 3931, 4092, 4253, 3609, 4092, 3770, 3930, 4091, 4252, 3609, 3769, 
                                 3930, 4091, 4252, 3608, 3769, 3930, 4091, 4252, 4413, 3286, 3608, 3769, 
                                 3930, 4252, 3769, 3930, 4091, 3608, 4253)
selected_y_coords_neuron_1b <- c(6310, 6151, 6151, 6151, 6152, 5992, 5834, 5833, 5833, 5832, 5834, 5675, 
                                 5674, 5674, 5673, 5675, 5515, 5515, 5515, 5514, 5514, 5517, 5516, 5356, 
                                 5356, 5355, 5197, 5197, 5197, 5198, 6469)

selected_x_coords_neuron_1c <- c(9241, 8436, 8758, 9080, 8114, 8275, 8597, 8919, 9240, 7953, 8275, 8436, 
                                 8597, 9080, 9240, 7953, 8114, 8436, 8597, 9079, 7953, 8114, 8435, 8596, 
                                 8757, 9079, 7952, 9079, 7952, 8274, 8596, 8757, 8918, 9079)
selected_y_coords_neuron_1c <- c(5662, 5664, 5663, 5662, 5505, 5505, 5504, 5503, 5185, 5188, 5187, 5186, 
                                 5186, 5185, 5025, 5029, 5028, 5027, 5027, 5026, 4869, 4869, 4868, 4868, 
                                 4868, 4867, 4710, 4708, 4551, 4551, 4550, 4549, 4549, 4549)


selected_x_coords_invasive_1a <- c(5385, 5546, 5707, 5868, 6029, 6190, 6351, 6512, 6833, 6994, 7155, 6673, 5224, 5385, 5546, 5707, 5868, 6029, 6190, 6351, 6512, 6834, 6995, 7156, 6673, 5867, 6028, 6189, 6350, 5867, 6028, 6189, 6350, 6511, 6833, 6672, 5707, 5868, 6028, 6189, 6350, 6511, 6833, 6672, 5385, 5546, 5707, 5868, 6029, 6190, 6350, 6511, 6833, 6994, 7155, 6672)
selected_y_coords_invasive_1a <- c(11078, 11078, 11077, 11077, 11077, 11076, 11076, 11076, 11075, 11074, 11074, 11075, 11238, 11237, 11237, 11237, 11236, 11236, 11235, 11235, 11235, 11234, 11234, 11233, 11234, 10441, 10441, 10440, 10440, 10600, 10600, 10599, 10599, 10598, 10598, 10598, 10759, 10759, 10759, 10758, 10758, 10758, 10757, 10757, 10919, 10919, 10918, 10918, 10918, 10917, 10917, 10917, 10916, 10915, 10915, 10916)


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

selected_coords_neuron_1c <- data.frame(
  selected_x_coords_neuron_1c,  # Your x coordinates
  selected_y_coords_neuron_1c   # Your y coordinates
)

selected_coords_gbm_1a <- data.frame(
  selected_x_coords_invasive_1a,  # Your x coordinates
  selected_y_coords_invasive_1a   # Your y coordinates
)

samp@meta.data$neuron_new <- "Other"  # Set all to "Other" initially

# Update for the first group: "Neuron-Healthy"
for (i in 1:nrow(selected_coords_neuron_healthy)) {
  selected_x <- selected_coords_neuron_healthy[i, "selected_x_coords_neuron_healthy"]
  selected_y <- selected_coords_neuron_healthy[i, "selected_y_coords_neuron_healthy"]
  match_idx <- which(samp@meta.data$x == selected_x & samp@meta.data$y == selected_y)
  samp@meta.data$neuron_new[match_idx] <- "Neuron-Healthy"
}

# Update for the second group: "Neuron-1A"
for (i in 1:nrow(selected_coords_neuron_1a)) {
  selected_x <- selected_coords_neuron_1a[i, "selected_x_coords_neuron_1a"]
  selected_y <- selected_coords_neuron_1a[i, "selected_y_coords_neuron_1a"]
  match_idx <- which(samp@meta.data$x == selected_x & samp@meta.data$y == selected_y)
  samp@meta.data$neuron_new[match_idx] <- "Neuron-1A"
}

for (i in 1:nrow(selected_coords_neuron_1b)) {
  selected_x <- selected_coords_neuron_1b[i, "selected_x_coords_neuron_1b"]
  selected_y <- selected_coords_neuron_1b[i, "selected_y_coords_neuron_1b"]
  match_idx <- which(samp@meta.data$x == selected_x & samp@meta.data$y == selected_y)
  samp@meta.data$neuron_new[match_idx] <- "Neuron-1B"
}

for (i in 1:nrow(selected_coords_neuron_1c)) {
  selected_x <- selected_coords_neuron_1c[i, "selected_x_coords_neuron_1c"]
  selected_y <- selected_coords_neuron_1c[i, "selected_y_coords_neuron_1c"]
  match_idx <- which(samp@meta.data$x == selected_x & samp@meta.data$y == selected_y)
  samp@meta.data$neuron_new[match_idx] <- "Neuron-1C"
}

# Update for the third group: "gbm"
for (i in 1:nrow(selected_coords_gbm_1a)) {
  selected_x <- selected_coords_gbm_1a[i, "selected_x_coords_invasive_1a"]
  selected_y <- selected_coords_gbm_1a[i, "selected_y_coords_invasive_1a"]
  match_idx <- which(samp@meta.data$x == selected_x & samp@meta.data$y == selected_y)
  samp@meta.data$neuron_new[match_idx] <- "GBM"
}


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

####P2 Coree, removing empty channel
# 1) Convert coords to data frame
orig_coords <- samp@images$GBM1031C@boundaries$centroids@coords
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
samp@images$GBM1031C@boundaries$centroids@coords <- coords_new

SpatialDimPlot(samp)


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

#Save metadata labels of selected regions 
write.csv(
  data.frame(neuron_new = samp@meta.data$neuron_new), 
  file = "neuron_new_metadata_P2_core.csv", 
  row.names = TRUE
)




#saveRDS(samp, file = "P2_Core_20250409.rds")




