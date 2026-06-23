library(tidyverse)
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

#P4 core
samp <- samples[[7]]
metadata <- samp@meta.data
# Define UI
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
      # Filter the data based on selection and Neuron type
      selected_data <- samp@meta.data %>%
        filter(
          x %in% event_data$x & 
            y %in% event_data$y & 
            CellType_GLSB_manual == "Invasive-front"
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

#Define selected regions
tumor_1x <- c(
  5693, 5774, 5854, 5934, 6015, 
  5773, 5854, 5934, 6014, 6095, 
  6175, 5693, 5773, 5853, 5934, 
  6094, 6175, 6255, 5692, 5773, 
  5853, 5933, 6014, 6094, 6174, 
  6255, 5933, 6013, 6093, 6174, 
  5772, 5853, 5933, 6013, 6094, 
  6174, 6254
)

tumor_1y <- c(
  4278, 4278, 4279, 4280, 4281, 
  4198, 4199, 4200, 4200, 4201, 
  4202, 4117, 4117, 4118, 4119, 
  4121, 4122, 4122, 4036, 4037, 
  4038, 4039, 4039, 4040, 4041, 
  4042, 3878, 3878, 3879, 3880, 
  3956, 3957, 3958, 3959, 3960, 
  3961, 3961
)

tumor_2x <- c(
  4331, 4412, 4492, 4572, 4653,
  4251, 4331, 4411, 4492, 4572, 
  4652, 4331, 4411, 4491, 4572, 
  4652, 4330, 4411, 4491, 4571
)

tumor_2y <- c(
  5068, 5069, 5070, 5071, 5072,
  4987, 4988, 4989, 4989, 4990, 
  4991, 4907, 4908, 4909, 4910, 
  4911, 4827, 4828, 4828, 4829
)

nearvasc_1x <- c(
  3854, 3935, 4015, 4095, 4336, 4416, 4497, 4577, 4657,
  3935, 4256, 4417, 4577, 4658,
  3934, 4015, 4095, 4175, 4255, 4336, 4416, 4496
)

nearvasc_1y <- c(
  6271, 6272, 6273, 6274, 6276, 6277, 6278, 6279, 6280,
  6352, 6356, 6357, 6359, 6360,
  6191, 6192, 6193, 6194, 6195, 6196, 6196, 6197
)

nearvasc_2x <- c(
  6579, 6659, 6739,
  6498, 6578, 6659,
  6498, 6578, 6658, 6739,
  6497, 6738,
  6497, 6577, 6658,
  6497, 6577, 6657, 6738,
  6657, 6737
)

nearvasc_2y <- c(
  4770, 4771, 4772,
  4689, 4690, 4690,
  4608, 4609, 4610, 4611,
  4528, 4530,
  4447, 4448, 4449,
  4367, 4367, 4368, 4369,
  4288, 4289
)

nearvasc_3x <- c(
  3530, 3208, 3690, 3289, 3610,
  3208, 3690, 3609,
  3689, 3368, 3449, 3609,
  3530, 3289, 3450,
  3209, 3289, 3370, 3450
)

nearvasc_3y <- c(
  5382, 5378, 5384, 5379, 5383,
  5298, 5303, 5302,
  5223, 5219, 5220, 5222,
  5462, 5460, 5462,
  5540, 5540, 5541, 5542
)

#Create data frame
selected_coords_tumor_1 <- data.frame(
  tumor_1x,  
  tumor_1y   
)

selected_coords_tumor_2 <- data.frame(
  tumor_2x, 
  tumor_2y   
)


selected_coords_nearvasc_1 <- data.frame(
  nearvasc_1x, 
  nearvasc_1y   
)

selected_coords_nearvasc_2 <- data.frame(
  nearvasc_2x,  
  nearvasc_2y   
)

selected_coords_nearvasc_3 <- data.frame(
  nearvasc_3x, 
  nearvasc_3y   
)

# Just adding metadata for new group 
samp@meta.data$vascular_tumor <- "Other"  # Set all to "Other" initially
for (i in 1:nrow(selected_coords_tumor_1)) {
  selected_x <- selected_coords_tumor_1[i, "tumor_1x"]
  selected_y <- selected_coords_tumor_1[i, "tumor_1y"]
  match_idx <- which(samp@meta.data$x == selected_x & samp@meta.data$y == selected_y)
  samp@meta.data$vascular_tumor[match_idx] <- "Invasive-front bulk 1"
}

for (i in 1:nrow(selected_coords_tumor_2)) {
  selected_x <- selected_coords_tumor_2[i, "tumor_2x"]
  selected_y <- selected_coords_tumor_2[i, "tumor_2y"]
  match_idx <- which(samp@meta.data$x == selected_x & samp@meta.data$y == selected_y)
  samp@meta.data$vascular_tumor[match_idx] <- "Invasive-front bulk 2"
}

for (i in 1:nrow(selected_coords_nearvasc_1)) {
  selected_x <- selected_coords_nearvasc_1[i, "nearvasc_1x"]
  selected_y <- selected_coords_nearvasc_1[i, "nearvasc_1y"]
  match_idx <- which(samp@meta.data$x == selected_x & samp@meta.data$y == selected_y)
  samp@meta.data$vascular_tumor[match_idx] <- "Invasive-front near vasculature 1"
}

for (i in 1:nrow(selected_coords_nearvasc_2)) {
  selected_x <- selected_coords_nearvasc_2[i, "nearvasc_2x"]
  selected_y <- selected_coords_nearvasc_2[i, "nearvasc_2y"]
  match_idx <- which(samp@meta.data$x == selected_x & samp@meta.data$y == selected_y)
  samp@meta.data$vascular_tumor[match_idx] <- "Invasive-front near vasculature 2"
}

for (i in 1:nrow(selected_coords_nearvasc_3)) {
  selected_x <- selected_coords_nearvasc_3[i, "nearvasc_3x"]
  selected_y <- selected_coords_nearvasc_3[i, "nearvasc_3y"]
  match_idx <- which(samp@meta.data$x == selected_x & samp@meta.data$y == selected_y)
  samp@meta.data$vascular_tumor[match_idx] <- "Invasive-front near vasculature 3"
}



# Visualize the spatial data
SpatialDimPlot(
  samp, 
  group.by = "vascular_tumor", 
  cols = c(
    "Other" = "grey", 
    "Invasive-front bulk 1" = "blue",
    "Invasive-front bulk 2" = "skyblue",
    "Invasive-front near vasculature 1" = "red4", 
    "Invasive-front near vasculature 2" = "red3", 
    "Invasive-front near vasculature 3" = "red"
  )
)

#Adjust boundary 
orig_coords <- samp@images$R1110C@boundaries$centroids@coords
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

bins_to_remove_row <- 43
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
samp@images$R1110C@boundaries$centroids@coords <- coords_new

SpatialDimPlot(samp)
SpatialDimPlot(
  samp, 
  group.by = "vascular_tumor", 
  cols = c(
    "Other" = "grey", 
    "Invasive-front bulk 1" = "deepskyblue2",
    "Invasive-front bulk 2" = "purple",
    "Invasive-front near vasculature 1" = "red", 
    "Invasive-front near vasculature 2" = "orange", 
    "Invasive-front near vasculature 3" = "yellow"
  )
)
write.csv(
  data.frame(vascular_tumor = samp@meta.data$vascular_tumor), 
  file = "P4Core_vasculartumor.csv", 
  row.names = TRUE
)
