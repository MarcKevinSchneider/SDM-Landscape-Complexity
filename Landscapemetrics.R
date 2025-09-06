#description: Script for calculating the landscape metrics
#author: Lukas Esselmann
#date: June 2025


# 1 - install and load packages  ####
#-----------------------------------#

# Define packages
list.of.packages <- c("terra", "landscapemetrics",
                      "dplyr", "raster", "tools")

# Check and install missing packages
for (pkg in list.of.packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    # Install the package if not already installed
    install.packages(pkg, dependencies = TRUE)
  }
  # Load the package
  library(pkg, character.only = TRUE)
}

# Load the packages
lapply(list.of.packages, library, character.only = TRUE)

set.seed(2962)

#for the functions
source("/src/Functions.R")

# 2 - Load the data  ####
#-----------------------------------#



# load landscape raster layers 
raster_files <- list.files("/Data/landscapes", pattern = "Landscape_.*\\.tif$", 
                           full.names = TRUE)

# class boundaries 
breaks <- seq(0, 1, length.out = 11)

# empty list
all_metrics <- list()
global_layer_counter <- 1

# loop
for (file in raster_files) {
  
  r_stack <- rast(file)
  file_name <- file_path_sans_ext(basename(file))
  
  classified_list <- list()
  
  # rescale and classify
  for (i in 1:nlyr(r_stack)) {
    r <- r_stack[[i]]
    
    r_min <- minmax(r)[1]
    r_max <- minmax(r)[2]
    
    if (r_min == r_max) {
      warning(paste("Layer", i, "in", file_name, "has no variance and was skipped."))
      next
    }
    
    r_rescaled <- (r - r_min) / (r_max - r_min)
    rcl <- cbind(breaks[-length(breaks)], breaks[-1], 1:10)
    r_classified <- classify(r_rescaled, rcl, right = FALSE)
    
    classified_list[[i]] <- r_classified
  }
  
  # stack classified layers
  classified_stack <- rast(classified_list)
  
  # terra to raster for landscapemetrics
  classified_raster_list <- lapply(1:nlyr(classified_stack), function(i) {
    raster(classified_stack[[i]])
  })
  
  # calc landscapemetrics
  for (i in seq_along(classified_raster_list)) {
    r <- classified_raster_list[[i]]
    
    patch_area <- lsm_p_area(r)
    mean_patch_area <- mean(patch_area$value, na.rm = TRUE)
    
    contiguity <- lsm_p_contig(r)
    mean_contig <- mean(contiguity$value, na.rm = TRUE)
    
    lsi <- lsm_l_lsi(r)$value
    
    all_metrics[[global_layer_counter]] <- data.frame(
      file = file_name,
      layer_in_file = i,
      global_layer = global_layer_counter,
      mean_patch_area = mean_patch_area,
      mean_contiguity = mean_contig,
      lsi = lsi
    )
    
    global_layer_counter <- global_layer_counter + 1
  }
}

# bind results
metrics_df <- bind_rows(all_metrics)

print(metrics_df)

# save as .csv
write.csv(metrics_df, "/Data/landscapes/landscape_metrics_all_files.csv", 
          row.names = FALSE)

