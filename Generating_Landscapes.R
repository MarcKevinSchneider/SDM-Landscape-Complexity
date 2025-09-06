#description: script for generating the neutral landscape models
#author: Lukas Esselmann
#date: June 2025



# 1 - install packages  ####
#-------------------------#

# Define packages
list.of.packages <- c("RandomFieldsUtils", "RandomFields", "NLMR", "landscapetools", 
                      "terra", "ggplot2")

# Check and install missing packages
for (pkg in list.of.packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    if (pkg %in% c("RandomFieldsUtils", "RandomFields")) {
      remotes::install_github(paste("cran", pkg, sep = "/"))
    } else {
      remotes::install_github(paste("ropensci", pkg, sep = "/"))
    }
  }
}

# Load the packages
lapply(list.of.packages, library, character.only = TRUE)


# 2 - Create 4 NLMs for 5 landscapes each ####
#-----------------------#

# Chosen country name
countryname <- "Landscape"

#seed for this script
set.seed(2962)

# 5 complexity steps from not complex to very
roughness_values <- c(0.1, 0.3, 0.5, 0.7, 0.9)
autocorr_range <- c(90,70,50,30,10)
mag <- c(1,3,5,7,9)
mosaic_germs <- c(3, 5, 7, 9, 11)
fbm_frac <- c(0.9,0.7,0.5,0.3,0.1)

for (idx in 1:5) {
  mpd <- rast(nlm_mpd(ncol = 500, nrow = 500, roughness = roughness_values[idx]))
  autocorr <- rast(nlm_gaussianfield(ncol=500, nrow=500,
                                     autocorr_range = autocorr_range[idx],
                                     mag_var = mag[idx], nug = mag[idx]))
  mosaic <- rast(nlm_mosaictess(ncol=500, nrow=500, germs = mosaic_germs[idx]))
  frac <- rast(nlm_fbm(ncol=500, nrow=500, fract_dim = fbm_frac[idx]))
  
  nlm_list <- list(mpd, autocorr, mosaic, frac)
  
  if (!all(sapply(nlm_list, inherits, "SpatRaster"))) {
    stop("One or more layers failed to generate correctly.")
  }
  
  ref <- nlm_list[[1]]
  for (j in 2:length(nlm_list)) {
    nlm_list[[j]] <- resample(nlm_list[[j]], ref, method = "bilinear")
  }
  
  r <- rast(nlm_list)
  crs(r) <- "EPSG:25832"
  names(r) <- c("mpd", "gaussian", "mosaic", "fbm")
  
  out_name <- paste0(countryname, "_", idx)
  writeRaster(r, paste0("/Data/landscapes/", out_name, ".tif"),
              overwrite = TRUE, filetype = "GTiff", datatype = "FLT4S")
  
  message(paste("Saved:", out_name))
}

