#description: script for creating and sampling the virtual species for this study
#author: Marc Kevin Schneider
#date: June 2025


# 1 - install and load packages  ####
#-----------------------------------#

# Define packages
list.of.packages <- c("terra", "sf", "predicts", "virtualspecies", "blockCV", "dismo", "dplyr")

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

source("D:/Universitaet/Master/SoSe_2025/Adv_Species_Distr/FinalProject/Scripts/Functions.R")

# 2 - creating three species with different niche breadths ####
#-----------------------------------------#
setwd("D:/Universitaet/Master/SoSe_2025/Adv_Species_Distr/FinalProject/Data")

#niche breadth
beta_vals <- c(0.25, 0.5, 0.75)

#iterating through the beta thresholds for all 5 landscapes
for (b in beta_vals) {
  for (j in 1:5) {
    create_vs(
      landscape_name = paste0("Landscape_", j, ".tif"),
      species_name = paste0("Species_beta", b, "_Landscape_", j, ".rds"),
      beta_val = b
    )
  }
}


# 3 - sampling the species ####
#-----------------------------------------#

species_l <- c("Species_beta0.25", "Species_beta0.5", "Species_beta0.75")
landscape_l <- c("Landscape_1", "Landscape_2", "Landscape_3", "Landscape_4", "Landscape_5")

for (i in 1:10){
  for (species in species_l){
    for (ls in landscape_l){
      #200 presence and 4000 background points
      sample_p_bkg(species, ls, i, 200, 4000)
      print(paste0("Finished landscape:", ls))
    }
    print(paste0("Finished species:", species))
  }
}



