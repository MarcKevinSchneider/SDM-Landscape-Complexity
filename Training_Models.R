#description: script for training and validating the models
#author: Marc Kevin Schneider
#date: June and July 2025


# 1 - install and load packages  ####
#-----------------------------------#

# Define packages
list.of.packages <- c("terra", "sf", "predicts", "virtualspecies", "blockCV", "dismo", "dplyr",
                      "dismo", #for MaxEnt
                      "randomForest", #for random forest
                      "mgcv", #for GAM
                      "gam", #for GLM
                      "xgboost", "caret", #for XGBoost
                      "biomod2", "maxnet", "mda", "earth", "Formula", #for biomod
                      "plotmo", "plotrix" #for biomod
                      )

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
source("/Scripts/Functions.R")

# 2 - training the six different model approaches using the sampled data from before ####
#-----------------------------------------#
setwd("/Data/")

species_l <- c("Species_beta0.25", "Species_beta0.5", "Species_beta0.75")
landscape_l <- c("Landscape_1", "Landscape_2", "Landscape_3", "Landscape_4", "Landscape_5")

#iterates through all 10 samples, all three species and all five landscapes
#and trains the models and saves them individually
for (i in 1:10){
  for (species in species_l){
    for (ls in landscape_l){
      gam_train(species, ls, i)
      glm_train(species, ls, i)
      rf_train(species, ls, i)
      xgboost_train(species, ls, i)
      biomod_train(species, ls, i)
      maxent_train(species, ls, i)
      print(paste0("Finished landscape:", ls))
    }
    print(paste0("Finished species:", species))
  }
  print(paste("Finished iteration", i))
}

model_l <- c("GAM", "GLM", "RF", "XGBoost", "Biomod", "MaxEnt")

for (species in species_l) {
  for (model_name in model_l) {
    results_list <- list()
    validation_results_overall <- data.frame()
    for (i in 1:10) {
      validation_results_iteration <- data.frame()
      for (ls in landscape_l) {
        predict_model(species, ls, model_name, i)
        #introducing a try-catch statement in case some of the predictions 
        #are empty (like for one prediction for MaxEnt)
        val <- tryCatch(
          read_val_data(species, ls, model_name, i),
          error = function(e) {
            warning(paste("Validation failed for", species, ls, model_name, "iter", i, ":", e$message))
            list(AUC = NA_real_, MSE = NA_real_, pearson_r = NA_real_)
          }
        )
        
        val_row <- data.frame(
          species = species,
          landscape = ls,
          iteration = i,
          model = model_name,
          AUC = round(val$AUC, 4),
          MSE = round(val$MSE, 4),
          pearson_r = round(val$pearson_r, 4)
        )
        print(val_row)
        validation_results_iteration <- rbind(validation_results_iteration, val_row)
      }
      validation_results_overall <- rbind(validation_results_overall, validation_results_iteration)
      print(paste("Finished iteration", i, "for model", model_name, "and species", species))
    }
    dir_val <- paste0("Data/results/",
                      model_name, "/")
    if (!dir.exists(dir_val)) dir.create(dir_val, recursive = TRUE)
    out_path <- paste0(dir_val, species, "validation_metrics.csv")
    write.csv(validation_results_overall, out_path, row.names = FALSE)
    print(paste("Saved results for species:", species, "and model:", model_name))
  }
}


