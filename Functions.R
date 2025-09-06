#description: Functions for the study
#author: Marc Kevin Schneider
#date: June & July 2025


# 1 - install and load packages  ####
#-----------------------------------#

# Define packages
list.of.packages <- c("terra", "sf", "predicts", "virtualspecies", "blockCV", 
                      "dismo", "dplyr", "geodata", "caret", "mgcv", "gam",
                      "randomForest", "xgboost", "biomod2", "caret", "doParallel",
                      "maxnet", "mda", "precrec", "ecospat", "ggplot2",
                      "ggridges", "tidyr", "purrr")

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

#setting the seed
set.seed(2962)

###########################################################
#### CHANGE ONLY THIS PATH FOR THE CODE TO WORK ###########
###########################################################
path <- "/Data/"

# 2 - Function for creating the virtual species ####
#-----------------------------------------#

create_vs <- function(landscape_name, species_name, beta_val){
  '
  Purpose: Takes a landscape name, a species name and a beta value to create a virtual species with
           a set response function
           
  Parameters:
  ---------------------------------------
  landscape_name: str
      Name of the neutral landscape model that should be used for the virtual species
  
  species_name: str
      Name of the virtual species that should get created
  
  beta_val: float
      Threshold for the presence-absence points
      
  
  Returns:
  ---------------------------------------
  A virtual species
  '
  
  wd <- paste0(path, "landscapes/")
  if(!dir.exists(wd)) dir.create(wd, recursive = TRUE)
  
  #Reading the landscape
  ls <- terra::rast(paste0(wd, landscape_name))
  
  #defining the response function
  params <- formatFunctions(mpd = c(fun = "dnorm", mean=0.4, sd=0.25),
                             gaussian = c(fun = "linearFun", a = 0.2, b = 0.3),
                             mosaic = c(fun = "quadraticFun", a = 0.1, b = 0.7, c = 0.3),
                             fbm = c(fun = "betaFun", p1 = 0.2, p2 = 0.4, alpha = 0.05, gamma = 0.2))
  
  #generating the species
  species <- generateSpFromFun(raster.stack = ls[[c("mpd", "gaussian", "mosaic", "fbm")]],
                               parameters = params,
                               formula = "mpd + gaussian + mosaic + fbm",
                               plot = TRUE)
  print("Generated the species...")
  
  #converting to presence-absence data
  speciesPA <- convertToPA(species, plot = TRUE, beta=beta_val)
  speciesList <-  list(species, speciesPA)
  species_path <- paste0(path, "virtual_species/")
  if(!dir.exists(species_path)) dir.create(species_path, recursive = TRUE)
  save_name <- paste0(species_path, species_name)
  saveRDS(speciesList, save_name)
  print("Saved the species...")
  return(speciesList)
}


# 3 - Function for sampling the presence and background points ####
#-----------------------------------------#

sample_p_bkg <- function(species_name, landscape_name, iter, sample_p, background_p){
  '
  Purpose: Samples presence and background points from a virtual species in a landscape,
           alongside a presence-absence dataset with 20,000 points
           
  Parameters:
  ---------------------------------------
  species_name: str
      Name of the virtual species
  
  landscape_name: str
      Name of the landscape
      
  iter: int
      Number of iteration
      
  sample_p: int
      Number of sampled presence points
      
  background_p: int
      Number of sampled background points
      
  
  Returns:
  ---------------------------------------
  A presence-absence dataset, a 80% training and 20% testing dataset and a background point dataset
  '
  
  wd <- paste0(path, "sample_points/")
  if(!dir.exists(wd)) dir.create(wd, recursive = TRUE)
  sp_path <- paste0(path, "virtual_species/")
  ls_path <- paste0(path, "landscapes/")
  
  #reading the species and landscape data
  species <- readRDS(paste0(sp_path, species_name, "_", landscape_name, ".rds"))
  
  landscape <- terra::rast(paste0(ls_path, landscape_name, ".tif"))
  
  
  #extracting the occurrence data
  presence <- species[[2]]
  #sampling the independent presence-absence dataset for validation later
  presence_absence <- sampleOccurrences(presence, n=20000, 
                                        type="presence-absence", sample.prevalence=0.5)
  
  #print("Sampled the presence-absence points...")
  species[[3]]<- presence_absence
  #saving the presence-absence data in a separate folder that is created if it doesnt exist yet
  dir_pres_abs <- paste0(wd, "Presence_Absence/", species_name, "/")
  if(!dir.exists(dir_pres_abs)) dir.create(dir_pres_abs, recursive = TRUE)
  saveRDS(species, paste0(dir_pres_abs, iter, "_", species_name, "_", landscape_name, "_Presence_Absence.RDS"))
  
  #sampling the 100 presence points
  presence_points <- sampleOccurrences(presence, n=sample_p, type="presence only", replacement=FALSE)
  presence_df <- as.data.frame(presence_points$sample.points)
  presence_sf <- sf::st_as_sf(presence_df,coords = c("x", "y"),crs = terra::crs(landscape),remove = F)
  
  #sampling the 2000 bkg points
  background_points <- sf::st_as_sf(as.data.frame(predicts::backgroundSample(mask=landscape, n=background_p)), 
                                    crs=terra::crs(landscape), coords=c("x","y"), remove=F)
  #print("Sampled all background points...")
  
  #extracting the data for the background points
  bg_extr=terra::extract(landscape, background_points)
  background_points =cbind(background_points,bg_extr);rm(bg_extr)
  
  #extracting the data for the presence-points
  species_data_extr <- terra::extract(landscape, presence_sf)
  species_data_compl <- cbind(presence_sf, species_data_extr)
  #print("Extracted the species and background data...")
  
  `%not_in%`<- Negate(`%in%`)
  
  #splitting into 20% testing 80% training data
  testData=dplyr::slice_sample(species_data_compl, prop=.2)
  trainData=dplyr::filter(species_data_compl, ID %not_in% testData$ID )
  
  #creating the directories for the test, training and background data
  #if they do not exist before this
  dir_test <- paste0(wd, "Testing_Data/", species_name, "/")
  if(!dir.exists(dir_test)) dir.create(dir_test, recursive = TRUE)
  
  dir_train <- paste0(wd, "Training_Data/", species_name, "/")
  if(!dir.exists(dir_train)) dir.create(dir_train, recursive = TRUE)
  
  dir_bkg <- paste0(wd, "Background_Data/", species_name, "/")
  if(!dir.exists(dir_bkg)) dir.create(dir_bkg, recursive = TRUE)
  
  #saving all three datasets separately and in separate folders
  sf::write_sf(testData, paste0(dir_test, iter, "_",
                                species_name, "_", landscape_name, "_testData.gpkg"))
  
  sf::write_sf(trainData, paste0(dir_train, iter, "_",
                                 species_name, "_", landscape_name, "_trainData.gpkg"))
  
  sf::write_sf(background_points, paste0(dir_bkg, iter, "_",
                                         species_name, "_", landscape_name, "bkg_points.gpkg"))
}


# 4 - Reading and formatting the training and background data ####
#-----------------------------------------#
read_format <- function(species_name, landscape_name, iter){
  '
  Helper function
  
  Purpose: Reads and formats the presence-only and background points for model training
           
  Parameters:
  ---------------------------------------
  species_name: str
      Name of the virtual species
  
  landscape_name: str
      Name of the landscape
      
  iter: int
      Number of iteration
      
  
  Returns:
  ---------------------------------------
  A DataFrame with the presence-only and background points data in a proper format
  as well as the model path for saving the model
  '
  #all paths needed
  train_path <- paste0(path, "sample_points/Training_Data/")
  bkg_path <- paste0(path, "sample_points/Background_Data/")
  model_path <- paste0(path, "models/")
  if(!dir.exists(model_path)) dir.create(model_path, recursive = TRUE)
  ls_path <- paste0(path, "landscapes/")
  
  #reading the landscape data
  landscape <- terra::rast(paste0(ls_path, landscape_name, ".tif"))
  vars <- names(landscape)
  
  #reading the training data
  po=sf::read_sf(paste0(train_path, species_name, "/", iter, "_",
                        species_name, "_", landscape_name, "_trainData.gpkg"))
  po$occ=1
  #setting Real and Observed to NULL so that the amounts of columns matches
  po$Real <- NULL
  po$Observed <- NULL
  bg=sf::read_sf(paste0(bkg_path, species_name, "/", iter, "_",
                        species_name, "_", landscape_name, "bkg_points.gpkg"))
  bg$occ=0
  
  #print(po)
  #print(bg)
  
  #binding the data
  data=rbind(po,bg)%>%as.data.frame()%>%dplyr::select(-geom)
  
  results <- list(data = data, model_path = model_path)
  return(results)
}



# 5 - Training the GAM ####
#-----------------------------------------#

gam_train <- function(species_name, landscape_name, iter){
  '
  Purpose: Trains a GAM model for a specific sample, landscape and iteration
           
  Parameters:
  ---------------------------------------
  species_name: str
      Name of the virtual species
  
  landscape_name: str
      Name of the landscape
      
  iter: int
      Number of iteration
      
  
  Returns:
  ---------------------------------------
  Saves a GAM model for the specific combination
  '
  #takes the results from the helper function
  result <- read_format(species_name, landscape_name, iter)
  data <- result$data
  model_path <- result$model_path
  #print(data)
  
  #defines the cross-validation
  ctrl = caret::trainControl(method="cv",
                             number=5)
  
  tuneGrid=expand.grid(mstop = c(50,100,150),
                       prune=c("no"))
  
  #create weights
  wt <- ifelse(data$occ == 1, 1, nrow(data[data$occ==1,])/nrow(data)) # down-weighting
  
  tmp <- Sys.time()
  
  #necessary for the gam function from mgcv
  #copied and adjusted from https://shorturl.at/lCw4n
  form <- occ ~ mpd + gaussian + mosaic + fbm
  
  #training the model
  gm <- mgcv::gam(
    formula = form,
    data = data,
    family = binomial(link = "logit"),
    weights = wt,
    method = "REML"
  )
  
  #for calculating the time it took to train the model
  time_train <- round(Sys.time() - tmp, 3)
  
  dir_mod <- paste0(model_path, "GAM/", species_name, "/")
  print(dir_mod)
  if(!dir.exists(dir_mod)) dir.create(dir_mod, recursive=TRUE)
  
  model_name <- paste0(dir_mod, iter, "_", species_name, "-",
                       landscape_name, "_", "GAM.RDS")
  
  #saving the model as a RDS file
  saveRDS(gm, model_name)
  print(paste0("Training the GAM took: ", time_train, " seconds"))
}


# 6 - Training the GLM ####
#-----------------------------------------#

glm_train <- function(species_name, landscape_name, iter){
  '
  Purpose: Trains a GLM model for a specific sample, landscape and iteration
           
  Parameters:
  ---------------------------------------
  species_name: str
      Name of the virtual species
  
  landscape_name: str
      Name of the landscape
      
  iter: int
      Number of iteration
      
  
  Returns:
  ---------------------------------------
  Saves a GLM model for the specific combination
  '
  #reads the data
  result <- read_format(species_name, landscape_name, iter)
  data <- result$data
  model_path <- result$model_path
  
  #calculating the weights
  pres_num <- as.numeric(table(data$occ)["1"])
  bkg_num <- as.numeric(table(data$occ)["0"])
  wt <- ifelse(data$occ == 1, 1, pres_num / bkg_num)
  
  tmp <- Sys.time()
  
  #simple GLM
  form <- occ ~ mpd + gaussian + mosaic + fbm
  
  glm <- glm(
    formula = form,
    data = data,
    weights = wt,
    family = binomial(link = "logit")
  )
  
  #model training time
  time_train <- round(Sys.time() - tmp, 3)
  
  dir_mod <- paste0(model_path, "GLM/", species_name, "/")
  #print(dir_mod)
  if(!dir.exists(dir_mod)) dir.create(dir_mod, recursive=TRUE)
  
  model_name <- paste0(dir_mod, iter, "_", species_name, "-",
                       landscape_name, "_", "GLM.RDS")
  
  #saving the model as a RDS file
  #saveRDS(lm_subset, model_name)
  saveRDS(glm, model_name)
  
  print(paste0("Training the GLM took: ", time_train, " seconds"))
}


# 7 - Training a downsampled RF ####
#-----------------------------------------#

rf_train <- function(species_name, landscape_name, iter){
  '
  Purpose: Trains a downsampled Random Forest model for a specific sample, landscape and iteration
           
  Parameters:
  ---------------------------------------
  species_name: str
      Name of the virtual species
  
  landscape_name: str
      Name of the landscape
      
  iter: int
      Number of iteration
      
  
  Returns:
  ---------------------------------------
  Saves a downsampled Random Forest model for the specific combination
  '
  #reads the data
  result <- read_format(species_name, landscape_name, iter)
  data <- result$data
  model_path <- result$model_path
  
  #converting the occurrence data to factor so that RF can use it
  data$occ <- factor(data$occ, levels = c("0", "1"))
  
  #for downsampling the occurrence and background classes
  #otherwise the data is imbalanced (which isnt great for RF)
  tbl <- table(data$occ)
  minority_size <- min(tbl)
  smpsize <- setNames(rep(minority_size, 2), names(tbl))
  
  #print(smpsize)
  
  tmp <- Sys.time()
  
  #training the downsampled RF
  rf_downsample <- randomForest(formula = occ ~ mpd + gaussian + mosaic + fbm,
                                data = data,
                                ntree = 1000,
                                sampsize = smpsize,
                                replace = TRUE)
  
  #model training time
  time_train <- round(Sys.time() - tmp, 3)
  
  dir_mod <- paste0(model_path, "RF/", species_name, "/")
  #print(dir_mod)
  if(!dir.exists(dir_mod)) dir.create(dir_mod, recursive=TRUE)
  
  model_name <- paste0(dir_mod, iter, "_", species_name, "-",
                       landscape_name, "_", "RF.RDS")
  
  #saving the model as a RDS file
  saveRDS(rf_downsample, model_name)
  
  print(paste0("Training the RF took: ", time_train, " seconds"))
}


# 8 - Training an XGBoost model ####
#-----------------------------------------#

xgboost_train <- function(species_name, landscape_name, iter){
  '
  Purpose: Trains a XGBoost model for a specific sample, landscape and iteration
           
  Parameters:
  ---------------------------------------
  species_name: str
      Name of the virtual species
  
  landscape_name: str
      Name of the landscape
      
  iter: int
      Number of iteration
      
  
  Returns:
  ---------------------------------------
  Saves a XGBoost model for the specific combination
  '
  #reads the data
  result <- read_format(species_name, landscape_name, iter)
  data <- result$data
  model_path <- result$model_path
  
  data$occ <- as.factor(data$occ)
  levels(data$occ) <- c("C0", "C1") #caret does not accept 0 and 1 as factor levels, so changed to this
  
  mytrControl <- trainControl(method = "cv",
                              number = 7, #7fold CV
                              summaryFunction = twoClassSummary,
                              classProbs = TRUE,
                              allowParallel = TRUE)
  
  
  #setting the range of parameters for grid search tuning
  mytuneGrid <- expand.grid(
    nrounds = seq(from = 500, to = 15000, by = 500),
    eta = 0.001,
    max_depth = 6,
    subsample = 0.75,
    gamma = 0,
    colsample_bytree = 0.8,
    min_child_weight = 1
  )
  
  tmp <- Sys.time()
  cluster <- makeCluster(6) #so that the code uses my 6 cores
  registerDoParallel(cluster)

  xgb_fit <- train(form = occ ~ mpd + gaussian + mosaic + fbm,
                   data = data,
                   method = "xgbTree",
                   metric = "Accuracy",
                   trControl = mytrControl,
                   tuneGrid = mytuneGrid,
                   verbose = TRUE)
  
  stopCluster(cluster)
  registerDoSEQ()
  
  #model training time
  time_train <- round(Sys.time() - tmp, 3)
  
  dir_mod <- paste0(model_path, "XGBoost/", species_name, "/")
  #print(dir_mod)
  if(!dir.exists(dir_mod)) dir.create(dir_mod, recursive=TRUE)
  
  model_name <- paste0(dir_mod, iter, "_", species_name, "-",
                       landscape_name, "_", "XGBoost.RDS")
  
  #saving the model as a RDS file
  saveRDS(xgb_fit, model_name)
  
  print(paste0("Training the XgBoost model took: ", time_train, " seconds"))
}

# 9 - Training a biomod model ####
#-----------------------------------------#

biomod_train <- function(species_name, landscape_name, iter){
  '
  Purpose: Trains a biomod ensemble model for a specific sample, landscape and iteration
           
  Parameters:
  ---------------------------------------
  species_name: str
      Name of the virtual species
  
  landscape_name: str
      Name of the landscape
      
  iter: int
      Number of iteration
      
  
  Returns:
  ---------------------------------------
  Saves a biomod model for the specific combination
  '
  result <- read_format(species_name, landscape_name, iter)
  data <- result$data
  model_path <- result$model_path
  
  covars <- c("mpd", "gaussian", "mosaic", "fbm")
  resp_var <- "occ"
  
  resp <- as.numeric(data[[resp_var]])
  
  #setting the background data to NA since biomod needs this for some reason
  resp[resp == 0] <- NA
  
  #print(table(resp))
  
  expl <- data[, covars]
  #biomod needs the x and y coordinates for some reason
  resp_xy <- data[, c("x", "y")]
  
  #format for the biomod models
  biomod_data <- BIOMOD_FormatingData(
    resp.var = resp,
    expl.var = expl,
    resp.name = resp_var,
    resp.xy = resp_xy,
    PA.nb.rep = 5, 
    PA.nb.absences = 2000, 
    PA.strategy = "random", 
    PA.dist.min = 0,
    na.rm = TRUE
  )
  
  #models used for the ensemble model
  #sadly CTA, FDA, GAM, GBM, GLM and MaxEnt are failing for some reason
  #models <- c("ANN", "CTA", "FDA", "GAM", "GBM", "GLM", "MARS",
            #  "MAXENT", "MAXNET", "RF", "RFd", "SRE", "XGBOOST")
  models <- c("ANN", "MARS", "MAXNET", "RF", "RFd", "SRE", "XGBOOST")
  
  #define modeling options
  biomod_options <- bm_ModelingOptions(
    data.type = "binary",
    models = models,
    strategy = "default"
  )
  #print("Hu")
  modeling_id <- paste0(iter, "_", species_name, "-",
                        landscape_name, "_")
  dir_mod <- paste0(model_path, "Biomod/", species_name, "/")
  #print(dir_mod)
  if(!dir.exists(dir_mod)) dir.create(dir_mod, recursive=TRUE)
  #print("He")
  tmp <- Sys.time()
  
  #actually training the models
  biomod_model_out <- BIOMOD_Modeling(
    bm.format = biomod_data,
    modeling.id = modeling_id,
    OPT.user = biomod_options,
    models = models,
    CV.strategy = "random",
    CV.nb.rep = 5, #has to match PA.nb.rep from BIOMOD_FormatingData
    CV.perc = 0.8,
    CV.do.full.models = TRUE,
    scale.models = FALSE,
    OPT.strategy = "default",
    metric.eval = c("TSS", "ROC"),
    var.import = 3,
    nb.cpu = 6
  )
  
  #print("Ha")
  
  #also using the ensemble model
  biomod_ensemble <- BIOMOD_EnsembleModeling(
    biomod_model_out,
    models.chosen = "all",
    em.by = "all",
    metric.eval = c("ROC"),
    metric.select.thresh = NULL,
    em.algo = c("EMmean")
  )
  
  model_name <- paste0(dir_mod, iter, "_", species_name, "-",
                       landscape_name, "_", "Biomod.RDS")
  
  ens_mod_name <- paste0(dir_mod, iter, "_", species_name, "-",
                       landscape_name, "_", "Ensemble_Biomod.RDS")
  
  saveRDS(biomod_model_out, model_name)
  saveRDS(biomod_ensemble, ens_mod_name)
  
  time_train <- round(Sys.time() - tmp, 3)
  print(paste0("Training the biomod2 models took: ", time_train, " seconds"))
}

# 9 - Tuning a MaxEnt model ####
#-----------------------------------------#

#from https://shorturl.at/f0BXw section MaxEnt and MaxNet
maxent_param <- function(data, y = "occ", k = 7, folds = NULL, filepath){
  '
  Helper function for tuning the MaxEnt model
  '
  if(is.null(folds)){
    # generate balanced CV folds
    folds <- caret::createFolds(y = as.factor(data$occ), k = k)
  }
  names(data)[which(names(data) == y)] <- "occ"
  covars <- c("mpd", "gaussian", "mosaic", "fbm")
  # regularization multipliers
  ms <- c(0.5, 1, 2, 3, 4)
  grid <- expand.grid(
    regmult = paste0("betamultiplier=", ms),
    features = list(
      c("noautofeature", "nothreshold"), # LQHP
      c("noautofeature", "nothreshold", "noproduct"), # LQH
      c("noautofeature", "nothreshold", "nohinge", "noproduct"), # LQ
      c("noautofeature", "nothreshold", "nolinear", "noquadratic", "noproduct"), # H
      c("noautofeature", "nothreshold", "noquadratic", "nohinge", "noproduct")), # L
    stringsAsFactors = FALSE
  )
  AUCs <- c()
  for(n in seq_along(grid[,1])){
    full_pred <- data.frame()
    for(i in seq_len(length(folds))){
      trainSet <- unlist(folds[-i])
      testSet <- unlist(folds[i])
      if(inherits(try(
        maxmod <- dismo::maxent(x = data[trainSet, covars],
                                p = data$occ[trainSet],
                                removeDuplicates = FALSE,
                                path = filepath,
                                args = as.character(unlist(grid[n, ]))
        )
      ), "try-error")){
        next
      }
      modpred <- predict(maxmod, data[testSet, covars], args = "outputformat=cloglog")
      pred_df <- data.frame(score = modpred, label = data$occ[testSet])
      full_pred <- rbind(full_pred, pred_df)
    }
    AUCs[n] <- precrec::auc(precrec::evalmod(scores = full_pred$score,
                                             labels = full_pred$label))[1,4]
  }
  best_param <- as.character(unlist(grid[which.max(AUCs), ]))
  return(best_param)
}



# 10 - Training a MaxEnt model ####
#-----------------------------------------#

maxent_train <- function(species_name, landscape_name, iter){
  '
  Purpose: Trains a MaxEnt model for a specific sample, landscape and iteration
           
  Parameters:
  ---------------------------------------
  species_name: str
      Name of the virtual species
  
  landscape_name: str
      Name of the landscape
      
  iter: int
      Number of iteration
      
  
  Returns:
  ---------------------------------------
  Saves a MaxEnt model for the specific combination
  '
  result <- read_format(species_name, landscape_name, iter)
  data <- result$data
  model_path <- result$model_path
  
  nfolds <- ifelse(sum(data$occ) < 10, 2, 7)
  
  covars <- c("mpd", "gaussian", "mosaic", "fbm")
  
  param_optim <- maxent_param(data = data,
                              k = nfolds,
                              filepath = paste0(model_path, 
                                                "output/maxent_files"))
  tmp <- Sys.time()
  
  #fit MaxEnt model with the tuned parameters
  maxmod <- dismo::maxent(x = data[, covars],
                          p = data$occ,
                          removeDuplicates = FALSE,
                          path = paste0(model_path,"output/maxent_files"),
                          args = param_optim)
  
  dir_mod <- paste0(model_path, "MaxEnt/", species_name, "/")
  if(!dir.exists(dir_mod)) dir.create(dir_mod, recursive=TRUE)
  
  model_name <- paste0(dir_mod, iter, "_", species_name, "-",
                       landscape_name, "_", "MaxEnt.RDS")
  
  #saving the model as a RDS file
  saveRDS(maxmod, model_name)
  
  time_train <- round(Sys.time() - tmp, 3)
  print(paste0("Training the MaxEnt model took: ", time_train, " seconds"))
}


# 11 - Predicting with the model ####
#-----------------------------------------#
predict_model <- function(species_name, landscape_name, model_name, iter) {
  '
  Purpose: Predicting for each landscape with every model
           
  Parameters:
  ---------------------------------------
  species_name: str
      Name of the virtual species
  
  landscape_name: str
      Name of the landscape
      
  model_name: str
      Name of the model
      
  iter: int
      Number of iteration
      
  
  Returns:
  ---------------------------------------
  A .tif file of the prediction of each model
  '
  model_path <- paste0(path, "models/")
  ls_path <- paste0(path, "landscapes/")
  pred_path <- paste0(path, "predictions/")
  if(!dir.exists(pred_path)) dir.create(pred_path, recursive=TRUE)
  test_path <- paste0(path, "sample_points/Testing_Data/")
  
  landscape <- terra::rast(paste0(ls_path, landscape_name, ".tif"))
  
  #checks for the model name, for biomod two models are needed
  #so it reads the simple and ensemble model
  if (model_name == "Biomod"){
    simple_model <- readRDS(
        paste0(model_path, "Biomod", "/", species_name, "/", iter, "_",
               species_name, "-", landscape_name, "_", model_name, ".RDS"))
    
    fitted_model <- readRDS(
      paste0(model_path, "Biomod", "/", species_name, "/", iter, "_",
             species_name, "-", landscape_name, "_", 
             "Ensemble_", model_name, ".RDS"))
    
    #print(class(simple_model))
    #print(class(fitted_model))

  } else {
    #reads the model the normal way it if isnt biomod
    fitted_model <- readRDS(
      paste0(model_path, model_name, "/", species_name, "/", iter, "_",
             species_name, "-", landscape_name, "_", model_name, ".RDS"))
  }
  print("Successfully loaded the model...")
  
  #extracts the model variables
  #here MaxEnt and Biomod have different formats compared to the other
  #four models, so they need to be extracted separately
  if (model_name == "MaxEnt"){
    lambda_info <- fitted_model@lambdas
    variable_names <- unique(sapply(strsplit(lambda_info, ","), "[", 1))
    expected_vars <- variable_names[1:4]
  } else if (model_name == "Biomod"){
    expected_vars <- fitted_model@expl.var.names
    proj_name <- paste0(iter, species_name, landscape_name, model_name)
  } else {
    expected_vars <- all.vars(stats::formula(fitted_model))[-1]
  }
          
  #just a check for the prediction
  names(landscape) <- expected_vars
  
  #biomod also predicts differently
  if (model_name == "Biomod"){
    #first a prediction with the simple model
    simple_proj <- BIOMOD_Projection(
      bm.mod = simple_model,
      new.env = landscape,
      proj.name = proj_name,
      selected.models = "all",
      binary.meth = "ROC",
      compress = TRUE,
      clamping.mask = TRUE
    )

    #then a prediction using the ensemble model
    ens_proj <- BIOMOD_EnsembleForecasting(bm.em = fitted_model,
                                       bm.proj = simple_proj,
                                       models.chosen = "all")
    
    #extracting the raster from the prediction
    pred_raster <- biomod2::get_predictions(ens_proj)
    #and using the ROC CV raster for the .tif file
    pred <- pred_raster$occ_EMmeanByROC_mergedData_mergedRun_mergedAlgo
    
    #(also normalizing the predictions)
    pred <- terra::app(pred, fun = function(x) {
      (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    })
  } else {
    #using terras predict for the other models
    pred <- terra::predict(object = landscape, model = fitted_model, na.rm = TRUE)
    #also normalizing here
    pred <- terra::app(pred, fun = function(x) {
      (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    })
  }
  
  dir_pred <- paste0(pred_path, model_name, "/", species_name, "/")
  if (!dir.exists(dir_pred)) dir.create(dir_pred, recursive = TRUE)
  
  pred_name <- paste0(dir_pred, iter, "_", species_name, "_",
                      landscape_name, "_", "prediction.tif")
  
  #and then just writing the data as a .tif prediction file
  terra::writeRaster(pred, pred_name, overwrite = TRUE)
  print("Successfully saved the prediction...")
}


# 12 - Reading the testing, background and presence-absence data for validation ####
#-----------------------------------------#
read_val_data <- function(species_name, landscape_name, model, iter){
  '
  Helper function
  
  Purpose: Reads and formats the testing, background and presence-absence data for validation
           
  Parameters:
  ---------------------------------------
  species_name: str
      Name of the virtual species
  
  landscape_name: str
      Name of the landscape
      
  iter: int
      Number of iteration
      
  
  Returns:
  ---------------------------------------
  A DataFrame with the testing, background and presence-absence data
  '
  #all paths needed
  test_path <- paste0(path, "sample_points/Testing_Data/")
  bkg_path <- paste0(path, "sample_points/Background_Data/")
  pr_ab_path <- paste0(path, "sample_points/Presence_Absence/")
  ls_path <- paste0(path, "landscapes/")
  pred_path <- paste0(path, "predictions/")
  
  #reading the landscape data
  landscape <- terra::rast(paste0(ls_path, landscape_name, ".tif"))
  vars <- names(landscape)
  
  #reading the prediction data
  pred <- terra::rast(paste0(pred_path, model, "/", species_name, "/", iter, "_",
                             species_name, "_", landscape_name, 
                             "_prediction.tif"))
  
  #reading the testing data
  testData <-  sf::read_sf(paste0(test_path, species_name, "/", iter, "_",
                               species_name, "_", landscape_name, 
                               "_testData.gpkg"))
  
  #reading the presence-absence data
  pres_abs_data <- readRDS(paste0(pr_ab_path, species_name, "/", iter, "_",
                                  species_name, "_", landscape_name, 
                                  "_Presence_Absence.RDS"))
  
  
  #reading the background data
  bg <- sf::read_sf(paste0(bkg_path, species_name, "/", iter, "_",
                        species_name, "_", landscape_name, "bkg_points.gpkg"))
  
  #boyceIndex <- ecospat::ecospat.boyce(fit=pred, 
  #                                     obs=sf::st_coordinates(testData),
  #                                     PEplot = FALSE) #false since true crashes

  
  #boyceIndex <- boyceIndex$cor
  
  #extracting the test and bck data from the prediction
  extrTest <- terra::extract(pred, testData)
  colnames(extrTest) <- c( "ID"  , "predicted")
  extrTest$observed <- 1
  
  extrBg <- terra::extract(pred, bg)
  colnames(extrBg) <- c( "ID"  , "predicted")
  extrBg$observed <- 0
  
  #bind both dataframes together:
  testData=rbind(extrTest, extrBg)
  rm(extrTest, extrBg, bg)
  
  #extract the presence-absence data
  pa_raw <- pres_abs_data[[3]]$sample.points
  pa_points <- sf::st_as_sf(pa_raw, coords = c("x", "y"), crs = terra::crs(pred))
  
  
  #extract the data from the prediction .tif file
  pa_extract <- terra::extract(pred, pa_points, method = "bilinear")
  
  colnames(pa_extract)[2] <- "predicted"
  
  pa_df <- cbind(pa_extract, observed = pa_raw$Observed)
  
  #print(pa_df)
  
  #pearson R
  pearson_r <- cor(
    as.numeric(pa_df$predicted), #without this the code would crash fsr
    as.numeric(pa_df$observed),
    method = "pearson",
    use = "complete.obs"
  )
  
  
  #calculate AUC and MSE
  AUC<-Metrics::auc(actual = testData$observed, predicted = testData$predicted)
  MSE<-Metrics::mse(actual = testData$observed, predicted = testData$predicted)
  #results <- list(boyceIndex = boyceIndex, AUC = AUC, MSE = MSE,
  #                pearson_r = pearson_r)
  
  results <- list(AUC = AUC, MSE = MSE, pearson_r = pearson_r)
  return(results)
}


# 13 - Function for plotting the results as boxplots ####
#-----------------------------------------#

plot_results_boxplot <- function(df, eval_metric) {
  '
  Purpose: plots the model results for each landscape in a 3x5 matrix
           
  Parameters:
  ---------------------------------------
  df: dataframe
      Dataframe of the model results
  
  eval_metric: str
      Name of evaluation metric
      possible are "Pearson correlation", "AUC" and "MSE"

      
  
  Returns:
  ---------------------------------------
  A 3x5 matrix of boxplots of the model results
  '
  #to lower key to avoid misspelling
  metric_key <- tolower(eval_metric)
  
  #match name to column name in the df
  metric <- if (metric_key == "pearson correlation") {
    "Pearson"
  } else if (metric_key == "auc") {
    "AUC"
  } else if (metric_key == "mse") {
    "MSE"
  } else {
    stop("Unsupported metric. Use 'Pearson correlation', 'AUC', or 'MSE'.")
  }
  
  #fit title
  plot_title <- paste0(eval_metric, " per model within each landscape")
  
  #plot nicely in a 3x5 matrix
  ggplot(df, aes(x = model, y = .data[[metric]], fill = model)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
    stat_summary(fun = median, geom = "point",
                 shape = 23, size = 2, fill = "black") +
    facet_grid(beta ~ landscape, scales = "free_y") +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
    labs(
      title = plot_title,
      x = "Model",
      y = eval_metric,
      fill = "Model"
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# 14 - Function for plotting the relationship between groups and eval metrics ####
#-----------------------------------------#

plot_group <- function(data, group_name) {
  '
  Purpose: Plots the relationship between eval and landscape metrics
           for the model groups in a 3x3 grids
           
  Parameters:
  ---------------------------------------
  data: dataframe
      Dataframe of the model results and landscape metrics
  
  group_name: str
      Name of model group

      
  
  Returns:
  ---------------------------------------
  A 3x3 matrix of plots of the relationship between eval and landscape metrics
  '
  if(group_name == "All") {
    df_plot <- data
  } else {
    df_plot <- filter(data, model_group == group_name)
  }
  
  ggplot(df_plot, aes(x = landscape_value, y = value)) +
    geom_point(alpha = 0.6, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "darkred") +
    stat_cor(
      aes(label = paste0("italic(R)^2 == ", ..r..^2, " ~ ',' ~ italic(p) == ", ..p..)),
      method = "pearson",
      label.x.npc = "left",
      label.y.npc = "top",
      size = 3
    ) +
    facet_grid(metric ~ landscape_metric, scales = "free_x") +
    labs(
      title = paste("Performance vs Landscape Heterogeneity –", group_name),
      x = "Landscape metric",
      y = "Performance metric"
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}


# 15 - Function for plotting the relationship between all models and one eval metric ####
#-----------------------------------------#

plot_metric_trends_by_model <- function(data, metric_name) {
  '
  Purpose: Plots the relationship between eval and landscape metrics
           for the all models in a 6x3 grid
           
  Parameters:
  ---------------------------------------
  data: dataframe
      Dataframe of the model results and landscape metrics
  
  metric_name: str
      Name of evaluation metric
      Possible are: "Pearson", "AUC", "MSE"

      
  
  Returns:
  ---------------------------------------
  A 3x3 matrix of plots of the relationship between eval and landscape metrics
  '
  #filter for the chosen evaluation metric
  df_plot <- filter(data, metric == metric_name)
  
  #compute r2 and p value per model and landscape
  r2_table <- df_plot %>%
    group_by(model, landscape_metric) %>%
    summarise(
      lm_fit = list(lm(value ~ landscape_value, data = cur_data())),
      .groups = "drop"
    ) %>%
    mutate(
      r2 = map_dbl(lm_fit, ~ summary(.x)$r.squared),
      p  = map_dbl(lm_fit, ~ summary(.x)$coefficients[2,4]),
      y_pos = 1.05 * max(df_plot$value, na.rm = TRUE)  # adjust vertical position
    )
  
  #plot the data
  ggplot(df_plot, aes(x = landscape_value, y = value)) +
    geom_point(alpha = 0.6, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "darkred") +
    geom_text(
      data = r2_table,
      aes(
        x = min(df_plot$landscape_value), 
        y = y_pos,
        label = paste0("R2=", round(r2, 2), ", p=", signif(p, 3))
      ),
      inherit.aes = FALSE,
      hjust = 0,
      size = 3
    ) +
    facet_grid(model ~ landscape_metric, scales = "free_x") +
    labs(
      x = "Landscape metric",
      y = metric_name
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}