#description: script for plotting the model results
#author: Lukas Esselmann & Marc Kevin Schneider
#date: July and August 2025


# 1 - install and load packages  ####
#-----------------------------------#

# Define packages
list.of.packages <- c("terra", "sf", "predicts", "virtualspecies", "blockCV",
                      "dplyr",
                      "dismo", #for MaxEnt
                      "randomForest", #for random forest
                      "mgcv", #for GAM
                      "gam", #for GLM
                      "xgboost", "caret", #for XGBoost
                      "biomod2", "maxnet", "mda", "earth", "Formula", #for biomod
                      "plotmo", "plotrix", #for biomod
                      "ggplot2", "ggridges", "ggpubr"
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



# 2 - Reading the model results ####
#-----------------------------------------#
path <- "/Data/results/"

# biomod
biomod_025 <- read.csv(paste0(path, "Biomod/Species_beta0.25validation_metrics.csv"))
biomod_050 <- read.csv(paste0(path, "Biomod/Species_beta0.5validation_metrics.csv"))
biomod_075 <- read.csv(paste0(path, "Biomod/Species_beta0.75validation_metrics.csv"))

# gam
gam_025 <- read.csv(paste0(path, "GAM/Species_beta0.25validation_metrics.csv"))
gam_050 <- read.csv(paste0(path, "GAM/Species_beta0.5validation_metrics.csv"))
gam_075 <- read.csv(paste0(path, "GAM/Species_beta0.75validation_metrics.csv"))

# glm
glm_025 <- read.csv(paste0(path, "GLM/Species_beta0.25validation_metrics.csv"))
glm_050 <- read.csv(paste0(path, "GLM/Species_beta0.5validation_metrics.csv"))
glm_075 <- read.csv(paste0(path, "GLM/Species_beta0.75validation_metrics.csv"))

# maxent
maxent_025 <- read.csv(paste0(path, "MaxEnt/Species_beta0.25validation_metrics.csv"))
maxent_050 <- read.csv(paste0(path, "MaxEnt/Species_beta0.5validation_metrics.csv"))
maxent_075 <- read.csv(paste0(path, "MaxEnt/Species_beta0.75validation_metrics.csv"))

# rf
rf_025 <- read.csv(paste0(path, "RF/Species_beta0.25validation_metrics.csv"))
rf_050 <- read.csv(paste0(path, "RF/Species_beta0.5validation_metrics.csv"))
rf_075 <- read.csv(paste0(path, "RF/Species_beta0.75validation_metrics.csv"))

# xgboost
xgboost_025 <- read.csv(paste0(path, "XGBoost/Species_beta0.25validation_metrics.csv"))
xgboost_050 <- read.csv(paste0(path, "XGBoost/Species_beta0.5validation_metrics.csv"))
xgboost_075 <- read.csv(paste0(path, "XGBoost/Species_beta0.75validation_metrics.csv"))


# 2 - loading the landscape metrics ####
#-----------------------------------------#

# read them
landscape_metrics <- read.csv("/Data/landscapes/landscape_metrics.csv")


# average landscape metric per file
landscape_means <- landscape_metrics %>%
  group_by(file) %>%  
  summarise(
    mean_patch_area = mean(mean_patch_area, na.rm = TRUE),
    mean_contiguity = mean(mean_contiguity, na.rm = TRUE),
    lsi = mean(lsi, na.rm = TRUE)
  ) %>%
  rename(landscape = file)  

# helper function to annotate the datasets with model name and beta value
load_model_data <- function(df, model, beta) {
  df$model <- model
  df$beta <- beta
  return(df)
}

# concat all datasets 
alle_modelle <- bind_rows(
  load_model_data(biomod_025, "Biomod", 0.25),
  load_model_data(biomod_050, "Biomod", 0.50),
  load_model_data(biomod_075, "Biomod", 0.75),
  load_model_data(gam_025, "GAM", 0.25),
  load_model_data(gam_050, "GAM", 0.50),
  load_model_data(gam_075, "GAM", 0.75),
  load_model_data(glm_025, "GLM", 0.25),
  load_model_data(glm_050, "GLM", 0.50),
  load_model_data(glm_075, "GLM", 0.75),
  load_model_data(maxent_025, "MaxEnt", 0.25),
  load_model_data(maxent_050, "MaxEnt", 0.50),
  load_model_data(maxent_075, "MaxEnt", 0.75),
  load_model_data(rf_025, "RF", 0.25),
  load_model_data(rf_050, "RF", 0.50),
  load_model_data(rf_075, "RF", 0.75),
  load_model_data(xgboost_025, "XGBoost", 0.25),
  load_model_data(xgboost_050, "XGBoost", 0.50),
  load_model_data(xgboost_075, "XGBoost", 0.75)
)

# merging model results and landscape metrics
model_results_ls_metrics <- alle_modelle %>%
  left_join(landscape_means, by = "landscape")

# save
write.csv(model_results_ls_metrics, 
          paste0(path, "Model_Landscape_results.csv"), 
          row.names = FALSE)

# reading the data again
model_results_ls_metrics <- read.csv(paste0(path, "Model_Landscape_results.csv"))

names(model_results_ls_metrics)[names(model_results_ls_metrics) == 'pearson_r'] <- 'Pearson'

# 3 - plotting the model results ####
#-----------------------------------------#

#pearson
plot_results_boxplot(model_results_ls_metrics, 'Pearson correlation')

#AUC
plot_results_boxplot(model_results_ls_metrics, 'AUC')

#mse
plot_results_boxplot(model_results_ls_metrics, 'MSE')


print(model_results_ls_metrics)



# 4 - model results vs. landscape metrics ####
#-----------------------------------------#

#grouping into the model categories
df_grouped <- model_results_ls_metrics |>
  mutate(model_group = case_when(
    model %in% c("Biomod", "MaxEnt") ~ "Other",
    model %in% c("XGBoost", "RF")    ~ "Tree",
    model %in% c("GAM", "GLM")       ~ "Regression",
    TRUE ~ model
  ))

#aggregating the model results within the categories
df_summary <- df_grouped |>
  group_by(model_group, landscape, beta) |>
  summarise(
    Pearson = mean(Pearson, na.rm = TRUE),
    AUC       = mean(AUC, na.rm = TRUE),
    MSE       = mean(MSE, na.rm = TRUE),
    mean_patch_area = mean(mean_patch_area, na.rm = TRUE),
    mean_contiguity = mean(mean_contiguity, na.rm = TRUE),
    lsi            = mean(lsi, na.rm = TRUE),
    .groups = "drop"
  )

#reshaping the model categories into long format
df_long <- df_summary |>
  pivot_longer(cols = c(Pearson, AUC, MSE),
               names_to = "metric", values_to = "value") |>
  pivot_longer(cols = c(mean_patch_area, mean_contiguity, lsi),
               names_to = "landscape_metric", values_to = "landscape_value")

#summary of all values
df_summary_all <- model_results_ls_metrics |>
  group_by(model, landscape, beta) |>
  summarise(
    Pearson = mean(Pearson, na.rm = TRUE),
    AUC       = mean(AUC, na.rm = TRUE),
    MSE       = mean(MSE, na.rm = TRUE),
    mean_patch_area = mean(mean_patch_area, na.rm = TRUE),
    mean_contiguity = mean(mean_contiguity, na.rm = TRUE),
    lsi            = mean(lsi, na.rm = TRUE),
    .groups = "drop"
  )

#also pivot table of all entries
df_long_all <- df_summary_all |>
  pivot_longer(
    cols = c(Pearson, AUC, MSE),
    names_to = "metric",
    values_to = "value"
  ) |>
  pivot_longer(
    cols = c(mean_patch_area, mean_contiguity, lsi),
    names_to = "landscape_metric",
    values_to = "landscape_value"
  )

#plotting them
plot_group(df_long_all, "All")
plot_group(df_long, "Other")
plot_group(df_long, "Tree")
plot_group(df_long, "Regression")

# rel. models vs. landscape metrics for Pearson
plot_metric_trends_by_model(df_long_all, "Pearson")

#for AUC
plot_metric_trends_by_model(df_long_all, "AUC")

#for MSE
plot_metric_trends_by_model(df_long_all, "MSE")

#grouping the data
df_long_group <- model_results_ls_metrics |>
  mutate(model_group = case_when(
    model %in% c("Biomod", "MaxEnt") ~ "Other",
    model %in% c("XGBoost", "RF")    ~ "Tree",
    model %in% c("GAM", "GLM")       ~ "Regression",
    TRUE ~ model
  )) |>
  pivot_longer(
    cols = c(Pearson, AUC, MSE),
    names_to = "metric",
    values_to = "value"
  )

#for position in the subplots
r2_table <- df_long_group %>%
  group_by(model_group, metric, model) %>%
  summarise(
    lm_fit = list(lm(value ~ as.numeric(factor(landscape)), data = cur_data())),
    .groups = "drop"
  ) %>%
  mutate(
    r2 = purrr::map_dbl(lm_fit, ~ summary(.x)$r.squared),
    p  = purrr::map_dbl(lm_fit, ~ summary(.x)$coefficients[2,4]),
    #vertical position for labels
    y_pos = case_when(model == "RF" ~ 0.95, 
                      model == "XGBoost" ~ 0.85,
                      model == "GAM" ~ 0.17,
                      model == "GLM" ~ 0.07,
                      model == "Biomod" ~ -0.05,
                      model == "MaxEnt" ~ -0.15)
  )

#plotting the data
ggplot(df_long_group, aes(x = landscape, y = value, color = model, group = model)) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_smooth(method = "lm", se = TRUE, aes(group = model)) +
  guides(color = guide_legend(override.aes = list(size = 0)))+
  geom_text(
    data = r2_table,
    aes(x = 1, y = y_pos, label = paste0("R²=", round(r2,2), ", p=", signif(p,3)), color = model),
    inherit.aes = FALSE,
    hjust = 0
  ) +
  facet_grid(model_group ~ metric, scales = "free_y") +
  theme_bw()

