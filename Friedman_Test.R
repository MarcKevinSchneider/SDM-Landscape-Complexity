#description: script for the Friedman-Aligned-Rank-Test
#author: Lukas Esselmann & Marc Kevin Schneider
#date: July and August 2025


# 1 - install and load packages  ####
#-----------------------------------#

# Define packages
list.of.packages <- c("terra", "sf", "predicts", "virtualspecies", "blockCV", 
                      "dismo", "dplyr", "tidyr",
                      "ggplot2", "ggridges"
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
source("D:/Universitaet/Master/SoSe_2025/Adv_Species_Distr/FinalProject/Scripts/Functions.R")

# 1 - friedman test  ####
#-----------------------------------#

path <- "D:/Universitaet/Master/SoSe_2025/Adv_Species_Distr/FinalProject/Data/results/"

df <- read.csv(paste0(path, "Model_Landscape_results.csv"))

# empty dataframe for saving the test results
friedman_results <- data.frame(
  landscape = numeric(),   
  chi_squared = numeric(), 
  df = integer(),          
  p_value = numeric()      
)

# extracting all landscapes
landscape_values <- unique(df$landscape)

# friedman test for each landscape
for (b in landscape_values) {
  # filter for the current landscape
  df_l <- df %>%
    filter(landscape == b) %>%
    dplyr::select(species, iteration, model, pearson_r) %>%  
    drop_na(pearson_r) %>%  #remove NAs
    pivot_wider(names_from = model, values_from = pearson_r)
  
  # check if enough iterations are available for the friedman test
  if (nrow(df_l) > 2) {
    # have to change to matrix since thats what the test needs
    friedman_matrix <- as.matrix(df_l[ , -c(1,2)])
    
    # execute test
    test_result <- friedman.test(friedman_matrix)
    
    # save in df
    friedman_results <- rbind(friedman_results, data.frame(
      landscape = b,
      chi_squared = test_result$statistic,
      df = test_result$parameter,
      p_value = test_result$p.value
    ))
  }
}

# test results
print("Friedman-test results per landscape:")
print(friedman_results)

str(friedman_results)

write.csv(friedman_results, 
          paste0(path, "FriedmanTest_Landscape_results.csv"), 
          row.names = FALSE)

# significance check
friedman_results <- friedman_results %>%
  mutate(significant = ifelse(p_value < 0.05, "Significant", "Not Significant"))

# plotting the p-value
ggplot(friedman_results, aes(x = landscape, y = p_value, group = 1, color = significant)) +
  geom_point(size = 4) +
  geom_line() +
  scale_y_log10() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 0.05, label = "0.05", color = "red", vjust = -0.5, hjust = 0) +
  labs(
    title = "Friedman test p-values per landscape (log scale)",
    x = "Landscape",
    y = "p-value (log scale)",
    color = "Significance"
  ) +
  theme_minimal()


