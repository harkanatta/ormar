# PERMANOVA Analysis

# Load required packages
library(vegan)
library(tidyr)

# Sediment data
sediment_env_data <- readRDS("data/raw/sediment_env_data.rds")


# Load cleaned data
#df_result <- result # just so this won't be accidentally changed
df <- df_result

# Prepare data for PERMANOVA
species_matrix <- df %>%
  select(station, year, species, adjusted_density) %>%
  pivot_wider(
    names_from = species,
    values_from = adjusted_density,
    values_fill = 0
  ) %>%
  select(-station, -year)

env_data <- df %>%
  select(station, year) %>%
  distinct()

# Perform PERMANOVA
permanova <- adonis2(species_matrix ~ station + year, 
                    data = env_data, 
                    permutations = 999, 
                    method = "bray")

# Print results
print(permanova)

# Save results
capture.output(print(permanova), file = "output/permanova_results.txt")





wide_species <- df_result %>% filter(year != 1999) %>% pivot_wider(names_from = species, values_from = adjusted_density, values_fill = 0)
sediment_env_data$year <- as.numeric(sediment_env_data$year)
merged_data <- merge(wide_species, sediment_env_data, by = c("station", "year"))
species_columns <- setdiff(names(merged_data), c("station", "year", "grain_20um", "grain_63um", "grain_125um", "grain_250um", "grain_1000um", "organic_content", "organic_sd", "depth", "rotten_herring"))
perma_global <- adonis2(as.matrix(merged_data[, species_columns]) ~ year, data = merged_data, method = "bray")


# Initialize pairwise_results as an empty list
pairwise_results <- list()
years <- unique(wide_species$year)

# Loop over unique year pairs
for (i in 1:(length(years) - 1)) {
  for (j in (i + 1):length(years)) {
    # Subset data for the two years being compared
    subset_data <- merged_data[merged_data$year %in% c(years[i], years[j]), ]
    
    # Run PERMANOVA if there are sufficient data points for each year
    if (nrow(subset_data) > 1) { # Ensure subset_data has multiple rows
      perma_pair <- adonis2(as.matrix(subset_data[, species_columns]) ~ year, 
                            data = subset_data, method = "bray")
      # Store result with labeled names
      pairwise_results[[paste0(years[i], "-", years[j])]] <- perma_pair
    }
  }
}

pairwise_results_df <- data.frame(
  Year_Comparison = names(pairwise_results),
  R2 = sapply(pairwise_results, function(x) x$R2[1]),
  F_value = sapply(pairwise_results, function(x) x$F[1]),
  p_value = sapply(pairwise_results, function(x) x$`Pr(>F)`[1])
) 








sediment_env_data$year <- as.factor(sediment_env_data$year)

library(car)
leveneTest(organic_content ~ year, data = sediment_env_data)

# Let's do all assumption checks:
# 1. Shapiro-Wilk test for normality
shapiro.test(sediment_env_data$organic_content)
