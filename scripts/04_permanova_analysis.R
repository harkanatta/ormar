# PERMANOVA Analysis

# Load required packages
library(vegan)
library(tidyr)

# Load cleaned data
df <- result

df <- result %>%
    filter(year != 1999) %>%
    dplyr::group_by(station, year, species) %>%
    dplyr::summarise(N = sum(adjusted_density), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = species, values_from = N, values_fill = 0) %>%
    mutate(
        Station = factor(station, levels = c("C4", "E3", "E4", "B5", "B8", "A7")),
        Year = factor(year)
    ) %>%
    dplyr::select(-station, -year)  # Remove the original stod and Artal columns
    
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
