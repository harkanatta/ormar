# Station Similarity Analysis
# Analyzes spatial patterns and temporal changes in community similarity between stations

# Load required packages
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)

# Function to calculate station similarities
calculate_station_similarities <- function(species_data) {
  # Create dissimilarity matrix
  species_matrix_wide <- species_data %>%
    select(-year_temp) %>%
    group_by(station_temp) %>%
    summarise(across(everything(), mean)) %>%
    column_to_rownames("station_temp")
  
  # Calculate Bray-Curtis dissimilarity
  station_dist_matrix <- vegdist(species_matrix_wide, method = "bray") %>%
    as.matrix()
  
  return(station_dist_matrix)
}

# Function to analyze temporal patterns
analyze_temporal_patterns <- function(species_data) {
  species_data %>%
    group_by(year_temp) %>%
    group_split() %>%
    lapply(function(year_data) {
      year_matrix <- year_data %>%
        select(-year_temp) %>%
        column_to_rownames("station_temp")
      
      dist_matrix <- vegdist(year_matrix, method = "bray") %>%
        as.matrix()
      
      # Calculate E3-E4 and mean other distances
      e3e4 <- dist_matrix["E3", "E4"]
      other_pairs <- dist_matrix[lower.tri(dist_matrix)]
      other_mean <- mean(other_pairs[other_pairs != e3e4])
      
      data.frame(
        year = unique(year_data$year_temp),
        e3e4_dist = e3e4,
        other_mean_dist = other_mean
      )
    }) %>%
    bind_rows()
}

# Calculate overall station similarities
station_dist_matrix <- calculate_station_similarities(species_matrix)

# Calculate temporal patterns
yearly_patterns <- analyze_temporal_patterns(species_matrix)

# Create summary of station pairs
station_pairs_summary <- data.frame(
  pair = c("E3-E4", "A7-B8", "B5-C4", "E3-A7", "E3-B8", "E4-B8"),
  dissimilarity = c(
    station_dist_matrix["E3", "E4"],
    station_dist_matrix["A7", "B8"],
    station_dist_matrix["B5", "C4"],
    station_dist_matrix["E3", "A7"],
    station_dist_matrix["E3", "B8"],
    station_dist_matrix["E4", "B8"]
  )
) %>%
  arrange(dissimilarity)

# Create visualization of temporal patterns
temporal_plot <- ggplot(yearly_patterns, aes(x = year)) +
  geom_line(aes(y = e3e4_dist, color = "E3-E4")) +
  geom_line(aes(y = other_mean_dist, color = "Other Pairs Mean")) +
  labs(
    title = "Station Similarity Patterns Over Time",
    x = "Year",
    y = "Bray-Curtis Dissimilarity",
    color = "Station Pairs"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom"
  )

# Save results
write.csv(station_pairs_summary, "output/station_pairs_summary.csv", row.names = FALSE)
write.csv(yearly_patterns, "output/yearly_similarity_patterns.csv", row.names = FALSE)
ggsave("output/temporal_similarity_patterns.png", temporal_plot, width = 8, height = 6, dpi = 300)

# Print summary statistics
cat("\nOverall Station Pair Summary:\n")
print(station_pairs_summary)

cat("\nTemporal Patterns Summary:\n")
print(yearly_patterns) 