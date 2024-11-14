# Load required packages
library(dplyr)
library(ggplot2)
library(vegan)
library(tidyr)
library(tibble)  # For column_to_rownames()

df <- readRDS("data/cleaned_data.rds") %>% filter(station %in% c("C4", "A7", "B5", "B8", "E4", "E3")) 

# Prepare data for diversity calculations
species_matrix <- df %>%
  select(species, sample_id, density) %>%
  pivot_wider(
    names_from = species,
    values_from = density,
    values_fill = 0,
    values_fn = sum
  ) %>%
  column_to_rownames(var = "sample_id") %>%
  mutate(across(everything(), as.numeric))

# Calculate species richness and diversity indices
richness <- specnumber(species_matrix)
shannon <- diversity(species_matrix, index = "shannon")
simpson <- diversity(species_matrix, index = "simpson")

# Create diversity dataframe
diversity_df <- data.frame(
  sample_id = rownames(species_matrix),
  Richness = richness,
  Shannon = shannon,
  Simpson = simpson
) %>%
  left_join(
    df %>%
      select(sample_id, year, station) %>%
      distinct(),
    by = "sample_id"
  )

# Create output directory if it doesn't exist
if(!dir.exists("output")) {
  dir.create("output")
}

# Plot diversity indices
richness_plot <- ggplot(diversity_df, aes(x = year, y = Richness, color = station, group = station)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Species Richness Over Time by Station",
       x = "Year",
       y = "Species Richness")

shannon_plot <- ggplot(diversity_df, aes(x = year, y = Shannon, color = station, group = station)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Shannon Diversity Over Time by Station",
       x = "Year",
       y = "Shannon Diversity Index")

# Save plots
ggsave("output/richness_plot.png", richness_plot, width = 10, height = 6, bg = "white")
ggsave("output/shannon_plot.png", shannon_plot, width = 10, height = 6, bg = "white")

# Summary statistics
summary_stats <- diversity_df %>%
  group_by(year, station) %>%
  summarise(
    Mean_Richness = mean(Richness),
    SD_Richness = sd(Richness),
    Mean_Shannon = mean(Shannon),
    SD_Shannon = sd(Shannon),
    .groups = 'drop'
  )

# Save summary statistics
write.csv(summary_stats, "output/diversity_summary_stats.csv", row.names = FALSE)

# Boxplots for richness and Shannon diversity by year
richness_boxplot <- ggplot(diversity_df, aes(x = as.factor(year), y = Richness)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Species Richness Distribution by Year",
       x = "Year",
       y = "Species Richness")

shannon_boxplot <- ggplot(diversity_df, aes(x = as.factor(year), y = Shannon)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Shannon Diversity Distribution by Year",
       x = "Year",
       y = "Shannon Diversity Index")

# Save boxplots
ggsave("output/richness_boxplot.png", richness_boxplot, width = 10, height = 6, bg = "white")
ggsave("output/shannon_boxplot.png", shannon_boxplot, width = 10, height = 6, bg = "white")

# Species accumulation curve
species_accum <- specaccum(species_matrix)

# Save species accumulation plot
png("output/species_accumulation.png", width = 800, height = 600)
plot(species_accum, ci.type = "poly", col = "blue", lwd = 2, ci.lty = 0, ci.col = "lightblue",
     xlab = "Number of Samples", ylab = "Species Richness",
     main = "Species Accumulation Curve")
dev.off()

# Top 10 most abundant species
top_species <- df %>%
  group_by(species) %>%
  summarise(total_abundance = sum(density)) %>%
  slice_max(order_by = total_abundance, n = 10)

# Bar plot of top 10 species
top_species_plot <- ggplot(top_species, aes(x = reorder(species, total_abundance), y = total_abundance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 Most Abundant Species",
       x = "Species",
       y = "Total Abundance")

# Save top species plot
ggsave("output/top_species_plot.png", top_species_plot, width = 10, height = 6, bg = "white")
