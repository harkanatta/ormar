# Load required packages
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(patchwork)

# Load the data with feeding guilds
henda_long_matched_with_guilds <- henda_long_matched %>%
  left_join(
    merged_data_allt %>% select(Flokkun, `Feeding guild`),
    by = c("species" = "Flokkun")
  )

# Prepare data for NMDS
species_matrix <- henda_long_matched_with_guilds %>%
  select(station, year, species, adjusted_density) %>%
  pivot_wider(
    names_from = species,
    values_from = adjusted_density,
    values_fill = 0
  ) %>%
  mutate(
    station_temp = station,
    year_temp = year
  ) %>%
  select(-station, -year)

# Perform NMDS
nmds <- metaMDS(select(species_matrix, -station_temp, -year_temp), distance = "bray")

# Extract site scores and add metadata
site_scores <- as.data.frame(scores(nmds, display = "sites"))
site_scores$Station <- species_matrix$station_temp
site_scores$Year <- as.factor(species_matrix$year_temp)

# Extract species scores and add feeding guild information
species_scores <- as.data.frame(scores(nmds, display = "species"))
species_scores$species <- rownames(species_scores)
species_scores <- species_scores %>%
  left_join(
    henda_long_matched_with_guilds %>% 
      select(species, `Feeding guild`) %>% 
      distinct(),
    by = "species"
  )

# Calculate centroids for feeding guilds
guild_centroids <- species_scores %>%
  group_by(`Feeding guild`) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2),
    n = n()
  ) %>%
  filter(!is.na(`Feeding guild`))

# Create enhanced NMDS plot with feeding guilds
nmds_guild_plot <- ggplot() +
  # Add points for sites
  geom_point(data = site_scores, 
             aes(x = NMDS1, y = NMDS2, color = Station, shape = Year),
             size = 3) +
  # Add feeding guild centroids
  geom_label(data = guild_centroids,
            aes(x = NMDS1, y = NMDS2, 
                label = paste0(`Feeding guild`, "\n(n=", n, ")")),
            alpha = 0.7,
            size = 3) +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  labs(title = "NMDS of Community Composition with Feeding Guilds")

# Save NMDS plot
ggsave("output/nmds_feeding_guilds.png", 
       nmds_guild_plot, 
       width = 12, 
       height = 10,
       dpi = 300,
       bg = "white")

# Prepare data for PCA and PERMANOVA by feeding guild
guild_abundances <- henda_long_matched_with_guilds %>%
  group_by(station, year, `Feeding guild`) %>%
  summarise(total_abundance = sum(adjusted_density), .groups = "drop") %>%
  pivot_wider(
    names_from = `Feeding guild`,
    values_from = total_abundance,
    values_fill = 0
  )

# Create abundance matrix for PERMANOVA (non-scaled)
guild_matrix <- guild_abundances %>%
  select(-station, -year)  # Keep raw abundances for PERMANOVA

# Perform PERMANOVA on raw abundances
guild_permanova <- adonis2(
  guild_matrix ~ station + year, 
  data = guild_abundances,
  method = "bray",
  permutations = 999
)

# For PCA, scale the data separately
pca_data <- scale(guild_matrix)
pca_result <- rda(pca_data)

# Extract PCA scores
pca_site_scores <- data.frame(
  scores(pca_result, display = "sites"),
  Station = guild_abundances$station,
  Year = factor(guild_abundances$year)
)

pca_guild_scores <- data.frame(
  scores(pca_result, display = "species"),
  Guild = colnames(select(guild_abundances, -station, -year))
)

# Calculate variance explained
var_explained <- summary(pca_result)$cont$importance[2,1:2] * 100

# Create PCA plot
pca_guild_plot <- ggplot() +
  geom_point(data = pca_site_scores,
             aes(x = PC1, y = PC2, color = Station, shape = Year),
             size = 3) +
  geom_segment(data = pca_guild_scores,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "darkgrey") +
  geom_text_repel(data = pca_guild_scores,
                  aes(x = PC1, y = PC2, label = Guild),
                  size = 4) +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  ) +
  labs(
    title = "PCA of Feeding Guild Composition",
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2])
  )

# Save PCA plot
ggsave("output/pca_feeding_guilds.png", 
       pca_guild_plot, 
       width = 12, 
       height = 10,
       dpi = 300,
       bg = "white")

# Save PERMANOVA results
sink("output/feeding_guild_analysis.txt")
cat("Feeding Guild Analysis Results\n")
cat("============================\n\n")

cat("1. PERMANOVA Results\n")
cat("-----------------\n")
print(guild_permanova)

cat("\n2. Feeding Guild Composition\n")
cat("-------------------------\n")
guild_summary <- henda_long_matched_with_guilds %>%
  group_by(`Feeding guild`) %>%
  summarise(
    n_species = n_distinct(species),
    total_abundance = sum(adjusted_density),
    relative_abundance = total_abundance / sum(henda_long_matched_with_guilds$adjusted_density)
  ) %>%
  arrange(desc(relative_abundance))

print(guild_summary)
sink() 
