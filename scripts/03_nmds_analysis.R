# NMDS Analysis

# Load required packages
library(vegan)
library(ggplot2)

# Load cleaned data
#df_result <- result # just so this won't be accidentally changed
df <- df_result

# Prepare data for NMDS - exclude 1999 data
species_matrix <- df %>%
  filter(year >= 2013) %>%  # Only include 2013-2017 data
  select(station, year, species, adjusted_density) %>%
  pivot_wider(
    names_from = species,
    values_from = adjusted_density,
    values_fill = 0
  ) %>%
  # Store station and year before removing them for NMDS
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
site_scores$Year <- as.factor(species_matrix$year_temp)  # Convert Year to factor

# Extract species scores
species_scores <- as.data.frame(scores(nmds, display = "species"))

# Calculate species contributions (using relative abundance)
species_rel_abundance <- colSums(select(species_matrix, -station_temp, -year_temp)) / 
    sum(select(species_matrix, -station_temp, -year_temp))
important_species <- names(sort(species_rel_abundance[species_rel_abundance > 0.01], decreasing = TRUE))

# Filter species scores to only show important species
species_scores_filtered <- species_scores[rownames(species_scores) %in% important_species, ]

env_data <- readRDS("data/raw/sediment_env_data.rds") %>%
  select(station, year, 
         organic_content, depth,
         grain_20um, grain_63um, grain_125um, 
         grain_250um, grain_1000um) %>%
  distinct()

 sig_env_vectors <- scores(env_fit, "vectors") %>%
 as.data.frame() %>%
 mutate(variable = rownames(.)) %>%
 filter(env_fit$vectors$pvals < 0.05)

# Create the plot
nmds_plot <- ggplot() +
    # Base NMDS plot elements
    geom_path(data = site_scores, 
             aes(x = NMDS1, y = NMDS2, color = Station, group = Station),
             arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
    geom_point(data = site_scores, 
               aes(x = NMDS1, y = NMDS2, color = Station, shape = Year),
               size = 3) +
    # Add significant environmental vectors
    geom_segment(data = sig_env_vectors,
                aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                arrow = arrow(length = unit(0.2, "cm")),
                color = "darkgrey") +
    geom_text(data = sig_env_vectors,
              aes(x = NMDS1*1.1, y = NMDS2*1.1, label = variable),
              size = 3, color = "darkgrey") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white")) +
    labs(title = "NMDS of Community Composition (2013-2017)",
         caption = "")

#"Figure X. Non-metric multidimensional scaling (NMDS) ordination of benthic community composition in Kolgrafafjörður from 2013-2017. Points represent sampling stations in different years, with arrows showing temporal trajectories for each station. Grey arrows show significant sediment variables (p < 0.05): organic content (r² = 0.54, p = 0.001), grain size <20μm (r² = 0.23, p = 0.036), grain size <63μm (r² = 0.31, p = 0.007), and grain size <125μm (r² = 0.21, p = 0.027). Stress = 0.126. Based on Bray-Curtis dissimilarity of species abundance data."


# Save plot
ggsave("output/nmds_plot.png", nmds_plot, width = 10, height = 8, bg = "white")

# 1. Find dominant species (>1% relative abundance)
dominant_species <- species_matrix %>%
  select(where(is.numeric)) %>%  
  select(-year_temp, -year) %>%    # Exclude ALL year columns
  summarise(across(everything(), sum)) %>%  
  pivot_longer(everything(), 
               names_to = "species", 
               values_to = "total_abundance") %>%
  arrange(desc(total_abundance)) %>%  
  mutate(relative_abundance = total_abundance/sum(total_abundance)) %>%
  filter(relative_abundance > 0.01)  

# 3. Analyze temporal patterns of dominant species
temporal_patterns <- species_matrix %>%
  select(where(is.numeric)) %>%     
  select(-year) %>%                 
  group_by(year_temp) %>%          # Group by year first
  summarise(across(everything(), sum)) %>%  # Sum across stations
  pivot_longer(-year_temp,          
               names_to = "species", 
               values_to = "abundance") %>%
  filter(species %in% dominant_species$species) %>%  
  arrange(year_temp, desc(abundance))  

# Print results
print("Dominant Species (>1% relative abundance):")
print(dominant_species, n = 20)

print("\nTemporal Patterns by Year:")
print(temporal_patterns, n = 50)




# Merge with NMDS scores
nmds_env <- site_scores %>%
  left_join(env_data, by = c("Station" = "station", "Year" = "year"))

# Calculate correlation with axes
axis_correlations <- data.frame(
  Variable = c("Organic Content", "Depth",
               "Grain Size <20μm", "Grain Size <63μm", 
               "Grain Size <125μm", "Grain Size <250μm",
               "Grain Size <1000μm"),
  NMDS1 = c(
    cor(nmds_env$NMDS1, nmds_env$organic_content),
    cor(nmds_env$NMDS1, nmds_env$depth),
    cor(nmds_env$NMDS1, nmds_env$grain_20um),
    cor(nmds_env$NMDS1, nmds_env$grain_63um),
    cor(nmds_env$NMDS1, nmds_env$grain_125um),
    cor(nmds_env$NMDS1, nmds_env$grain_250um),
    cor(nmds_env$NMDS1, nmds_env$grain_1000um)
  ),
  NMDS2 = c(
    cor(nmds_env$NMDS2, nmds_env$organic_content),
    cor(nmds_env$NMDS2, nmds_env$depth),
    cor(nmds_env$NMDS2, nmds_env$grain_20um),
    cor(nmds_env$NMDS2, nmds_env$grain_63um),
    cor(nmds_env$NMDS2, nmds_env$grain_125um),
    cor(nmds_env$NMDS2, nmds_env$grain_250um),
    cor(nmds_env$NMDS2, nmds_env$grain_1000um)
  )
)

# Format correlations for printing
correlations_formatted <- axis_correlations %>%
  mutate(across(c(NMDS1, NMDS2), round, 3))

# Print formatted table
print(correlations_formatted)
