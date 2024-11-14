library(readr)

# Read in the henda dataset, specifying column types
henda <- read_csv(
  "data/processed/henda.csv",
  col_types = cols(
    Year = col_integer(),
    Flokkun = col_character(),  # Ensure Flokkun is read as character
    A7 = col_integer(),
    B5 = col_integer(),
    B8 = col_integer(),
    C4 = col_integer(),
    E3 = col_integer(),
    E4 = col_integer()
  ),
  na = c("", "NA")
) %>%
  mutate(across(where(is.numeric), ~ replace_na(., 0)))  %>% # Replace NA with 0 for numeric columns
  select(Flokkun:last_col()) 

# View the dataset to confirm changes
print(henda, n=1000)



henda_long <- henda %>%
  pivot_longer(
    cols = A7:E4,               # Specify the columns to pivot
    names_to = "Station",       # New column for the station names
    values_to = "Count"         # New column for the count values
  ) %>%
  rename(
    station = Station,          # Rename Station to station
    year = Year,                # Rename Year to year
    species = Flokkun,          # Rename Flokkun to species
    adjusted_density = Count    # Rename Count to adjusted_density
  )

species_matrix_all <- henda_long %>%
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

# Find dominant species (>1% relative abundance) including 1999
dominant_species <- species_matrix_all %>%
  select(where(is.numeric)) %>%  
  select(-year_temp) %>%
  summarise(across(everything(), sum)) %>%  
  pivot_longer(everything(), 
               names_to = "species", 
               values_to = "total_abundance") %>%
  arrange(desc(total_abundance)) %>%  
  mutate(relative_abundance = total_abundance/sum(total_abundance)) %>%
  filter(relative_abundance > 0.01)









# PCA plot function with scaling
create_year_pca <- function(year_data, year, important_species) {
  # Scale the data
  scaled_data <- scale(year_data)
  
  # Perform PCA for this year
  species_pca <- rda(scaled_data)  # Use scaled data
  
  # Extract scores
  site_scores <- data.frame(
    scores(species_pca, display = "sites")
  )
  
  species_loadings <- data.frame(
    scores(species_pca, display = "species"),
    Species = rownames(scores(species_pca, display = "species"))
  ) %>%
    filter(Species %in% important_species)
  
  # Calculate variance explained
  var_explained <- round(summary(species_pca)$cont$importance[2,1:2] * 100, 1)
  
  # Create PCA biplot
  pca_plot <- ggplot() +
    # Add site scores
    geom_point(data = site_scores, 
               aes(x = PC1, y = PC2),
               size = 3) +
    # Add species arrows
    geom_segment(data = species_loadings,
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "cm"))) +
    # Add species labels
    geom_text_repel(data = species_loadings,
                    aes(x = PC1, y = PC2, label = Species),
                    size = 3) +
    theme_bw() +
    labs(
      title = paste("PCA of benthic community composition -", year),
      x = paste0("PC1 (", var_explained[1], "%)"),
      y = paste0("PC2 (", var_explained[2], "%)"),
      subtitle = "Arrows show species with >1% relative abundance"
    )
  
  return(pca_plot)
}

# Create PCA plots for each year
years <- sort(unique(species_matrix_all$year_temp))
pca_plots <- list()

for(year in years) {
  year_data <- species_matrix_all %>%
    filter(year_temp == year)
  
  pca_plots[[as.character(year)]] <- create_year_pca(
    year_data, 
    year,
    dominant_species$species
  )
  
  ggsave(
    paste0("output/pca_biplot_", year, ".png"),
    pca_plots[[as.character(year)]], 
    width = 10, height = 8, 
    dpi = 300, 
    bg = "white"
  )
}

# Combine all plots
combined_plot <- wrap_plots(pca_plots, ncol = 2) +
  plot_annotation(
    title = "Temporal changes in community composition (1999, 2013-2017)",
    caption = "Based on species abundance data. Only species with >1% relative abundance shown."
  )

ggsave("output/pca_biplots_all_years.png", 
       combined_plot, 
       width = 20, height = 15, 
       dpi = 300, 
       bg = "white")

# Combine all plots
combined_plot <- wrap_plots(pca_plots, ncol = 2) +
  plot_annotation(
    title = "Temporal changes in community composition (1999, 2013-2017)",
    caption = "Based on species abundance data. Only species with >1% relative abundance shown."
  )

ggsave("output/pca_biplots_all_years.png", 
       combined_plot, 
       width = 20, height = 15, 
       dpi = 300, 
       bg = "white")










library(ggrepel)
library(patchwork)

# PCA plot function with reversed roles
# Function to create reversed PCA biplot with scaled arrows
create_year_pca_reversed <- function(year_data, year, important_species) {
  # Scale the data
  scaled_data <- scale(year_data)
  
  # Perform PCA for this year
  species_pca <- rda(scaled_data)  # Use scaled data
  
  # Extract scores - note the reversed roles
  species_scores <- data.frame(
    scores(species_pca, display = "species"),
    Species = rownames(scores(species_pca, display = "species"))
  )
  
  station_loadings <- data.frame(
    scores(species_pca, display = "sites"),
    Station = rownames(scores(species_pca, display = "sites"))
  )
  
  # Calculate variance explained
  var_explained <- round(summary(species_pca)$cont$importance[2,1:2] * 100, 1)
  
  # Create reversed PCA biplot
  pca_plot <- ggplot() +
    # Add species points
    geom_point(data = species_scores, 
               aes(x = PC1, y = PC2),
               size = 3) +
    # Add station arrows with scaling
    geom_segment(data = station_loadings,
                 aes(x = 0, y = 0, xend = PC1 * 0.1, yend = PC2 * 0.1, color = Station),  # Scale arrows
                 arrow = arrow(length = unit(0.2, "cm"))) +
    # Add species labels
    geom_text_repel(data = species_scores,
                    aes(x = PC1, y = PC2, label = Species),
                    size = 3) +
    theme_bw() +
    labs(
      title = paste("PCA of benthic community composition -", year),
      subtitle = "Species shown as points, stations as arrows",
      x = paste0("PC1 (", var_explained[1], "%)"),
      y = paste0("PC2 (", var_explained[2], "%)")
    ) +
    coord_fixed()  # Ensure equal scaling of axes
  
  return(pca_plot)
}

# Create reversed PCA plots for each year
years <- sort(unique(species_matrix_all$year_temp))
pca_plots_reversed <- list()

for(year in years) {
  year_data <- species_matrix_all %>%
    filter(year_temp == year)
  
  # Call the reversed PCA function
  pca_plots_reversed[[as.character(year)]] <- create_year_pca_reversed(
    year_data, 
    year,
    dominant_species$species
  )
  
  # Save the reversed PCA plot with a unique filename
  ggsave(
    paste0("output/pca_biplot_reversed_", year, ".png"),
    pca_plots_reversed[[as.character(year)]], 
    width = 10, height = 8, 
    dpi = 300, 
    bg = "white"
  )
}

# Combine all reversed plots
combined_plot_reversed <- wrap_plots(pca_plots_reversed, ncol = 2) +
  plot_annotation(
    title = "Temporal changes in community composition (1999, 2013-2017)",
    subtitle = "Species shown as points, stations as arrows",
    caption = "Based on species abundance data. Only species with >1% relative abundance shown."
  )

# Save the combined reversed plots with a unique filename
ggsave("output/pca_biplots_reversed_all_years.png", 
       combined_plot_reversed, 
       width = 20, height = 15, 
       dpi = 300, 
       bg = "white")

