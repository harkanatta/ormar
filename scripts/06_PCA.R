library(ggrepel)
library(patchwork)

# Load cleaned data
#df_result <- result # just so this won't be accidentally changed
df <- df_result

# Prepare data including 1999
species_matrix_all <- df %>%
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


format_species_names <- function(species_names) {
  formatted_names <- sapply(species_names, function(name) {
    # If there's a space, it's already formatted (species name or family), leave it as is
    if (grepl(" ", name)) {
      return(name)
    }
    # If it's a single word, add "sp."
    paste0(name, " sp.")
  })
  return(formatted_names)
}


# PCA plot function
create_year_pca <- function(year_data, year, important_species) {
  # Perform PCA for this year
  species_pca <- rda(select(year_data, -station_temp, -year_temp))
  
  # Extract scores
  site_scores <- data.frame(
    scores(species_pca, display = "sites"),
    Station = year_data$station_temp
  )
  
  species_loadings <- data.frame(
    scores(species_pca, display = "species"),
    Species = rownames(scores(species_pca, display = "species"))
  ) %>%
    filter(Species %in% important_species) %>%
    mutate(Species = format_species_names(Species))  # Format species names
  
  # Calculate variance explained
  var_explained <- round(summary(species_pca)$cont$importance[2,1:2] * 100, 1)
  
  # Set consistent axis limits across all plots
  axis_limits <- max(abs(c(
    site_scores$PC1, site_scores$PC2,
    species_loadings$PC1, species_loadings$PC2
  ))) * 1.2  # Add 20% padding
  
  # Create PCA biplot
  pca_plot <- ggplot() +
    # Add site scores
    geom_point(data = site_scores,
               aes(x = PC1, y = PC2),
               size = 2,
               shape = 19,
               color = "black") +
    # Add site labels
    geom_text_repel(
      data = site_scores, 
      aes(x = PC1, y = PC2, label = Station),
      size = 3.5,
      fontface = "bold",
      box.padding = 0.5,
      point.padding = 0.3
    ) +
    # Add species arrows
    geom_segment(data = species_loadings,
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = "darkgrey") +
    # Add species labels
    geom_text_repel(data = species_loadings,
                    aes(x = PC1, y = PC2, label = Species),
                    size = 2.8,
                    fontface = "italic",
                    color = "grey30",
                    box.padding = 0.5,
                    point.padding = 0.3) +
    # Set consistent axis limits
    coord_fixed(xlim = c(-axis_limits, axis_limits),
               ylim = c(-axis_limits, axis_limits)) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey95"),
      axis.text = element_text(color = "black", size = 8),
      axis.title = element_text(size = 10),
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    ) +
    labs(
      title = paste("PCA of benthic community composition -", year),
      subtitle = "Stations shown as points, species as arrows",
      x = paste0("PC1 (", var_explained[1], "%)"),
      y = paste0("PC2 (", var_explained[2], "%)")
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
    subtitle = "Stations shown as points, species as arrows",
    caption = "Based on species abundance data. Only species with >1% relative abundance shown.",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      plot.caption = element_text(size = 10, hjust = 0)
    )
  )

ggsave("output/pca_biplots_all_years.png", 
       combined_plot, 
       width = 20, height = 15, 
       dpi = 300, 
       bg = "white")










library(ggrepel)
library(patchwork)

# Function to create reversed PCA biplot
create_year_pca_reversed <- function(year_data, year, important_species) {
    # Perform PCA
    species_pca <- rda(select(year_data, -station_temp, -year_temp))
    
    # Extract scores - note the reversed roles
    species_scores <- data.frame(
        scores(species_pca, display = "species"),
        Species = rownames(scores(species_pca, display = "species"))
    ) %>%
        mutate(Species = format_species_names(Species))  # Format species names
    
    station_loadings <- data.frame(
        scores(species_pca, display = "sites"),
        Station = year_data$station_temp
    )
    
    # Calculate variance explained
    var_explained <- round(summary(species_pca)$cont$importance[2,1:2] * 100, 1)
    
    # Set consistent axis limits across all plots
    axis_limits <- max(abs(c(
        species_scores$PC1, species_scores$PC2,
        station_loadings$PC1, station_loadings$PC2
    ))) * 1.2  # Add 20% padding
    
    # Create reversed PCA biplot
    pca_plot <- ggplot() +
        # Add species points
        geom_point(data = species_scores,
                   aes(x = PC1, y = PC2),
                   size = 2,
                   shape = 19,
                   color = "black") +
        # Add station arrows
        geom_segment(data = station_loadings,
                     aes(x = 0, y = 0, xend = PC1, yend = PC2),
                     arrow = arrow(length = unit(0.2, "cm")),
                     color = "darkgrey") +
        # Add station labels
        geom_text_repel(data = station_loadings,
                        aes(x = PC1, y = PC2, label = Station),
                        size = 3.5,
                        fontface = "bold",
                        box.padding = 0.5,
                        point.padding = 0.3) +
        # Add species labels
        geom_text_repel(data = species_scores,
                        aes(x = PC1, y = PC2, label = Species),
                        size = 2.8,
                        fontface = "italic",
                        color = "grey30",
                        box.padding = 0.5,
                        point.padding = 0.3) +
        # Set consistent axis limits
        coord_fixed(xlim = c(-axis_limits, axis_limits),
                    ylim = c(-axis_limits, axis_limits)) +
        theme_bw() +
        theme(
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "grey95"),
            axis.text = element_text(color = "black", size = 8),
            axis.title = element_text(size = 10),
            plot.title = element_text(size = 11, face = "bold"),
            plot.subtitle = element_text(size = 9),
            plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
        ) +
        labs(
            title = paste("PCA of benthic community composition -", year),
            subtitle = "Species shown as points, stations as arrows",
            x = paste0("PC1 (", var_explained[1], "%)"),
            y = paste0("PC2 (", var_explained[2], "%)")
        )
    
    return(pca_plot)
}

# Create reversed PCA plots for each year
years <- sort(unique(species_matrix_all$year_temp))
pca_plots_reversed <- list()

for(year in years) {
  year_data <- species_matrix_all %>%
    filter(year_temp == year)
  
  pca_plots_reversed[[as.character(year)]] <- create_year_pca_reversed(
    year_data, 
    year,
    dominant_species$species
  )
  
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

ggsave("output/pca_biplots_reversed_all_years.png", 
       combined_plot_reversed, 
       width = 20, height = 15, 
       dpi = 300, 
       bg = "white")
