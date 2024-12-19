library(ggrepel)
library(vegan)
library(readxl)
library(patchwork)

# Find dominant species (>1% relative abundance)
dominant_species <- henda %>%
  pivot_longer(cols = c(A7, B5, B8, C4, E3, E4), 
              names_to = "Station", 
              values_to = "abundance") %>%
  group_by(Flokkun) %>%
  summarise(total_abundance = sum(abundance)) %>%
  ungroup() %>%
  mutate(relative_abundance = total_abundance/sum(total_abundance)) %>%
  filter(relative_abundance > 0.01) %>%  # Keep species with >1% relative abundance
  rename(Species = Flokkun) %>%
  arrange(desc(relative_abundance))

# Print to verify
print("Dominant species identified:")
print(dominant_species)

# Function to process and create PCA for one year
create_year_pca <- function(year, data = henda) {
  # Calculate dominant species for this specific year
  year_dominant_species <- data %>%
    filter(Year == year) %>%
    pivot_longer(cols = c(A7, B5, B8, C4, E3, E4), 
                names_to = "Station", 
                values_to = "abundance") %>%
    group_by(Flokkun) %>%
    summarise(total_abundance = sum(abundance)) %>%
    ungroup() %>%
    mutate(relative_abundance = total_abundance/sum(total_abundance)) %>%
    filter(relative_abundance > 0.01) %>%  # Keep species with >1% relative abundance
    rename(Species = Flokkun)
  
  # Filter data for specific year and prepare
  year_data <- data %>%
    filter(Year == year) %>%
    select(-Year) %>%
    as.data.frame()
  
  # Prepare data for PCA
  rownames(year_data) <- year_data$Flokkun
  year_data <- year_data %>% select(-Flokkun)
  
  # Handle NAs
  year_data[is.na(year_data)] <- 0
  
  # Run PCA
  species_pca <- rda(scale(year_data))
  
  # Extract scores
  site_scores <- data.frame(
    scores(species_pca, display = "sites"),
    Species = rownames(scores(species_pca, display = "sites"))
  ) %>%
    # Add indicator for dominant species
    mutate(is_dominant = Species %in% year_dominant_species$Species)
  
  species_loadings <- data.frame(
    scores(species_pca, display = "species"),
    Station = colnames(year_data)
  )
  
  # Calculate variance explained
  var_explained <- round(summary(species_pca)$cont$importance[2,1:2] * 100, 1)
  
  # Define title and subtitle
  plot_title <- paste("PCA of benthic community composition -", year)
  plot_subtitle <- paste0("Species shown as points (>1% abundance labeled, n=", 
                          sum(site_scores$is_dominant), "), stations as arrows")
  
  # Print title and subtitle
  print(plot_title)
  print(plot_subtitle)
  
  # Create PCA plot
  pca_plot <- ggplot() +
    # Add points for all species
    geom_point(data = site_scores,
               aes(x = PC1, y = PC2),
               size = 2,
               shape = 19,
               color = "black") +
    # Add arrows for stations
    geom_segment(data = species_loadings,
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = "darkgrey") +
    # Add station labels
    geom_text_repel(data = species_loadings,
                    aes(x = PC1, y = PC2, label = Station),
                    size = 3.5,
                    fontface = "bold",
                    box.padding = 0.5,
                    point.padding = 0.3) +
    # Add labels only for dominant species
    geom_text_repel(data = filter(site_scores, is_dominant),
                    aes(x = PC1, y = PC2, label = Species),
                    size = 2.8,
                    fontface = "italic",
                    color = "grey30",
                    box.padding = 0.5,
                    point.padding = 0.3) +
    coord_fixed() +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey95"),
      axis.text = element_text(size = 14, color = "black"),     # Increased from 12
      axis.title = element_text(size = 16),
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10)
    ) +
    # Comment out title and subtitle in the plot
     labs(
    #   title = plot_title,
    #   subtitle = plot_subtitle,
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")
     )
  
  return(pca_plot)
}

# Create PCA plots for all years
years <- c(1999, 2013:2017)
pca_plots <- list()

for(year in years) {
  pca_plots[[as.character(year)]] <- create_year_pca(year)
  
  # Save individual plots
  ggsave(
    paste0("output/pca_biplot_", year, ".png"),
    pca_plots[[as.character(year)]], 
    width = 10, height = 8, 
    dpi = 300, 
    bg = "white"
  )
}

# Combine all plots with a layout that stretches them evenly
combined_plot <- wrap_plots(pca_plots, ncol = 2) +
  plot_layout(guides = "collect") +  # Ensure legends are collected and aligned
  plot_annotation(
    title = "Temporal changes in community composition (1999, 2013-2017)",
    subtitle = "Species shown as points, stations as arrows",
    caption = "Based on species abundance data",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 16, color = "black"),
      axis.title = element_text(size = 18),    
      plot.subtitle = element_text(size = 14),
      plot.caption = element_text(size = 10, hjust = 0)
    )
  )

# Save combined plot with A4 dimensions
ggsave("output/pca_biplots_all_years.png", 
       combined_plot, 
       width = 8.27, height = 11.69,  # A4 dimensions in inches
       dpi = 300, 
       bg = "white")

# After the PCA plots creation, add this section for results reporting

# Modified extract_pca_results function to calculate dominant species for each year
extract_pca_results <- function(year, species_pca) {
  # Calculate dominant species for this specific year
  year_dominant_species <- henda %>%
    filter(Year == year) %>%
    pivot_longer(cols = c(A7, B5, B8, C4, E3, E4), 
                names_to = "Station", 
                values_to = "abundance") %>%
    group_by(Flokkun) %>%
    summarise(total_abundance = sum(abundance)) %>%
    ungroup() %>%
    mutate(relative_abundance = total_abundance/sum(total_abundance)) %>%
    filter(relative_abundance > 0.01) %>%  # Keep species with >1% relative abundance
    rename(Species = Flokkun) %>%
    arrange(desc(relative_abundance))
  
  # Get variance explained
  var_explained <- round(summary(species_pca)$cont$importance[2,] * 100, 1)
  
  # Get species scores and sort by absolute contribution to PC1 and PC2
  species_scores <- scores(species_pca, display = "sites")
  top_contributors <- data.frame(
    Species = rownames(species_scores),
    PC1_loading = species_scores[,1],
    PC2_loading = species_scores[,2]
  ) %>%
    mutate(
      PC1_abs = abs(PC1_loading),
      PC2_abs = abs(PC2_loading)
    ) %>%
    # Join with year-specific dominant species information
    left_join(
      year_dominant_species %>% 
        select(Species, relative_abundance),
      by = "Species"
    )
  
  station_scores <- data.frame(
    scores(species_pca, display = "species"),
    Station = colnames(year_data)
  )
  
  return(list(
    year = year,
    variance = var_explained,
    top_species = top_contributors,
    dominant_species = year_dominant_species,  # Include the year's dominant species list
    stations = station_scores
  ))
}

# Store results for all years
pca_results <- list()
for(year in years) {
  # Run PCA again to get results
  year_data <- henda %>%
    filter(Year == year) %>%
    select(-Year) %>%
    as.data.frame()
  
  rownames(year_data) <- year_data$Flokkun
  year_data <- year_data %>% select(-Flokkun)
  
  species_pca <- rda(scale(year_data))
  pca_results[[as.character(year)]] <- extract_pca_results(year, species_pca)
}

# Create summary tables for the manuscript
# 1. Variance explained table
variance_table <- data.frame(
  Year = years,
  PC1 = sapply(pca_results, function(x) x$variance[1]),
  PC2 = sapply(pca_results, function(x) x$variance[2]),
  Total = sapply(pca_results, function(x) sum(x$variance[1:2]))
)

# Print variance table in publication format
cat("\nTable 1. Percentage of variance explained by the first two principal components\n")
print(knitr::kable(variance_table, 
                   digits = 1,
                   caption = "Percentage of variance explained by PCA axes"))

# 2. Top contributing species for each year
cat("\nTable 2. Top contributing species to community composition by year\n")
for(year in years) {
  cat(sprintf("\nYear %d:\n", year))
  print(knitr::kable(pca_results[[as.character(year)]]$top_species[1:5,],
                     digits = 3,
                     caption = sprintf("Top 5 contributing species in %d", year)))
}

# 3. Station patterns summary
cat("\nStation patterns across years:\n")
for(year in years) {
  cat(sprintf("\nYear %d:\n", year))
  stations <- pca_results[[as.character(year)]]$stations
  print(knitr::kable(stations, 
                     digits = 3,
                     caption = sprintf("Station coordinates in %d", year)))
}

# Save results to file for manuscript
sink("output/pca_results_summary.txt")
cat("PCA Results Summary for Marine Biology Manuscript\n")
cat("================================================\n\n")

cat("1. Variance Explained\n")
cat("-------------------\n")
print(knitr::kable(variance_table, digits = 1))

cat("\n\n2. Species Contributions by Year\n")
cat("------------------------------\n")
for(year in years) {
  cat(sprintf("\nYear %d:\n", year))
  
  # Print year-specific dominant species first
  cat("\nDominant Species for this year (>1% relative abundance):\n")
  print(knitr::kable(
    pca_results[[as.character(year)]]$dominant_species %>%
      mutate(relative_abundance = relative_abundance * 100),
    col.names = c("Species", "Total Abundance", "Relative Abundance (%)"),
    digits = 2
  ))
  
  # Print PCA results for dominant species
  cat("\nPCA loadings for dominant species:\n")
  dominant_results <- pca_results[[as.character(year)]]$top_species %>%
    filter(!is.na(relative_abundance)) %>%
    arrange(desc(PC1_abs + PC2_abs)) %>%
    select(Species, PC1_loading, PC2_loading, relative_abundance) %>%
    mutate(relative_abundance = relative_abundance * 100)
  print(knitr::kable(dominant_results, digits = 3))
  
  # Print all species contributions
  cat("\nTop Contributing Species (including rare species):\n")
  all_results <- pca_results[[as.character(year)]]$top_species %>%
    arrange(desc(PC1_abs + PC2_abs)) %>%
    select(Species, PC1_loading, PC2_loading) %>%
    slice_head(n = 10)
  print(knitr::kable(all_results, digits = 3))
}

cat("\n\n3. Station Patterns\n")
cat("-----------------\n")
for(year in years) {
  cat(sprintf("\nYear %d:\n", year))
  print(knitr::kable(pca_results[[as.character(year)]]$stations, digits = 3))
}

# Add summary of overall dominant species
cat("\n\n4. Overall Community Dominance Patterns\n")
cat("-----------------------------------\n")
print(knitr::kable(
  dominant_species %>%
    arrange(desc(relative_abundance)) %>%
    mutate(
      relative_abundance = relative_abundance * 100,
      cumulative_abundance = cumsum(relative_abundance)
    ),
  col.names = c("Species", "Total Abundance", "Relative Abundance (%)", "Cumulative Abundance (%)"),
  digits = 2
))

sink()
