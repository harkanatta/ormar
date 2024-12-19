#===============================================================================
# SIMPER Analysis for Trait-Based Community Changes
# Author: [Your Name]
# Date: [Current Date]
# Description: Analyzes temporal changes in benthic community traits using SIMPER
#===============================================================================

#-------------------------------------------------------------------------------
# Load Required Libraries
#-------------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(vegan)
library(kableExtra)
library(ggplot2)
library(stringr)

#-------------------------------------------------------------------------------
# Helper Functions
#-------------------------------------------------------------------------------
# Read and parse trait descriptions from text.csv
parse_trait_descriptions <- function() {
  text <- read.csv("data/traits_data/text.csv")
  
  # Extract feeding guild descriptions (contains all components)
  feeding_guild_text <- text$text[text$find_name == "feeding_guild"]
  
  # Parse position components from feeding guild text
  position_matches <- stringr::str_extract_all(
    feeding_guild_text,
    "\\((\\w+)\\s*\\((\\w+)\\)(?=[,;])"
  )[[1]]
  position_lookup <- setNames(
    c("Epibenthic", "Surface", "Subsurface"),
    c("EP", "SR", "SS")
  )
  
  # Parse feeding type components from feeding guild text
  feeding_type_matches <- stringr::str_extract_all(
    feeding_guild_text,
    "([^,;(]+?)\\s*\\(([^)]+?)\\s*;\\s*([^)]+?)\\)"
  )[[1]]
  feeding_type_lookup <- setNames(
    sapply(stringr::str_match_all(feeding_type_matches, "([^(]+)\\s*\\(([^;]+);\\s*([^)]+)\\)"), 
           function(x) trimws(x[1])),
    sapply(stringr::str_match_all(feeding_type_matches, "\\(([^;)]+)\\)"), 
           function(x) trimws(x[2]))
  )
  
  # Parse size/type components from feeding guild text
  size_matches <- stringr::str_extract_all(
    feeding_guild_text,
    "([^,(]+?)\\s*\\((?:e\\.g\\.)?\\s*[^,(]*?([^,)]+?)\\)"
  )[[1]]
  size_lookup <- setNames(
    sapply(stringr::str_match_all(size_matches, "([^(]+)\\s*\\("), 
           function(x) trimws(x[2])),
    sapply(stringr::str_match_all(size_matches, "\\(([^)]+)\\)"), 
           function(x) trimws(x[2]))
  )
  
  # Print parsed values for verification
  cat("\nParsed from text.csv:\n")
  cat("\nPosition components:\n")
  str(position_lookup)
  cat("\nFeeding type components:\n")
  str(feeding_type_lookup)
  cat("\nSize/type components:\n")
  str(size_lookup)
  
  list(
    position = position_lookup,
    feeding_type = feeding_type_lookup,
    size = size_lookup
  )
}

# Get lookup tables from text.csv
trait_lookups <- parse_trait_descriptions()

# Decode feeding guild abbreviations
decode_feeding_guild <- function(guild) {
  # Use parsed lookup tables
  position_lookup <- trait_lookups$position
  feeding_type_lookup <- trait_lookups$feeding_type
  size_lookup <- trait_lookups$size
  
  # Function to decode a single guild
  decode_single_guild <- function(single_guild) {
    if (is.na(single_guild)) return(NA)
    
    # Handle compound feeding modes (e.g., "Br/Gr")
    if (!grepl("-", single_guild)) {
      # For feed_mode column that uses "/" separator
      components <- strsplit(single_guild, "/")[[1]]
      decoded <- sapply(components, function(x) feeding_type_lookup[[x]] %||% x)
      return(paste(decoded, collapse = "/"))
    }
    
    # For feeding_guild column that uses "-" separator
    components <- strsplit(single_guild, "-")[[1]]
    
    # Decode each component based on its position
    decoded <- character(length(components))
    
    for (i in seq_along(components)) {
      component <- components[i]
      # Strip any whitespace
      component <- trimws(component)
      if (i == 1) {
        # First component is usually position
        decoded[i] <- position_lookup[[component]] %||% component
      } else if (i == length(components)) {
        # Last component is usually size
        decoded[i] <- size_lookup[[component]] %||% component
      } else {
        # Middle components are usually feeding type
        decoded[i] <- feeding_type_lookup[[component]] %||% component
      }
    }
    
    # Join with spaces
    paste(decoded, collapse = " ")
  }
  
  # Apply the decoding function to each guild in the input vector
  sapply(guild, decode_single_guild)
}

#-------------------------------------------------------------------------------
# Data Import and Preprocessing
#-------------------------------------------------------------------------------
# Read the merged data with feeding guilds
merged_data_allt <- read.csv("data/merged_data_alltB.csv")

# Define taxa to remove
taxa_to_remove <- c(
  "Foraminifera",
  "Nematoda",
  "Cirripedia",
  "Porifera",
  "Cnidaria",
  "Bryozoa",
  "Sipuncula",
  "Platyhelminthes",
  "Nemertea",
  #"Oligochaeta",
  "Ostracoda"
)

# Filter out unwanted taxa and then filter for specific stations
merged_data <- merged_data_allt %>%
  filter(!matched_taxon_name %in% taxa_to_remove) %>%
  subset(stod %in% c("C4", "A7", "B5", "B8", "E4", "E3")) %>%
  janitor::clean_names()

#-------------------------------------------------------------------------------
# SIMPER Analysis Setup
#-------------------------------------------------------------------------------
# Aggregate data by year and feeding guilds (combining all stations)
yearly_guild_data <- merged_data %>%
  group_by(artal, feeding_guild) %>%
  summarise(total_abundance = sum(nfm), .groups = "drop") %>%
  pivot_wider(
    names_from = feeding_guild,
    values_from = total_abundance,
    values_fill = 0
  )

# Create guild matrix (only numeric columns)
guild_matrix <- yearly_guild_data %>% 
  select(-artal) %>% 
  as.matrix()

# Check data structure before SIMPER
print("Checking data structure:")
print(str(yearly_guild_data))
print("Number of years:")
print(years)
print("Guild matrix dimensions:")
print(dim(guild_matrix))

# Perform SIMPER analysis for consecutive year comparisons
years <- sort(unique(yearly_guild_data$artal))
sequential_simper <- list()

for(i in 1:(length(years)-1)) {
  year_pair <- years[i:(i+1)]
  
  cat("\nProcessing years:", year_pair[1], "vs", year_pair[2], "\n")
  
  # Get indices for the current year pair
  year1_idx <- which(yearly_guild_data$artal == year_pair[1])
  year2_idx <- which(yearly_guild_data$artal == year_pair[2])
  
  # Extract data for these years
  year_data <- guild_matrix[c(year1_idx, year2_idx), ]
  
  # Only proceed if we have valid data
  if(all(rowSums(year_data) > 0)) {
    # Create a proper community data matrix
    comm_data <- as.data.frame(year_data)
    
    # Create a proper group factor
    groups <- factor(c(year_pair[1], year_pair[2]))
    
    # Remove any columns that sum to zero (species not present in either year)
    comm_data <- comm_data[, colSums(comm_data) > 0]
    
    cat("Dimensions after removing zero columns:", dim(comm_data), "\n")
    cat("Number of non-zero columns:", sum(colSums(comm_data) > 0), "\n")
    
    comparison_name <- paste(year_pair[1], year_pair[2], sep="_vs_")
    
    tryCatch({
      # Perform SIMPER on raw data
      sim_result <- simper(comm_data, 
                           group = groups,
                           permutations = 999,
                           trace = FALSE)
      
      sequential_simper[[comparison_name]] <- sim_result
      
      # Print summary of results
      cat("\nSIMPER results for", comparison_name, ":\n")
      sim_summary <- summary(sim_result)
      print(sim_summary)
      
      # Save detailed results to a CSV file
      results_df <- data.frame(
        Comparison = comparison_name,
        Species = rownames(sim_summary$species),
        Average_Contribution = sim_summary$species[, "average"],
        Overall_Contribution = sim_summary$species[, "cumsum"],
        stringsAsFactors = FALSE
      )
      
      write.csv(results_df,
                file = paste0("output/simper_results_", comparison_name, ".csv"),
                row.names = FALSE)
      
    }, error = function(e) {
      cat("\nError in SIMPER analysis for", comparison_name, ":\n")
      cat("Error message:", e$message, "\n")
      cat("Data summary:\n")
      print(summary(comm_data))
    })
  } else {
    warning(paste("Zero abundance rows found in years:", 
                  year_pair[1], "vs", year_pair[2]))
  }
}

# After the loop, create a summary of all comparisons
if(length(sequential_simper) > 0) {
  cat("\nOverall SIMPER Analysis Summary:\n")
  cat("Number of successful comparisons:", length(sequential_simper), "\n")
  cat("Completed comparisons:", paste(names(sequential_simper), collapse=", "), "\n")
  
  # Create a summary plot of the top contributing species across all comparisons
  all_results <- data.frame()
  
  for(comp_name in names(sequential_simper)) {
    sim_summary <- summary(sequential_simper[[comp_name]])
    temp_df <- data.frame(
      Comparison = comp_name,
      Species = rownames(sim_summary$species),
      Contribution = sim_summary$species[, "average"],
      stringsAsFactors = FALSE
    )
    all_results <- rbind(all_results, temp_df)
  }
  
  # Save overall summary
  write.csv(all_results,
            file = "output/simper_results_summary.csv",
            row.names = FALSE)
  
  # Create summary plot
  contrib_plot <- ggplot(all_results, 
                         aes(x = reorder(Species, Contribution), 
                             y = Contribution)) +
    geom_boxplot() +
    coord_flip() +
    theme_minimal() +
    labs(title = "Species Contributions to Community Dissimilarity",
         x = "Species",
         y = "Average Contribution")
  
  print(contrib_plot)
  ggsave("output/simper_contributions_plot.png", 
         contrib_plot, 
         width = 12, 
         height = 8,
         dpi = 300,
         bg = "white")
}

#-------------------------------------------------------------------------------
# Visualization Functions
#-------------------------------------------------------------------------------
# Helper functions for creating plots
create_temporal_plot <- function(data, trait_column) {
  ggplot(data, 
         aes(x = factor(artal), 
             y = relative_abundance, 
             fill = !!sym(trait_column))) +
    geom_bar(stat = "identity", position = "stack") +
    theme_minimal() +
    labs(
      title = paste("Temporal Changes in", tools::toTitleCase(gsub("_", " ", trait_column)), "Composition"),
      x = "Year",
      y = "Relative Abundance",
      fill = tools::toTitleCase(gsub("_", " ", trait_column))
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
}

create_contribution_plot <- function(simper_results, trait_column) {
  if (length(simper_results) > 0) {
    contribution_data <- lapply(simper_results, function(x) {
      if (!is.null(x)) {
        summary_data <- summary(x)
        data.frame(
          trait = rownames(summary_data$species),
          contribution = summary_data$species[, "average"],
          comparison = rep(names(x), nrow(summary_data$species)),
          stringsAsFactors = FALSE
        )
      }
    }) %>% bind_rows()
    
    ggplot(contribution_data, 
           aes(x = reorder(trait, contribution), y = contribution)) +
      geom_boxplot(fill = "lightblue") +
      coord_flip() +
      theme_minimal() +
      labs(
        title = paste(tools::toTitleCase(gsub("_", " ", trait_column)), 
                      "Contributions to Community Dissimilarity"),
        x = tools::toTitleCase(gsub("_", " ", trait_column)),
        y = "Average Contribution"
      )
  } else {
    warning("No valid SIMPER results available for plotting")
    NULL
  }
}

create_change_plot <- function(data, trait_column) {
  yearly_changes <- data %>%
    group_by(!!sym(trait_column)) %>%
    mutate(
      change = relative_abundance - lag(relative_abundance)
    ) %>%
    filter(!is.na(change))
  
  ggplot(yearly_changes, 
         aes(x = artal, y = change, fill = !!sym(trait_column))) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(
      title = paste("Year-to-Year Changes in", 
                    tools::toTitleCase(gsub("_", " ", trait_column)), 
                    "Composition"),
      x = "Year",
      y = "Change in Relative Abundance",
      fill = tools::toTitleCase(gsub("_", " ", trait_column))
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

#-------------------------------------------------------------------------------
# Create and Save Visualizations
#-------------------------------------------------------------------------------
# 1. Temporal changes in feeding guild composition
yearly_guild_summary <- merged_data %>%
  group_by(artal, feeding_guild) %>%
  summarise(
    total_abundance = sum(nfm),
    .groups = "drop"
  ) %>%
  group_by(artal) %>%
  mutate(relative_abundance = total_abundance / sum(total_abundance))

# Create temporal plot (this should work regardless of SIMPER results)
temporal_plot <- ggplot(yearly_guild_summary, 
                       aes(x = factor(artal), 
                           y = relative_abundance, 
                           fill = decode_feeding_guild(`feeding_guild`))) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(
    title = "Temporal Changes in Feeding Guild Composition",
    x = "Year",
    y = "Relative Abundance",
    fill = "Feeding Guild"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.text = element_text(size = 8)  # Smaller text to fit decoded names
  )

# Display the temporal plot
print(temporal_plot)

# Save the temporal plot
ggsave("output/feeding_guild_temporal_change.png", 
       temporal_plot, 
       width = 12, 
       height = 8,
       dpi = 300,
       bg = "white")

# 2. Contribution of each feeding guild to dissimilarity
if(length(sequential_simper) > 0) {
  contribution_data <- lapply(sequential_simper, function(x) {
    if (!is.null(x)) {
      summary_data <- summary(x)
      data.frame(
        guild = rownames(summary_data$species),
        contribution = summary_data$species[, "average"],
        comparison = rep(names(x), nrow(summary_data$species)),
        stringsAsFactors = FALSE
      )
    }
  }) %>% bind_rows()
  
  # Only create plot if we have valid contribution data
  if(nrow(contribution_data) > 0) {
    contribution_plot <- ggplot(contribution_data, 
                                aes(x = reorder(guild, contribution), 
                                    y = contribution)) +
      geom_boxplot(fill = "lightblue") +
      coord_flip() +
      theme_minimal() +
      labs(
        title = "Feeding Guild Contributions to Community Dissimilarity",
        x = "Feeding Guild",
        y = "Average Contribution"
      )
    
    # Display the plot in the viewer
    print(contribution_plot)
    
    # Save the plot as a PNG file
    ggsave("output/feeding_guild_contribution.png", 
           contribution_plot, 
           width = 12, 
           height = 8,
           dpi = 300,
           bg = "white")
  } else {
    warning("No valid contribution data available for plotting")
  }
} else {
  warning("No valid SIMPER results available for plotting")
}

# 3. Year-to-year changes in guild composition
yearly_changes <- yearly_guild_summary %>%
  group_by(feeding_guild) %>%
  mutate(
    change = relative_abundance - lag(relative_abundance)
  ) %>%
  filter(!is.na(change))

change_plot <- ggplot(yearly_changes, 
                      aes(x = artal, y = change, fill = feeding_guild)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(
    title = "Year-to-Year Changes in Feeding Guild Composition",
    x = "Year",
    y = "Change in Relative Abundance",
    fill = "Feeding Guild"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Display the plot in the viewer
print(change_plot)

# Save the plot as a PNG file
ggsave("output/feeding_guild_changes.png", 
       change_plot, 
       width = 12, 
       height = 8,
       dpi = 300,
       bg = "white")

# Create a summary table of feeding guild changes
guild_summary_table <- yearly_guild_summary %>%
  pivot_wider(
    names_from = artal,
    values_from = relative_abundance,
    values_fill = 0
  ) %>%
  arrange(desc(total_abundance))

# Save outputs
ggsave("output/feeding_guild_temporal_changes.png", 
       temporal_plot, 
       width = 12, 
       height = 8,
       dpi = 300,
       bg = "white")

write.csv(guild_summary_table, 
          "output/feeding_guild_summary.csv", 
          row.names = FALSE)

# Compare feeding guilds between 1999 and 2017
guild_comparison <- yearly_guild_summary %>%
  filter(artal %in% c(1999, 2017)) %>%
  pivot_wider(
    names_from = artal,
    values_from = relative_abundance,
    names_prefix = "year_"
  ) %>%
  mutate(
    year_1999 = replace_na(year_1999, 0),
    year_2017 = replace_na(year_2017, 0),
    difference = year_2017 - year_1999
  ) %>%
  arrange(difference)

# Create comparison plot
comparison_plot <- ggplot(guild_comparison, 
                          aes(x = reorder(feeding_guild, -difference), 
                              y = difference)) +
  geom_bar(stat = "identity", 
           aes(fill = difference < 0)) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Changes in Feeding Guild Composition: 2017 vs 1999",
    x = "Feeding Guild",
    y = "Change in Relative Abundance (2017 - 1999)",
    fill = "Decreased"
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  )

print(comparison_plot)

# Print summary of missing guilds
missing_guilds <- guild_comparison %>%
  filter(year_1999 > 0, year_2017 == 0) %>%
  select(feeding_guild, year_1999)

if(nrow(missing_guilds) > 0) {
  cat("\nFeeding guilds present in 1999 but absent in 2017:\n")
  print(missing_guilds)
}

# Create detailed feeding guild summary
guild_species_summary <- merged_data %>%
  group_by(feeding_guild, matched_taxon_name) %>%
  summarise(
    total_abundance = sum(nfm),
    years_present = n_distinct(artal),
    first_year = min(artal),
    last_year = max(artal),
    .groups = "drop"
  ) %>%
  arrange(feeding_guild, desc(total_abundance)) %>%
  group_by(feeding_guild) %>%
  mutate(
    rank_in_guild = row_number(),
    prop_of_guild = total_abundance / sum(total_abundance)
  ) %>%
  filter(rank_in_guild <= 5)  # Top 5 species per guild

# Create year-specific summaries
year_specific_changes <- merged_data %>%
  group_by(artal, feeding_guild, matched_taxon_name) %>%
  summarise(
    abundance = sum(nfm),
    .groups = "drop"
  ) %>%
  group_by(feeding_guild, matched_taxon_name) %>%
  summarise(
    years_found = paste(sort(unique(artal)), collapse = ", "),
    total_abundance = sum(abundance),
    .groups = "drop"
  ) %>%
  arrange(feeding_guild, desc(total_abundance))

# Save summaries to files
write.csv(guild_species_summary, 
          "output/feeding_guild_species_summary.csv", 
          row.names = FALSE)

write.csv(year_specific_changes,
          "output/feeding_guild_temporal_changes_by_species.csv",
          row.names = FALSE)

# Print summary to console
cat("\nTop species by feeding guild:\n")
print(guild_species_summary, n = 50)

# Calculate and print guild-level statistics
guild_stats <- merged_data %>%
  group_by(feeding_guild) %>%
  summarise(
    total_species = n_distinct(matched_taxon_name),
    total_abundance = sum(nfm),
    mean_abundance_per_year = total_abundance / n_distinct(artal),
    years_present = n_distinct(artal),
    .groups = "drop"
  ) %>%
  arrange(desc(total_abundance))

cat("\nFeeding guild summary statistics:\n")
print(guild_stats)

#-------------------------------------------------------------------------------
# Trait Analysis Functions
#-------------------------------------------------------------------------------
perform_trait_analysis <- function(data, trait_column) {
  # Validate input trait column
  valid_traits <- c("food_source", "motility", "habit", 
                    "om_ca_he", "food_size_type", "feed_mode", "feeding_guild")
  
  if (!trait_column %in% valid_traits) {
    stop("Invalid trait column. Must be one of: ", paste(valid_traits, collapse = ", "))
  }
  
  # Create output directory
  dir.create("output", showWarnings = FALSE)
  
  # Check if trait column has data
  if (all(is.na(data[[trait_column]]))) {
    stop("No data found in column: ", trait_column)
  }
  
  # Print unique values in trait column for debugging
  cat("\nUnique values in", trait_column, ":\n")
  print(unique(data[[trait_column]]))
  
  # Aggregate data by year and selected trait
  yearly_trait_data <- data %>%
    filter(!is.na(!!sym(trait_column))) %>%  # Remove NA values
    group_by(artal, !!sym(trait_column)) %>%
    summarise(total_abundance = sum(nfm), .groups = "drop") %>%
    pivot_wider(
      names_from = !!sym(trait_column),
      values_from = total_abundance,
      values_fill = 0
    )
  
  # Print structure of aggregated data for debugging
  cat("\nStructure of yearly trait data:\n")
  print(str(yearly_trait_data))
  
  # Create trait matrix
  trait_matrix <- yearly_trait_data %>% 
    select(-artal) %>% 
    as.matrix()
  
  # Check if matrix has data
  if (ncol(trait_matrix) == 0 || nrow(trait_matrix) == 0) {
    stop("No valid data for SIMPER analysis after aggregation")
  }
  
  # Perform SIMPER analysis
  years <- sort(unique(yearly_trait_data$artal))
  sequential_simper <- list()
  
  for(i in 1:(length(years)-1)) {
    year_pair <- years[i:(i+1)]
    year1_idx <- which(yearly_trait_data$artal == year_pair[1])
    year2_idx <- which(yearly_trait_data$artal == year_pair[2])
    
    year_data <- trait_matrix[c(year1_idx, year2_idx), ]
    year_groups <- factor(c(year_pair[1], year_pair[2]))
    
    year_data_norm <- decostand(year_data, method = "total")
    
    comparison_name <- paste(year_pair[1], year_pair[2], sep="_vs_")
    tryCatch({
      sequential_simper[[comparison_name]] <- simper(year_data_norm, 
                                                     group = year_groups, 
                                                     permutations = 999)
    }, error = function(e) {
      warning(paste("Error in SIMPER analysis for years", year_pair[1], 
                    "vs", year_pair[2], ":", e$message))
    })
  }
  
  # Create visualizations
  yearly_trait_summary <- data %>%
    group_by(artal, !!sym(trait_column)) %>%
    summarise(
      total_abundance = sum(nfm),
      .groups = "drop"
    ) %>%
    group_by(artal) %>%
    mutate(relative_abundance = total_abundance / sum(total_abundance))
  
  # Temporal changes plot
  temporal_plot <- create_temporal_plot(yearly_trait_summary, trait_column)
  
  # Contribution plot
  contribution_plot <- create_contribution_plot(sequential_simper, trait_column)
  
  # Year-to-year changes plot
  change_plot <- create_change_plot(yearly_trait_summary, trait_column)
  
  # Save outputs
  ggsave(paste0("output/", trait_column, "_temporal_changes.png"),
         temporal_plot,
         width = 12,
         height = 8,
         dpi = 300,
         bg = "white")
  
  # Return results as a list
  return(list(
    simper_results = sequential_simper,
    temporal_plot = temporal_plot,
    contribution_plot = contribution_plot,
    change_plot = change_plot,
    summary_data = yearly_trait_summary
  ))
}

#-------------------------------------------------------------------------------
# Execute Trait Analysis
#-------------------------------------------------------------------------------
# Run analysis for all traits
traits <- c("food_source", "motility", "habit", 
            "om_ca_he", "food_size_type", "feed_mode", "feeding_guild")

# To run analysis for all traits
trait_results <- lapply(traits, function(trait) {
  perform_trait_analysis(merged_data, trait)
})

# Access specific results
# For example, to view the temporal plot for food source:
print(results$temporal_plot)

# Save plots for a specific trait
ggsave("output/food_source_temporal.png", results$temporal_plot, 
       width = 12, height = 8, dpi = 300, bg="white")
ggsave("output/food_source_contribution.png", results$contribution_plot, 
       width = 12, height = 8, dpi = 300, bg="white")
ggsave("output/food_source_changes.png", results$change_plot, 
       width = 12, height = 8, dpi = 300, bg="white")

# Extract data for 1999 and 2017
years_compare <- c(1999, 2017)
comparison_data <- merged_data %>%
  filter(artal %in% years_compare) %>%  # 'artal' is the year column
  select(artal, taxon_name, nfm, food_size_type) %>%  # correct column names
  group_by(artal, taxon_name, food_size_type) %>%
  summarise(abundance = sum(nfm), .groups = "drop")

# Reshape to wide format for easier comparison
comparison_wide <- comparison_data %>%
  pivot_wider(
    names_from = artal,
    values_from = abundance,
    names_prefix = "year_"
  ) %>%
  mutate(
    year_1999 = replace_na(year_1999, 0),
    year_2017 = replace_na(year_2017, 0),
    abs_difference = abs(year_2017 - year_1999),
    relative_change = case_when(
      year_1999 == 0 & year_2017 == 0 ~ 0,
      year_1999 == 0 ~ Inf,
      TRUE ~ (year_2017 - year_1999) / year_1999
    )
  )

# 1. Species present in only one year
cat("\nSpecies present in only one year:\n")
unique_species <- comparison_wide %>%
  filter(year_1999 == 0 | year_2017 == 0) %>%
  arrange(food_size_type, desc(abs_difference))
print(unique_species, n = Inf)

# 2. Species with large proportional changes
cat("\nSpecies with major changes in abundance:\n")
major_changes <- comparison_wide %>%
  filter(year_1999 > 0 & year_2017 > 0) %>%
  filter(abs(relative_change) >= 0.5) %>%  # 50% or more change
  arrange(food_size_type, desc(abs_difference))
print(major_changes, n = Inf)

# Calculate proportions by feeding size type for each year
feeding_size_props <- comparison_data %>%
  group_by(artal, food_size_type) %>%
  summarise(
    total_abundance = sum(abundance),
    .groups = "drop"
  ) %>%
  group_by(artal) %>%
  mutate(proportion = total_abundance / sum(total_abundance))

# Print proportions
cat("\nFeeding size type proportions by year:\n")
print(feeding_size_props, n = Inf)

# Create visualization
size_plot <- ggplot(feeding_size_props, 
                    aes(x = factor(artal), y = proportion, fill = food_size_type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(
    title = "Changes in Feeding Size Type Proportions: 1999 vs 2017",
    x = "Year",
    y = "Proportion",
    fill = "Feeding Size Type"
  )

# Display the plot
print(size_plot)

# Save the plot
ggsave("output/feeding_size_changes.png", 
       size_plot, 
       width = 12, 
       height = 8,
       dpi = 300,
       bg = "white")

# Open sink to capture output
sink("output/feeding_size_analysis_summary.txt")

cat("=============================================================\n")
cat("FEEDING SIZE ANALYSIS SUMMARY\n")
cat("=============================================================\n\n")

# 1. Species present in only one year
cat("SPECIES PRESENT IN ONLY ONE YEAR:\n")
cat("-------------------------------------------------------------\n")
print(unique_species, n = Inf)
cat("\n\n")

# 2. Species with major changes
cat("SPECIES WITH MAJOR CHANGES IN ABUNDANCE (>50% change):\n")
cat("-------------------------------------------------------------\n")
print(major_changes, n = Inf)
cat("\n\n")

# 3. Feeding size type proportions
cat("FEEDING SIZE TYPE PROPORTIONS BY YEAR:\n")
cat("-------------------------------------------------------------\n")
print(feeding_size_props, n = Inf)
cat("\n\n")

# 4. Summary statistics
cat("SUMMARY STATISTICS:\n")
cat("-------------------------------------------------------------\n")
summary_stats <- comparison_data %>%
  group_by(food_size_type) %>%
  summarise(
    n_species = n_distinct(taxon_name),
    total_abundance_1999 = sum(abundance[artal == 1999]),
    total_abundance_2017 = sum(abundance[artal == 2017]),
    percent_change = ((total_abundance_2017 - total_abundance_1999) / total_abundance_1999) * 100
  )
print(summary_stats)
cat("\n\n")

# 5. Most dramatic changes
cat("MOST DRAMATIC CHANGES:\n")
cat("-------------------------------------------------------------\n")
cat("Top decreases (species present in 1999 but absent or rare in 2017):\n")
print(comparison_wide %>% 
        filter(year_1999 > 0) %>% 
        arrange(relative_change) %>% 
        head(10))
cat("\n")
cat("Top increases (species absent or rare in 1999 but present in 2017):\n")
print(comparison_wide %>% 
        filter(year_2017 > 0) %>% 
        arrange(desc(relative_change)) %>% 
        head(10))

# Close sink
sink()

cat("Analysis summary has been written to 'output/feeding_size_analysis_summary.txt'\n")
