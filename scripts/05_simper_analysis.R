# Load necessary libraries
library(dplyr)
library(tidyr)
library(vegan)
library(kableExtra)
library(BBI)

# Add this function after the libraries and before the data processing
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

# Check if df_result is loaded
if (!exists("df_result")) {
  stop("df_result is not defined. Please load or create this data frame before running the script.")
}

# Ensure necessary columns are present
required_columns <- c("station", "year", "species", "adjusted_density")
missing_columns <- setdiff(required_columns, names(df_result))
if (length(missing_columns) > 0) {
  stop("The following required columns are missing from df_result: ", paste(missing_columns, collapse = ", "))
}

# Summarize and reshape data
species_data <- df_result %>%
  group_by(station, year, species) %>%
  summarise(N = sum(adjusted_density), .groups = "drop") %>%
  pivot_wider(names_from = species, values_from = N, values_fill = 0)


# Separate environmental and species data
env_data <- data.frame(
  Station = factor(species_data$station, levels = c("C4", "E3", "E4", "B5", "B8", "A7")),
  Year = factor(species_data$year)
)

# Create species matrix (only numeric columns)
species_matrix <- species_data %>% 
  select(-station, -year) %>% 
  as.matrix()


# Use cleaned data from Jörundur

species_data <- henda_long_matched %>% # see 03_nmds_analysis.R
  group_by(station, year, species) %>%
  summarise(N = sum(adjusted_density), .groups = "drop") %>%
  pivot_wider(names_from = species, values_from = N, values_fill = 0)

# Separate environmental and species data
env_data <- data.frame(
  Station = factor(species_data$station, levels = c("C4", "E3", "E4", "B5", "B8", "A7")),
  Year = factor(species_data$year)
)

# Create species matrix (only numeric columns)
species_matrix <- species_data %>% 
  select(-station, -year) %>% 
  as.matrix()


# Perform SIMPER analysis for year-to-year comparisons
years <- sort(unique(as.character(env_data$Year)))
sequential_simper <- list()

for(i in 1:(length(years)-1)) {
  # Select data for consecutive years
  year_pair <- years[i:(i+1)]
  year_data <- species_matrix[env_data$Year %in% year_pair, ]
  year_groups <- env_data$Year[env_data$Year %in% year_pair]
  
  # Perform SIMPER analysis for this pair of years
  comparison_name <- paste(year_pair[1], year_pair[2], sep="_")
  sequential_simper[[comparison_name]] <- simper(year_data, group = year_groups, permutations = 999)
}

# Baseline comparison (if 1999 data exists)
if ("1999" %in% years) {
  baseline_data <- species_matrix[env_data$Year %in% c("1999", years), ]
  baseline_groups <- env_data$Year[env_data$Year %in% c("1999", years)]
  simper_baseline <- simper(baseline_data, group = baseline_groups, permutations = 999)
}

# Calculate Benthic Biotic Indices (BBI) for each year and station
calculate_biotic_indices <- function(data) {
  required_cols <- c("year", "station", "species", "adjusted_density")
  stopifnot(
    "Missing required columns" = all(required_cols %in% colnames(data)),
    "No data provided" = nrow(data) > 0,
    "Negative density values found" = all(data$adjusted_density >= 0),
    "Missing values found" = !any(is.na(data$adjusted_density))
  )
  
  target_stations <- c("C4", "A7", "B5", "B8", "E4", "E3")
  indices_list <- list()
  species_ambi_mapping <- data.frame(species = character(), AMBI = numeric())
  
  for (current_year in sort(unique(data$year))) {
    year_data <- data %>%
      filter(year == current_year, station %in% target_stations) %>%
      group_by(station, species) %>%
      summarise(density = sum(adjusted_density), .groups = 'drop') %>%
      pivot_wider(names_from = station, values_from = density, values_fill = 0)
    
    bbi_results <- BBI(year_data)
    
    # Extract and accumulate species-AMBI mappings
    current_species_groups <- data.frame(
      species = bbi_results$taxa[, 1],
      AMBI = bbi_results$table$AMBI
    )
    
    # Combine with existing mappings, keeping the first non-NA AMBI value for each species
    species_ambi_mapping <- bind_rows(species_ambi_mapping, current_species_groups) %>%
      group_by(species) %>%
      summarise(AMBI = first(na.omit(AMBI))) %>%
      ungroup()
    
    indices_list[[as.character(current_year)]] <- list(
      indices = as.data.frame(cbind(bbi_results$BBI, year = current_year)),
      classification = bbi_results$BBIclass,
      normalized = as.data.frame(nEQR(bbi_results$BBI)[1])
    )
  }
  
  # Add the species-AMBI mapping to the return list
  indices_list$species_groups <- species_ambi_mapping
  
  return(indices_list)
}

# Use the result dataframe that already has the correct density calculations
biotic_indices <- calculate_biotic_indices(df_result)

# Extract species groups from the accumulated mappings
species_groups <- biotic_indices$species_groups

# Function to summarize one comparison and add ecological group
summarize_one_comparison <- function(comparison_data, df_result, year1, year2, species_groups) {
  # Construct the dynamic key
  comparison_key <- paste(year1, year2, sep = "_")
  
  # Access cusum and other components dynamically
  species_names <- names(comparison_data[[comparison_key]]$cusum)
  cumulative_contributions <- unname(comparison_data[[comparison_key]]$cusum)
  individual_contributions <- c(cumulative_contributions[1], diff(cumulative_contributions))
  
  
  abundance_changes <- sapply(species_names, function(species) {
    abundance_year1 <- sum(df_result$adjusted_density[df_result$species == species & df_result$year == year1])
    abundance_year2 <- sum(df_result$adjusted_density[df_result$species == species & df_result$year == year2])
    change <- abundance_year2 - abundance_year1
    return(change)
  })
  
  ecological_groups <- sapply(species_names, function(species) {
    match_index <- match(species, species_groups$species)
    if (!is.na(match_index)) {
      return(species_groups$AMBI[match_index])
    } else {
      return(NA)
    }
  })
  
  df <- data.frame(
    Species = species_names,
    Cumulative_Percent = cumulative_contributions * 100,
    Individual_Percent = individual_contributions * 100,
    Abundance_Change = abundance_changes,
    Ecological_Group = ecological_groups
  )
  
  return(df)
}

# Process sequential comparisons
print("Sequential Year Comparisons:")
for(i in seq_along(sequential_simper)) {
  comparison_name <- names(sequential_simper)[i]
  years <- unlist(strsplit(comparison_name, "_"))
  year1 <- years[1]
  year2 <- years[2]
  
  cat(sprintf("\nComparison: %s\n", comparison_name))
  
  # Process this comparison
  summary_data <- summarize_one_comparison(sequential_simper[[i]], df_result, year1, year2, species_groups) %>%
    filter(Cumulative_Percent <= 70) %>%
    mutate(
      Species = format_species_names(Species)
    ) %>%
    mutate_if(is.numeric, ~ round(.x, 2))  # Round all numeric columns to 2 decimal places
  rownames(summary_data) <- NULL
  # Apply line breaks to column names
  colnames(summary_data) <- gsub("_", "<br>", colnames(summary_data))
  
  if(nrow(summary_data) > 0) {
    # Export the summary data to a CSV file with a descriptive filename
    #file_name <- sprintf("species_contributions_sequential_%s_%s.csv", year1, year2)
    #write.csv(summary_data, file = file_name, row.names = FALSE)
    
    # Print the table in console for reference
    print(
      kableExtra::kbl(
        summary_data,
        caption = sprintf("Species contributions for %s", comparison_name),
        digits = 2,
        format = "html",
        escape = F
      ) %>%
        kableExtra::kable_styling(
          full_width = FALSE,
          bootstrap_options = c("striped", "hover", "condensed", "responsive")
        ) %>%
        kableExtra::row_spec(0, bold = TRUE, color = "black", background = "#f0f0f0") %>%  # Header styling
        kableExtra::column_spec(1, bold = TRUE, width = "5em") %>%  # First column bold
        kableExtra::add_header_above(c(" " = 1, "Species Contributions" = ncol(summary_data) - 1))  # Custom header
    )
    
    cat(sprintf("Exported summary data for %s to %s\n", comparison_name, file_name))
  } else {
    cat("No species met the contribution threshold for this comparison.\n")
  }
}


# Assuming 'years' is your list of all years and 'sequential_simper' has already been populated

# Generate all possible pairs of years for arbitrary comparison
arbitrary_year_pairs <- combn(c("1999","2017"), 2, simplify = FALSE)

# Initialize a new list to store results for arbitrary year pairs
sequential_simper_arbitrary <- list()

# New loop for arbitrary year pairs
for (pair in arbitrary_year_pairs) {
  # Extract the two years from the current pair
  year1 <- pair[1]
  year2 <- pair[2]
  
  # Select data for these two specific years
  year_data <- species_matrix[env_data$Year %in% c(year1, year2), ]
  year_groups <- env_data$Year[env_data$Year %in% c(year1, year2)]
  
  # Generate a unique name for this comparison
  comparison_name <- paste(year1, year2, sep = "_")
  
  # Perform SIMPER analysis for this arbitrary pair of years
  sequential_simper_arbitrary[[comparison_name]] <- simper(year_data, group = year_groups, permutations = 999)
}

# Now, `sequential_simper_arbitrary` contains SIMPER results for all arbitrary pairs

print("Arbitrary Year Comparisons:")
for(i in seq_along(sequential_simper_arbitrary)) {
  # Retrieve the comparison name for the arbitrary pair
  comparison_name <- names(sequential_simper_arbitrary)[i]
  years <- unlist(strsplit(comparison_name, "_"))
  year1 <- years[1]
  year2 <- years[2]
  
  cat(sprintf("\nComparison: %s\n", comparison_name))
  
  # Process this arbitrary comparison
  summary_data <- summarize_one_comparison(sequential_simper_arbitrary[[i]], df_result, year1, year2, species_groups) %>%
    filter(Cumulative_Percent <= 70) %>%
    mutate(
      Species = format_species_names(Species)
    ) %>%
    mutate_if(is.numeric, ~ round(.x, 2))
  rownames(summary_data) <- NULL
  summary_data <- summary_data %>% mutate_if(is.numeric, ~ round(.x, 2))
  # Apply line breaks to column names
  colnames(summary_data) <- gsub("_", "<br>", colnames(summary_data))
  
  if(nrow(summary_data) > 0) {
    # Export the summary table to a CSV file
    #file_name <- sprintf("species_contributions_%s.csv", comparison_name)
    #write.csv(summary_data, file = file_name, row.names = FALSE)
    print(
      kableExtra::kbl(
        summary_data,
        caption = sprintf("Species contributions for %s", comparison_name),
        digits = 2,
        format = "html",
        escape = F
      ) %>%
        kableExtra::kable_styling(
          full_width = FALSE,
          bootstrap_options = c("striped", "hover", "condensed", "responsive")
        ) %>%
        kableExtra::row_spec(0, bold = TRUE, color = "black", background = "#f0f0f0") %>%  # Header styling
        kableExtra::column_spec(1, bold = TRUE, width = "5em") %>%  # First column bold
        kableExtra::add_header_above(c(" " = 1, "Species Contributions" = ncol(summary_data) - 1))  # Custom header
    )
  } else {
    cat("No species met the contribution threshold for this comparison.\n")
  }
}

# Create a function to prepare data for visualization
prepare_simper_viz_data <- function(sequential_simper, df_result, species_groups) {
  viz_data <- data.frame()
  
  for(i in seq_along(sequential_simper)) {
    comparison_name <- names(sequential_simper)[i]
    years <- unlist(strsplit(comparison_name, "_"))
    year1 <- years[1]
    year2 <- years[2]
    
    # Get summary data for this comparison
    summary_data <- summarize_one_comparison(sequential_simper[[i]], df_result, year1, year2, species_groups) %>%
      filter(Cumulative_Percent <= 70)
    
    # Calculate mean abundance for each species in each year
    year_data <- bind_rows(
      df_result %>%
        filter(year == year1, species %in% summary_data$Species) %>%
        group_by(species) %>%
        summarise(abundance = mean(adjusted_density), year = year1),
      df_result %>%
        filter(year == year2, species %in% summary_data$Species) %>%
        group_by(species) %>%
        summarise(abundance = mean(adjusted_density), year = year2)
    )
    
    # Add contribution information
    year_data <- year_data %>%
      left_join(summary_data %>% 
                select(Species, Individual_Percent, Ecological_Group),
                by = c("species" = "Species"))
    
    viz_data <- bind_rows(viz_data, year_data)
  }
  
  return(viz_data)
}

# Prepare visualization data
viz_data <- prepare_simper_viz_data(sequential_simper, df_result, species_groups)

# Create main SIMPER visualization
simper_plot <- ggplot(viz_data %>% 
                     filter(Individual_Percent >= 5), # Focus on major contributors
       aes(x = year, y = abundance, fill = factor(Ecological_Group))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~species, scales = "free_y") +
  scale_fill_brewer(palette = "Set2", name = "AMBI Group") +
  theme_bw() +
  labs(
    x = "Year",
    y = expression(paste("Abundance (individuals/m"^2, ")")),
    title = "Temporal changes in key species abundance (1999 and 2013-2017)",
    subtitle = "Species contributing >5% to between-year dissimilarity"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "italic"),
    legend.position = "bottom"
  )

# Create summary table
summary_table <- viz_data %>%
  group_by(species) %>%
  summarise(
    AMBI_Group = first(Ecological_Group),
    contribution = first(Individual_Percent),
    `2013` = mean(abundance[year == "2013"], na.rm = TRUE),
    `2014` = mean(abundance[year == "2014"], na.rm = TRUE),
    `2015` = mean(abundance[year == "2015"], na.rm = TRUE),
    `2016` = mean(abundance[year == "2016"], na.rm = TRUE),
    `2017` = mean(abundance[year == "2017"], na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    species = format_species_names(species)
  ) %>%
  filter(rowSums(!is.na(across(starts_with("20")))) > 0) %>%
  mutate(across(starts_with("20"), ~replace_na(., 0))) %>%
  mutate(
    across(starts_with("20"), ~round(., 0)),
    contribution = round(contribution, 1)
  ) %>%
  arrange(desc(contribution))

# Now create the formatted table
kableExtra::kbl(
  summary_table,
  caption = "Species contributions to community dissimilarity (2013-2017)",
  format = "html",
  col.names = c("Species", "AMBI Group", "% Contribution", "2013", "2014", "2015", "2016", "2017"),
  align = c("l", "c", "r", "r", "r", "r", "r", "r")
) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE) %>%
  row_spec(0, bold = TRUE) %>%
  column_spec(1, italic = TRUE) %>%
  add_header_above(c(" " = 3, "Abundance (ind./m²)" = 5))

# Export plot and table
ggsave("output/simper_temporal_changes.png", simper_plot, 
       width = 12, height = 8, dpi = 300)

write.csv(summary_table, "output/simper_summary_table.csv", row.names = FALSE)

# New SIMPER analysis with updated data structure
# Prepare data for SIMPER analysis
species_data_new <- merged_data_allt %>% 
  mutate(Density = Nfm, Year = Artal) %>%
  select(Year, Density,
         `Major Group`,
         `Taxon name`,
         `Food Source`,
         Motility,
         Habit,
         `Om/Ca/He`,
         `Food size/type`,
         FeedMode,
         `Feeding guild`
  ) %>%
  group_by(Year, `Taxon name`) %>%
  summarise(N = sum(Density), .groups = "drop") %>%
  pivot_wider(names_from = `Taxon name`, values_from = N, values_fill = 0)

# Inspect the data
print("First few rows of species_data_new:")
print(head(species_data_new))

# Look at data for 2016 and 2017 specifically
print("\nData for 2016 and 2017:")
print(species_data_new %>% filter(Year %in% c(2016, 2017)))

# Check for non-zero values in each year
print("\nNumber of non-zero species in each year:")
species_counts <- species_data_new %>%
  gather(species, abundance, -Year) %>%
  group_by(Year) %>%
  summarise(
    total_species = n(),
    non_zero_species = sum(abundance > 0),
    total_abundance = sum(abundance)
  )
print(species_counts)

# Separate environmental and species data
env_data_new <- data.frame(
  Year = factor(species_data_new$Year)
)

# Create species matrix (only numeric columns)
species_matrix_new <- species_data_new %>% 
  select(-Year) %>% 
  as.matrix()

# Perform SIMPER analysis for year-to-year comparisons
years_new <- sort(unique(as.character(env_data_new$Year)))
sequential_simper_new <- list()

# Add error checking and debugging for the year pairs
for(i in 1:(length(years_new)-1)) {
  # Select data for consecutive years
  year_pair <- years_new[i:(i+1)]
  
  # Debug print statements
  cat("\nProcessing year pair:", year_pair[1], "-", year_pair[2], "\n")
  
  # Get data for each year
  year1_data <- species_matrix_new[env_data_new$Year == year_pair[1], , drop = FALSE]
  year2_data <- species_matrix_new[env_data_new$Year == year_pair[2], , drop = FALSE]
  
  # Combine data for SIMPER
  year_data <- rbind(year1_data, year2_data)
  year_groups <- factor(c(year_pair[1], year_pair[2]))
  
  # Check dimensions and data before SIMPER analysis
  cat("Dimensions of year_data:", dim(year_data)[1], "x", dim(year_data)[2], "\n")
  cat("Year groups:", paste(year_groups, collapse=", "), "\n")
  
  # Find species with non-zero values in either year
  species_present <- colSums(year_data) > 0
  
  if(sum(species_present) > 0) {
    # Keep only species that are present
    year_data <- year_data[, species_present, drop = FALSE]
    cat("Number of species included:", sum(species_present), "\n")
    
    # Perform SIMPER analysis
    comparison_name <- paste(year_pair[1], year_pair[2], sep="_")
    tryCatch({
      sequential_simper_new[[comparison_name]] <- simper(year_data, group = year_groups, permutations = 999)
      cat("Successfully completed SIMPER analysis for", comparison_name, "\n")
      
      # Print top contributing species
      top_species <- summary(sequential_simper_new[[comparison_name]])
      cat("\nTop contributing species:\n")
      print(head(top_species))
      
    }, error = function(e) {
      cat("Error in SIMPER analysis for", comparison_name, ":", conditionMessage(e), "\n")
    })
  } else {
    cat("No species with non-zero values found between years", year_pair[1], "and", year_pair[2], "\n")
  }
}

# Modified visualization code to work with the original data
viz_data_new <- merged_data_allt %>%
  filter(Artal %in% years_new) %>%
  group_by(Artal, `Taxon name`, `Major Group`, `Feeding guild`) %>%
  summarise(Density = sum(Nfm), .groups = 'drop') %>%
  rename(Year = Artal)

# Get the species that contributed to SIMPER analysis
contributing_species <- unique(unlist(lapply(sequential_simper_new, function(x) {
  if(!is.null(x)) {
    # Get species that contribute to 70% cumulative dissimilarity
    summary_data <- summary(x)
    species_contrib <- rownames(summary_data)[summary_data$cumsum <= 70]
    return(species_contrib)
  }
  return(NULL)
})))

# Filter visualization data for contributing species
if(length(contributing_species) > 0) {
  viz_data_new <- viz_data_new %>%
    filter(`Taxon name` %in% contributing_species)
}

# Create enhanced SIMPER visualization
if(nrow(viz_data_new) > 0) {
  simper_plot_new <- ggplot(viz_data_new %>% 
                             group_by(`Taxon name`) %>%
                             filter(sum(Density) > 0), 
                           aes(x = Year, y = Density, fill = `Feeding guild`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~`Taxon name`, scales = "free_y") +
    scale_fill_brewer(palette = "Set3", name = "Feeding Guild") +
    theme_bw() +
    labs(
      x = "Year",
      y = expression(paste("Density (individuals/m"^2, ")")),
      title = "Temporal changes in species density with feeding guilds",
      subtitle = "Species contributing to dissimilarity between years"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "italic"),
      legend.position = "bottom"
    )
  
  # Save the new visualization
  ggsave("output/simper_temporal_changes_ecological.png", simper_plot_new, 
         width = 14, height = 10, dpi = 300)
} else {
  warning("No data available for visualization after SIMPER analysis")
}
