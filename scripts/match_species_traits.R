# 01_SPECIES_MATCHING.R
library(dplyr)
library(readr)
library(tidyr)
library(stringr)

# Load necessary libraries
library(dplyr)
library(tidyr)

# Load your datasets
#df <- read_csv("path_to_df.csv")
df <- df_result
merged_data_allt <- read_csv("data/processed/merged_data_allt.csv")
trait_matrix <- read_csv("./data/traits_data/TraitMatrix.csv")




# 1. Find dominant species (>1% relative abundance)
dominant_species <- species_matrix %>% #from 03_nmds
  select(year_temp, where(is.numeric)) %>%  # Keep year_temp and numeric columns
  group_by(year_temp) %>%                   # Group by year first
  summarise(across(everything(), sum)) %>%    # Sum everything except year_temp
  pivot_longer(-year_temp,                  # Exclude year_temp from pivot
               names_to = "species", 
               values_to = "total_abundance") %>%
  group_by(year_temp) %>%                   # Group by year again for relative abundance
  mutate(relative_abundance = total_abundance/sum(total_abundance)) %>%
  arrange(year_temp, desc(total_abundance)) %>%
  filter(relative_abundance > 0.01)

str(dominant_species)
str(merged_data_allt)
str(trait_matrix)


# Replace "Cistenides granulata" with "Sabella granulata" in the species column
dominant_species <- dominant_species %>%
  mutate(species = ifelse(species == "Cistenides granulata", "Sabella granulata", species))

trait_matrix %>% 
  filter(tolower(Genus) %in% tolower(extract_first_word(unique(dominant_species$species)))) %>%
  DT::datatable()

# Create a new column for the full species name
not_in_traits_matrix <- dominant_species %>% 
  mutate(
    species_low = tolower(extract_first_word(species)),
    full_species_name = tolower(species)  # Store the full species name in lowercase
  ) %>% 
  filter(!species_low %in% tolower(trait_matrix$Genus))

# Function to filter, select taxonomic levels, and merge with abundance data
filter_taxonomic_levels <- function(data, levels, not_in_traits, abundance_data) {
  # Filter by taxonomic levels and get distinct matches
  filtered_data <- lapply(levels, function(level) {
    data %>%
      filter(tolower(!!sym(level)) %in% not_in_traits$species_low) %>%
      distinct(matched_taxon_name, .keep_all = TRUE)
  })
  
  combined_filtered <- bind_rows(filtered_data) %>%
    distinct(matched_taxon_name, .keep_all = TRUE)
  
  # Merge with abundance data
  merged_data <- combined_filtered %>%
    left_join(
      abundance_data %>% 
        select(year_temp, species, total_abundance, relative_abundance),
      by = c("Genus" = "species") # Adjust key columns as needed
    )
  
  return(merged_data)
}


# Define the taxonomic levels to filter by
taxonomic_levels <- c("Genus", "Family", "Order")


# Apply the function
merged_traits <- filter_taxonomic_levels(
  data = merged_data_allt,
  levels = taxonomic_levels,
  not_in_traits = not_in_traits_matrix,
  abundance_data = dominant_species
)




# Apply the function to filter by multiple taxonomic levels
all_filtered <- filter_taxonomic_levels(merged_data_allt, taxonomic_levels, not_in_traits_matrix)

# Extract species from combined filtered results
matched_species <- tolower(all_filtered$matched_taxon_name)

# Identify species in not_in_traits_matrix that are not in matched_species using full species names
unmatched_species <- not_in_traits_matrix %>%
  filter(!full_species_name %in% matched_species)

# Display unmatched species
print("Unmatched species:")
print(unmatched_species)

# Filter trait_matrix to include rows that meet at least one of the conditions
filtered_traits <- trait_matrix %>%
  filter(
    tolower(Genus) %in% tolower(extract_first_word(unique(dominant_species$species))) |
    tolower(Genus) %in% tolower(all_filtered$Genus) |
    tolower(Family) %in% tolower(all_filtered$Family) |
    tolower(Order) %in% tolower(all_filtered$Order)
  )

# Display the filtered traits
print(filtered_traits)

# Save the filtered_traits data frame as an RDS file
saveRDS(filtered_traits, file = "data/traits_data/filtered_traits.rds")







# nákvæm pörun:
find_exact_match <- function(x, y) {
  sapply(x, function(Flokkun) {
    matched <- y[Flokkun == y]
    if (length(matched) > 0) {
      return(matched[1])
    } else {
      return(NA_character_)
    }
  }, simplify = "character")
}


# Define a function to find the taxon match based on the exact match and taxonomic hierarchy
find_taxon_match <- function(x, kolgr_taxa, y) {
  # Initialize a vector to store the matched taxon names
  matched_taxon <- rep(NA_character_, length(x))
  
  # Define the taxonomic hierarchy order
  taxon_levels <- c("Genus", "Family", "Superfamily", "Order", "Superorder", "Subterclass", "Infraclass", "Subclass", "Class", "Phylum", "Kingdom")
  
  # Iterate through each element in x (df$Flokkun)
  for (i in 1:length(x)) {
    found_match <- FALSE
    
    # Try to find an exact match between the current element and the cleaned taxon names
    taxon_name <- find_exact_match(x[i], y$`matched_taxon_name`)
    if (!is.na(taxon_name)) {
      matched_taxon[i] <- taxon_name
      found_match <- TRUE
    }
    
    # If no exact match was found, search for a match in the taxonomic hierarchy
    if (!found_match) {
      for (level in taxon_levels) {
        taxon_name <- find_exact_match(kolgr_taxa[[level]][i], y$`matched_taxon_name`)
        if (!is.na(taxon_name)) {
          matched_taxon[i] <- taxon_name
          break
        }
      }
    }
  }
  
  # Return the vector of matched taxon names
  return(matched_taxon)
}


# Find the taxon match for each element in df$Flokkun
df$matched_taxon_name <- find_taxon_match(df$species, df, merged_data_allt)

# Identify rows with unmatched taxon names (empty strings)
unmatched_rows <-  df$matched_taxon_name == ""

# Define a function to find the closest string match between two sets of strings
find_closest_match <- function(x, y) {
  string_distances <- stringdist::stringdistmatrix(tolower(x), tolower(y), method = "jw")
  y[apply(string_distances, 1, which.min)]
}

# For the unmatched rows, find the closest string match between Flokkun and Cleaned Taxon name
df$matched_taxon_name[unmatched_rows] <- find_closest_match(df$species[unmatched_rows], merged_data_allt$`matched_taxon_name`)

# Merge the data frames based on the matched taxon names
merged_data_kolgr <- merge(df, merged_data_allt, by.x = "matched_taxon_name", by.y = "matched_taxon_name")


