#' Get Biological Trait Descriptions for a Taxon
#' 
#' @param taxon_name Character string of the taxon name to search for
#' @param filtered_traits Dataframe containing trait data
#' @param trait_descriptions Dataframe containing trait descriptions
#' @return Dataframe with trait information and descriptions, or NULL if taxon not found
#' 
#' 

filtered_traits <- readRDS("data/traits_data/filtered_traits.rds")

get_trait_description <- function(taxon_name, filtered_traits, trait_descriptions) {
  # Find taxon at any taxonomic level
  taxon_row <- filtered_traits %>%
    filter(
      Genus == taxon_name |
        Order == taxon_name |
        Family == taxon_name
    ) %>%
    select(-c(Genus, Order, Family, AphiaID))  # Exclude taxonomic and ID columns
  
  # Return NULL if taxon not found
  if (nrow(taxon_row) == 0) return(NULL)
  
  # Get non-zero traits
  trait_info <- data.frame(
    trait = names(taxon_row),
    value = unlist(taxon_row[1,])
  ) %>%
    filter(value > 0) %>%
    mutate(
      trait_category = sub("^([^_]+)_.*$", "\\1", trait),
      trait_subcategory = sub("^[^_]+_(.*)$", "\\1", trait)
    )
  
  # Match with descriptions and arrange results
  results <- trait_info %>%
    left_join(
      trait_descriptions,
      by = c(
        "trait_category" = "Trait",
        "trait_subcategory" = "Category"
      )
    ) %>%
    select(trait_category, trait_subcategory, value, Description) %>%
    arrange(trait_category, desc(value))
  
  return(results)
}


# Read trait descriptions
trait_descriptions <- read.csv("data/traits_data/Biological trait descriptions.txt", 
                               stringsAsFactors = FALSE) %>%
  mutate(
    Trait = case_when(
      Trait == "Maximum size (sr)" ~ "sr",
      Trait == "Morphology (m)" ~ "m",
      Trait == "Lifespan (l)" ~ "l",
      Trait == "Egg development location (ed)" ~ "ed",
      Trait == "Larva development location (ld)" ~ "ld",
      Trait == "Living habit (lh)" ~ "lh",
      Trait == "Sediment position (sp)" ~ "sp",
      Trait == "Feeding mode (f)" ~ "f",
      Trait == "Mobility (mob)" ~ "mob",
      Trait == "Bioturbation mode (b)" ~ "b",
      TRUE ~ Trait
    )
  )

# Get trait descriptions for a taxon
#spionidae_traits <- get_trait_description("Spionidae", filtered_traits, trait_descriptions)

get_top_species_traits <- function(year, dominant_species, filtered_traits, trait_descriptions, n = 5) {
  # Get top n species for the given year
  top_species <- dominant_species %>%
    filter(year_temp == year) %>%
    head(n) %>%
    pull(species)
  
  # Create a named list to store results
  trait_list <- list()
  
  # Get trait descriptions for each species
  for(sp in top_species) {
    # Extract genus from species name (first word)
    genus <- word(sp, 1)
    
    # Get trait description
    traits <- get_trait_description(genus, filtered_traits, trait_descriptions)
    
    # Add to list with full species name as identifier
    if(!is.null(traits)) {
      trait_list[[sp]] <- traits
    }
  }
  
  return(trait_list)
}

# Example usage for 1999
traits_1999 <- get_top_species_traits(1999, dominant_species, filtered_traits, trait_descriptions)

# Print results for each species
for(species in names(traits_1999)) {
  cat("\nTrait descriptions for", species, ":\n")
  print(traits_1999[[species]])
  cat("\n-------------------\n")
}

#' Write Trait Descriptions for Top Species to File
#' 
#' @param year Numeric year to analyze
#' @param output_file Character string specifying the output file path
write_top_species_traits <- function(year, dominant_species, filtered_traits, trait_descriptions, 
                                   output_file, n = 5) {
  # Open the sink file
  sink(output_file)
  
  # Get traits for top species
  traits_list <- get_top_species_traits(year, dominant_species, filtered_traits, trait_descriptions, n)
  
  # Write header
  cat(sprintf("\nTop %d Species Trait Descriptions for Year %d\n", n, year))
  cat("=================================================\n\n")
  
  # Print results for each species
  for(species in names(traits_list)) {
    cat(sprintf("\n%s\n", species))
    cat(paste(rep("-", nchar(species)), collapse = ""), "\n")
    print(traits_list[[species]])
    cat("\n\n")
  }
  
  # Close the sink
  sink()
}

# Example usage: Write traits for 1999 to a file
# write_top_species_traits(
#   year = 1999,
#   dominant_species = dominant_species,
#   filtered_traits = filtered_traits,
#   trait_descriptions = trait_descriptions,
#   output_file = "output/trait_descriptions/top_species_traits_1999.txt"
# )

#' Write Trait Descriptions for Multiple Years
#' 
#' @param years Vector of years to analyze
write_multiple_years_traits <- function(years, dominant_species, filtered_traits, trait_descriptions, n = 5) {
  for(year in years) {
    output_file <- sprintf("output/trait_descriptions/top_species_traits_%d.txt", year)
    write_top_species_traits(year, dominant_species, filtered_traits, trait_descriptions, 
                           output_file, n)
  }
}

# Example usage for multiple years
years_to_analyze <- c(1999,2013, 2014, 2015, 2016, 2017,)
write_multiple_years_traits(years_to_analyze, dominant_species, filtered_traits, trait_descriptions)
