trait_matrix <- read_csv("./data/traits_data/TraitMatrix.csv")
df <- read_rds("data/cleaned_data.rds") %>%
  filter(station %in% c("C4", "A7", "B5", "B8", "E4", "E3")) %>% 
  # First, apply the subdivision to the count
  mutate(adjusted_count = count * subdivision) %>%
  
  # Calculate the correct number of samples and total area for each station-year combination
  group_by(station, year) %>%
  mutate(
    n_samples = n_distinct(sample_id),
    grab_area = case_when(
      year == 1999 ~ 0.0225,
      TRUE ~ 0.04
    ),
    total_area = grab_area * n_samples
  ) %>%
  ungroup() %>%
  
  # Now group by station, year, and species to get total counts
  group_by(station, year, species) %>%
  summarise(
    original_total_count = sum(count),
    adjusted_total_count = sum(adjusted_count),
    n_samples = first(n_samples),  # Use the previously calculated n_samples
    total_area = first(total_area),  # Use the previously calculated total_area
    .groups = 'keep'
  ) %>%
  
  # Calculate densities
  mutate(
    original_density = original_total_count / total_area,
    adjusted_density = adjusted_total_count / total_area,
    density_difference = adjusted_density - original_density,
    density_ratio = adjusted_density / original_density
  ) %>%
  
  # Ungroup
  ungroup()

# Get vectors of valid taxonomic levels from trait matrix
valid_genera <- trait_matrix$Genus
valid_families <- trait_matrix$Family
valid_orders <- trait_matrix$Order
valid_classes <- trait_matrix$Class
valid_phyla <- trait_matrix$Phylum

# Try matching at different taxonomic levels
df <- df %>%
  mutate(
    Genus = word(species, 1),
    Genus = case_when(
      str_count(species, "\\s+") >= 1 ~ word(species, 1),
      TRUE ~ species
    ),
    matched_level = case_when(
      Genus %in% valid_genera ~ "Genus",
      Genus %in% valid_families ~ "Family",
      Genus %in% valid_orders ~ "Order",
      Genus %in% valid_classes ~ "Class",
      Genus %in% valid_phyla ~ "Phylum",
      TRUE ~ "Unmatched"
    )
  )

# Print summary of matches at each level
match_summary <- df %>%
  group_by(matched_level) %>%
  summarise(count = n_distinct(Genus)) %>%
  arrange(matched_level)

print("Matches by taxonomic level:")
print(match_summary)

# Get unmatched taxa
unmatched <- df %>%
  filter(matched_level == "Unmatched") %>%
  distinct(Genus) %>%
  pull(Genus)

print("\nUnmatched taxa:")
print(unmatched)





df <- df %>% 
  filter(!species %in% tolower(c("harpacticoida", "Campanulariidae sp"))) %>% 
  mutate(
    species = case_when(
      # Handle previously unmatched taxa
      tolower(species) == "campanulariidae sp" ~ NA_character_,  # Remove as it's too general
      tolower(species) == "tunicata" ~ NA_character_,  # Remove as it's too general
      tolower(species) %in% c("nephthys", "nepthys") ~ "Nephtys",  # Fix spelling
      tolower(species) %in% c("pleurogoniuminosissimum", "plerugoniuminosissum") ~ "Pleurogonium spinosissimum",  # Fix spelling
      tolower(species) == "leda" ~ "Nuculana pernula",  # Update to current taxonomy
      
      # General taxonomic standardization (no year-specific)
      str_detect(tolower(species), "^sipuncul") ~ NA_character_,
      tolower(species) == "ostracoda" ~ "Ostracoda",
      tolower(species) == "terebellides stroemi" ~ "Terebellides stroemii",
      tolower(species) == "bivalvia" ~ NA_character_,
      tolower(species) == "pectinaria koreni" ~ "Lagis koreni",
      tolower(species) == "sphaerosyllis" ~ "Sphaerosyllis erinaceus",
      tolower(species) %in% c("maldane sarsi", "maldanidae", "praxillella praetermissa") ~ "Praxillella",
      #tolower(species) == "oligochaeta" ~ NA_character_,
      #tolower(species) == "Tubificidae" ~ NA_character_,
      
      # 1999 specific matches
      year == "1999" & tolower(species) %in% c("ampharetinae", "ampharete acutifrons") ~ "Ampharetidae",
      year == "1999" & tolower(species) == "bivalvia" ~ NA_character_,
      year == "1999" & tolower(species) == "leucon acutirostris" ~ "Cumacea",
      year == "1999" & tolower(species) == "eudorella emarginata" ~ "Cumacea",
      year == "1999" & tolower(species) == "musculus discors" ~ "Musculus",
      year == "1999" & tolower(species) == "aricidea (acmira) cerrutii" ~ "Aricidea suecica",
      year == "1999" & tolower(species) == "mediomastus filiformis" ~ "Heteromastus filiformis",
      
      # 2013 specific matches
      year == "2013" & tolower(species) %in% c("harmothoe", "polynoidae") ~ "Harmothoe extenuata",
      year == "2013" & tolower(species) == "spionidae" ~ "Spio",
      
      # 2014 specific matches
      year == "2014" & tolower(species) == "ampharete acutifrons" ~ "Ampharetidae",
      year == "2014" & tolower(species) == "spionidae" ~ "Spio",
      year == "2014" & tolower(species) == "terebellidae" ~ "Terebellides stroemii",
      year == "2014" & tolower(species) == "caprellidae" ~ "Caprella septemtrionalis",
      year == "2014" & tolower(species) == "praxillella" ~ "Praxillella praetermissa",
      year == "2014" & tolower(species) == "syllidae" ~ "Syllis cornuta",
      year == "2014" & tolower(species) == "mya" ~ "Mya arenaria",
      year == "2014" & tolower(species) == "mytilus edulis" ~ "Mytilidae",
      year == "2014" & tolower(species) == "nephtyidae" ~ "Nephthys",
      year == "2014" & tolower(species) == "amphipoda" ~ "Protomedeia fasciata",
      year == "2014" & tolower(species) == "tubificidae" ~ "Tubificoides kozloffi",
      year == "2014" & tolower(species) == "harmothoe" ~ "Harmothoe extenuata",
      
      # 2015 specific matches
      year == "2015" & tolower(species) == "spionidae" ~ "Spio",
      year == "2015" & tolower(species) == "terebellidae" ~ "Terebellides stroemii",
      year == "2015" & tolower(species) == "crenella" ~ "Crenella decussata",
      year == "2015" & tolower(species) == "syllidae" ~ "Syllis cornuta",
      year == "2015" & tolower(species) == "mya" ~ "Mya arenaria",
      year == "2015" & tolower(species) == "mytilidae" ~ "Mytilus edulis",
      year == "2015" & tolower(species) == "ampharetidae" ~ NA_character_,
      year == "2015" & tolower(species) == "pholoe" ~ "Pholoe minuta",
      year == "2015" & tolower(species) == "harmothoe" ~ "Harmothoe extenuata",
      year == "2015" & tolower(species) == "amphipoda" ~ "Protomedeia fasciata",
      
      # 2016 specific matches
      year == "2016" & tolower(species) == "leaena ebranchiata" ~ "Terebellides stroemii",
      year == "2016" & tolower(species) == "pseudopolydora antennata" ~ "Polydora",
      year == "2016" & tolower(species) == "cistenides hyperborea" ~ "Cistenides granulata",
      year == "2016" & tolower(species) == "paraonidae" ~ "Aricidea suecica",
      #year == "2016" & tolower(species) == "dorvilleidae" ~ "Capitella capitata",
      year == "2016" & tolower(species) == "orbiniidae" ~ "Scoloplos armiger",
      year == "2016" & tolower(species) == "polychaeta" ~ NA_character_,
      year == "2016" & tolower(species) == "polynoidae" ~ "Harmothoe imbricata",
      year == "2016" & tolower(species) %in% c("spio filicornis", "spionidae") ~ "Spionidae",
      year == "2016" & tolower(species) %in% c("syllis", "syllidae") ~ "Syllis cornuta",
      year == "2016" & tolower(species) %in% c("terebellidae", "terebelliformia") ~ "Terebellides stroemii",
      year == "2016" & tolower(species) == "amphipoda" ~ NA_character_,
      year == "2016" & tolower(species) == "bivalvia" ~ NA_character_,
      year == "2016" & tolower(species) == "cardium" ~ NA_character_,
      year == "2016" & tolower(species) == "aricidea" ~ "Amphitrite cirrata",
      year == "2016" & tolower(species) == "capitellidae" ~ "Capitella capitata",
      year == "2016" & tolower(species) == "cirratulidae" ~ "Cirratulus cirratus",
      year == "2016" & tolower(species) == "cossuridae" ~ "Cossura longocirrata",
      year == "2016" & tolower(species) == "pectinariidae" ~ "Pectinaria koreni",
      year == "2016" & tolower(species) == "phyllodocida" ~ "Phyllodoce maculata",
      year == "2016" & tolower(species) == "lumbrineridae" ~ "Lumbrineris",
      year == "2016" & tolower(species) == "pholoe" ~ "Pholoe minuta",
      year == "2016" & tolower(species) == "balanus balanus" ~ "Balanus",
      year == "2016" & tolower(species) %in% c("priapulidae", "priapulidae") ~ "Priapulus caudatus",
      
      # 2017 specific matches
      #year == "2017" & tolower(species) == "oligochaeta" ~ NA_character_,
      #year == "2017" & tolower(species) == "tubificidae" ~ NA_character_,
      year == "2017" & tolower(species) == "mediomastus filiformis" ~ "Heteromastus filiformis",
      year == "2017" & tolower(species) == "spionidae" ~ "Spio",
      year == "2017" & tolower(species) == "terebellidae" ~ "Terebellides stroemii",
      year == "2017" & tolower(species) == "caprellidae" ~ "Caprella septentrionalis",
      year == "2017" & tolower(species) == "opisthobranchia" ~ "Retusa pertenuis",
      year == "2017" & tolower(species) == "mya" ~ "Mya arenaria",
      year == "2017" & tolower(species) == "maldanidae" ~ "Praxillella praetermissa",
      year == "2017" & tolower(species) %in% c("ampharete", "ampharete acutifrons") ~ "Ampharetidae",
      year == "2017" & tolower(species) == "pholoe" ~ "Pholoe minuta",
      year == "2017" & tolower(species) == "harmothoe" ~ "Harmothoe extenuata",
      
      # Default case
      TRUE ~ species
    )
  ) %>% 
  drop_na(species)    


df <- df %>%
  mutate(
    Genus = word(species, 1),  # Extract first word
    # Only keep first word if species name has multiple words
    Genus = case_when(
      str_count(species, "\\s+") >= 1 ~ word(species, 1),
      TRUE ~ species
    ),
    is_valid_genus = Genus %in% valid_genera
  )

# Check results
n_matches <- sum(df$is_valid_genus)
print(paste("Number of valid genus matches:", n_matches))

# View unmatched genera if needed
unmatched <- unique(df$Genus[!df$is_valid_genus])
print("Unmatched genera:")
print(unmatched)

# Function to search for a taxon across all taxonomic columns
find_taxon <- function(taxon, data) {
  # List of taxonomic columns to search
  tax_cols <- c("matched_taxon_name", "Kingdom", "Phylum", "Class", "Subclass", 
                "Order", "Family", "Genus", "Species", "Subfamily", "Tribe")
  
  # Search each column
  matches <- lapply(tax_cols, function(col) {
    matches <- grep(taxon, data[[col]], value = TRUE, ignore.case = TRUE)
    if(length(matches) > 0) {
      return(data.frame(search_term = taxon, 
                       found_in = col, 
                       matches = paste(unique(matches), collapse = "; ")))
    }
    return(NULL)
  })
  
  # Combine results
  do.call(rbind, matches)
}

# Search for each unmatched taxon
unmatched_results <- lapply(unmatched, find_taxon, merged_data)
results_df <- do.call(rbind, unmatched_results)

# Print results in a readable format
if(!is.null(results_df)) {
  print("Found matches:")
  for(i in 1:nrow(results_df)) {
    cat("\n", results_df$search_term[i], "found in", results_df$found_in[i], "as:", results_df$matches[i])
  }
} else {
  print("No matches found")
}

# Print taxa with no matches at all
no_matches <- unmatched[!unmatched %in% results_df$search_term]
if(length(no_matches) > 0) {
  cat("\n\nTaxa with no matches in any column:\n")
  print(no_matches)
}

# First create a taxonomic matching dataframe from df
df_taxa <- df %>%
  mutate(
    match_taxon = case_when(
      # If it's found as a genus, use that
      Genus %in% trait_matrix$Genus ~ Genus,
      # If it's found as a family, use that
      Genus %in% trait_matrix$Family ~ Genus,
      # If it's found as an order, use that
      Genus %in% trait_matrix$Order ~ Genus,
      # If it's found as a class, use that
      Genus %in% trait_matrix$Class ~ Genus,
      # If it's found as a phylum, use that
      Genus %in% trait_matrix$Phylum ~ Genus,
      TRUE ~ NA_character_
    ),
    taxon_rank = case_when(
      Genus %in% trait_matrix$Genus ~ "Genus",
      Genus %in% trait_matrix$Family ~ "Family",
      Genus %in% trait_matrix$Order ~ "Order",
      Genus %in% trait_matrix$Class ~ "Class",
      Genus %in% trait_matrix$Phylum ~ "Phylum",
      TRUE ~ NA_character_
    )
  )

# First, let's clean up trait_matrix by removing rows with NA genus
trait_matrix_clean <- trait_matrix %>%
  filter(!is.na(Genus))

# Calculate average traits for each taxonomic level
family_traits <- trait_matrix_clean %>%
  group_by(Family) %>%
  summarise(across(starts_with(c("f_", "sr_", "m_", "l_", "ed_", "ld_", "lh_", "sp_", "mob_", "b_")), 
                   ~mean(., na.rm = TRUE))) %>%
  mutate(taxon_rank = "Family")

order_traits <- trait_matrix_clean %>%
  group_by(Order) %>%
  summarise(across(starts_with(c("f_", "sr_", "m_", "l_", "ed_", "ld_", "lh_", "sp_", "mob_", "b_")), 
                   ~mean(., na.rm = TRUE))) %>%
  mutate(taxon_rank = "Order")

class_traits <- trait_matrix_clean %>%
  group_by(Class) %>%
  summarise(across(starts_with(c("f_", "sr_", "m_", "l_", "ed_", "ld_", "lh_", "sp_", "mob_", "b_")), 
                   ~mean(., na.rm = TRUE))) %>%
  mutate(taxon_rank = "Class")

# Join with higher taxonomic levels
df_with_traits_all <- df_taxa %>%
  # First try genus match
  left_join(
    trait_matrix_clean %>% select(Genus, starts_with(c("f_", "sr_", "m_", "l_", "ed_", "ld_", "lh_", "sp_", "mob_", "b_"))),
    by = c("match_taxon" = "Genus")
  ) %>%
  # For those without genus traits, try family
  left_join(
    family_traits %>% select(Family, starts_with(c("f_", "sr_", "m_", "l_", "ed_", "ld_", "lh_", "sp_", "mob_", "b_"))),
    by = c("match_taxon" = "Family")
  ) %>%
  # Then order
  left_join(
    order_traits %>% select(Order, starts_with(c("f_", "sr_", "m_", "l_", "ed_", "ld_", "lh_", "sp_", "mob_", "b_"))),
    by = c("match_taxon" = "Order")
  ) %>%
  # Then class
  left_join(
    class_traits %>% select(Class, starts_with(c("f_", "sr_", "m_", "l_", "ed_", "ld_", "lh_", "sp_", "mob_", "b_"))),
    by = c("match_taxon" = "Class")
  ) %>%
  # Combine traits using coalesce for each trait type
  mutate(
    # Feeding traits
    f_Suspension = coalesce(f_Suspension.x, f_Suspension.y),
    f_Surface_deposit = coalesce(f_Surface_deposit.x, f_Surface_deposit.y),
    f_Subsurface_deposit = coalesce(f_Subsurface_deposit.x, f_Subsurface_deposit.y),
    f_Scavenger = coalesce(f_Scavenger.x, f_Scavenger.y),
    f_Predator = coalesce(f_Predator.x, f_Predator.y),
    f_Parasite = coalesce(f_Parasite.x, f_Parasite.y),
    
    # Size traits
    sr_Less_than_10 = coalesce(sr_Less_than_10.x, sr_Less_than_10.y),
    sr_11_to_20 = coalesce(sr_11_to_20.x, sr_11_to_20.y),
    sr_21_to_100 = coalesce(sr_21_to_100.x, sr_21_to_100.y),
    sr_101_to_200 = coalesce(sr_101_to_200.x, sr_101_to_200.y),
    sr_201_to_500 = coalesce(sr_201_to_500.x, sr_201_to_500.y),
    sr_More_than_500 = coalesce(sr_More_than_500.x, sr_More_than_500.y),
    
    # Morphology traits
    m_Soft = coalesce(m_Soft.x, m_Soft.y),
    m_Tunic = coalesce(m_Tunic.x, m_Tunic.y),
    m_Exoskeleton = coalesce(m_Exoskeleton.x, m_Exoskeleton.y),
    m_Crustose = coalesce(m_Crustose.x, m_Crustose.y),
    m_Cushion = coalesce(m_Cushion.x, m_Cushion.y),
    m_Stalked = coalesce(m_Stalked.x, m_Stalked.y),
    
    # Lifespan traits
    l_Less_than_1 = coalesce(l_Less_than_1.x, l_Less_than_1.y),
    l_1_to_3 = coalesce(l_1_to_3.x, l_1_to_3.y),
    l_3_to_10 = coalesce(l_3_to_10.x, l_3_to_10.y),
    l_More_than_10 = coalesce(l_More_than_10.x, l_More_than_10.y),
    
    # Egg development traits
    ed_Asexual = coalesce(ed_Asexual.x, ed_Asexual.y),
    ed_Sexual_pelagic = coalesce(ed_Sexual_pelagic.x, ed_Sexual_pelagic.y),
    ed_Sexual_benthic = coalesce(ed_Sexual_benthic.x, ed_Sexual_benthic.y),
    ed_Sexual_brooded = coalesce(ed_Sexual_brooded.x, ed_Sexual_brooded.y),
    
    # Larval development traits
    ld_Pelagic_planktotrophic = coalesce(ld_Pelagic_planktotrophic.x, ld_Pelagic_planktotrophic.y),
    ld_Pelagic_lecithotrophic = coalesce(ld_Pelagic_lecithotrophic.x, ld_Pelagic_lecithotrophic.y),
    ld_Benthic_direct = coalesce(ld_Benthic_direct.x, ld_Benthic_direct.y),
    
    # Living habit traits
    lh_Tube_dwelling = coalesce(lh_Tube_dwelling.x, lh_Tube_dwelling.y),
    lh_Burrow_dwelling = coalesce(lh_Burrow_dwelling.x, lh_Burrow_dwelling.y),
    lh_Free_living = coalesce(lh_Free_living.x, lh_Free_living.y),
    lh_Crevice_hole_under_stones = coalesce(lh_Crevice_hole_under_stones.x, lh_Crevice_hole_under_stones.y),
    lh_Epi_endo_biotic = coalesce(lh_Epi_endo_biotic.x, lh_Epi_endo_biotic.y),
    lh_Attached_to_substratum = coalesce(lh_Attached_to_substratum.x, lh_Attached_to_substratum.y),
    
    # Sediment position traits
    sp_Surface = coalesce(sp_Surface.x, sp_Surface.y),
    sp_Shallow_infauna_0_to_5cm = coalesce(sp_Shallow_infauna_0_to_5cm.x, sp_Shallow_infauna_0_to_5cm.y),
    sp_Mid_depth_infauna_5_to_10cm = coalesce(sp_Mid_depth_infauna_5_to_10cm.x, sp_Mid_depth_infauna_5_to_10cm.y),
    sp_Deep_infauna_more_than_10cm = coalesce(sp_Deep_infauna_more_than_10cm.x, sp_Deep_infauna_more_than_10cm.y),
    
    # Mobility traits
    mob_Sessile = coalesce(mob_Sessile.x, mob_Sessile.y),
    mob_Swim = coalesce(mob_Swim.x, mob_Swim.y),
    mob_Crawl_creep_climb = coalesce(mob_Crawl_creep_climb.x, mob_Crawl_creep_climb.y),
    mob_Burrower = coalesce(mob_Burrower.x, mob_Burrower.y),
    
    # Bioturbation traits
    b_Diffusive_mixing = coalesce(b_Diffusive_mixing.x, b_Diffusive_mixing.y),
    b_Surface_deposition = coalesce(b_Surface_deposition.x, b_Surface_deposition.y),
    b_Upward_conveyor = coalesce(b_Upward_conveyor.x, b_Upward_conveyor.y),
    b_Downward_conveyer = coalesce(b_Downward_conveyer.x, b_Downward_conveyer.y),
    b_None = coalesce(b_None.x, b_None.y)
  ) %>%
  # Remove the duplicate columns
  select(-ends_with(".x"), -ends_with(".y"))

# Check updated trait coverage
trait_coverage_new <- df_with_traits_all %>%
  group_by(taxon_rank) %>%
  summarise(
    n_taxa = n_distinct(match_taxon),
    n_with_traits = n_distinct(match_taxon[!is.na(f_Suspension)])
  )

print("Updated trait coverage by taxonomic rank:")
print(trait_coverage_new)

# Save the final dataset
saveRDS(df_with_traits_all, "data/df_with_traits_all.rds")

# Quick look at feeding guild distribution
feeding_summary <- df_with_traits_all %>%
  group_by(match_taxon) %>%
  slice(1) %>%  # Take first occurrence of each taxon
  summarise(
    suspension = mean(f_Suspension, na.rm = TRUE),
    deposit = mean(f_Surface_deposit, na.rm = TRUE),
    subsurface = mean(f_Subsurface_deposit, na.rm = TRUE),
    predator = mean(f_Predator, na.rm = TRUE)
  ) %>%
  filter(if_any(everything(), ~!is.na(.)))  # Remove rows where all feeding values are NA

print("\nFeeding guild summary (non-NA values only):")
print(summary(feeding_summary))

# Check columns in final dataset
print("Columns in final dataset:")
print(names(df_with_traits_all))

# Check columns in trait matrix for comparison
print("\nColumns in trait matrix:")
print(names(trait_matrix))

# Check for missing trait columns
trait_cols <- names(trait_matrix)[str_detect(names(trait_matrix), "^(f|sr|m|l|ed|ld|lh|sp|mob|b)_")]
final_cols <- names(df_with_traits_all)[str_detect(names(df_with_traits_all), "^(f|sr|m|l|ed|ld|lh|sp|mob|b)_")]

missing_cols <- setdiff(trait_cols, final_cols)
print("\nMissing trait columns in final dataset:")
print(missing_cols)

# Calculate feeding guild proportions with mutually exclusive categories
feeding_props_complete <- df_with_traits_all %>%
  filter(tolower(species) != "oligochaeta") %>% 
  # First replace NAs with 0 for feeding traits
  mutate(
    f_Surface_deposit = replace_na(f_Surface_deposit, 0),
    f_Subsurface_deposit = replace_na(f_Subsurface_deposit, 0),
    f_Suspension = replace_na(f_Suspension, 0),
    f_Predator = replace_na(f_Predator, 0)
  ) %>%
  group_by(year) %>%
  summarise(
    total_density = sum(adjusted_density),
    # Strictly mutually exclusive categories
    deposit = sum(adjusted_density * ((f_Surface_deposit > 0 | f_Subsurface_deposit > 0) & f_Suspension == 0 & f_Predator == 0)) / total_density,
    deposit_suspension = sum(adjusted_density * ((f_Surface_deposit > 0 | f_Subsurface_deposit > 0) & f_Suspension > 0 & f_Predator == 0)) / total_density,
    suspension = sum(adjusted_density * (f_Surface_deposit == 0 & f_Subsurface_deposit == 0 & f_Suspension > 0 & f_Predator == 0)) / total_density,
    predator = sum(adjusted_density * (f_Predator > 0 & f_Surface_deposit == 0 & f_Subsurface_deposit == 0 & f_Suspension == 0)) / total_density,
    predator_deposit = sum(adjusted_density * (f_Predator > 0 & (f_Surface_deposit > 0 | f_Subsurface_deposit > 0) & f_Suspension == 0)) / total_density,
    predator_suspension = sum(adjusted_density * (f_Predator > 0 & f_Suspension > 0)) / total_density,
    other = sum(adjusted_density * (f_Surface_deposit == 0 & f_Subsurface_deposit == 0 & f_Suspension == 0 & f_Predator == 0)) / total_density
  ) %>%
  mutate(
    # Verify total equals 100%
    total = deposit + deposit_suspension + suspension + predator + predator_deposit + predator_suspension + other
  )

# Print the complete proportions and totals
print("\nFeeding guild proportions (complete, mutually exclusive) by year:")
print(feeding_props_complete %>%
        mutate(across(where(is.numeric), ~round(., 3))))

# Prepare data for plotting
feeding_props_plot <- feeding_props_complete %>%
  select(-total, -total_density) %>%
  pivot_longer(-year, names_to = "feeding_guild", values_to = "proportion") %>%
  # Reorder factors for plotting
  mutate(feeding_guild = factor(feeding_guild, 
                              levels = c("deposit", "deposit_suspension", "suspension", 
                                       "predator", "predator_deposit", "predator_suspension", "other")))

# Create plot with complete data
feeding_guild_plot_complete <- ggplot(feeding_props_plot, 
                                    aes(x = factor(year), 
                                        y = proportion, 
                                        fill = feeding_guild)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Set2", 
                   labels = c("Deposit feeders", 
                            "Deposit/suspension feeders",
                            "Suspension feeders",
                            "Predators",
                            "Predator/deposit feeders",
                            "Predator/suspension feeders",
                            "Other")) +
  labs(title = "Temporal Changes in Feeding Guild Composition",
       x = "Year",
       y = "Relative Abundance",
       fill = "Feeding Guild") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_y_continuous(labels = scales::percent_format())

# Save the complete plot
ggsave("output/feeding_guild_temporal_changes.png",
       feeding_guild_plot_complete,
       width = 12,
       height = 8,
       dpi = 300,
       bg = "white")

# Count taxa in each feeding guild by year
taxa_counts <- df_with_traits_all %>%
  # First replace NAs with 0 for feeding traits
  mutate(
    f_Surface_deposit = replace_na(f_Surface_deposit, 0),
    f_Subsurface_deposit = replace_na(f_Subsurface_deposit, 0),
    f_Suspension = replace_na(f_Suspension, 0),
    f_Predator = replace_na(f_Predator, 0)
  ) %>%
  mutate(
    feeding_guild = case_when(
      (f_Surface_deposit > 0 | f_Subsurface_deposit > 0) & f_Suspension == 0 & f_Predator == 0 ~ "deposit",
      (f_Surface_deposit > 0 | f_Subsurface_deposit > 0) & f_Suspension > 0 & f_Predator == 0 ~ "deposit_suspension",
      f_Surface_deposit == 0 & f_Subsurface_deposit == 0 & f_Suspension > 0 & f_Predator == 0 ~ "suspension",
      f_Predator > 0 & f_Surface_deposit == 0 & f_Subsurface_deposit == 0 & f_Suspension == 0 ~ "predator",
      f_Predator > 0 & (f_Surface_deposit > 0 | f_Subsurface_deposit > 0) & f_Suspension == 0 ~ "predator_deposit",
      f_Predator > 0 & f_Suspension > 0 ~ "predator_suspension",
      TRUE ~ "other"
    )
  ) %>%
  group_by(year) %>%
  summarise(
    total_taxa = n_distinct(match_taxon),
    deposit_taxa = n_distinct(match_taxon[feeding_guild == "deposit"]),
    deposit_suspension_taxa = n_distinct(match_taxon[feeding_guild == "deposit_suspension"]),
    suspension_taxa = n_distinct(match_taxon[feeding_guild == "suspension"]),
    predator_taxa = n_distinct(match_taxon[feeding_guild == "predator"]),
    predator_deposit_taxa = n_distinct(match_taxon[feeding_guild == "predator_deposit"]),
    predator_suspension_taxa = n_distinct(match_taxon[feeding_guild == "predator_suspension"]),
    other_taxa = n_distinct(match_taxon[feeding_guild == "other"])
  )

print("\nNumber of taxa in each feeding guild by year:")
print(taxa_counts)

# Verify feeding guild claims
print("\nVerifying feeding guild claims:")

# 1. Check 1999 baseline proportions
baseline_check <- df_with_traits_all %>%
  filter(year == 1999) %>%
  summarise(
    total_density = sum(adjusted_density),
    # Calculate proportions for comparison with claims
    deposit_prop = sum(adjusted_density * ((f_Surface_deposit > 0 | f_Subsurface_deposit > 0) & f_Suspension == 0 & f_Predator == 0)) / total_density,
    deposit_suspension_prop = sum(adjusted_density * ((f_Surface_deposit > 0 | f_Subsurface_deposit > 0) & f_Suspension > 0)) / total_density,
    suspension_prop = sum(adjusted_density * (f_Surface_deposit == 0 & f_Subsurface_deposit == 0 & f_Suspension > 0)) / total_density
  )

print("\n1999 Baseline proportions (claimed vs. calculated):")
print("Claimed: Deposit (40.5%), Deposit/suspension (26.2%), Suspension (10.0%)")
print("Calculated:")
print(baseline_check %>% mutate(across(where(is.numeric), ~round(., 3))))

# 2. Check predator patterns in 2014
predator_check_2014 <- df_with_traits_all %>%
  filter(year == 2014) %>%
  filter(f_Predator > 0) %>%
  group_by(species) %>%
  summarise(
    total_abundance = sum(adjusted_total_count),
    density = sum(adjusted_density)
  ) %>%
  arrange(desc(total_abundance))

print("\nPredator species abundance in 2014:")
print(predator_check_2014)

# 3. Check Capitella capitata decline 2013-2014
capitella_check <- df_with_traits_all %>%
  filter(year %in% c(2013, 2014), 
         species == "Capitella capitata") %>%
  group_by(year) %>%
  summarise(
    abundance = sum(adjusted_total_count),
    density = sum(adjusted_density)
  )

print("\nCapitella capitata abundance 2013-2014:")
print(capitella_check)

# 4. Check meiofauna predator proportions across years
meiofauna_pred_check <- df_with_traits_all %>%
  filter(species %in% c("Syllis", "Monoculodes")) %>%
  group_by(year) %>%
  summarise(
    total_density = sum(adjusted_density),
    community_total = sum(df_with_traits_all$adjusted_density[df_with_traits_all$year == first(year)]),
    proportion = total_density / community_total
  )

print("\nMeiofauna predator proportions by year:")
print(meiofauna_pred_check %>% mutate(across(where(is.numeric), ~round(., 3))))

# 5. Compare 2017 to 1999 baseline
guild_comparison <- df_with_traits_all %>%
  filter(year %in% c(1999, 2017)) %>%
  group_by(year) %>%
  summarise(
    total_density = sum(adjusted_density),
    deposit = sum(adjusted_density * ((f_Surface_deposit > 0 | f_Subsurface_deposit > 0) & f_Suspension == 0 & f_Predator == 0)) / total_density,
    deposit_suspension = sum(adjusted_density * ((f_Surface_deposit > 0 | f_Subsurface_deposit > 0) & f_Suspension > 0)) / total_density,
    suspension = sum(adjusted_density * (f_Surface_deposit == 0 & f_Subsurface_deposit == 0 & f_Suspension > 0)) / total_density
  )

print("\nComparison of 1999 and 2017 feeding guild proportions:")
print(guild_comparison %>% mutate(across(where(is.numeric), ~round(., 3))))

# Check specific species mentioned in claims
print("\nChecking specific species claims:")

# Check Eteone longa and Harmothoe sp. in 2014
predator_species_2014 <- df_with_traits_all %>%
  filter(year == 2014,
         species %in% c("Eteone longa", "Harmothoe extenuata", "Asteroidea", "Anthozoa")) %>%
  group_by(species) %>%
  summarise(
    abundance = sum(adjusted_total_count),
    density = sum(adjusted_density)
  )

print("\nKey predator species in 2014:")
print(predator_species_2014)

# Check subsurface predators
subsurface_pred_2014 <- df_with_traits_all %>%
  filter(year == 2014,
         species %in% c("Nephthys", "Priapulus caudatus")) %>%
  group_by(species) %>%
  summarise(
    abundance = sum(adjusted_total_count),
    density = sum(adjusted_density)
  )

print("\nSubsurface predators in 2014:")
print(subsurface_pred_2014)

# Check meiofauna predator species
meiofauna_species <- df_with_traits_all %>%
  filter(species %in% c("Syllis cornuta", "Monoculodes")) %>%
  group_by(year, species) %>%
  summarise(
    abundance = sum(adjusted_total_count),
    density = sum(adjusted_density),
    total_community = sum(df_with_traits_all$adjusted_density[df_with_traits_all$year == first(year)]),
    proportion = density / total_community,
    .groups = 'drop'
  ) %>%
  arrange(year, species)

print("\nMeiofauna predator species by year:")
print(meiofauna_species %>% mutate(across(where(is.numeric), ~round(., 3))))


