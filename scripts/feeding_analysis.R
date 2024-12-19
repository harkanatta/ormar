# Improved feeding guild analysis
library(tidyverse)

# Function to calculate feeding guild proportions
calculate_guild_proportions <- function(data) {
  data %>%
    # Replace NAs with 0 for feeding traits
    mutate(across(starts_with("f_"), ~replace_na(., 0))) %>%
    group_by(year) %>%
    summarise(
      total_density = sum(adjusted_density),
      # Calculate mutually exclusive categories
      deposit = sum(adjusted_density * ((f_Surface_deposit > 0 | f_Subsurface_deposit > 0) & 
                                    f_Suspension == 0 & f_Predator == 0)) / total_density,
      deposit_suspension = sum(adjusted_density * ((f_Surface_deposit > 0 | f_Subsurface_deposit > 0) & 
                                               f_Suspension > 0 & f_Predator == 0)) / total_density,
      suspension = sum(adjusted_density * (f_Surface_deposit == 0 & f_Subsurface_deposit == 0 & 
                                      f_Suspension > 0 & f_Predator == 0)) / total_density,
      predator = sum(adjusted_density * (f_Predator > 0 & f_Surface_deposit == 0 & 
                                     f_Subsurface_deposit == 0 & f_Suspension == 0)) / total_density,
      predator_deposit = sum(adjusted_density * (f_Predator > 0 & 
                                             (f_Surface_deposit > 0 | f_Subsurface_deposit > 0) & 
                                             f_Suspension == 0)) / total_density,
      predator_suspension = sum(adjusted_density * (f_Predator > 0 & f_Suspension > 0)) / total_density,
      other = sum(adjusted_density * (f_Surface_deposit == 0 & f_Subsurface_deposit == 0 & 
                                  f_Suspension == 0 & f_Predator == 0)) / total_density
    ) %>%
    mutate(
      total = deposit + deposit_suspension + suspension + predator + 
              predator_deposit + predator_suspension + other
    )
}

# Function to create the stacked bar plot
create_guild_plot <- function(data) {
  # Prepare data for plotting
  plot_data <- data %>%
    select(-total) %>%
    pivot_longer(-c(year, total_density), 
                 names_to = "feeding_guild", 
                 values_to = "proportion") %>%
    mutate(
      year = as.factor(year),
      # Ensure no NAs
      proportion = replace_na(proportion, 0),
      feeding_guild = factor(feeding_guild, 
                           levels = c("deposit", "deposit_suspension", "suspension",
                                    "predator", "predator_deposit", 
                                    "predator_suspension", "other"),
                           labels = c("Deposit feeders",
                                    "Deposit/suspension feeders",
                                    "Suspension feeders",
                                    "Predators",
                                    "Predator/deposit feeders",
                                    "Predator/suspension feeders",
                                    "Other"))
    )
  
  # Print data before plotting for debugging
  cat("\nPlot data check:\n")
  print(plot_data %>% 
    group_by(year) %>% 
    summarise(total = sum(proportion, na.rm = TRUE)))
  
  # Create plot
  ggplot(plot_data, 
         aes(x = year, 
             y = proportion, 
             fill = feeding_guild)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_brewer(palette = "Set2") +
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
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                      limits = c(0, 1))  # Set y-axis limits from 0 to 1 (0-100%)
}

# Function to check missing traits and find potential matches
check_missing_traits <- function(data, trait_matrix) {
  # Get species missing traits
  missing_traits <- data %>%
    filter(year == 1999,
           is.na(f_Surface_deposit) & is.na(f_Subsurface_deposit) & 
           is.na(f_Suspension) & is.na(f_Predator)) %>%
    select(species, match_taxon, taxon_rank, adjusted_density) %>%
    arrange(desc(adjusted_density))
  
  cat("\nSpecies missing feeding traits in 1999:\n")
  print(missing_traits, n = Inf)
  
  # Extract genus names
  genera <- missing_traits %>%
    mutate(genus = word(species, 1)) %>%
    pull(genus) %>%
    unique()
  
  # Check for matches in trait matrix
  cat("\nPotential matches in trait matrix:\n")
  for(g in genera) {
    matches <- trait_matrix %>%
      filter(str_detect(tolower(Genus), tolower(g)) |
             str_detect(tolower(Family), tolower(g)) |
             str_detect(tolower(Order), tolower(g)) |
             str_detect(tolower(Class), tolower(g)))
    
    if(nrow(matches) > 0) {
      cat("\nMatches for", g, ":\n")
      print(matches %>% select(Genus, Family, Order, Class))
    }
  }
}

# Main analysis function
analyze_feeding_guilds <- function() {
  # Load data
  df_with_traits_all <- readRDS("data/df_with_traits_all.rds")
  trait_matrix <- read_csv("data/traits_data/TraitMatrix.csv", show_col_types = FALSE)
  
  # Check year format
  cat("\nYear format check:\n")
  print(unique(df_with_traits_all$year))
  
  # First standardize taxonomy and ensure year is numeric
  df_with_traits_fixed <- df_with_traits_all %>%
    mutate(
      year = as.numeric(year),  # Ensure year is numeric
      species = case_when(
        # Update taxonomy first
        tolower(species) == "leda" ~ "Nuculana pernula",
        year == 1999 & tolower(species) == "leucon acutirostris" ~ "Cumacea",
        year == 1999 & tolower(species) == "eudorella emarginata" ~ "Cumacea",
        TRUE ~ species
      )
    )
  
  # Then fix traits - replace NA with 0 first
  df_with_traits_fixed <- df_with_traits_fixed %>%
    mutate(
      across(starts_with("f_"), ~replace_na(., 0))
    ) %>%
    mutate(
      # Fix Oligochaeta traits (keep them in dataset)
      f_Surface_deposit = case_when(
        species == "Oligochaeta" ~ 1,
        TRUE ~ f_Surface_deposit
      ),
      f_Subsurface_deposit = case_when(
        species == "Oligochaeta" ~ 1,
        TRUE ~ f_Subsurface_deposit
      ),
      # Fix Cumacea traits
      f_Surface_deposit = case_when(
        species == "Cumacea" ~ 1,
        TRUE ~ f_Surface_deposit
      ),
      f_Suspension = case_when(
        species == "Cumacea" ~ 0.5,
        TRUE ~ f_Suspension
      ),
      # Fix Nuculana pernula (formerly Leda) traits
      f_Surface_deposit = case_when(
        species == "Nuculana pernula" ~ 1,
        TRUE ~ f_Surface_deposit
      ),
      # Fix Lanassa venusta traits (Terebellid deposit feeder)
      f_Surface_deposit = case_when(
        species == "Lanassa venusta" ~ 1,
        TRUE ~ f_Surface_deposit
      ),
      # Fix Yoldia hyperborea (deposit feeding bivalve)
      f_Surface_deposit = case_when(
        species == "Yoldia hyperborea" ~ 1,
        TRUE ~ f_Surface_deposit
      ),
      # Fix Paroediceros lynceus (predatory amphipod)
      f_Predator = case_when(
        species == "Paroediceros lynceus" ~ 1,
        TRUE ~ f_Predator
      )
    )
  
  # Calculate detailed feeding guild proportions
  feeding_props_detailed <- df_with_traits_fixed %>%
    group_by(year) %>%
    summarise(
      total_density = sum(adjusted_density),
      # Pure deposit feeders
      deposit = sum(adjusted_density * (
        (f_Surface_deposit > 0 | f_Subsurface_deposit > 0) & 
        f_Suspension == 0 & f_Predator == 0
      ), na.rm = TRUE),
      # Deposit/suspension feeders
      deposit_suspension = sum(adjusted_density * (
        (f_Surface_deposit > 0 | f_Subsurface_deposit > 0) & 
        f_Suspension > 0 & f_Predator == 0
      ), na.rm = TRUE),
      # Pure suspension feeders
      suspension = sum(adjusted_density * (
        f_Surface_deposit == 0 & f_Subsurface_deposit == 0 & 
        f_Suspension > 0 & f_Predator == 0
      ), na.rm = TRUE),
      # Pure predators
      predator = sum(adjusted_density * (
        f_Predator > 0 & f_Surface_deposit == 0 & 
        f_Subsurface_deposit == 0 & f_Suspension == 0
      ), na.rm = TRUE),
      # Predator combinations
      predator_deposit = sum(adjusted_density * (
        f_Predator > 0 & (f_Surface_deposit > 0 | f_Subsurface_deposit > 0) & 
        f_Suspension == 0
      ), na.rm = TRUE),
      predator_suspension = sum(adjusted_density * (
        f_Predator > 0 & f_Suspension > 0
      ), na.rm = TRUE),
      # Other
      other = sum(adjusted_density * (
        f_Surface_deposit == 0 & f_Subsurface_deposit == 0 & 
        f_Suspension == 0 & f_Predator == 0
      ), na.rm = TRUE)
    ) %>%
    mutate(
      # Calculate proportions
      across(c(deposit, deposit_suspension, suspension,
               predator, predator_deposit, predator_suspension, other),
             ~. / total_density),
      # Calculate total to verify
      total = deposit + deposit_suspension + suspension + 
              predator + predator_deposit + predator_suspension + other
    )
  
  # Print raw values before proportions for debugging
  cat("\nRaw values before proportions:\n")
  print(feeding_props_detailed %>%
    select(year, total_density, deposit, deposit_suspension, suspension) %>%
    mutate(across(where(is.numeric), ~round(., 1))))
  
  # Calculate species-specific densities for verification
  species_densities <- df_with_traits_fixed %>%
    filter(year %in% c(2013, 2014)) %>%
    group_by(year, species) %>%
    summarise(
      density = sum(adjusted_density),
      .groups = 'drop'
    ) %>%
    arrange(year, desc(density))
  
  # Calculate predator proportions over time
  predator_trends <- df_with_traits_fixed %>%
    group_by(year) %>%
    summarise(
      total_density = sum(adjusted_density),
      predator_density = sum(adjusted_density * (f_Predator > 0), na.rm = TRUE),
      predator_proportion = predator_density / total_density
    )
  
  # Create and save plot
  guild_plot <- create_guild_plot(feeding_props_detailed)
  ggsave("output/feeding_guild_temporal_changes_improved.png",
         guild_plot,
         width = 12,
         height = 8,
         dpi = 300,
         bg = "white")
  
  # Return all results
  list(
    proportions = feeding_props_detailed,
    species_densities = species_densities,
    predator_trends = predator_trends,
    plot = guild_plot,
    data = df_with_traits_fixed
  )
}

# Run the analysis
results <- analyze_feeding_guilds()

# Print detailed results
cat("\nFeeding guild proportions by year (%):\n")
print(results$proportions %>%
  mutate(across(where(is.numeric), ~round(. * 100, 1))))

cat("\nKey species densities 2013-2014 (ind./m²):\n")
print(results$species_densities %>%
  filter(species %in% c("Capitella capitata", "Eteone longa", 
                       "Harmothoe extenuata", "Syllis cornuta",
                       "Priapulus caudatus", "Nephthys")) %>%
  arrange(year, desc(density)))

cat("\nPredator proportion trends:\n")
print(results$predator_trends %>%
  mutate(across(where(is.numeric), ~round(., 3))))

# Verify feeding guild claims
cat("\nVerifying feeding guild claims:\n")

# 1. Check 1999 baseline proportions
baseline_check <- results$data %>%
  filter(year == 1999) %>%
  summarise(
    total_density = sum(adjusted_density),
    # Calculate proportions for comparison with claims
    deposit_prop = sum(adjusted_density * ((f_Surface_deposit > 0 | f_Subsurface_deposit > 0) & 
                                         f_Suspension == 0 & f_Predator == 0)) / total_density,
    deposit_suspension_prop = sum(adjusted_density * ((f_Surface_deposit > 0 | f_Subsurface_deposit > 0) & 
                                                    f_Suspension > 0)) / total_density,
    suspension_prop = sum(adjusted_density * (f_Surface_deposit == 0 & f_Subsurface_deposit == 0 & 
                                            f_Suspension > 0)) / total_density
  )

cat("\n1999 Baseline proportions (claimed vs. calculated):\n")
cat("Claimed: Deposit (40.5%), Deposit/suspension (26.2%), Suspension (10.0%)\n")
cat("Calculated:\n")
print(baseline_check %>% mutate(across(where(is.numeric), ~round(. * 100, 1))))

# 2. Check predator patterns in 2014
predator_check_2014 <- results$data %>%
  filter(year == 2014, f_Predator > 0) %>%
  group_by(species) %>%
  summarise(
    abundance = sum(adjusted_total_count),
    density = sum(adjusted_density)
  ) %>%
  arrange(desc(density))

cat("\nPredator species in 2014:\n")
print(predator_check_2014)

# 3. Check Capitella capitata decline
capitella_check <- results$data %>%
  filter(year %in% c(2013, 2014), 
         species == "Capitella capitata") %>%
  group_by(year) %>%
  summarise(
    abundance = sum(adjusted_total_count),
    density = sum(adjusted_density)
  )

cat("\nCapitella capitata abundance 2013-2014:\n")
print(capitella_check)

# 4. Check meiofauna predator trends
meiofauna_check <- results$data %>%
  filter(species %in% c("Syllis cornuta", "Monoculodes")) %>%
  group_by(year, species) %>%
  summarise(
    density = sum(adjusted_density),
    total_community = sum(results$data$adjusted_density[results$data$year == first(year)]),
    proportion = density / total_community,
    .groups = 'drop'
  ) %>%
  arrange(year, species)

cat("\nMeiofauna predator proportions:\n")
print(meiofauna_check %>% mutate(across(where(is.numeric), ~round(., 3)))) 