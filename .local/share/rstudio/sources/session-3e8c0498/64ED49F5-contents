# NMDS Analysis

# Load required packages
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)


# Load cleaned data
#df_result <- result # just so this won't be accidentally changed
df <- df_result


# skráin hans Jörundar
henda <- read_csv(
  "data/processed/hendaB.csv",
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

# breyta á long format

# Transform henda from wide to long format
henda_long <- henda %>%
  pivot_longer(cols = c("A7","B5", "B8",  "C4",  "E3", "E4") , # Select station columns
               names_to = "station",                     # New column for station names
               values_to = "adjusted_density") %>%       # New column for density values
  rename(species = Flokkun,                              # Rename 'Flokkun' to 'species'
         year = Year) %>%                                # Rename 'Year' to 'year'
  select(station, year, species, adjusted_density)       # Select and order columns

# henda út flokkunareiningum
# Anthozoa               
# Balanus                
# Bryozoa                
# Ciliatocardium ciliatum
# Cystenides granulata   
# Maldane sarsi          
# Ostracoda              
# Paraonidae             
# Philomedes globosus    
# Platyhelminthes        
# Sipuncula              
# Verruca stroemia       


# Required packages
library(dplyr)
library(stringr)

# Step 1: Remove " sp." from species names
henda_long_cleaned <- henda_long %>%
  mutate(species = str_replace(species, " sp\\.", ""))

# Step 2: Match with reference dataset and remove non-matching species
henda_long_matched <- henda_long_cleaned %>%
  semi_join(df, by = "species")

# Step 3: Quick validation check
cat("Data Cleaning Summary:\n",
    "Original rows:", nrow(henda_long), "\n",
    "Final rows:", nrow(henda_long_matched), "\n",
    "Rows removed:", nrow(henda_long) - nrow(henda_long_matched), "\n")

# Step 4: View removed species for validation
removed_species <- henda_long_cleaned %>%
  anti_join(df, by = "species") %>%
  group_by(species) %>%
  summarise(
    observations = n(),
    total_abundance = sum(adjusted_density),
    years = paste(sort(unique(year)), collapse = ", "),
    stations = paste(sort(unique(station)), collapse = ", ")
  ) %>%
  arrange(desc(observations))

# Print removed species summary
print(removed_species, n = Inf, width = Inf)


# Prepare data for NMDS - exclude 1999 data
#species_matrix <- df %>% #Mín gögn
species_matrix <- henda_long_matched %>% #Gögn Jörundar
  #filter(year >= 2013) %>%  # Only include 2013-2017 data
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
   # Add points for sites
   geom_point(data = site_scores, 
              aes(x = NMDS1, y = NMDS2, color = Station, shape = Year),
              size = 4) +
   # Add smaller ellipses around years with solid lines
    # stat_ellipse(data = site_scores,
    #              aes(x = NMDS1, y = NMDS2, group = Year),
    #              type = "norm",
    #              level = 0.7,
    #              linetype = 1,
    #              size = 0.5) +
   # Add significant environmental vectors
   # geom_segment(data = sig_env_vectors,
   #              aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
   #              arrow = arrow(length = unit(0.2, "cm")),
   #              color = "darkgrey") +
   # geom_text(data = sig_env_vectors,
   #           aes(x = NMDS1*1.1, y = NMDS2*1.1, label = variable),
   #           size = 3, color = "darkgrey") +
   theme_bw() +
   theme(
     panel.background = element_rect(fill = "white"),
     axis.text = element_text(size = 16, color = "black"),     # Increased from 12
     axis.title = element_text(size = 18),                     # Increased from 14
     legend.text = element_text(size = 16),                    # Increased from 12
     legend.title = element_text(size = 16, face = "bold"),    # Increased from 12
     plot.title = element_text(size = 20, face = "bold"),      # Increased from 16
     legend.position = "right"
   ) +
   labs(title = "NMDS of Community Composition (1999 and 2013-2017)",
        caption = "") +
   guides(color = guide_legend("Station"),
          shape = guide_legend("Year"))
 
 # Save with higher resolution and size
 ggsave("output/nmds_plot_1999_an-hringja.png", 
        nmds_plot, 
        width = 12,      # Increased width
        height = 10,     # Increased height
        dpi = 300,       # Higher resolution
        bg = "white")

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
