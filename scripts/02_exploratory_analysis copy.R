# Load required packages
library(readr)
library(dplyr)
library(ggplot2)
library(vegan)
library(tidyr)
library(tibble)  # For column_to_rownames()

KolgrTaxa <- read_csv(file("data/raw/KolgrTaxa.csv", encoding = "UTF-8"), na = "empty")
# Convert taxa to remove to lowercase
taxa_to_remove <- tolower(c(
  "Foraminifera", "Nematoda", "Cirripedia", "Porifera", 
  "Cnidaria", "Bryozoa", "Sipuncula", "Platyhelminthes", 
  "Nemertea", "Ostracoda", "harpacticoida"#, "Oligochaeta",
))

# First ensure all character columns are properly encoded
KolgrTaxa <- KolgrTaxa %>%
  mutate(across(where(is.character), ~ iconv(., from = "UTF-8", to = "UTF-8", sub = "")))

# Then filter using if_any
KolgrTaxa <- KolgrTaxa %>%
  filter(!if_any(where(is.character), ~ tolower(.) %in% taxa_to_remove))

# Define patterns to remove and create regex pattern
remove_patterns <- c("nýsestir", "ungviði", "ungv", "ungv\\.", "juv")
pattern <- paste(remove_patterns, collapse = "|")

# Filter out rows with matching patterns in 'gamalt'
KolgrTaxa <- KolgrTaxa[!grepl(pattern, tolower(KolgrTaxa$gamalt)), ]

# Check final dimensions
print(paste("Final dimensions:", nrow(KolgrTaxa), "rows,", ncol(KolgrTaxa), "columns"))

   KolgrTaxa <- KolgrTaxa %>% 
     mutate(
       Flokkun = case_when(
         # General taxonomic standardization (no year-specific)
         str_detect(tolower(Flokkun), "^sipuncul") ~ NA_character_,
         tolower(Flokkun) == "terebellides stroemi" ~ "Terebellides stroemii",
         tolower(Flokkun) == "bivalvia" ~ NA_character_,
         # 2016 specific matches
         Artal == "2016" & tolower(Flokkun) == "leaena ebranchiata" ~ "Terebellides stroemii",
         Artal == "2016" & tolower(Flokkun) == "pseudopolydora antennata" ~ "Polydora",
         Artal == "2016" & tolower(Flokkun) == "cistenides hyperborea" ~ "Cistenides granulata",
         Artal == "2016" & tolower(Flokkun) == "amphipoda" ~ NA_character_,
         Artal == "2016" & tolower(Flokkun) == "bivalvia" ~ NA_character_,
         Artal == "2016" & tolower(Flokkun) == "cardium" ~ NA_character_,
         
         # Default case
         TRUE ~ Flokkun
       )
     ) %>% 
     drop_na(Flokkun)    

KolgrTaxa <- KolgrTaxa  %>%
   rename(
    species = Flokkun,
    sample_id = id,
    year = Artal,
    station = stod,
    subdivision = skipting,
    count = N,
    density = Nu
  ) %>% 
   filter(station %in% c("C4", "A7", "B5", "B8", "E4", "E3")) %>%
  # Rename stations
  mutate(legacy_stations = station) %>%
  mutate(station = case_when(
      legacy_stations == "E3" ~ "St-1",  # start near shore?
      legacy_stations == "E4" ~ "St-2",
      legacy_stations == "C4" ~ "St-3",
      legacy_stations == "B5" ~ "St-4",
      legacy_stations == "A7" ~ "St-5",
      legacy_stations == "B8" ~ "St-6",
      TRUE ~ legacy_stations # Keep original name if not in the list (optional safeguard)
  ))


  # First, apply the subdivision to the count
  result <- KolgrTaxa %>%
  mutate(adjusted_count = count * subdivision) %>%
  
  # Calculate the correct number of samples and total area for each station-year combination
  group_by(legacy_stations,station, year) %>%
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
  group_by(legacy_stations,station, year, species) %>%
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


### summary statistic

# various Benthic Biotic Indices (BBI) for different years and stations. It then creates heatmaps to visualize these indices. This provides a good overview of the ecological status of the benthic communities over time and across different sampling sites.
library(BBI)
library(gplots)

# Calculate Benthic Biotic Indices (BBI) for each year and station
calculate_biotic_indices <- function(data) {
  # Validate input data
  required_cols <- c("year", "station", "species", "adjusted_density")  # Changed from count to adjusted_density
  stopifnot(
    "Missing required columns" = all(required_cols %in% colnames(data)),
    "No data provided" = nrow(data) > 0,
    "Negative density values found" = all(data$adjusted_density >= 0),
    "Missing values found" = !any(is.na(data$adjusted_density))
  )
  
  # Define target stations
  #target_stations <- c("C4", "A7", "B5", "B8", "E4", "E3")
  
  # Initialize lists to store results
  indices_list <- list()
  
  # Calculate indices for each year
  for (current_year in sort(unique(data$year))) {
    # Prepare data for BBI calculation using density values
    year_data <- data %>%
      filter(year == current_year #,
             #station %in% target_stations) 
             ) %>%
      group_by(station, species) %>%
      summarise(density = sum(adjusted_density), .groups = 'drop') %>%  # Using adjusted_density
      pivot_wider(names_from = station, 
                  values_from = density, 
                  values_fill = 0)
    
    # Calculate BBI indices
    bbi_results <- BBI(year_data)
    
    # Store results
    indices_list[[as.character(current_year)]] <- list(
      indices = as.data.frame(cbind(bbi_results$BBI, year = current_year)),
      classification = bbi_results$BBIclass,
      normalized = as.data.frame(nEQR(bbi_results$BBI)[1])
    )
  }
  
  return(indices_list)
}

# Use the result dataframe that already has the correct density calculations
biotic_indices <- calculate_biotic_indices(result)

# Extract normalized EQR values for further analysis
normalized_eqr <- do.call(rbind, lapply(biotic_indices, function(x) x$normalized))

# Prepare data for heatmap
names(normalized_eqr) <- c("nAMBI", "nISI", "nNSI", "nNQI1", "nShannon", "nEQR")
#stations <- rep(c("C4", "A7", "B5", "B8", "E4", "E3"), length(unique(df$year)))
stations <- rep(c("St-3", "St-5", "St-4", "St-6", "St-2", "St-1"), length(unique(df$year)))
years <- rep(sort(unique(df$year)), each = 6)
row_labels <- paste(stations, years, sep = " - ")

normalized_eqr_subset <- normalized_eqr[, c(1, 4, 5)]  # Select columns 1 (nAMBI), 4 (nNQI1), 5 (nShannon)
# Create a matrix for heatmap
heatmap_matrix <- as.matrix(normalized_eqr_subset)
rownames(heatmap_matrix) <- row_labels

# Define color palette
color_palette <- colorRampPalette(c("red", "yellow", "green"))(100)

# Create heatmap
png("output/BBI_heatmap.png", width = 1000, height = 1000, res = 100)
heatmap.2(heatmap_matrix,
          main = "Benthic Biotic Indices\nby Station and Year",
          xlab = "Indices",
          ylab = "",
          col = color_palette,
          scale = "none",
          key = TRUE,
          keysize = 1.5,
          density.info = "none",
          trace = "none",
          Rowv = T,
          Colv = FALSE,
          dendrogram = "row",
          margins = c(5, 8),
          cexRow = 0.9,
          cexCol = 1.2,
          srtCol = 45,
          adjCol = c(1, 1),
          cellnote = round(heatmap_matrix, 2),
          notecol = "black",
          notecex = 0.7,
          labRow = row_labels,
          labCol = colnames(heatmap_matrix))
dev.off()



# Function to extract normalized values and calculate summary statistics
get_yearly_stats <- function(year_data) {
    df <- year_data$normalized
    
    # Calculate means and SDs for required indices
    means <- c(
        mean(df$nEQR.nAMBI),
        mean(df$nEQR.nNQI1),
        mean(df$nEQR.nShannon)
    )
    
    sds <- c(
        sd(df$nEQR.nAMBI),
        sd(df$nEQR.nNQI1),
        sd(df$nEQR.nShannon)
    )
    
    # Calculate mean status
    mean_status <- mean(c(
        mean(df$nEQR.nAMBI),
        mean(df$nEQR.nNQI1),
        mean(df$nEQR.nShannon)
    ))
    
    status <- case_when(
        mean_status < 0.2 ~ "Bad",
        mean_status < 0.4 ~ "Poor",
        mean_status < 0.6 ~ "Moderate",
        mean_status < 0.8 ~ "Good",
        TRUE ~ "High"
    )
    
    return(c(means, sds, status))
}

# Create empty dataframe for results
years <- names(biotic_indices)
summary_table <- data.frame(
    Year = as.numeric(years),
    matrix(NA, nrow = length(years), ncol = 7)
)

# Fill the dataframe
for(i in seq_along(years)) {
    stats <- get_yearly_stats(biotic_indices[[years[i]]])
    summary_table[i, 2:8] <- stats
}

# Set column names
colnames(summary_table) <- c("Year", 
                             "nAMBI_mean", "nNQI1_mean", "nShannon_mean",
                             "nAMBI_sd", "nNQI1_sd", "nShannon_sd",
                             "Status")

# Format the table for presentation
# Format the table for presentation
final_table <- summary_table %>%
    mutate(
        across(ends_with("_mean"), as.numeric),
        across(ends_with("_sd"), as.numeric),
        nAMBI = sprintf("%.2f ± %.2f", as.numeric(nAMBI_mean), as.numeric(nAMBI_sd)),
        nNQI1 = sprintf("%.2f ± %.2f", as.numeric(nNQI1_mean), as.numeric(nNQI1_sd)),
        nShannon = sprintf("%.2f ± %.2f", as.numeric(nShannon_mean), as.numeric(nShannon_sd))
    ) %>%
    select(Year, nAMBI, nNQI1, nShannon, Status)

# Alternative approach if the above still gives errors
final_table <- summary_table %>%
    mutate(
        nAMBI = paste0(round(as.numeric(nAMBI_mean), 2), " ± ", 
                       round(as.numeric(nAMBI_sd), 2)),
        nNQI1 = paste0(round(as.numeric(nNQI1_mean), 2), " ± ", 
                       round(as.numeric(nNQI1_sd), 2)),
        nShannon = paste0(round(as.numeric(nShannon_mean), 2), " ± ", 
                          round(as.numeric(nShannon_sd), 2))
    ) %>%
    select(Year, nAMBI, nNQI1, nShannon, Status)

# Print using kable
library(knitr)
kable(final_table,
      caption = "Mean normalized benthic indices (± SD) across all stations showing temporal changes in Kolgrafafjörður fjord from 1999 to 2017.",
      align = c("c", "c", "c", "c", "l"))

library(flextable)
library(officer)  # <- this is required for fp_border()

ft <- flextable(final_table) %>%
    # set_caption(
    #   caption = "Mean normalized benthic indices (± SD) across all stations showing temporal changes in Kolgrafafjörður fjord from 1999 to 2017. *Ecological Status is assigned based on classification from the Marine Research Institute (MRI).",
    #   autonum = NULL,
    #   fp_p = fp_par(text.align = "center"),
    #   caption_pos = "bottom"
    # ) %>%
    align(j = 1:4, align = "center") %>%
    align(j = 5, align = "left") %>%
    set_header_labels(
        Year = "Year",
        nAMBI = "nAMBI",
        nNQI1 = "nNQI1",
        nShannon = "nShannon",
        Status = "Ecological Status*"
    ) %>%
    theme_box() %>%
    fontsize(size = 10) %>%
    padding(padding = 4) %>%
    border_outer(border = fp_border(width = 1.5)) %>%
    border_inner_h(border = fp_border(width = 0.5)) %>%
    autofit() %>%
    bold(part = "header") %>%
    bg(bg = "white", part = "all") %>%
    colformat_double(j = "Year", digits = 0, big.mark = "")

save_as_docx(
  ft,
  path = "output/benthic_indices_table.docx"
)



# Calculate the statistics (sjá hvaða tegundir eru mest eða minnst)
summarise_station_stats <- function(data, slice_n = 1) {
    stats <- data %>%
        #filter(species == "Capitella capitata") %>%
        group_by(station, legacy_stations, year) %>%
        summarise(
            n_taxa = n_distinct(species),
            total_density = sum(adjusted_density),
            .groups = 'drop'
        )
    
    get_top_n <- function(df, col, decreasing = FALSE) {
        df %>%
            arrange(if (decreasing) desc(.data[[col]]) else .data[[col]]) %>%
            slice_head(n = slice_n)
    }
    
    list(
        stats = stats,
        min_taxa = get_top_n(stats, "n_taxa", decreasing = FALSE),
        max_taxa = get_top_n(stats, "n_taxa", decreasing = TRUE),
        min_density = get_top_n(stats, "total_density", decreasing = FALSE),
        max_density = get_top_n(stats, "total_density", decreasing = TRUE)
    )
}

summary_results <- summarise_station_stats(result, slice_n = 3)
summary_results$min_taxa      # Will return up to 3 rows with min n_taxa
summary_results$max_density   # Up to 3 rows with max total_density



# Myndir fig. 3 Mín gögn
# Define the desired order of years
year_levels <- c("1999", "2013", "2014", "2015", "2016", "2017")

# Process the data and convert 'year' to a factor with specified levels
df_processed <- KolgrTaxa %>% 
  # Uncomment the line below if you need to exclude zero densities
  # filter(adjusted_density > 0) %>%  # Exclude zero densities if needed
  group_by(station, year) %>%
  summarise(n_species = n_distinct(species), .groups = "drop") %>%  # Count unique species
  mutate(
    year = as.character(year),
    year_factor = factor(year, levels = year_levels)
  )

# Data for solid lines (excluding 1999)
df_solid <- df_processed %>%
  filter(year != "1999") %>%
  mutate(year_factor = factor(year, levels = year_levels))

# Data for dotted lines connecting 1999 to the next year
df_dotted <- df_processed %>%
  group_by(station) %>%
  arrange(as.numeric(year)) %>%
  mutate(year_numeric = as.numeric(year)) %>%
  filter(year %in% c("1999", as.character(min(year_numeric[year_numeric > 1999])))) %>%
  ungroup() %>%
  mutate(year_factor = factor(year, levels = year_levels))

# Create the plot
p <- ggplot() +
  # Solid lines excluding 1999
  geom_line(
    data = df_solid,
    aes(x = year_factor, y = n_species, color = station, group = station),
    size = 1
  ) +
  # Dotted lines connecting 1999 to the next year
  geom_line(
    data = df_dotted,
    aes(x = year_factor, y = n_species, color = station, group = station),
    linetype = "dotted", size = 1
  ) +
  # Points for all years
  geom_point(
    data = df_processed,
    aes(x = year_factor, y = n_species, color = station),
    size = 2
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Number of Species Across Stations (1999-2017)",
    x = "",
    y = "Number of Species",
    color = "Station"
  ) +
  scale_x_discrete(limits = year_levels) +
  theme(
    axis.text.x = element_text(hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

# Display the plot
print(p)


ggsave("number_of_species_across_stations_2013-2017_min-gogn.png", plot = p, width = 8, height = 6, units = "in", dpi = 600, bg="white")


###########################
# Myndir fig3 Jörundar gögn
###########################

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
  semi_join(KolgrTaxa, by = "species")

# Step 3: Quick validation check
cat("Data Cleaning Summary:\n",
    "Original rows:", nrow(henda_long), "\n",
    "Final rows:", nrow(henda_long_matched), "\n",
    "Rows removed:", nrow(henda_long) - nrow(henda_long_matched), "\n")

# Step 4: View removed species for validation
removed_species <- henda_long_cleaned %>%
  anti_join(KolgrTaxa, by = "species") %>%
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


year_levels <- c("1999", "2013", "2014", "2015", "2016", "2017")

# Process the data and convert 'year' to a factor with specified levels
df_processed <- henda_long_cleaned %>% 
  # Uncomment the line below if you need to exclude zero densities
  filter(adjusted_density > 0) %>%  # Exclude zero densities if needed
  group_by(station, year) %>%
  summarise(n_species = n_distinct(species), .groups = "drop") %>%  # Count unique species
  mutate(
    year = as.character(year),
    year_factor = factor(year, levels = year_levels)
  )

# Data for solid lines (excluding 1999)
df_solid <- df_processed %>%
  filter(year != "1999") %>%
  mutate(year_factor = factor(year, levels = year_levels))

# Data for dotted lines connecting 1999 to the next year
df_dotted <- df_processed %>%
  group_by(station) %>%
  arrange(as.numeric(year)) %>%
  mutate(year_numeric = as.numeric(year)) %>%
  filter(year %in% c("1999", as.character(min(year_numeric[year_numeric > 1999])))) %>%
  ungroup() %>%
  mutate(year_factor = factor(year, levels = year_levels))

# Create the plot
p <- ggplot() +
  # Solid lines excluding 1999
  geom_line(
    data = df_solid,
    aes(x = year_factor, y = n_species, color = station, group = station),
    size = 1
  ) +
  # Dotted lines connecting 1999 to the next year
  geom_line(
    data = df_dotted,
    aes(x = year_factor, y = n_species, color = station, group = station),
    linetype = "dotted", size = 1
  ) +
  # Points for all years
  geom_point(
    data = df_processed,
    aes(x = year_factor, y = n_species, color = station),
    size = 2
  ) +
  theme_minimal(base_size = 12) +
  labs(
    #title = "Number of Species Across Stations (1999-2017)",
    x = "",
    y = "Number of Species",
    color = "Station"
  ) +
  scale_x_discrete(limits = year_levels) +
#theme( # miðað við fulla breidd á blaðsíðu
#  axis.text.x = element_text(hjust = 1, size = 8),
#  axis.text.y = element_text(size = 8),
#  axis.title = element_text(size = 10, face = "bold"),
#  legend.title = element_text(size = 10),
#  legend.text = element_text(size = 8),
#  plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
#)

theme(# fyrir hálfa blaðsíðu
  axis.text.x = element_text(hjust = 1, size = 9),
  axis.text.y = element_text(size = 9),
  axis.title = element_text(size = 10, face = "bold"),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),
  plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
)



# Display the plot
print(p)

ggsave( # miðað við fulla breidd á blaðsíðu
  "supplementary/number_of_species_across_stations_2013-2017.png",
  plot = p,
  width = 174,
  height = 117,  # or whatever height you want
  units = "mm",  # key: use mm to match journal
  dpi = 600,
  bg = "white"
)
ggsave( # miðað við hálfa breidd á blaðsíðu
  "supplementary/number_of_species_across_stations_2013-2017_half-column.png",
  plot = p,
  width = 87,
  height = 117,  # e.g. 60–100 mm
  units = "mm",
  dpi = 600,
  bg = "white"
)




##########
#n$hannon#
##########

# Use the result dataframe that already has the correct density calculations
biotic_indices <- calculate_biotic_indices(henda_long_cleaned) #hendalong cleaned er gögnin hans Jörundar.

# Extract normalized EQR values for further analysis
normalized_eqr <- do.call(rbind, lapply(biotic_indices, function(x) x$normalized))

# Prepare data for heatmap
names(normalized_eqr) <- c("nAMBI", "nISI", "nNSI", "nNQI1", "nShannon", "nEQR")

normalized_eqr$Year <- as.numeric(sub("\\..*", "", rownames(normalized_eqr)))
normalized_eqr$Station <- sub(".*\\.", "", rownames(normalized_eqr))
# Generate the plot

library(ggplot2)
library(dplyr)

# Define the desired order of years
year_levels <- c("1999", "2013", "2014", "2015", "2016", "2017")

# Convert Year to character and then to factor with specified levels
normalized_eqr <- normalized_eqr %>%
  mutate(
    Year = as.character(Year),
    Year_Factor = factor(Year, levels = year_levels)
  )

# Main plot
p <- ggplot(normalized_eqr, aes(x = Year_Factor, y = nShannon, color = Station, group = Station)) +
  
  # Main solid lines (excluding 1999)
  geom_line(
    data = normalized_eqr %>%
      filter(Year != "1999") %>%
      mutate(Year_Factor = factor(Year, levels = year_levels)),
    size = 1
  ) +
  
  # Dotted lines connecting 1999 to the next year
  geom_line(
    data = normalized_eqr %>%
      group_by(Station) %>%
      mutate(Year_Num = as.numeric(Year)) %>%
      arrange(Year_Num) %>%
      filter(Year %in% c("1999", as.character(min(Year_Num[Year_Num > 1999])))) %>%  # Only 1999 and the next year
      ungroup() %>%
      mutate(Year_Factor = factor(Year, levels = year_levels)),
    aes(x = Year_Factor, y = nShannon, color = Station, group = Station),
    linetype = "dotted", size = 1
  ) +
  
  # Add points for all years
  geom_point(aes(x = Year_Factor), size = 2) +
  
  theme_minimal(base_size = 12) +
  labs(
    #title = "nShannon Trends Across Stations (1999-2017)",
    x = "",
    y = "nShannon",
    color = "Station"
  ) +
#theme( # miðað við heila blaðsíðu
#  axis.text.x = element_text(hjust = 1, size = 8),
#  axis.text.y = element_text(size = 8),
#  axis.title = element_text(size = 10, face = "bold"),
#  legend.title = element_text(size = 10),
#  legend.text = element_text(size = 8),
#  plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
#) +
theme( # miðað við hálfa blaðsíðu
  axis.text.x = element_text(hjust = 1, size = 9),
  axis.text.y = element_text(size = 9),
  axis.title = element_text(size = 10, face = "bold"),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),
  plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
) +
  scale_x_discrete(limits = year_levels)



# Display the plot
#print(p)

ggsave( # miðað við heila blaðsíðu
  "supplementary/nshannon_fullwidth.png",
  plot = p,
  width = 174,
  height = 117,  # or adjust height as needed
  units = "mm",
  dpi = 600,
  bg = "white"
)

ggsave( # miðað við hálfa blaðsíðu
  "supplementary/nshannon_halfwidth.png",
  plot = p,
  width = 87,
  height = 100,  # Adjust height for layout balance
  units = "mm",
  dpi = 600,
  bg = "white"
)






# Ensure ggpubr is installed: install.packages("ggpubr")

create_combined_plot_shared_legend <- function(henda_long_cleaned, normalized_eqr, year_levels) {
    # Load necessary libraries
    library(ggplot2)
    library(dplyr)
    library(cowplot)
    library(ggpubr)
    library(grid)
    
    # --- Process data (same as before) ---
    df_species <- henda_long_cleaned %>%
        filter(adjusted_density > 0) %>%
        rename(Station = station) %>%
        group_by(Station, year) %>%
        summarise(n_species = n_distinct(species), .groups = "drop") %>%
        mutate(
            year = as.character(year),
            year_factor = factor(year, levels = year_levels)
        )
    df_eqr <- normalized_eqr %>%
        mutate(
            Year = as.character(Year),
            Year_Factor = factor(Year, levels = year_levels)
        )
    unique_stations <- sort(unique(c(df_species$Station, df_eqr$Station)))
    df_species$Station <- factor(df_species$Station, levels = unique_stations)
    df_eqr$Station <- factor(df_eqr$Station, levels = unique_stations)
    if (length(unique_stations) < 2) {
        warning("Only one unique station found. Legend might not be necessary.")
    }
    message("Unique stations for legend (sorted): ", paste(unique_stations, collapse=", "))
    
    # --- Create minimal plot for legend (same as before) ---
    legend_data <- data.frame(Station = factor(unique_stations, levels = unique_stations), dummy_y = 1)
    p_minimal_for_legend <- ggplot(legend_data, aes(x = 1, y = dummy_y, color = Station)) +
        geom_point(size = 2, show.legend = TRUE) +
        scale_color_discrete(name = "Station", limits = unique_stations, drop = FALSE) +
        theme_minimal() +
        theme(
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.title = element_text(hjust = 0.5),
            legend.box = "horizontal"
        ) +
        guides(color = guide_legend(nrow = 1, title.position = "top"))
    
    # --- Extract the legend using ggpubr::get_legend ---
    shared_legend_raw <- ggpubr::get_legend(p_minimal_for_legend) # Store raw output
    
    # --- Refined Check and Extraction ---
    shared_legend <- NULL # Initialize as NULL
    is_valid_legend <- FALSE
    
    if (inherits(shared_legend_raw, "gtable")) {
        shared_legend <- shared_legend_raw
        is_valid_legend <- TRUE
        message("ggpubr::get_legend returned a gtable directly.")
    } else if (is.list(shared_legend_raw) && length(shared_legend_raw) > 0 && inherits(shared_legend_raw[[1]], "gtable")) {
        shared_legend <- shared_legend_raw[[1]] # *** Explicitly extract the gtable ***
        is_valid_legend <- TRUE
        message("ggpubr::get_legend returned a list; extracted gtable from first element.")
    } else {
        # Fallback to cowplot (same as before)
        message("ggpubr::get_legend did not return a recognizable gtable or list containing one (got ", class(shared_legend_raw)[1], "). Trying cowplot::get_legend...")
        shared_legend_cowplot <- cowplot::get_legend(p_minimal_for_legend)
        if(inherits(shared_legend_cowplot, "gtable")) {
            message("cowplot::get_legend succeeded as fallback.")
            shared_legend <- shared_legend_cowplot
            is_valid_legend <- TRUE
        } else {
            warning("Both get_legend methods failed (cowplot got ", class(shared_legend_cowplot)[1], "). Proceeding without legend.")
            # shared_legend remains NULL
            is_valid_legend <- FALSE
        }
    }
    
    # --- Create p1 and p2 (same as before, ensure legend.position='none') ---
    p1 <- ggplot(df_species, aes(x = year_factor, y = n_species, color = Station, group = Station)) +
        geom_line(data = df_species %>% filter(year != "1999"), linewidth = 1) +
        geom_line(
            data = df_species %>%
                group_by(Station) %>%
                arrange(as.numeric(year)) %>%
                filter(year %in% c("1999", as.character(min(as.numeric(year)[as.numeric(year) > 1999])))) %>%
                ungroup(),
            aes(x = year_factor, y = n_species),
            linetype = "dotted", linewidth = 1
        ) +
        geom_point(size = 2) +
        scale_x_discrete(limits = year_levels) +
        scale_color_discrete(limits = unique_stations, drop = FALSE) +
        labs(x = "", y = "Number of Species") +
        theme_minimal(base_size = 12) +
        theme(
            axis.text.x = element_text(size = 8, hjust = 1),
            axis.text.y = element_text(size = 8),
            axis.title = element_text(size = 10, face = "bold"),
            legend.position = "none", # REMOVE LEGEND
            plot.title = element_blank()
        )
    
    p2 <- ggplot(df_eqr, aes(x = Year_Factor, y = nShannon, color = Station, group = Station)) +
        geom_line(data = df_eqr %>% filter(Year != "1999"), linewidth = 1) +
        geom_line(
            data = df_eqr %>%
                group_by(Station) %>%
                arrange(as.numeric(Year)) %>%
                filter(Year %in% c("1999", as.character(min(as.numeric(Year)[as.numeric(Year) > 1999])))) %>%
                ungroup(),
            aes(x = Year_Factor, y = nShannon),
            linetype = "dotted", linewidth = 1
        ) +
        geom_point(size = 2) +
        scale_x_discrete(limits = year_levels) +
        scale_color_discrete(limits = unique_stations, drop = FALSE) +
        labs(x = "", y = "nShannon") +
        theme_minimal(base_size = 12) +
        theme(
            axis.text.x = element_text(size = 8, hjust = 1),
            axis.text.y = element_text(size = 8),
            axis.title = element_text(size = 10, face = "bold"),
            legend.position = "none", # ENSURE LEGEND IS OFF
            plot.title = element_blank()
        )
    
    
    # --- Combine Plots without legends ---
    plots_row <- plot_grid(p1, p2, nrow = 1, align = "hv", rel_widths = c(1, 1))
    message("Class of plots_row: ", paste(class(plots_row), collapse=", ")) # Should be ggplot or gtable
    
    # --- Combine Legend and Plot Row (only if legend is valid) ---
    if (is_valid_legend && !is.null(shared_legend)) {
        # *** Add a final check on the legend object ***
        if (!inherits(shared_legend, "gtable")) {
            warning("Object 'shared_legend' is NOT a gtable right before final plot_grid. Class: ", paste(class(shared_legend), collapse=", "))
            # Fallback: don't add legend if it's somehow wrong type here
            combined_plot_final <- plots_row
            message("Final plot created WITHOUT shared legend due to unexpected type.")
        } else {
            message("Attempting final plot_grid. Class of shared_legend: ", paste(class(shared_legend), collapse=", ")) # Should be gtable
            combined_plot_final <- plot_grid(
                shared_legend, # This MUST be a gtable
                plots_row,     # This should be ok
                ncol = 1,
                rel_heights = c(0.15, 1)
            )
            message("Final plot created with shared legend.")
        }
        
    } else {
        # If legend extraction failed initially
        message("Final plot created WITHOUT shared legend (extraction failed or object invalid).")
        combined_plot_final <- plots_row
    }
    
    return(combined_plot_final)
}

# --- Example Usage ---
 final_plot <- create_combined_plot_shared_legend(henda_long_cleaned, normalized_eqr, year_levels)
 print(final_plot)

ggsave( # miðað við heila blaðsíðu
  "supplementary/nshannon_and_species_fullwidth_combined.png",
  plot = final_plot,
  width = 174,
  height = 117,  # Adjust height for layout balance
  units = "mm",
  dpi = 600,
  bg = "white"
)




#################
#fig3 organic matter#
#################


kornast_organic_tap_dypi <- structure(list(stod = c("A7", "A7", "A7", "A7", "A7", "B5", "B5", 
"B5", "B5", "B5", "B8", "B8", "B8", "B8", "B8", "C4", "C4", "C4", 
"C4", "C4", "E3", "E3", "E3", "E3", "E3", "E4", "E4", "E4", "E4", 
"E4"), Artal = c("2013", "2014", "2015", "2016", "2017", "2013", 
"2014", "2015", "2016", "2017", "2013", "2014", "2015", "2016", 
"2017", "2013", "2014", "2015", "2016", "2017", "2013", "2014", 
"2015", "2016", "2017", "2013", "2014", "2015", "2016", "2017"
), `20um` = c(0.628496042216359, 0.343339587242026, 0.453474676089517, 
0.390734265734266, 0.643465909090909, 0.386002120890774, 0.321766561514196, 
0.386792452830189, 0.19559585492228, 0.546637744034707, 0.592529025744573, 
0.515981735159817, 0.369913686806412, 0.424997544719429, 0.221565731166913, 
0.332928794917986, 0.153116531165312, 0.234791889007471, 0.190773067331671, 
0.152190051967335, 0.385025817555938, 0.211453744493392, 0.157303370786517, 
0.0827048768225239, 0.247269116186693, 0.426682936750051, 0.406354515050167, 
0.214331413947537, 0.0125582337451894, 0.286830357142857), `63um` = c(0.350395778364116, 
0.24577861163227, 0.215547703180212, 0.215909090909091, 0.151988636363636, 
0.52661717921527, 0.342271293375394, 0.304245283018868, 0.348445595854922, 
0.31236442516269, 0.317491166077739, 0.23972602739726, 0.276202219482121, 
0.278886610994081, 0.282127031019202, 0.306115483075756, 0.276422764227642, 
0.257203842049093, 0.17643391521197, 0.197475872308834, 0.379173838209983, 
0.177679882525698, 0.143392188336009, 0.0935143288084465, 0.200595829195631, 
0.438071995118975, 0.29933110367893, 0.211772232885477, 0.0256734859226251, 
0.266741071428571), `125um` = c(0.0100263852242744, 0.315196998123827, 
0.267373380447585, 0.284965034965035, 0.177556818181818, 0.0752916224814422, 
0.32807570977918, 0.295990566037736, 0.399611398963731, 0.134490238611714, 
0.0261736496718829, 0.237442922374429, 0.327990135635019, 0.240590007791825, 
0.370753323485968, 0.319370437091116, 0.485094850948509, 0.44076840981857, 
0.350997506234414, 0.528582034149963, 0.207194492254733, 0.281938325991189, 
0.173354735152488, 0.283810960281549, 0.201588877855015, 0.104535285743339, 
0.249163879598662, 0.19065898912348, 0.0339781243670245, 0.280133928571429
), `250um` = c(0.0058047493403694, 0.0956848030018762, 0.0459363957597173, 
0.0734265734265734, 0.00852272727272727, 0.00742311770943797, 
0.00788643533123028, 0.0129716981132075, 0.0518134715025907, 
0.00650759219088937, 0.0553003533568905, 0.00684931506849315, 
0.0135635018495684, 0.0384999321551487, 0.0782865583456425, 0.0337347112923107, 
0.0569105691056911, 0.0448239060832444, 0.180174563591022, 0.0972531551596139, 
0.0236660929432013, 0.0558002936857562, 0.0358480470840021, 0.0834590246354952, 
0.041708043694141, 0.0213544844417328, 0.0133779264214047, 0.0300703774792067, 
0.00212679765039498, 0.0223214285714286), `1000um` = c(0.00527704485488127, 
0, 0.0176678445229682, 0.034965034965035, 0.0184659090909091, 
0.00466595970307529, 0, 0, 0.00453367875647668, 0, 0.00850580514891469, 
0, 0.0123304562268804, 0.0170259043395175, 0.0472673559822747, 
0.00785057362283114, 0.0284552845528455, 0.0224119530416222, 
0.101620947630923, 0.0244988864142539, 0.00493975903614458, 0.273127753303965, 
0.490101658640984, 0.456510809451986, 0.30883813306852, 0.00935529794590197, 
0.0317725752508361, 0.353166986564299, 0.925663358314766, 0.143973214285714
), tap = c(0.181515414138913, 0.117573294968505, 0.118865354458575, 
0.10855175801997, 0.144215480661367, 0.174146595044186, 0.119127447894571, 
0.124407020872865, 0.105198658268468, 0.115848718019756, 0.178503108560497, 
0.12940849965294, 0.117179453375411, 0.111245191819955, 0.116370341491164, 
0.142787808192584, 0.113006923837784, 0.132290387475618, 0.0969106601341345, 
0.10220582768636, 0.215412621359223, 0.101006035443262, 0.090402518790167, 
0.0885715142546809, 0.0796504565039956, 0.135196812364163, 0.126232974506719, 
0.103923910552975, 0.092065868263473, 0.101525813669742), sd = c(0.00544534378975074, 
0.000976059208807937, 0.00105675176697187, 0.00129099190584701, 
0.010280862872346, 0.00245055200549895, 0.000348361831630225, 
0.00083859912379809, 0.001823914492673, 0.000495255850064288, 
0.010990421037094, 0.000600243017338078, 0.00277626938387639, 
0.000420021748488332, 0.0014744503375705, 0.0199049666612412, 
0.00034970661779751, 0.00314117917812511, 0.000215745910008022, 
0.000435189983246506, 0.127862755335917, 0.000119010118739825, 
0.00226462935582809, 0.00132600786385895, 0.00243540659124004, 
0.0118334943096431, 0.000449569977393725, 0.000285425329099763, 
0.00246993386642007, 0.000814994916882806), dypi = c(16.6, 16.6, 
16.6, 16.6, 16.6, 24.5, 24.5, 24.5, 24.5, 24.5, 16.6, 16.6, 16.6, 
16.6, 16.6, 34, 34, 34, 34, 34, 11.2, 11.2, 11.2, 11.2, 11.2, 
11, 11, 11, 11, 11), `Síld 2013` = c(0L, 0L, 0L, 0L, 0L, 1L, 
1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L)), row.names = c(NA, -30L), class = "data.frame")


# Save figure as PNG at 300 dpi, 18 cm wide
png("supplementary/fig3.png", width = 17.4, height = 7, units = "cm", res = 300)

#png("organic_matter_boxplots.png", width = 10, height = 5, units = "in", res = 300)

# Set up the plotting area
par(mfrow=c(1,2), mar = c(5, 4, 4, 2) + 0.1)

# Create the boxplots
boxplot(tap ~ Artal, data = kornast_organic_tap_dypi, main="Organic Matter by Year", 
        xlab="Year", ylab="Organic Matter Content",
        cex.axis = 0.8, cex.lab = 0.9, cex.main = 1)
boxplot(tap ~ stod, data = kornast_organic_tap_dypi, main="Organic Matter by Station", 
        xlab="Station", ylab="Organic Matter Content",
        cex.axis = 0.8, cex.lab = 0.9, cex.main = 1)

# Close the PNG device
#dev.off()

# Close device
dev.off()



























# Prepare data for diversity calculations
species_matrix <- KolgrTaxa %>%
  select(species, sample_id, density) %>%
  pivot_wider(
    names_from = species,
    values_from = density,
    values_fill = 0,
    values_fn = sum
  ) %>%
  column_to_rownames(var = "sample_id") %>%
  mutate(across(everything(), as.numeric))

# Calculate species richness and diversity indices
richness <- specnumber(species_matrix)
shannon <- diversity(species_matrix, index = "shannon")
simpson <- diversity(species_matrix, index = "simpson")

# Create diversity dataframe
diversity_df <- data.frame(
  sample_id = rownames(species_matrix),
  Richness = richness,
  Shannon = shannon,
  Simpson = simpson
) %>%
  left_join(
    KolgrTaxa %>%
      select(sample_id, year, station) %>%
      distinct(),
    by = "sample_id"
  )

# Create output directory if it doesn't exist
#if(!dir.exists("output")) {
#  dir.create("output")
#}

# Plot diversity indices
richness_plot <- diversity_df %>% 
  filter(!year == 1999 & station %in% c("St-3", "St-5", "St-4", "St-6", "St-2", "St-1")) %>% 
  ggplot(aes(x = year, y = Richness, color = station, group = station)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Species Richness Over Time by Station",
       x = "Year",
       y = "Species Richness")

shannon_plot <- diversity_df %>% 
  filter(!year == 1999 & station %in% c("C4", "A7", "B5", "B8", "E4", "E3")) %>% 
  ggplot(aes(x = year, y = Shannon, color = station, group = station)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Shannon Diversity Over Time by Station",
       x = "Year",
       y = "Shannon Diversity Index")

# Save plots
ggsave("output/richness_plot.png", richness_plot, width = 10, height = 6, bg = "white")
ggsave("output/shannon_plot.png", shannon_plot, width = 10, height = 6, bg = "white")

# Summary statistics
summary_stats <- diversity_df %>%
  group_by(year, station) %>%
  summarise(
    Mean_Richness = mean(Richness),
    SD_Richness = sd(Richness),
    Mean_Shannon = mean(Shannon),
    SD_Shannon = sd(Shannon),
    .groups = 'drop'
  )

# Save summary statistics
#write.csv(summary_stats, "output/diversity_summary_stats.csv", row.names = FALSE)

# Boxplots for richness and Shannon diversity by year
richness_boxplot <- ggplot(diversity_df, aes(x = as.factor(year), y = Richness)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Species Richness Distribution by Year",
       x = "Year",
       y = "Species Richness")

shannon_boxplot <- ggplot(diversity_df, aes(x = as.factor(year), y = Shannon)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Shannon Diversity Distribution by Year",
       x = "Year",
       y = "Shannon Diversity Index")

# Save boxplots
ggsave("output/richness_boxplot.png", richness_boxplot, width = 10, height = 6, bg = "white")
ggsave("output/shannon_boxplot.png", shannon_boxplot, width = 10, height = 6, bg = "white")

# Species accumulation curve
species_accum <- specaccum(species_matrix)

# Save species accumulation plot
png("output/species_accumulation.png", width = 800, height = 600)
plot(species_accum, ci.type = "poly", col = "blue", lwd = 2, ci.lty = 0, ci.col = "lightblue",
     xlab = "Number of Samples", ylab = "Species Richness",
     main = "Species Accumulation Curve")
dev.off()

# Top 10 most abundant species
top_species <- df %>%
  group_by(species) %>%
  summarise(total_abundance = sum(density)) %>%
  slice_max(order_by = total_abundance, n = 10)

# Bar plot of top 10 species
top_species_plot <- ggplot(top_species, aes(x = reorder(species, total_abundance), y = total_abundance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 Most Abundant Species",
       x = "Species",
       y = "Total Abundance")

# Save top species plot
ggsave("output/top_species_plot.png", top_species_plot, width = 10, height = 6, bg = "white")
