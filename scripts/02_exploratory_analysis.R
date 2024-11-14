# Load required packages
library(readr)
library(dplyr)
library(ggplot2)
library(vegan)
library(tidyr)
library(tibble)  # For column_to_rownames()

KolgrTaxa <- read_csv(file("data/raw/KolgrTaxa.csv", encoding = "UTF-8"), na = "empty")
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

# Filter out the specified taxa using if_all()
KolgrTaxa <- KolgrTaxa %>%
  filter(if_all(where(is.character), ~ !. %in% taxa_to_remove))

KolgrTaxa$gamalt <- iconv(KolgrTaxa$gamalt, from = "UTF-8", to = "UTF-8", sub = NA)
   
   remove_list <- paste(c(
     "nýsestir",
     "ungviði",
     "ungv",
     "ungv.",
     "juv",
     "harpacticoida"
   ), collapse = '|') 
   
   remove_ind <- lapply(strsplit(remove_list , "\\|")[[1]] , \(x) grep(x , KolgrTaxa$gamalt , fixed = T)) |> 
     unlist() |> 
     unique()
   
   ekkiungvidi <- KolgrTaxa[-remove_ind,] 

   jorundur <- ekkiungvidi %>% 
     filter(!Flokkun %in% c("harpacticoida", "Campanulariidae")) %>% 
     mutate(
       Flokkun = case_when(
         # General taxonomic standardization (no year-specific)
         str_detect(tolower(Flokkun), "^sipuncul") ~ NA_character_,
         tolower(Flokkun) == "ostracoda" ~ "Ostracoda",
         tolower(Flokkun) == "terebellides stroemi" ~ "Terebellides stroemii",
         tolower(Flokkun) == "bivalvia" ~ NA_character_,
         tolower(Flokkun) == "pectinaria koreni" ~ "Lagis koreni",
         tolower(Flokkun) == "sphaerosyllis" ~ "Sphaerosyllis erinaceus",
         tolower(Flokkun) %in% c("maldane sarsi", "maldanidae", "praxillella praetermissa") ~ "Praxillella",
         tolower(Flokkun) == "oligochaeta" ~ NA_character_,
         tolower(Flokkun) == "Tubificidae" ~ NA_character_,
         
         # 1999 specific matches
         Artal == "1999" & tolower(Flokkun) %in% c("ampharetinae", "ampharete acutifrons") ~ "Ampharetidae",
         Artal == "1999" & tolower(Flokkun) == "bivalvia" ~ NA_character_,
         Artal == "1999" & tolower(Flokkun) == "leucon acutirostris" ~ "Cumacea",
         Artal == "1999" & tolower(Flokkun) == "eudorella emarginata" ~ "Cumacea",
         Artal == "1999" & tolower(Flokkun) == "musculus discors" ~ "Musculus",
         Artal == "1999" & tolower(Flokkun) == "aricidea (acmira) cerrutii" ~ "Aricidea suecica",
         Artal == "1999" & tolower(Flokkun) == "mediomastus filiformis" ~ "Heteromastus filiformis",
         
         # 2013 specific matches
         Artal == "2013" & tolower(Flokkun) %in% c("harmothoe", "polynoidae") ~ "Harmothoe extenuata",
         Artal == "2013" & tolower(Flokkun) == "spionidae" ~ "Spio",
         
         # 2014 specific matches
         Artal == "2014" & tolower(Flokkun) == "ampharete acutifrons" ~ "Ampharetidae",
         Artal == "2014" & tolower(Flokkun) == "spionidae" ~ "Spio",
         Artal == "2014" & tolower(Flokkun) == "terebellidae" ~ "Terebellides stroemii",
         Artal == "2014" & tolower(Flokkun) == "caprellidae" ~ "Caprella septemtrionalis",
         Artal == "2014" & tolower(Flokkun) == "praxillella" ~ "Praxillella praetermissa",
         Artal == "2014" & tolower(Flokkun) == "syllidae" ~ "Syllis cornuta",
         Artal == "2014" & tolower(Flokkun) == "mya" ~ "Mya arenaria",
         Artal == "2014" & tolower(Flokkun) == "mytilus edulis" ~ "Mytilidae",
         Artal == "2014" & tolower(Flokkun) == "nephtyidae" ~ "Nephthys",
         Artal == "2014" & tolower(Flokkun) == "amphipoda" ~ "Protomedeia fasciata",
         Artal == "2014" & tolower(Flokkun) == "tubificidae" ~ "Tubificoides kozloffi",
         Artal == "2014" & tolower(Flokkun) == "harmothoe" ~ "Harmothoe extenuata",
         
         # 2015 specific matches
         Artal == "2015" & tolower(Flokkun) == "spionidae" ~ "Spio",
         Artal == "2015" & tolower(Flokkun) == "terebellidae" ~ "Terebellides stroemii",
         Artal == "2015" & tolower(Flokkun) == "crenella" ~ "Crenella decussata",
         Artal == "2015" & tolower(Flokkun) == "syllidae" ~ "Syllis cornuta",
         Artal == "2015" & tolower(Flokkun) == "mya" ~ "Mya arenaria",
         Artal == "2015" & tolower(Flokkun) == "mytilidae" ~ "Mytilus edulis",
         Artal == "2015" & tolower(Flokkun) == "ampharetidae" ~ NA_character_,
         Artal == "2015" & tolower(Flokkun) == "pholoe" ~ "Pholoe minuta",
         Artal == "2015" & tolower(Flokkun) == "harmothoe" ~ "Harmothoe extenuata",
         Artal == "2015" & tolower(Flokkun) == "amphipoda" ~ "Protomedeia fasciata",
         
         # 2016 specific matches
         Artal == "2016" & tolower(Flokkun) == "leaena ebranchiata" ~ "Terebellides stroemii",
         Artal == "2016" & tolower(Flokkun) == "pseudopolydora antennata" ~ "Polydora",
         Artal == "2016" & tolower(Flokkun) == "cistenides hyperborea" ~ "Cistenides granulata",
         Artal == "2016" & tolower(Flokkun) == "paraonidae" ~ "Aricidea suecica",
         Artal == "2016" & tolower(Flokkun) == "dorvilleidae" ~ "Capitella capitata",
         Artal == "2016" & tolower(Flokkun) == "orbiniidae" ~ "Scoloplos armiger",
         Artal == "2016" & tolower(Flokkun) == "polychaeta" ~ NA_character_,
         Artal == "2016" & tolower(Flokkun) == "polynoidae" ~ "Harmothoe imbricata",
         Artal == "2016" & tolower(Flokkun) %in% c("spio filicornis", "spionidae") ~ "Spionidae",
         Artal == "2017" & tolower(Flokkun) == "spionidae" ~ "Spio",
         Artal == "2016" & tolower(Flokkun) %in% c("syllis", "syllidae") ~ "Syllis cornuta",
         Artal == "2016" & tolower(Flokkun) %in% c("terebellidae", "terebelliformia") ~ "Terebellides stroemii",
         Artal == "2016" & tolower(Flokkun) == "amphipoda" ~ NA_character_,
         Artal == "2016" & tolower(Flokkun) == "bivalvia" ~ NA_character_,
         Artal == "2016" & tolower(Flokkun) == "cardium" ~ NA_character_,
         Artal == "2016" & tolower(Flokkun) == "aricidea" ~ "Amphitrite cirrata",
         Artal == "2016" & tolower(Flokkun) == "capitellidae" ~ "Capitella capitata",
         Artal == "2016" & tolower(Flokkun) == "cirratulidae" ~ "Cirratulus cirratus",
         Artal == "2016" & tolower(Flokkun) == "cossuridae" ~ "Cossura longocirrata",
         Artal == "2016" & tolower(Flokkun) == "pectinariidae" ~ "Pectinaria koreni",
         Artal == "2016" & tolower(Flokkun) == "phyllodocida" ~ "Phyllodoce maculata",
         Artal == "2016" & tolower(Flokkun) == "lumbrineridae" ~ "Lumbrineris",
         Artal == "2016" & tolower(Flokkun) == "pholoe" ~ "Pholoe minuta",
         Artal == "2016" & tolower(Flokkun) == "balanus balanus" ~ "Balanus",
         Artal == "2016" & tolower(Flokkun) %in% c("priapulidae", "priapulidae") ~ "Priapulus caudatus",
         
         # 2017 specific matches
         Artal == "2017" & tolower(Flokkun) == "oligochaeta" ~ NA_character_,
         Artal == "2017" & tolower(Flokkun) == "tubificidae" ~ NA_character_,
         Artal == "2017" & tolower(Flokkun) == "mediomastus filiformis" ~ "Heteromastus filiformis",
         Artal == "2017" & tolower(Flokkun) == "spionidae" ~ "Spio",
         Artal == "2017" & tolower(Flokkun) == "terebellidae" ~ "Terebellides stroemii",
         Artal == "2017" & tolower(Flokkun) == "caprellidae" ~ "Caprella septentrionalis",
         Artal == "2017" & tolower(Flokkun) == "opisthobranchia" ~ "Retusa pertenuis",
         Artal == "2017" & tolower(Flokkun) == "mya" ~ "Mya arenaria",
         Artal == "2017" & tolower(Flokkun) == "maldanidae" ~ "Praxillella praetermissa",
         Artal == "2017" & tolower(Flokkun) %in% c("ampharete", "ampharete acutifrons") ~ "Ampharetidae",
         Artal == "2017" & tolower(Flokkun) == "pholoe" ~ "Pholoe minuta",
         Artal == "2017" & tolower(Flokkun) == "harmothoe" ~ "Harmothoe extenuata",
         
         # Default case
         TRUE ~ Flokkun
       )
     ) %>% 
     drop_na(Flokkun)    

df <- jorundur  %>%
  rename(
    species = Flokkun,
    sample_id = id,
    year = Artal,
    station = stod,
    subdivision = skipting,
    count = N,
    density = Nu
  )

result <- df %>%
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
  target_stations <- c("C4", "A7", "B5", "B8", "E4", "E3")
  
  # Initialize lists to store results
  indices_list <- list()
  
  # Calculate indices for each year
  for (current_year in sort(unique(data$year))) {
    # Prepare data for BBI calculation using density values
    year_data <- data %>%
      filter(year == current_year,
             station %in% target_stations) %>%
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
stations <- rep(c("C4", "A7", "B5", "B8", "E4", "E3"), length(unique(df$year)))
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


ft <- flextable(final_table) %>%
  set_caption(caption = "Mean normalized benthic indices (± SD) across all stations showing temporal changes in Kolgrafafjörður fjord from 1999 to 2017.") %>%
  align(j = 1:4, align = "center") %>%
  align(j = 5, align = "left") %>%
  set_header_labels(
    Year = "Year",
    nAMBI = "nAMBI",
    nNQI1 = "nNQI1",
    nShannon = "nShannon",
    Status = "Ecological Status"
  ) %>%
  theme_box() %>%
  fontsize(size = 10) %>%
  padding(padding = 4) %>%
  border_outer(border = fp_border(width = 1.5)) %>%
  border_inner_h(border = fp_border(width = 0.5)) %>%
  autofit() %>%
  bold(part = "header") %>%
  bg(bg = "white", part = "all")

save_as_docx(
  ft,
  path = "output/benthic_indices_table.docx"
)



# Calculate the statistics
stats <- result %>%
    filter(!year == 1999 & station %in% c("C4", "A7", "B5", "B8", "E4", "E3")) %>% 
    group_by(station, year) %>%
    summarise(
        n_taxa = n_distinct(species),
        total_density = sum(adjusted_density),
        .groups = 'drop'
    )

min_taxa_info <- stats %>% 
    filter(n_taxa == min(n_taxa)) %>% 
    slice(1)  # In case there are multiple stations with the same minimum

max_taxa_info <- stats %>% 
    filter(n_taxa == max(n_taxa)) %>% 
    slice(1)  # In case there are multiple stations with the same maximum

min_density_info <- stats %>% 
    filter(total_density == min(total_density)) %>% 
    slice(1)

max_density_info <- stats %>% 
    filter(total_density == max(total_density)) %>% 
    slice(1)

# Print the results
print(min_taxa_info)
print(max_taxa_info)
print(min_density_info)
print(max_density_info)

# Prepare data for diversity calculations
species_matrix <- df %>%
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
    df %>%
      select(sample_id, year, station) %>%
      distinct(),
    by = "sample_id"
  )

# Create output directory if it doesn't exist
if(!dir.exists("output")) {
  dir.create("output")
}

# Plot diversity indices
richness_plot <- diversity_df %>% 
  filter(!year == 1999 & station %in% c("C4", "A7", "B5", "B8", "E4", "E3")) %>% 
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
write.csv(summary_stats, "output/diversity_summary_stats.csv", row.names = FALSE)

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
