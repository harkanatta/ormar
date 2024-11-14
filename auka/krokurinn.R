
# Example data for two stations
station1 <- data.frame(
  species = c("Scoloplos armiger", "Eteone longa", "ostracoda", "harpacticoida", "nematoda", "krabbalirfa", "Capitella capitata", "Malacoceros fuliginosus", "amphipoda", "Spionidae sp."),
  count = c(190, 1, 1, 4, 2, 1, 1, 24, 2, 1)
) %>% filter(!species %in% c("ostracoda", "harpacticoida", "nematoda", "krabbalirfa","amphipoda"))

station2 <- data.frame(
  species = c("Scoloplos armiger", "Microphthalmus aberrans", "Malacoceros fuliginosus", "Eteone longa"),
  count = c(42, 2, 18, 1)
)

total_area <- 0.0225 * 3

# Add station labels
station1$station <- "st1"
station2$station <- "st2"

# Add total_area and adjusted count
station1$total_area <- total_area
station1$adjusted_count <- station1$count / total_area
station2$total_area <- total_area
station2$adjusted_count <- station2$count / total_area

# Combine both stations into a single long-format data frame
long_format_df <- rbind(station1, station2)

# Reorder columns for clarity
long_format_df <- long_format_df[, c("station", "species", "total_area", "count", "adjusted_count")]

# Display the long-format data frame
print(long_format_df)


library(dplyr)
library(tidyr)

# Custom function for calculating biotic indices, adjusted to your data
calculate_biotic_indices_krokurinn <- function(data) {
  # Check for required columns in the provided data
  required_cols <- c("station", "species", "adjusted_count")
  stopifnot(
    "Missing required columns" = all(required_cols %in% colnames(data)),
    "No data provided" = nrow(data) > 0,
    "Negative density values found" = all(data$adjusted_count >= 0),
    "Missing values found" = !any(is.na(data$adjusted_count))
  )
  
  # Initialize list to store results
  indices_list <- list()
  
  # Prepare data by grouping by station and species
  biotic_data <- data %>%
    group_by(station, species) %>%
    summarise(adjusted_density = sum(adjusted_count), .groups = 'drop') %>%
    pivot_wider(names_from = station, values_from = adjusted_density, values_fill = 0)
  
  # Placeholder for biotic index calculations (replace BBI with your actual formula)
  # For example purposes, we’ll assume a simple calculation.
  # Replace this with the actual BBI calculation logic
  bbi_results <- list(
    BBI = rowMeans(biotic_data[, -1]),  # Simple example: taking the mean adjusted density across stations
    BBIclass = ifelse(rowMeans(biotic_data[, -1]) > 10, "High", "Low"),  # Example classification based on arbitrary threshold
    nEQR = rowMeans(biotic_data[, -1]) / max(rowMeans(biotic_data[, -1]))  # Normalized example based on max value
  )
  
  # Store results
  indices_list <- list(
    indices = as.data.frame(cbind(biotic_data$species, BBI = bbi_results$BBI)),
    classification = as.data.frame(cbind(biotic_data$species, BBIclass = bbi_results$BBIclass)),
    normalized = as.data.frame(cbind(biotic_data$species, nEQR = bbi_results$nEQR))
  )
  
  return(indices_list)
}

# Example of running the function
biotic_indices <- calculate_biotic_indices_krokurinn(long_format_df)

# Display the result
print(biotic_indices)


