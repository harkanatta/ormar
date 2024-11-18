# Load necessary libraries
library(dplyr)
library(tidyr)
library(BBI)  # Ensure the BBI package is installed


# Example data for two stations
station1 <- data.frame(
  species = c("Scoloplos armiger", "Eteone longa", "ostracoda", "harpacticoida", "nematoda", "krabbalirfa", "Capitella capitata", "Malacoceros fuliginosus", "amphipoda", "Spionidae sp."),
  count = c(190, 1, 1, 4, 2, 1, 1, 24, 2, 1)
) %>% filter(!species %in% c("ostracoda", "harpacticoida", "nematoda", "krabbalirfa"))

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

# Prepare your data: species as rows, stations as columns
biotic_data <- long_format_df %>%
  select(station, species, adjusted_count) %>%
  pivot_wider(names_from = station, values_from = adjusted_count, values_fill = 0)

# Run BBI to calculate indices
bbi_results <- BBI(biotic_data)


# Calculate nEQR using the BBI results
eqr_results <- nEQR(bbi_results$BBI)


AMBI <- as.data.frame(bbi_results$BBI)["AMBI"]
NQI1 <- as.data.frame(bbi_results$BBI)["NQI1"]
NQI1_Class <- as.data.frame(bbi_results$BBIclass)["NQI1"]

# Búa til töflu með stöðvarheitum, NQI1 og NQI1 Class
filtered_data <- data.frame(
  station = rownames(as.data.frame(bbi_results$BBI)),
  AMBI = AMBI$AMBI,
  NQI1 = NQI1$NQI1,
  NQI1_Class = NQI1_Class$NQI1
)
print(filtered_data)

## Fyrir báðar stöðvar saman
krokurinn_data <- long_format_df %>%
  mutate(station = "Krokurinn") %>%
  group_by(station, species) %>%
  summarise(
    total_area = sum(total_area),
    count = sum(count),
    adjusted_count = sum(adjusted_count),
    .groups = 'drop'
  )

library(dplyr)

# Original data for "Krokurinn"
krokurinn_data <- long_format_df %>%
  mutate(station = "Krokurinn") %>%
  group_by(station, species) %>%
  summarise(
    total_area = 0.0225 * 6,
    count = sum(count),
    adjusted_count = sum(adjusted_count),
    .groups = 'drop'
  )

# Duplicate the data frame and change the station name
krokurinn2_data <- krokurinn_data %>%
  mutate(station = "Krokurinn2")

# Combine both data frames
combined_data <- bind_rows(krokurinn_data, krokurinn2_data)


# Prepare your data: species as rows, stations as columns
biotic_data <- combined_data %>%
  select(station, species, adjusted_count) %>%
  pivot_wider(names_from = station, values_from = adjusted_count, values_fill = 0)

# Run BBI to calculate indices
bbi_results <- BBI(biotic_data)
AMBI <- as.data.frame(bbi_results$BBI)["AMBI"]
NQI1 <- as.data.frame(bbi_results$BBI)["NQI1"]
NQI1_Class <- as.data.frame(bbi_results$BBIclass)["NQI1"]

# Búa til töflu með stöðvarheitum, NQI1 og NQI1 Class
filtered_data <- data.frame(
  station = rownames(as.data.frame(bbi_results$BBI)),
  AMBI = AMBI$AMBI,
  NQI1 = NQI1$NQI1,
  NQI1_Class = NQI1_Class$NQI1
)

# Prenta töfluna
print(filtered_data)
