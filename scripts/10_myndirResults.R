# nokkrar myndir með gögnum frá Jörundi

#p <- henda_long_cleaned %>%
library(ggplot2)
library(dplyr)

# Define the desired order of years
year_levels <- c("1999", "2013", "2014", "2015", "2016", "2017")

# Process the data and convert 'year' to a factor with specified levels
df_processed <- df %>% 
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


ggsave("number_of_species_across_stations_2013-2017.png", plot = p, width = 8, height = 6, units = "in", dpi = 600, bg="white")


##########
#n$hannon#
##########

# Gögn koma úr 02_explor...
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
p <- ggplot(normalized_eqr, aes(x = Year_Factor, y = nEQR.nShannon, color = Station, group = Station)) +
  
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
    aes(x = Year_Factor, y = nEQR.nShannon, color = Station, group = Station),
    linetype = "dotted", size = 1
  ) +
  
  # Add points for all years
  geom_point(aes(x = Year_Factor), size = 2) +
  
  theme_minimal(base_size = 12) +
  labs(
    title = "nShannon Trends Across Stations (1999-2017)",
    x = "",
    y = "nShannon",
    color = "Station"
  ) +
  theme(
    axis.text.x = element_text(hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  scale_x_discrete(limits = year_levels)

# Display the plot
print(p)
ggsave("nshannon_trends_across_stations_2013_2017.png", plot = p, width = 8, height = 6, units = "in", dpi = 600, bg="white")
