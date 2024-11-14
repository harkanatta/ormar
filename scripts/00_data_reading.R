# Data Reading

# Load required packages
library(dplyr)
library(tidyr)
library(readr)
library(here)

# Function to load 2016 data
load_data <- function(dataPath, classes) { 
  datafiles <- list.files(path=dataPath, pattern = "csv", full.names = TRUE, recursive = TRUE)
  datafiles <- datafiles[-19] # Exclude specific file if necessary
  tables <- lapply(datafiles, read.csv, colClasses=classes, na.strings=c("NA", ""))
  names(tables) <- substr(dirname(datafiles), nchar(dirname(datafiles)) - 2, nchar(dirname(datafiles)))
  data.table::rbindlist(tables, idcol = "id")
}

# Read data for each year
t2013 <- read.csv(here("Kolgr2013", list.files(here("Kolgr2013"), pattern = "csv")), encoding="latin1")[ ,-c(1,2)] %>%
  pivot_longer(!Flokkun, names_to = "dolla", values_to = "N") %>% 
  na.omit(N) %>% 
  mutate(Artal = 2013)

t2014 <- read.csv(here("Kolgr2014", list.files(here("Kolgr2014"))), encoding="latin1") %>%
  pivot_longer(!Flokkun, names_to = "dolla", values_to = "N") %>% 
  na.omit(N) %>% 
  mutate(Artal = 2014)

t2015 <- read.csv(here("Kolgr2015", list.files(here("Kolgr2015"), pattern = "csv")), encoding="latin1") %>%
  pivot_longer(!Flokkun, names_to = "dolla", values_to = "N") %>% 
  na.omit(N) %>% 
  mutate(Artal = 2015)

t2017 <- read.csv(here("Kolgr2017", list.files(here("Kolgr2017"), pattern = "csv")), encoding="latin1") %>%
  pivot_longer(!Flokkun, names_to = "dolla", values_to = "N") %>% 
  na.omit(N) %>% 
  mutate(Artal = 2017)

GVH <- rbind(t2013, t2014, t2015, t2017) %>% 
  mutate(skipting = substr(dolla, nchar(dolla), nchar(dolla))) %>% 
  mutate(skipting = case_when(grepl("[[:alpha:]]", skipting) ~ 1,
                              TRUE ~ as.numeric(as.character(skipting))),
         dolla = substr(dolla, 1, 3),
         stod = substr(dolla, 1, 2),
         Nu = N * skipting) %>%
  rename(id = dolla)

data2016 <- load_data(here("Kolgr2016"), c("character", "integer", "factor", "character"))
data2016$stod <- substr(data2016$id, 1, 2)
data2016$skipting <- gsub("¼|1/4|0.25", 4, data2016$skipting) %>% as.numeric()

gogn <- data2016 %>%
  as.data.frame() %>%
  group_by(Flokkun, id, N, skipting, stod) %>%
  summarize(Artal = 2016, Nu = sum(N * skipting)) %>%
  ungroup() %>%
  select(Flokkun, id, N, Artal, skipting, stod, Nu)

agnar <- read_csv("./Kolgr1999/taflaAgnarStytt.csv") %>% 
  filter(rowSums(.[, 2:7], na.rm = TRUE) != 0) %>% 
  replace(is.na(.), 0) %>% 
  pivot_longer(cols = !Flokkun, names_to = 'stod', values_to = 'N') %>% 
  filter(N != 0) %>% 
  mutate(Artal = 1999, skipting = 1, id = stod, Nu = N) %>% 
  select(Flokkun, id, N, Artal, skipting, stod, Nu)

allartalningar <- rbind(GVH, gogn, agnar)

# Save the combined data
# saveRDS(allartalningar, "data/allartalningar.rds")
