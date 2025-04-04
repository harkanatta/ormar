# Data Cleaning and Preparation

# Load required packages
library(dplyr)
library(tidyr)
library(stringr)

# Load the combined data
allartalningar <- readRDS("data/allartalningar.rds")
df <- allartalningar

# Capitalize the first letter of each word in the 'Flokkun' column
df <- df %>%
  mutate(Flokkun = str_to_sentence(Flokkun))

# Create a 'gamalt' column with the original 'Flokkun' values
df <- df %>%
  mutate(gamalt = Flokkun)


# Standardize species names
df <- df %>% mutate(Flokkun = case_when(
  str_detect(Flokkun, "Campanulariidae sp") ~ "Campanulariidae",
  str_detect(Flokkun, "Cycloterus lumpus") ~ "Cyclopterus lumpus",
  str_detect(Flokkun, "Götungar") ~ "Foraminifera",
  str_detect(Flokkun, "Möttuldýr?" ) ~ "Tunicata",
  str_detect(Flokkun, "Nephtiydae" ) ~ "Nephtyidae",
  str_detect(Flokkun, "Phyllodoce maculata" ) ~ "Phyllodoce maculata",
  str_detect(Flokkun, "Polynoida" ) ~ "Polynoidae",
  str_detect(Flokkun, "Sternapsis scutata") ~ "Sternaspis scutata",
  str_detect(Flokkun, "Terribellides stroemi") ~ "Terebellides stroemii",
  str_detect(Flokkun, "Tubificoides benedict") ~ "Tubificoides benedii",
  str_detect(Flokkun, "Clinocardium cillaturn" ) ~ "Ciliatocardium ciliatum ciliatum",
  str_detect(Flokkun, "Gattyana cirrosa" ) ~ "Gattyana cirrhosa",
  str_detect(Flokkun, "Möttuldýr" ) ~ "Tunicata",
  str_detect(Flokkun, "Nemertea=nemertina" ) ~ "Nemertea",
  str_detect(Flokkun, "Nudibranch"  ) ~ "Nudibranchia",
  str_detect(Flokkun, "Priapulus candatus" ) ~ "Priapulus caudatus",
  str_detect(Flokkun, "Priapulus camelus" ) ~ "Priapulus caudatus",
  str_detect(Flokkun, "Sipunculidea/Sipunculidae") ~ "Sipunculidae",
  str_detect(Flokkun, "Terribellides kozloffi" ) ~ "Tubificoides kozloffi",
  str_detect(Flokkun, "Tubicoides kozloffi" ) ~ "Tubificoides kozloffi",
  str_detect(Flokkun, "Tubificoides benedict" ) ~ "Tubificoides benedii",
  str_detect(Flokkun, "Tubificoides benedi" ) ~ "Tubificoides benedii",
  str_detect(Flokkun, "Terebellides benedi" ) ~ "Tubificoides benedii",
  str_detect(Flokkun, "Astartidae borealis" ) ~ "Astarte borealis",
  str_detect(Flokkun, "Cerastoderma ovale") ~ "Parvicardium pinnulatum",
  str_detect(Flokkun, "Exogone verrugera") ~ "Exogone verugera",
  str_detect(Flokkun, "Nuculana tenuis" ) ~ "Ennucula tenuis",
  str_detect(Flokkun, "Opistobranchia") ~ "Opisthobranchia",
  str_detect(Flokkun, "Serripes groenlandica") ~ "Serripes groenlandicus",
  str_detect(Flokkun, "Cardidae"   ) ~ "Cardiidae",
  str_detect(Flokkun, "Henricia sanguinolenta") ~ "Henricia sanguinolenta",
  str_detect(Flokkun, "Macoma calcaria" ) ~ "Macoma calcarea",
  str_detect(Flokkun, "Ophryotrocha cf Cosmetandra") ~ "Ophryotrocha cosmetandra",
  str_detect(Flokkun, "Praxillella m"  ) ~ "Praxillella",
  str_detect(Flokkun, "Ranaormur"  ) ~ "Nemertea",
  str_detect(Flokkun, "Terebellidae"  ) ~ "Terebellidae",
  str_detect(Flokkun, "Terebellides benedi"  ) ~ "Terebellides", 
  str_detect(Flokkun, "Tubicoides benedi" ) ~ "Tubificoides benedii",
  str_detect(Flokkun, "Cardidae" ) ~ "Cardiidae", 
  str_detect(Flokkun, "Corophium bonellii" ) ~ "Crassicorophium bonellii",
  str_detect(Flokkun, "Corophium bonnellii" ) ~ "Crassicorophium bonellii",
  str_detect(Flokkun, "Plergonium spinosum" ) ~ "Pleurogonium spinosissimum",
  str_detect(Flokkun, "spio" ) ~ "Spio",
  str_detect(Flokkun, "Skeljar" ) ~ "Bivalvia",
  TRUE ~ Flokkun
))




# Clean the 'Flokkun' column by removing specific patterns
df <- df %>%
  mutate(Flokkun = str_remove_all(Flokkun, 
                                  " kemur líka Heteromastus filiformis| Sipunculidea/|\\.| TUNICATA EÐA FLEIRI?| nýsestir| ungviði| ungv\\.| ungv| juv| sp\\.| sp| ath"
  ))


# Remove unwanted species
remove_list <- paste(c("Campanulariidae sp", "—Ssekkt", "—ßekkt", "Foraminifera", "Götungar", "Harpacticoida", "Hydrozoa", "Lirfur", "lirfur", "lirfa", "Skordýr", "Ranaormur", "Lirfurogdrasl", "Bandormur", "Lirfa", "egg", "Óþekkt", "Nematoda", "Nemertea", "Möttuldýr?", "Porifera", "Plerugoniumnosissum", "Terebellides benedi", "Plergoniumnosum", "egg/lirfur", "Ostracoda", "Copepoda", "Collembola", "Cyclopterus lumpus"), collapse = '|')

remove_ind <- lapply(strsplit(remove_list, "\\|")[[1]], \(x) grep(x, df$Flokkun, fixed = TRUE)) |> 
  unlist() |> 
  unique()
df <- df[-remove_ind,]

# Remove specific patterns from 'gamalt' column
remove_list <- paste(c("nýsestir", "ungviði", "ungv", "ungv.", "juv", "harpacticoida"), collapse = '|') 

remove_ind <- lapply(strsplit(remove_list, "\\|")[[1]], \(x) grep(x, df$gamalt, fixed = TRUE)) |> 
  unlist() |> 
  unique()
df <- df[-remove_ind,]

# Rename columns for clarity
allartalningar <- df %>%
  rename(
    species = Flokkun,
    sample_id = id,
    year = Artal,
    station = stod,
    subdivision = skipting,
    count = N,
    density = Nu
  )

# Save the cleaned data
saveRDS(allartalningar, "data/cleaned_data.rds")
