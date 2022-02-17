library(dplyr)
library(sp)
library(raster)
library(rgdal)

# Set directory structure
relPath <- "F:/Chapter1/"
occPath <- "Data/Occurrences/Raw/"
occFile <- "occurrences.csv"
predPath <- "Data/CHELSA/Processed/Current/current.tif"
outPath <- "Data/Occurrences/Subsets"

# Read in occurrence data (n = 42728)
occurrences <- read.csv(paste0(relPath, occPath, occFile), header = T, stringsAsFactors = F)

# Subset occurrences with longitude/latitude (n = 42728)
occurrences <- subset(occurrences, !is.na(latitude) & !is.na(longitude))

# Subset species based on number of occurrences (10 or more)
occurrences <- occurrences %>% group_by(scrubbed_species_binomial) %>% mutate(speciesCount = n()) %>% filter(speciesCount >= 10)



# Read in predictor data
preds <- stack(paste0(relPath, predPath))

# Extract predictor data
occExtr <- extract(preds, occSp, df = T)[, 2:20]
occurrences <- bind_cols(occurrences, occExtr)

# Loop over species and remove species with less than 10 unique climate observations
species <- unique(occurrences$scrubbed_species_binomial)
tmpAll <- occurrences
for (i in species) {
  tmp <- tmpAll %>% filter(scrubbed_species_binomial == i) %>% dplyr::select(current.1:current.19) %>% distinct()
  if (nrow(tmp) < 10) occurrences <- occurrences %>% filter(scrubbed_species_binomial != i)
}

# Make data spatial
occSp <- SpatialPointsDataFrame(data.frame(LON = occurrences$longitude, LAT = occurrences$latitude),
                                occurrences, # Remove latitude and longitude
                                proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# Project to equal-area CRS (WGS 84/NSIDC EASE-Grid 2.0 Global)
occSp <- spTransform(occSp, CRS("+init=epsg:6933"))

# Export occurrences
writeOGR(occSp, dsn = paste0(relPath, outPath), "occurrences10", driver = "ESRI Shapefile", overwrite_layer = T)
write.csv(occurrences, paste0(relPath, outPath, "/occurrences10.csv"), row.names = F)