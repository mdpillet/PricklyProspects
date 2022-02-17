library(raster)
library(dplyr)
library(sp)

# Set directories
relPath <- "D:/Chapter1/"
predPath <- "Data/CHELSA/Raw/Current/"
outPath <- "Data/CHELSA/Processed/Current/"
occPath <- "Data/Occurrences/Raw/"
occFile <- "occurrences.csv"
tmpDir <- "Temp/"

# Set temporary directory
rasterOptions(tmpdir = paste0(relPath, tmpDir))

# Read in predictor data
predFiles <- list.files(paste0(relPath, predPath), "tif$")
preds <- stack(paste0(relPath, predPath, predFiles))
NAvalue(preds) <- -32768

# Aggregate predictors to 300 arc seconds (~10km)
preds <- aggregate(preds, 10)

# Read in occurrence data and make spatial
occurrences <- read.csv(paste0(relPath, occPath, occFile), header = T, stringsAsFactors = F)
occurrences <- subset(occurrences, !is.na(latitude) & !is.na(longitude))
occurrences <- occurrences %>% group_by(scrubbed_species_binomial) %>% mutate(speciesCount = n()) %>% filter(speciesCount >= 10)
occSp <- SpatialPointsDataFrame(data.frame(LON = occurrences$longitude, LAT = occurrences$latitude),
                                occurrences, # Remove latitude and longitude
                                proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
cropExtent <- floor(extent(occSp) + 10) # Add 5 degrees in each direction
# Crop predictors
preds <- crop(preds, cropExtent)

# Convert units
div10 <- function(x) x / 10
div1000 <- function(x) x / 1000
beginCluster()
for (i in c(1, 2, 5:11)) preds[[i]] <- clusterR(preds[[i]], calc, args = list(fun = div10))
for (i in c(3, 4)) preds[[i]] <- clusterR(preds[[i]], calc, args = list(fun = div1000))
endCluster()
preds <- stack(preds)

# Name layers
names(preds) <- paste0("BIO", 1:19)

# Project predictors to equal-area CRS (WGS 84/NSIDC EASE-Grid 2.0 Global)
predsProj <- projectRaster(preds, crs = CRS("+init=epsg:6933"))

# Export aggregated predictors
writeRaster(predsProj, paste0(relPath, outPath, "current.tif"), overwrite = T)