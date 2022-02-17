library(raster)
library(dplyr)
library(sp)
library(foreach)
library(doParallel)

# Set directories
relPath <- "F:/Chapter1/"
predPath <- "Data/CHELSA/Raw/Future/"
outPath <- "Data/CHELSA/Processed/Future/"
occPath <- "Data/Occurrences/Raw/"
occFile <- "occurrences.csv"
tmpDir <- "Temp/"
outFile <- "Out/prepFutureCHELSAData.txt"

# Set temporary directory
rasterOptions(tmpdir = paste0(relPath, tmpDir))

# Read in predictor data
predFiles <- list.files(paste0(relPath, predPath), "tif$")

# Create function to process predictors
processPreds <- function(x, y) {
  # Aggregate predictors to 300 arc seconds (~10km)
  print(paste0("Aggregating...", y))
  x <- aggregate(x, 10)
  
  # Read in occurrence data and make spatial
  occurrences <- read.csv(paste0(relPath, occPath, occFile), header = T, stringsAsFactors = F)
  occurrences <- subset(occurrences, !is.na(latitude) & !is.na(longitude))
  occurrences <- occurrences %>% group_by(scrubbed_species_binomial) %>% mutate(speciesCount = n()) %>% filter(speciesCount >= 10)
  occSp <- SpatialPointsDataFrame(data.frame(LON = occurrences$longitude, LAT = occurrences$latitude),
                                  occurrences, # Remove latitude and longitude
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  cropExtent <- floor(extent(occSp) + 10) # Add 5 degrees in each direction
  # Crop predictors
  print(paste0("Cropping...", y))
  x <- crop(x, cropExtent)
  
  # Convert units
  print(paste0("Unit conversion...", y))
  div10 <- function(x) x / 10
  div1000 <- function(x) x / 1000
  beginCluster()
  for (i in c(1, 2, 5:11)) x[[i]] <- clusterR(x[[i]], calc, args = list(fun = div10))
  for (i in c(3, 4)) x[[i]] <- clusterR(x[[i]], calc, args = list(fun = div1000))
  endCluster()
  x <- stack(x)
  
  # Name layers
  names(x) <- paste0("BIO", 1:19)
  
  # Project predictors to equal-area CRS (WGS 84/NSIDC EASE-Grid 2.0 Global)
  print(paste0("Projecting...", y))
  predsProj <- projectRaster(x, crs = CRS("+init=epsg:6933"))
  
  # Export aggregated predictors
  writeRaster(predsProj, paste0(relPath, outPath, y, ".tif"), overwrite = T)
}

# Set up cluster
cl <- makeCluster(detectCores() - 1, outfile = paste0(relPath, outFile))
registerDoParallel(cl)

# Loop over year/GCM/RCP combinations
years <- c("2041-2060", "2061-2080")
RCPs <- c("rcp26", "rcp45", "rcp85")
GCMs <- c("NorESM1-M", "GFDL-ESM2M", "HadGEM2-AO", "CCSM4", "BNU-ESM")
combs <- expand.grid(years, RCPs, GCMs, stringsAsFactors = F)
foreach(i = 1:nrow(combs), .inorder = F, .packages = c("raster", "dplyr", "sp")) %dopar% {
  # Subset climate layers
  tmp <- predFiles[grepl(combs[i, 1], predFiles) & grepl(combs[i, 2], predFiles) & grepl(combs[i, 3], predFiles)]
  # Sort and stack climate layers
  preds <- stack(paste0(relPath, predPath, tmp[c(1, 12:19, 2:11)]))
  NAvalue(preds) <- -32768
  # Set projection
  crs(preds) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
  # Process future climate stack
  if (!file.exists(paste0(relPath, outPath, paste(as.vector(combs[i, ]), collapse = "_"), ".tif"))) processPreds(preds, paste(as.vector(combs[i, ]), collapse = "_"))
  return(0)
}

# Close cluster connection
stopCluster(cl)