library(dplyr)
library(raster)
library(doParallel)
library(foreach)

# Set directory structure
relPath <- "F:/Chapter1/"
statsPath <- "Data/ModelStats/modelStats.csv"
suitPath <- "Data/SuitabilityMaps/"
binaryPath <- "Data/BinaryMaps/"
outPath <- "Out/"
tempPath <- "Temp/"

# Set temporary directory
rasterOptions(tmpdir = paste0(relPath, tempPath))

# Read in model statistics
modStats <- read.csv(paste0(relPath, statsPath), header = T, stringsAsFactors = F)

# Extract taxon
modStats$Taxon <- unlist(lapply(lapply(strsplit(modStats$Model, "_"), "[", 1:2), paste, collapse = "_"))

# Subset taxa with AUC >= 0.5 and complete thresholds
subList <- unique(subset(modStats, AUC < 0.5 | !is.finite(MaxTSS) | !is.finite(Omiss5))$Taxon)
modStats <- filter(modStats, !(Taxon %in% subList))

# Get list of suitability maps
maps <- list.files(paste0(relPath, suitPath), "tif$")

# Get list of binary map
binaryMaps <- list.files(paste0(relPath, binaryPath), "tif$")

# Set up cluster
cl <- makeCluster(detectCores() - 1, outfile = paste0(relPath, outPath, "thresholdMaps.txt"))
registerDoParallel(cl)

# Loop over models
foreach(i = 1:nrow(modStats), .packages = c("raster")) %dopar% {
  print(i)
  modTmp <- strsplit(modStats[i, "Model"], ".", fixed = T)[[1]][1]
  mapTmp <- maps[grepl(modTmp, maps)]
  mapTmp <- mapTmp[grepl(paste0("model", toupper(modStats[i, "Features"]), ".tif"), mapTmp)]
  # Check if all binary maps have been created
  bMapTmp <- binaryMaps[grepl(modTmp, binaryMaps)]
  bMapTmp <- bMapTmp[grepl(paste0("model", toupper(modStats[i, "Features"]), "_"), bMapTmp)]
  if (length(bMapTmp) != 186) {
    for (j in 1:length(mapTmp)) {
      map <- raster(paste0(relPath, suitPath, mapTmp[j])) / 10000
      reclassMatrix <- matrix(data = c(0, modStats[i, "MaxTSS"], 0, modStats[i, "MaxTSS"], Inf, 1), nrow = 2, ncol = 3, byrow = T)  
      binaryMap <- reclassify(map, reclassMatrix, right = F)
      writeRaster(binaryMap, paste0(relPath, binaryPath, map@data@names, "_MaxTSS.tif"), format = "GTiff", datatype = "INT2S", options = "COMPRESS=LZW", overwrite = T)
      reclassMatrix <- matrix(data = c(0, modStats[i, "Omiss5"], 0, modStats[i, "Omiss5"], Inf, 1), nrow = 2, ncol = 3, byrow = T)  
      binaryMap <- reclassify(map, reclassMatrix, right = F)
      writeRaster(binaryMap, paste0(relPath, binaryPath, map@data@names, "_Omiss5.tif"), format = "GTiff", datatype = "INT2S", options = "COMPRESS=LZW", overwrite = T)  
    }
  }
  removeTmpFiles(h = 1)
}

# Close cluster connection
stopCluster(cl)