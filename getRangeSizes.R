library(raster)
library(doParallel)
library(foreach)

# Set directory structure
relPath <- "/xdisk/benquist/mich/"
binaryPath <- "Data/BinaryMaps/BinaryMaps/"
rangeSizePath <- "Data/RangeSizes/BySpecies/"
outPath <- "Out/"
tempPath <- "Temp/"

# Set temporary directory
rasterOptions(tmpdir = paste0(relPath, tempPath))

# Get list of binary map
binaryMaps <- list.files(paste0(relPath, binaryPath), "tif$")

# Get list of taxa
taxa <- unlist(lapply(lapply(strsplit(binaryMaps, "_"), "[", 1:2), paste, collapse = "_"))
taxa <- unique(taxa)

# Set up cluster
cl <- makeCluster(detectCores(), outfile = paste0(relPath, outPath, "getRangeSizes.txt"), type = "FORK")
registerDoParallel(cl)

# Calculate range sizes
foreach(i = 1:length(taxa), .packages = c("raster")) %dopar% {
  print(taxa[i])
  # Subset maps
  mapsTmp <- binaryMaps[grepl(taxa[i], binaryMaps)]
  # Create range size data frame
  rangeSizes <- data.frame(Map = mapsTmp, Count0 = NA, Count1 = NA, CountNA = NA)
  if (!file.exists(paste0(relPath, rangeSizePath, taxa[i], ".csv"))) {
    for (j in 1:nrow(rangeSizes)) {
      # Read raster
      tmpRaster <- raster(paste0(relPath, binaryPath, rangeSizes[j, "Map"]))
      # Get frequency table
      tmpFreq <- as.data.frame(freq(tmpRaster, useNA = "always"))
      # Set sizes
      count0 <- tmpFreq[!is.na(tmpFreq$value) & tmpFreq$value == 0, "count"]
      if (length(count0) == 0) rangeSizes[j, "Count0"] <- 0
      else rangeSizes[j, "Count0"] <- count0
      count1 <- tmpFreq[!is.na(tmpFreq$value) & tmpFreq$value == 1, "count"]
      if (length(count1) == 0) rangeSizes[j, "Count1"] <- 0
      else rangeSizes[j, "Count1"] <- count1
      countNA <- tmpFreq[is.na(tmpFreq$value), "count"]
      if (length(countNA) == 0) rangeSizes[j, "countNA"] <- 0
      else rangeSizes[j, "CountNA"] <- countNA
    }
    write.csv(rangeSizes, paste0(relPath, rangeSizePath, taxa[i], ".csv"), row.names = F)
  }
  removeTmpFiles(h = 1)
  return(0)
} 

# Close cluster connection
stopCluster(cl)