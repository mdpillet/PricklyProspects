library(raster)
library(doParallel)
library(foreach)

# Set directory structure
relPath <- "/xdisk/benquist/mich/"
binaryPath <- "Data/BinaryMaps/BinaryMaps/"
centroidPath <- "Data/Centroids/BySpecies/"
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
cl <- makeCluster(detectCores(), outfile = paste0(relPath, outPath, "getCentroids.txt"), type = "FORK")
registerDoParallel(cl)

# Calculate centroids
foreach(i = 1:length(taxa), .packages = c("raster")) %dopar% {
  print(taxa[i])
  # Subset maps
  mapsTmp <- binaryMaps[grepl(taxa[i], binaryMaps)]
  # Create centroid data frame
  centroids <- data.frame(Map = mapsTmp, Lon = NA, Lat = NA)
  if (!file.exists(paste0(relPath, centroidPath, taxa[i], ".csv"))) {
    for (j in 1:nrow(centroids)) {
      # Read raster
      tmpRaster <- raster(paste0(relPath, binaryPath, centroids[j, "Map"]))
      # Calculate centroid
      tmpCoord <- xyFromCell(tmpRaster, which(tmpRaster[] == 1))
      centroid <- colMeans(tmpCoord)
      centroids[j, 2:3] <- centroid
    }
    write.csv(centroids, paste0(relPath, centroidPath, taxa[i], ".csv"), row.names = F)
  }
  removeTmpFiles(h = 1)
  return(0)
}

# Close cluster connection
stopCluster(cl)