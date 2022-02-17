library(raster)
library(doParallel)
library(foreach)

# Set directory structure
relPath <- "F:/Chapter1/"
binaryPath <- "Data/BinaryMaps/"
richnessPath <- "Data/DiversityMaps/Raw/"
extendMap <- "Data/CHELSA/Processed/Current/current.tif"
outPath <- "Out/"
tempPath <- "Temp/"

# Set temporary directory
rasterOptions(tmpdir = paste0(relPath, tempPath))

# Get list of binary map
binaryMaps <- list.files(paste0(relPath, binaryPath), "tif$")

# Get list of richness parameters
richnessPam <- unlist(lapply(lapply(strsplit(binaryMaps, "_"), "[", 3:11), paste, collapse = "_"))
richnessPam <- unique(richnessPam)

# Exclude species
subsetMaps <- binaryMaps[!grepl("Opuntia_ficus.indica|Hylocereus_undatus", binaryMaps)]

# Load map for extending/masking
extendMap <- stack(paste0(relPath, extendMap))

# Set up cluster
cl <- makeCluster(detectCores(), outfile = paste0(relPath, outPath, "createDiversityMaps.txt"), type = "FORK")
registerDoParallel(cl)

# Loop over maps
foreach(i = 1:length(richnessPam), .packages = c("raster")) %dopar% {
  if (!file.exists(paste0(relPath, richnessPath, richnessPam[i]))) {
    print(richnessPam[i])
    # Subset maps
    richnessMaps <- subsetMaps[grepl(richnessPam[i], subsetMaps)]
    # Read in rasters and create a list
    rasterList <- lapply(paste0(relPath, binaryPath, richnessMaps), raster)
    # Add function call arguments to list
    rasterList$fun <- sum
    rasterList$na.rm <- TRUE
    # Mosaic rasters
    richness <- do.call(mosaic, rasterList)
    # Extend and mask map
    richness <- extend(richness, extent(extendMap))
    richness <- mask(richness, extendMap[[1]], maskvalue = NA, updateValue = NA)
    # Export richness map
    writeRaster(richness, paste0(relPath, richnessPath, richnessPam[i]))
  }
  return(0)
}

# Close cluster connection
stopCluster(cl)