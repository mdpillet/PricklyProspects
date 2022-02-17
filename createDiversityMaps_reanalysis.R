### THIS REANALYSIS IS BASED ON PROBABILISTIC ESTIMATION OF RICHNESS

library(raster)
library(doParallel)
library(foreach)

# Set directory structure
relPath <- "F:/Chapter1/"
suitPath <- "Data/SuitabilityMaps/"
richnessPath <- "Data/DiversityMaps/RawReanalysis/"
extendMap <- "Data/CHELSA/Processed/Current/current.tif"
speciesPath <- "Data/EcologicalCovariates/covariates.csv"
outPath <- "Out/"
tempPath <- "Temp/"

# Set temporary directory
rasterOptions(tmpdir = paste0(relPath, tempPath))

# Get list of suitability maps
suitMaps <- list.files(paste0(relPath, suitPath), "tif$")

# Get list of richness parameters
richnessPam <- unlist(lapply(lapply(strsplit(suitMaps, "_"), "[", 3:10), paste, collapse = "_"))
richnessPam <- unique(richnessPam)

# Exclude species
speciesList <- read.csv(paste0(relPath, speciesPath), header = T, stringsAsFactors = F)$Taxon
speciesList <- paste(speciesList, collapse = "|")
subsetMaps <- suitMaps[grepl(speciesList, suitMaps)]
subsetMaps <- subsetMaps[!grepl("Opuntia_ficus.indica|Hylocereus_undatus", subsetMaps)]

# Load map for extending/masking
extendMap <- stack(paste0(relPath, extendMap))

# Set up cluster
# cl <- makeCluster(detectCores(), outfile = paste0(relPath, outPath, "createDiversityMaps.txt"), type = "FORK")
cl <- makeCluster(detectCores(), outfile = paste0(relPath, outPath, "createDiversityMaps.txt"))
registerDoParallel(cl)

# Loop over maps
# for (i in 1:length(richnessPam)) {
foreach(i = 1:length(richnessPam), .packages = c("raster")) %dopar% {
  if (!file.exists(paste0(relPath, richnessPath, richnessPam[i]))) {
    print(richnessPam[i])
    # Subset maps
    richnessMaps <- subsetMaps[grepl(richnessPam[i], subsetMaps)]
    # Read in rasters and create a list
    rasterList <- lapply(paste0(relPath, suitPath, richnessMaps), raster)
    # Add function call arguments to list
    rasterList$fun <- sum
    rasterList$na.rm <- TRUE
    # Mosaic rasters
    richness <- do.call(mosaic, rasterList)
    # Extend and mask map
    richness <- extend(richness, extent(extendMap))
    richness <- mask(richness, extendMap[[1]], maskvalue = NA, updateValue = NA)
    # Divide by storage constant
    richness <- richness / 10000
    # Export richness map
    writeRaster(richness, paste0(relPath, richnessPath, richnessPam[i]))
  }
  return(0)
}

# Close cluster connection
stopCluster(cl)