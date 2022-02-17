library(raster)
library(sp)
library(rgdal)
library(foreach)
library(doParallel)

# Set directories
relPath <- "F:/Chapter1/"
predPath <- "Data/CHELSA/Processed/Current/current.tif"
bufferPath <- "Data/Buffers"
PCAPath <- "Data/PCA/"
samplePath <- "Data/PredSamples/"
outPath <- "Out/runPCA.txt"
tmpPath <- "Temp/"

# List buffers
obs_files <- list.files(path = paste0(relPath, bufferPath), pattern = "shp$")
obs_files <- obs_files[grepl("_100km|_500km", obs_files)]

# Read in predictor stack
preds <- stack(paste0(relPath, predPath))

# Set up cluster
cl <- makeCluster(detectCores() - 1, outfile = paste0(relPath, outPath))
registerDoParallel(cl)

# Set seed
set.seed(2020)

# Loop through species
sampleSize <- 1000
foreach(i = 1:length(obs_files), .packages = c("raster", "sp", "rgdal"), .inorder = T) %dopar% {
  # Get taxon name
  taxon <- strsplit(obs_files[i], ".", fixed = T)[[1]][1]
  if (file.exists(paste0(relPath, PCAPath, taxon, "_corrF.rda")) &
      file.exists(paste0(relPath, PCAPath, taxon, "_corrT.rda"))) return(0)
  # Load buffer
  buffer <- readOGR(paste0(relPath, bufferPath), taxon)
  # Sample predictors
  print(paste0("Sampling ", taxon))
  tmp <- crop(preds, buffer)
  tmp <- raster::mask(tmp, buffer)
  predSample <- as.data.frame(sampleRandom(tmp, min(ncell(tmp), sampleSize), na.rm = T, cells = T))
  names(predSample) <- c("CellNumber", paste0("BIO", 1:19))
  write.csv(predSample, paste0(relPath, samplePath, taxon, ".csv"))
  # Remove zero-variance variables
  predSampleProc <- predSample[, which(apply(predSample[, 2:20], 2, var) != 0) + 1]
  # Run PCA without correlation analysis
  PCAmodel_corrF <- prcomp(predSampleProc, center = T, scale. = T)
  # Iteratively remove variables with highest number of correlations over threshold, until no correlations over threshold remain
  threshold <- 0.7
  exclName <- "init"
  tmpSample <- predSampleProc
  while (!is.na(exclName)) {
    tmpCorr <- cor(tmpSample, use = "complete.obs")
    exclName <- names(sort(colSums(abs(tmpCorr) > threshold)[(colSums(abs(tmpCorr) > threshold)) > 1], decreasing = T)[1])
    if (!is.na(exclName)) tmpSample <- tmpSample[, colnames(tmpSample) != exclName]
  }
  lowCorrVars <- names(tmpSample)
  # Run PCA with correlation analysis
  PCAmodel_corrT <- prcomp(predSampleProc[, lowCorrVars], center = T, scale. = T)
  # Export PCA results
  save(PCAmodel_corrF, file = paste0(relPath, PCAPath, taxon, "_corrF.rda"))
  save(PCAmodel_corrT, file = paste0(relPath, PCAPath, taxon, "_corrT.rda"))
  # Free memory
  removeTmpFiles(h = 0)
}

# Close cluster connection
stopCluster(cl)