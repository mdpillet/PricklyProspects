library(sp)
library(raster)
library(rgdal)
library(maxnet)
library(foreach)
library(doParallel)

# Set directory structure
relPath <- "F:/Chapter1/"
modelPath <- "Data/Models/"
PCAPath <- "Data/PCA/"
occPath <- "Data/Occurrences/Species/"
currPath <- "Data/CHELSA/Processed/Current/current.tif"
futPath <- "Data/CHELSA/Processed/Future/"
bufferPath <- "Data/Buffers"
mapPath <- "Data/SuitabilityMaps/"
outPath <- "Out/"
tempPath <- "Temp/"

# Set temporary directory
rasterOptions(tmpdir = paste0(relPath, tempPath))

# Get model list
models <- list.files(paste0(relPath, modelPath), "rda$")

# Get species list
species <- unique(sapply(lapply(strsplit(models, "_"), "[", 1:2), paste, collapse = "_"))

# Get predictor list
preds <- c(currPath, paste0(futPath, list.files(paste0(relPath, futPath), "tif$")))

# Set up cluster
cl <- makeCluster(detectCores() - 1, outfile = paste0(relPath, outPath, "predictSuitability.txt"))
registerDoParallel(cl)

# Predict suitability
foreach(i = 1:length(species), .inorder = F, .packages = c("sp", "raster", "rgdal", "maxnet")) %dopar% {
  mapList <- list.files(paste0(relPath, mapPath), species[i])
  if (length(mapList) != 3348) {
    print(paste0("Processing ", species[i]))
    # Subset models
    modelsTmp <- models[grepl(species[i], models)]
    # Load buffer
    buffer <- readOGR(paste0(relPath, bufferPath), paste0(species[i], "_500km"), verbose = F)
    for (j in preds) {
      # Set predictor name for export
      predName <- strsplit(j, "\\.|/")[[1]][5]
      if (predName == "current") predName <- "current_NA_NA"
      if (length(mapList[grepl(predName, mapList)]) != 108) {
        # Load predictors
        predsTmp <- stack(paste0(relPath, j))
        names(predsTmp) <- paste0("BIO", 1:19)
        # Crop and mask predictors
        predsTmp <- crop(predsTmp, buffer)
        predsTmp <- raster::mask(predsTmp, buffer)
        for (k in modelsTmp) {
          # Load model
          load(paste0(relPath, modelPath, k))
          # Predict suitability
          varSel <- strsplit(k, "_", fixed = T)[[1]][4]
          predsTmp_df <- as.data.frame(getValues(predsTmp))
          if (varSel == "varPCAtrans") {
            # Load PCA model
            load(paste0(relPath, PCAPath, gsub(paste0("_", varSel), "", k, fixed = T)))
            if (strsplit(strsplit(k, "_", fixed = T)[[1]][5], ".", fixed = T)[[1]][1] == "corrF") PCAmodel <- PCAmodel_corrF
            else PCAmodel <- PCAmodel_corrT
            # Transform predictors
            predsTmp_df <- as.data.frame(predict(PCAmodel, predsTmp_df))
          }
          # Get predictions
          predsTmp_df[, c("suit_l", "suit_lq", "suit_lqh")] <- NA
          predictCases_l <- complete.cases(predsTmp_df[, names(model_l$levels)])
          predictCases_lq <- complete.cases(predsTmp_df[, names(model_lq$levels)])
          predictCases_lqh <- complete.cases(predsTmp_df[, names(model_lqh$levels)])
          predsTmp_df[predictCases_l, "suit_l"] <- predict(model_l, predsTmp_df[predictCases_l, ], clamp = F, type = "cloglog")
          predsTmp_df[predictCases_lq, "suit_lq"] <- predict(model_lq, predsTmp_df[predictCases_lq,], clamp = F, type = "cloglog")
          predsTmp_df[predictCases_lqh, "suit_lqh"] <- predict(model_lqh, predsTmp_df[predictCases_lqh,], clamp = F, type = "cloglog")
          # Convert to integer
          predsTmp_df[, c("suit_l", "suit_lq", "suit_lqh")] <- round(predsTmp_df[, c("suit_l", "suit_lq", "suit_lqh")] * 10000)
          # Export predictions to raster
          suit_l <- setValues(predsTmp[[1]], predsTmp_df[, "suit_l"])
          suit_lq <- setValues(predsTmp[[1]], predsTmp_df[, "suit_lq"])
          suit_lqh <- setValues(predsTmp[[1]], predsTmp_df[, "suit_lqh"])
          # Crop to different buffers and export
          for (l in c("500", "100", "0")) {
            if (l != "500") {
              bufferTmp <- readOGR(paste0(relPath, bufferPath), paste0(species[i], "_", l, "km"), verbose = F)
              suit_l <- crop(suit_l, bufferTmp)
              suit_l <- raster::mask(suit_l, bufferTmp)
              suit_lq <- crop(suit_lq, bufferTmp)
              suit_lq <- raster::mask(suit_lq, bufferTmp)
              suit_lqh <- crop(suit_lqh, bufferTmp)
              suit_lqh <- raster::mask(suit_lqh, bufferTmp)
            }
            # Export suitability layers
            expName <- strsplit(k, ".", fixed = T)[[1]][1]
            expName <- paste(expName, predName, paste0("d", l), sep = "_")
            writeRaster(suit_l, paste0(relPath, mapPath, expName, "_modelL.tif"), format = "GTiff", datatype = "INT2S", options = "COMPRESS=LZW", overwrite = T)
            writeRaster(suit_lq, paste0(relPath, mapPath, expName, "_modelLQ.tif"), format = "GTiff", datatype = "INT2S", options = "COMPRESS=LZW", overwrite = T)
            writeRaster(suit_lqh, paste0(relPath, mapPath, expName, "_modelLQH.tif"), format = "GTiff", datatype = "INT2S", options = "COMPRESS=LZW", overwrite = T)
          }
        }
      }  
    }
  }
  removeTmpFiles(h = 3)
}

# Close cluster connection
stopCluster(cl)