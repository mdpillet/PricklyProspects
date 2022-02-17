library(sp)
library(raster)
library(rgdal)
library(maxnet)
library(foreach)
library(doParallel)
library(caret)

# Set directory structure
relPath <- "F:/BIEN/"
occPath <- "Data/Occurrences/PerSpecies/MaxEntSpecies"
predRadius <- "500km"
predPath <- paste0("Data/CurrentPredComb/Cropped/", predRadius, "/")
modelPath <- "Data/ModelsCV/PCAraw/"
PCAPath <- "Data/PCA/"
outPath <- "OutFiles/"

# Read in predictors
pred_files <- list.files(path = paste0(relPath, predPath), pattern = "tif$", recursive = T)

# List occurrence files
occ_files <- list.files(path = paste0(relPath, occPath), pattern = "shp$")

# # Set up cluster
cl <- makeCluster(8, outfile = paste0(relPath, outPath, "buildMaxEntModels_PCAraw_CV.txt"))
registerDoParallel(cl)

# Set seed
seedNo <- 2020
set.seed(seedNo)

# Loop through species and build MaxEnt models
foreach(i = 1:length(pred_files), .inorder = F, .packages = c("sp", "raster", "rgdal", "maxnet", "caret")) %dopar% {
  # for (i in 1:length(pred_files)) {
  # Get taxon
  predFileName <- strsplit(pred_files[i], ".", fixed = T)[[1]][1]
  taxonName <- predFileName
  print(paste0("Modeling ", predFileName))
  # Loop through PCA iterations
  PCAmodels <- list.files(path = paste0(relPath, PCAPath), pattern = taxonName, recursive = T)
  PCAmodels <- PCAmodels[grepl(predRadius, PCAmodels)]
  for (j in PCAmodels) {
    if (!file.exists(paste0(relPath, modelPath, strsplit(j, split = ".", fixed = T)[[1]][1], "_regOver.rda")) |
        !file.exists(paste0(relPath, modelPath, strsplit(j, split = ".", fixed = T)[[1]][1], "_regUnder.rda"))) {
      print(paste0("Modeling ", j))
      # Load PCA model
      load(paste0(relPath, PCAPath, j))
      # Extract names of 6 variables with highest absolute loadings in PC1
      # If less than 6 variables available, all are used
      loadings <- PCAmodel$rotation[, 1]
      varNames <- names(sort(abs(loadings), decreasing = T, index.return = T)$x[1:6])
      # Read predictor file and subset variables based on PCA loadings
      preds <- stack(paste0(relPath, predPath, pred_files[i]))
      names(preds) <- c(paste0("BIO", 1:19), "pH")
      preds <- subset(preds, varNames)
      # Read occurrence data
      occ <- readOGR(dsn = paste0(relPath, occPath), layer = taxonName, verbose = F)
      # Extract presence data
      presData <- extract(preds, occ)
      # Subset unique records
      presData <- unique(presData)
      # Split data frame into 5 folds
      folds <- createFolds(1:nrow(presData), k = 5, list = TRUE, returnTrain = FALSE)
      if ((nrow(presData) / 5) < 2) folds <- createFolds(1:nrow(presData), k = 4, list = TRUE, returnTrain = FALSE)  
      if ((nrow(presData) / 4) < 2) folds <- createFolds(1:nrow(presData), k = 3, list = TRUE, returnTrain = FALSE)  
      if ((nrow(presData) / 3) < 2) folds <- createFolds(1:nrow(presData), k = 2, list = TRUE, returnTrain = FALSE) 
      # Extract background data
      bgData <- sampleRandom(preds, min(10000, ncell(preds)), na.rm = T)
      # Loop over folds and fit model
      if (!file.exists(paste0(relPath, modelPath, strsplit(j, split = ".", fixed = T)[[1]][1], "_regOver.rda"))) {
        model <- list()
        for (k in 1:length(folds)) {
          foldData <- presData[folds[[k]],]
          predData <- as.data.frame(rbind(foldData, bgData))
          rowsFoldData <- nrow(foldData)
          model[[k]] <- maxnet(c(rep(1, rowsFoldData), rep(0, nrow(bgData))),
                               predData, addsamplestobackground = T)
          predData$presAbs <- c(rep(1, rowsFoldData), rep(0, nrow(bgData)))
        }
        save(model, presData, folds, bgData, file = paste0(relPath, modelPath, strsplit(j, split = ".", fixed = T)[[1]][1], "_regOver.rda"))
      }
      if (!file.exists(paste0(relPath, modelPath, strsplit(j, split = ".", fixed = T)[[1]][1], "_regUnder.rda"))) {
        model <- list()
        for (k in 1:length(folds)) {
          foldData <- presData[folds[[k]],]
          predData <- as.data.frame(rbind(foldData, bgData))
          rowsFoldData <- nrow(foldData)
          model[[k]] <- maxnet(c(rep(1, rowsFoldData), rep(0, nrow(bgData))),
                               predData, maxnet.formula(c(rep(1, rowsFoldData), rep(0, nrow(bgData))),
                                                        predData, classes = "lq"),
                               addsamplestobackground = T)
          predData$presAbs <- c(rep(1, rowsFoldData), rep(0, nrow(bgData)))
        }
        save(model, presData, folds, bgData, file = paste0(relPath, modelPath, strsplit(j, split = ".", fixed = T)[[1]][1], "_regUnder.rda"))
      }
      # Remove temporary files
      rm(PCAmodel, preds, occ, presData, bgData, predData, model)
      removeTmpFiles()
    }
  }
}

# Close cluster connection
stopCluster(cl)