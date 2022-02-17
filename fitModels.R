library(sp)
library(rgdal)
library(maxnet)
library(paran)
library(foreach)
library(doParallel)

# Set directory structure
relPath <- "F:/Chapter1/"
occPath <- "Data/Occurrences/Species"
modelPath <- "Data/Models/"
predSamplePath <- "Data/PredSamples/"
predDataPath <- "Data/PredData/"
PCAPath <- "Data/PCA/"
outPath <- "Out/"

# List occurrence files
occ_files <- list.files(path = paste0(relPath, occPath), pattern = "shp$")

# Set up cluster
cl <- makeCluster(detectCores() - 1, outfile = paste0(relPath, outPath, "fitModels.txt"))
registerDoParallel(cl)

# Loop through species and build MaxEnt models
foreach(i = 1:length(occ_files), .inorder = F, .packages = c("sp", "rgdal", "maxnet", "paran")) %dopar% {
  # Set seed
  seedNo <- 2020
  set.seed(seedNo)
  # Get taxon
  taxon <- strsplit(occ_files[i], ".", fixed = T)[[1]][1]
  print(paste0("Modeling ", taxon))
  # Read occurrence data
  occ <- readOGR(dsn = paste0(relPath, occPath), layer = taxon, verbose = F)
  # Extract presence data
  presData <- as.data.frame(occ[,42:60])[,1:19]
  colnames(presData) <- paste0("BIO", 1:19)
  presData <- presData[complete.cases(presData),]
  # Loop over 100km and 500km predictor samples
  for (j in c("100km", "500km")) {
    # Read background sample
    bgData <- read.csv(paste0(relPath, predSamplePath, taxon, "_", j, ".csv"), header = T, stringsAsFactors = F)[,3:21]
    # Combine presence and background data
    predData <- as.data.frame(rbind(presData, bgData))
    
    if (!file.exists(paste0(relPath, modelPath, taxon, "_", j, "_varBIEN_corrNA.rda"))) {
      # Subset BIEN variables
      predDataSubset <- predData[, c("BIO1", "BIO2", "BIO12", "BIO15", "BIO18", "BIO19")]
      # Create models
      model_l <- maxnet(c(rep(1, nrow(presData)), rep(0, nrow(bgData))),
                        predDataSubset,
                        maxnet.formula(c(rep(1, nrow(presData)), rep(0, nrow(bgData))), predDataSubset, classes = "l"),
                        addsamplestobackground = T, maxit = 100000)
      model_lq <- maxnet(c(rep(1, nrow(presData)), rep(0, nrow(bgData))),
                        predDataSubset,
                        maxnet.formula(c(rep(1, nrow(presData)), rep(0, nrow(bgData))), predDataSubset, classes = "lq"),
                        addsamplestobackground = T, maxit = 100000)
      model_lqh <- maxnet(c(rep(1, nrow(presData)), rep(0, nrow(bgData))),
                        predDataSubset,
                        maxnet.formula(c(rep(1, nrow(presData)), rep(0, nrow(bgData))), predDataSubset, classes = "lqh"),
                        addsamplestobackground = T, maxit = 100000)
      # Save models
      save(model_l, model_lq, model_lqh, file = paste0(relPath, modelPath, taxon, "_", j, "_varBIEN_corrNA.rda"))
      # Save predictor data
      predDataSubset$Presence <- c(rep(1, nrow(presData)), rep(0, nrow(bgData)))
      save(predDataSubset, file = paste0(relPath, predDataPath, taxon, "_", j, "_varBIEN_corrNA.rda"))
    }
    
    if (!file.exists(paste0(relPath, modelPath, taxon, "_", j, "_varRandom_corrNA.rda"))) {
      # Subset random variables
      predDataSubset <- predData[, sample(names(predData), size = 6, replace = F)]
      # Create models
      model_l <- maxnet(c(rep(1, nrow(presData)), rep(0, nrow(bgData))),
                        predDataSubset,
                        maxnet.formula(c(rep(1, nrow(presData)), rep(0, nrow(bgData))), predDataSubset, classes = "l"),
                        addsamplestobackground = T, maxit = 100000)
      model_lq <- maxnet(c(rep(1, nrow(presData)), rep(0, nrow(bgData))),
                         predDataSubset,
                         maxnet.formula(c(rep(1, nrow(presData)), rep(0, nrow(bgData))), predDataSubset, classes = "lq"),
                         addsamplestobackground = T, maxit = 100000)
      model_lqh <- maxnet(c(rep(1, nrow(presData)), rep(0, nrow(bgData))),
                          predDataSubset,
                          maxnet.formula(c(rep(1, nrow(presData)), rep(0, nrow(bgData))), predDataSubset, classes = "lqh"),
                          addsamplestobackground = T, maxit = 100000)
      # Save models
      save(model_l, model_lq, model_lqh, file = paste0(relPath, modelPath, taxon, "_", j, "_varRandom_corrNA.rda"))
      # Save predictor data
      predDataSubset$Presence <- c(rep(1, nrow(presData)), rep(0, nrow(bgData)))
      save(predDataSubset, file = paste0(relPath, predDataPath, taxon, "_", j, "_varRandom_corrNA.rda"))
    }
    
    # Fit PCA models
    for (k in c("corrF", "corrT")) {
      if (!file.exists(paste0(relPath, modelPath, taxon, "_", j, "_varPCAraw_", k, ".rda"))) {
        # Load PCA
        load(paste0(relPath, PCAPath, taxon, "_", j, "_", k, ".rda"))
        if (k == "corrF") PCAmodel <- PCAmodel_corrF
        else PCAmodel <- PCAmodel_corrT
        # Extract names of 6 variables with highest absolute loadings in PC1
        # If less than 6 variables available, all are used
        loadings <- PCAmodel$rotation[, 1]
        varNames <- names(sort(abs(loadings), decreasing = T, index.return = T)$x[1:min(6, length(loadings))])
        # Subset variables
        predDataSubset <- predData[, varNames]
        # Create models
        model_l <- maxnet(c(rep(1, nrow(presData)), rep(0, nrow(bgData))),
                          predDataSubset,
                          maxnet.formula(c(rep(1, nrow(presData)), rep(0, nrow(bgData))), predDataSubset, classes = "l"),
                          addsamplestobackground = T, maxit = 100000)
        model_lq <- maxnet(c(rep(1, nrow(presData)), rep(0, nrow(bgData))),
                           predDataSubset,
                           maxnet.formula(c(rep(1, nrow(presData)), rep(0, nrow(bgData))), predDataSubset, classes = "lq"),
                           addsamplestobackground = T, maxit = 100000)
        model_lqh <- maxnet(c(rep(1, nrow(presData)), rep(0, nrow(bgData))),
                            predDataSubset,
                            maxnet.formula(c(rep(1, nrow(presData)), rep(0, nrow(bgData))), predDataSubset, classes = "lqh"),
                            addsamplestobackground = T, maxit = 100000)
        # Save models
        save(model_l, model_lq, model_lqh, file = paste0(relPath, modelPath, taxon, "_", j, "_varPCAraw_", k, ".rda"))
        # Save predictor data
        predDataSubset$Presence <- c(rep(1, nrow(presData)), rep(0, nrow(bgData)))
        save(predDataSubset, file = paste0(relPath, predDataPath, taxon, "_", j, "_varPCAraw_", k, ".rda"))
      }
      
      if (!file.exists(paste0(relPath, modelPath, taxon, "_", j, "_varPCAtrans_", k, ".rda"))) {
        # Perform parallel analysis
        if (k == "corrF") parAnalysis <- paran(bgData, centile = 95, seed = seedNo, iterations = 1000, quietly = T, status = F)$Retained
        else parAnalysis <- paran(bgData[, row.names(PCAmodel$rotation)], centile = 95, seed = seedNo, iterations = 1000, quietly = T, status = F)$Retained
        # Transform predictor data
        predDataSubset <- as.data.frame(predict(PCAmodel, predData)[, 1:max(parAnalysis, 2)])
        # Create models
        model_l <- maxnet(c(rep(1, nrow(presData)), rep(0, nrow(bgData))),
                          predDataSubset,
                          maxnet.formula(c(rep(1, nrow(presData)), rep(0, nrow(bgData))), predDataSubset, classes = "l"),
                          addsamplestobackground = T, maxit = 100000)
        model_lq <- maxnet(c(rep(1, nrow(presData)), rep(0, nrow(bgData))),
                           predDataSubset,
                           maxnet.formula(c(rep(1, nrow(presData)), rep(0, nrow(bgData))), predDataSubset, classes = "lq"),
                           addsamplestobackground = T, maxit = 100000)
        model_lqh <- maxnet(c(rep(1, nrow(presData)), rep(0, nrow(bgData))),
                            predDataSubset,
                            maxnet.formula(c(rep(1, nrow(presData)), rep(0, nrow(bgData))), predDataSubset, classes = "lqh"),
                            addsamplestobackground = T, maxit = 100000)
        # Save models
        save(model_l, model_lq, model_lqh, file = paste0(relPath, modelPath, taxon, "_", j, "_varPCAtrans_", k, ".rda"))
        # Save predictor data
        predDataSubset$Presence <- c(rep(1, nrow(presData)), rep(0, nrow(bgData)))
        save(predDataSubset, file = paste0(relPath, predDataPath, taxon, "_", j, "_varPCAtrans_", k, ".rda"))
      }
    }
  }
}

# Close cluster connection
stopCluster(cl)