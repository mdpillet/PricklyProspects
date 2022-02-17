library(pROC)

# Set directory structure
relPath <- "F:/Chapter1/"
modelPath <- "Data/Models/"
predDataPath <- "Data/PredData/"
statsPath <- "Data/ModelStats/"

# Load models
models <- list.files(path = paste0(relPath, modelPath), "rda$")

# Create function to calculate thresholds and AUC
calcStats <- function(model, preds) {
  # Get model predictions
  prediction <- as.numeric(predict(model, preds, clamp = F, type = "cloglog"))
  # Calculate AUROC
  auroc <- roc(response = preds$Presence,
               predictor = prediction,
               smoothed = T, ci = T, quiet = T)
  auc <- as.numeric(auroc$auc)
  # Set up data frame with thresholds
  tmpThresholds <- data.frame(Threshold = auroc$thresholds, 
                              Specificity = auroc$specificities,
                              Sensitivity = auroc$sensitivities)
  # Calculate TSS at each threshold
  tmpThresholds$TSS <- tmpThresholds$Specificity + tmpThresholds$Sensitivity - 1
  # Find threshold with maximum TSS
  maxtss <- tmpThresholds[tmpThresholds$TSS == max(tmpThresholds$TSS), "Threshold"][1]
  # Get suitability predictions at presence locations
  binary <- subset(preds, Presence == 1)
  binary$Suitability <- predict(model, binary, clamp = F, type = "cloglog")
  # Calculate omission rate at each threshold
  for (j in 1:nrow(tmpThresholds)) {
    tmpThreshold <- tmpThresholds[j, "Threshold"]
    tmpThresholds[j, "OmissionRate5"] <- nrow(subset(binary, binary$Suitability < tmpThreshold)) / nrow(binary)
  }
  # Get first threshold with 5% omission rate
  index <- min(which(tmpThresholds[order(tmpThresholds$OmissionRate5),]$OmissionRate5 >= 0.05))
  omiss <- tmpThresholds[order(tmpThresholds$OmissionRate5),][index, "Threshold"]
  # Return statistics
  return(c(auc, maxtss, omiss))
}

# Set up data frames for thresholds
thresholds_l <- data.frame(Model = models, Features = "l", AUC = NA, MaxTSS = NA, Omiss5 = NA)
thresholds_lq <- data.frame(Model = models, Features = "lq", AUC = NA, MaxTSS = NA, Omiss5 = NA)
thresholds_lqh <- data.frame(Model = models, Features = "lqh", AUC = NA, MaxTSS = NA, Omiss5 = NA)

# Loop over models
for (i in 1:nrow(thresholds_l)) {
  print(i)
  # Load model and predictor data
  load(paste0(relPath, modelPath, models[i]))
  load(paste0(relPath, predDataPath, models[i]))
  # Get thresholds and AUC
  thresholds_l[i, 3:5] <- calcStats(model_l, predDataSubset)
  thresholds_lq[i, 3:5] <- calcStats(model_lq, predDataSubset)
  thresholds_lqh[i, 3:5] <- calcStats(model_lqh, predDataSubset)
}

# Combine thresholds
thresholds <- rbind(thresholds_l, thresholds_lq, thresholds_lqh)

# Export thresholds
write.csv(thresholds, paste0(relPath, statsPath, "modelStats.csv"), row.names = F)