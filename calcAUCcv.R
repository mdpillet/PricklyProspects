library(pROC)

# Set directory structure
relPath <- "F:/Chapter1/"
modelPath <- "Data/ModelsCV/"
predDataPath <- "Data/PredDataCV/"
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
  # Return AUC
  return(auc)
}

# Set up data frames for thresholds
auc_l <- data.frame(Model = models, Features = "l", AUC = NA)
auc_lq <- data.frame(Model = models, Features = "lq", AUC = NA)
auc_lqh <- data.frame(Model = models, Features = "lqh", AUC = NA)

# Loop over models
for (i in 1:nrow(auc_l)) {
  print(i)
  # Load model and predictor data
  load(paste0(relPath, modelPath, models[i]))
  load(paste0(relPath, predDataPath, models[i]))
  # Get thresholds and AUC
  auc_l[i, "AUC"] <- calcStats(model_l, predDataSubset)
  auc_lq[i, "AUC"] <- calcStats(model_lq, predDataSubset)
  auc_lqh[i, "AUC"] <- calcStats(model_lqh, predDataSubset)
}

# Combine thresholds
aucs <- rbind(auc_l, auc_lq, auc_lqh)

# Export thresholds
write.csv(aucs, paste0(relPath, statsPath, "modelStatsCV.csv"), row.names = F)