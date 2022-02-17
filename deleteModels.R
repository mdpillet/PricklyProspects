# Set directory structure
relPath <- "F:/Chapter1/"
modelPath <- "Data/Models/"
occPath <- "Data/Occurrences/Species/"

# List model files
models <- list.files(paste0(relPath, modelPath), "rda$")

# Get species list
species <- list.files(paste0(relPath, occPath), "shp$")
species <- sapply(strsplit(species, ".", fixed = T), "[", 1)

# Find species with missing models
for (i in species) {
  tmp <- models[grepl(i, models)]
  if (length(tmp) < 12 & length(tmp) != 0) {
    print(i)
    file.remove(paste0(relPath, modelPath, tmp))
  }
}