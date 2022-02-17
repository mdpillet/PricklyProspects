library(rgdal)

# Set directory structure
relPath <- "F:/Chapter1/"
occPath <- "Data/Occurrences/Subsets"
outPath <- "Data/Occurrences/Species"

# Read occurrence data
occurrences <- readOGR(paste0(relPath, occPath), "occurrences10")

# Split by species
species <- as.character(unique(occurrences$scrbbd_s_))
for (i in species) {
  print(i)
  tmp <- subset(occurrences, scrbbd_s_ == i)
  writeOGR(tmp, dsn = paste0(relPath, outPath), gsub(" ", "_", i), driver = "ESRI Shapefile", overwrite_layer = T)
}