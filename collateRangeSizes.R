relPath <- "F:/Chapter1/"
rangeSizePath <- "Data/RangeSizes/BySpecies/"
outPath <- "Data/RangeSizes/All/rangeSizes.csv"

# List range size files
rangeSizes <- list.files(paste0(relPath, rangeSizePath), "csv$", full.names = T)

# Collate range sizes
out <- do.call(rbind, lapply(rangeSizes, read.csv))

# Export range sizes
write.csv(out, paste0(relPath, outPath), row.names = F)