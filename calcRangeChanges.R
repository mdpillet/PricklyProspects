relPath <- "F:/Chapter1/"
rangeSizePath <- "Data/RangeSizes/All/rangeSizes.csv"
rangeChangePath <- "Data/RangeChanges/rangeChanges.csv"

# Read range sizes
rangeSizes <- read.csv(paste0(relPath, rangeSizePath), header = T, stringsAsFactors = F)

# Parse out parameters
rangeSizesSplit <- strsplit(rangeSizes[, "Map"], "_")
rangeSizes$Taxon <- unlist(lapply(lapply(rangeSizesSplit, "[", 1:2), paste, collapse = "_"))
rangeSizes$SamplingDistance <- sapply(rangeSizesSplit, "[", 3)
rangeSizes$VariableSelection <- sapply(rangeSizesSplit, "[", 4)
rangeSizes$CorrelationFilter <- sapply(rangeSizesSplit, "[", 5)
rangeSizes$Time <- sapply(rangeSizesSplit, "[", 6)
rangeSizes$RCP <- sapply(rangeSizesSplit, "[", 7)
rangeSizes$GCM <- sapply(rangeSizesSplit, "[", 8)
rangeSizes$ProjectionDistance <- sapply(rangeSizesSplit, "[", 9)
rangeSizes$Model <- sapply(rangeSizesSplit, "[", 10)
rangeSizes$Threshold <- sapply(rangeSizesSplit, "[", 11)

# Calculate total size
rangeSizes$CountTotal <- rangeSizes$Count0 + rangeSizes$Count1

# Loop over data frame and add current range size
currentMaps <- subset(rangeSizes, Time == "current")
for (i in 1:nrow(rangeSizes)) {
  print(i)
  if (rangeSizes[i, "Time"] == "current") {
    rangeSizes[i, "CurrentCount1"] <- NA
    rangeSizes[i, "CurrentCountTotal"] <- NA
  }
  else {
    tmp <- subset(currentMaps, Taxon == rangeSizes[i, "Taxon"] & 
                    SamplingDistance == rangeSizes[i, "SamplingDistance"] &
                    VariableSelection == rangeSizes[i, "VariableSelection"] &
                    CorrelationFilter == rangeSizes[i, "CorrelationFilter"] &
                    Time == "current" & 
                    ProjectionDistance == rangeSizes[i, "ProjectionDistance"] &
                    Model == rangeSizes[i, "Model"] &
                    Threshold == rangeSizes[i, "Threshold"])
    rangeSizes[i, "CurrentCount1"] <- tmp[1, "Count1"]
    rangeSizes[i, "CurrentCountTotal"] <- tmp[1, "CountTotal"]
  }
}

# Calculate range changes (will result in NA if CurrentCount1 is 0)
rangeSizes$RangeChange <- (rangeSizes$Count1 / rangeSizes$CountTotal) / (rangeSizes$CurrentCount1 / rangeSizes$CurrentCountTotal)

# Export range changes
write.csv(rangeSizes, paste0(relPath, rangeChangePath), row.names = F)