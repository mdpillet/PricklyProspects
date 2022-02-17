library(raster)
library(rasterVis)

# Set directory structure
relPath <- "F:/Chapter1/"
diversityPath <- "Data/DiversityMaps/Raw/"
uncPath <- "Data/UncertaintyRichness/"
figPath <- "Figures/"

# List diversity maps
diversityMaps <- list.files(paste0(relPath, diversityPath), "tif$")

# Stack diversity maps and convert to data frame
richnessStack <- stack(paste0(relPath, diversityPath, diversityMaps))
aggrStack <- aggregate(richnessStack[[1]], 10)
for (i in 2:nlayers(richnessStack)) {
  print(i)
  aggrStack <- stack(aggrStack, aggregate(richnessStack[[i]], 10))
}
richnessVector <- getValues(aggrStack)
richnessVector <- t(richnessVector)
richnessVector <- as.data.frame(richnessVector)

# Add columns for model settings
richnessVector$SamplingDistance <- sapply(strsplit(row.names(richnessVector), "_"), "[", 1)
richnessVector$VariableSelection <- sapply(strsplit(row.names(richnessVector), "_"), "[", 2)
richnessVector$CorrelationFilter <- sapply(strsplit(row.names(richnessVector), "_"), "[", 3)
richnessVector$Time <- sapply(strsplit(row.names(richnessVector), "_"), "[", 4)
richnessVector$RCP <- sapply(strsplit(row.names(richnessVector), "_"), "[", 5)
richnessVector$GCM <- sapply(strsplit(row.names(richnessVector), "_"), "[", 6)
richnessVector$ProjectionDistance <- sapply(strsplit(row.names(richnessVector), "_"), "[", 7)
richnessVector$Model <- sapply(strsplit(row.names(richnessVector), "_"), "[", 8)
richnessVector$Threshold <- sapply(strsplit(row.names(richnessVector), "_"), "[", 9)

# Remove current maps
richnessVector <- subset(richnessVector, Time != "current")

# Run ANOVA
anovList <- list()
for (i in 1:ncell(aggrStack)) {
  skipVector <- unique(richnessVector[, i])
  skipVector <- skipVector[!is.nan(skipVector)]
  if (length(skipVector) > 1) {
    tmp <- richnessVector[!is.nan(richnessVector[, i]), ]
    if (length(unique(tmp$ProjectionDistance)) > 1) {
      print(i)
      anovList[i] <- summary(aov(tmp[, i] ~ Time + RCP + GCM + Threshold + VariableSelection + Model + ProjectionDistance + SamplingDistance + CorrelationFilter, data = tmp))    
    }
    else anovList[i] <- NA
  }
  else anovList[i] <- NA
}

# Export ANOVAs
saveRDS(anovList, paste0(relPath, uncPath, "ANOVA.rds"))
anovList <- readRDS(paste0(relPath, uncPath, "ANOVA.rds"))

# Calculate uncertainties
uncComponents <- data.frame(Driver = c("Time", "RCP", "GCM", "Threshold", "VariableSelection", "Model", "ProjectionDistance", "SamplingDistance", "CorrelationFilter"))
for (i in 1:length(anovList)) {
  print(i)
  if (is.na(anovList[i])) uncComponents[, i+1] <- NA
  else uncComponents[i+1] <- (anovList[[i]]$`Sum Sq` / sum(anovList[[i]]$`Sum Sq`) * 100)[1:9]
}

# Find driver with highest uncertainty explained for each cell
uncComponents$Driver <- as.character(uncComponents$Driver)
uncDriver <- character(length(anovList))
for (i in 2:ncol(uncComponents)) {
  if (is.na(unique(uncComponents[1, i]))) tmp <- NA
  else tmp <- uncComponents[uncComponents[, i] == max(uncComponents[, i]), 1]
  print(tmp)
  uncDriver[i - 1] <- tmp
}

# Fill raster
uncMap <- aggrStack[[1]]
uncMap <- setValues(uncMap, as.factor(uncDriver))
writeRaster(uncMap, paste0(relPath, uncPath, "uncertaintyMap.tif"))
uncMap <- raster(paste0(relPath, uncPath, "uncertaintyMap.tif"))

# Get top uncertainty source distribution
uncMap@data@attributes
table(values(uncMap)) / sum(!is.na(values(uncMap)))

# Plot uncertainty raster
plot(uncMap, col = rainbow(5))
legend("topright", legend = c("CorrelationFilter", "Model", "ProjectionDistance", "Threshold", "VariableSelection"),
       fill = rainbow(5))
uncPlot <- gplot(uncMap) + geom_tile(aes(fill = factor(value,
                                                       labels = c("CorrelationFilter", "Model", "ProjectionDistance", "Threshold", "VariableSelection")))) +
  coord_equal() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 42, face = "bold")) +
  scale_fill_discrete(name = "Driver", labels = c("CorrelationFilter", "Model", "ProjectionDistance", "Threshold", "VariableSelection"), na.value = NA)
ggsave(paste0(relPath, figPath, "topUncertaintySource.png"))