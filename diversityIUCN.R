library(raster)
library(rgdal)
### THIS REANALYSIS IS BASED ON PROBABILISTIC ESTIMATION OF RICHNESS

library(spatialEco)

# Set directories
relPath <- "F:/Chapter1/"
rangePath <- "Data/IUCNRanges"
templatePath <- "Data/DiversityMaps/SummaryReanalysis/mean_current.tif"
divMapPath <- "Data/DiversityMaps/IUCN/"
exportPath <- "Data/DiversityMaps/IUCN/richness.tif"
extentPath <- "Data/CHELSA/Processed/Current/current.tif"

# Read in template
div <- raster(paste0(relPath, templatePath))
div <- setValues(div, 0)

# Loop through non-MaxEnt files
rangeFiles <- list.files(path = paste0(relPath, rangePath), pattern = glob2rx("*shp"), full.names = TRUE)
for (i in 1:length(rangeFiles)) {
  taxon <- strsplit(strsplit(rangeFiles[i], "IUCNRanges/")[[1]][2], ".shp")[[1]][1]
  print(taxon)
  tmp <- readOGR(dsn = paste0(relPath, rangePath), layer = taxon)
  extr <- extract(div, tmp, cellnumbers = T)
  extr <- do.call("rbind", extr)
  extr <- unique(extr)
  if (nrow(extr) == 0) next
  extr <- as.data.frame(extr)
  div[extr$cell] <- div[extr$cell] + 1
  print(i / length(rangeFiles)*100)
}

plot(div)

# Read total extent
extentMap <- stack(paste0(relPath, extentPath))

# Mask diversity map with a CHELSA layer
div <- mask(div, extentMap[[1]], maskvalue = NA, updateValue = NA)

# Export expert map
writeRaster(div, paste0(relPath, exportPath), overwrite = T)

# Import expert map
div <- raster(paste0(relPath, exportPath))

# Load MaxEnt map
maxent <- raster(paste0(relPath, templatePath))

# Calculate overall correlation
cor.test(values(maxent), values(div))

# Calculate spatial correlation
corRaster <- rasterCorrelation(maxent, div, s = 9)
writeRaster(corRaster, paste0(relPath, divMapPath, "correlation_windowSize9.tif"))

# Correlation between spatial correlation and expert map
cor.test(values(div), values(corRaster))

# Calculate correlation by quantile
richnessValues <- data.frame(MaxEnt = values(maxent), IUCN = values(div))
richnessValues <- richnessValues[complete.cases(richnessValues),]
tmp <- subset(richnessValues, MaxEnt <= 1)
cor.test(tmp$MaxEnt, tmp$IUCN)
tmp <- subset(richnessValues, MaxEnt <= 3 & MaxEnt > 1)
cor.test(tmp$MaxEnt, tmp$IUCN)
tmp <- subset(richnessValues, MaxEnt <= 7 & MaxEnt > 3)
cor.test(tmp$MaxEnt, tmp$IUCN)
tmp <- subset(richnessValues, MaxEnt > 7)
cor.test(tmp$MaxEnt, tmp$IUCN)
cor(richnessValues$MaxEnt, richnessValues$IUCN)

# Recalculate correlation at different scales
cor.test(values(aggregate(div, 5)), values(aggregate(maxent, 5)))
cor.test(values(aggregate(div, 10)), values(aggregate(maxent, 10)))