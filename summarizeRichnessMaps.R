library(raster)
library(rasterVis)
library(parallel)

# Set directory structure
relPath <- "F:/Chapter1/"
diversityPath <- "Data/DiversityMaps/Raw/"
outPath <- "Data/DiversityMaps/Summary/"

# List diversity maps
diversityMaps <- list.files(paste0(relPath, diversityPath), "tif$")

# Split maps by future and current
current <- diversityMaps[grepl("current", diversityMaps)]
cStack <- stack(paste0(relPath, diversityPath, current))
future50 <- diversityMaps[grepl("2041", diversityMaps)]
f50Stack <- stack(paste0(relPath, diversityPath, future50))
future70 <- diversityMaps[grepl("2061", diversityMaps)] 
f70Stack <- stack(paste0(relPath, diversityPath, future70))

# Calculate mean richness maps
beginCluster(detectCores() - 1)
cMean <- clusterR(cStack, calc, args = list(mean, na.rm = T))
writeRaster(cMean, paste0(relPath, outPath, "mean_current.tif"))
f50Mean <- clusterR(f50Stack, calc, args = list(mean, na.rm = T))
writeRaster(f50Mean, paste0(relPath, outPath, "mean_2041.2060.tif"))
f70Mean <- clusterR(f70Stack, calc, args = list(mean, na.rm = T))
writeRaster(f70Mean, paste0(relPath, outPath, "mean_2061.2080.tif"))
endCluster()

# Calculate SD richness maps
beginCluster(detectCores() - 1)
cSD <- clusterR(cStack, calc, args = list(sd, na.rm = T))
writeRaster(cSD, paste0(relPath, outPath, "sd_current.tif"))
f50SD <- clusterR(f50Stack, calc, args = list(sd, na.rm = T))
writeRaster(f50SD, paste0(relPath, outPath, "sd_2041.2060.tif"))
f70SD <- clusterR(f70Stack, calc, args = list(sd, na.rm = T))
writeRaster(f70SD, paste0(relPath, outPath, "sd_2061.2080.tif"))
endCluster()

# Calculate CV richness maps
cMean <- raster(paste0(relPath, outPath, "mean_current.tif"))
f50Mean <- raster(paste0(relPath, outPath, "mean_2041.2060.tif"))
f70Mean <- raster(paste0(relPath, outPath, "mean_2061.2080.tif"))
cSD <- raster(paste0(relPath, outPath, "sd_current.tif"))
f50SD <- raster(paste0(relPath, outPath, "sd_2041.2060.tif"))
f70SD <- raster(paste0(relPath, outPath, "sd_2061.2080.tif"))
cCV <- cSD / cMean
f50CV <- f50SD / f50Mean
f70CV <- f70SD / f70Mean
writeRaster(cCV, paste0(relPath, outPath, "cv_current.tif"))
writeRaster(f50CV, paste0(relPath, outPath, "cv_2041.2060.tif"))
writeRaster(f70CV, paste0(relPath, outPath, "cv_2061.2080.tif"))

# Split maps by time and RCP
future50rcp26 <- diversityMaps[grepl("2041", diversityMaps) & grepl("rcp26", diversityMaps)]
future50rcp45 <- diversityMaps[grepl("2041", diversityMaps) & grepl("rcp45", diversityMaps)]
future50rcp85 <- diversityMaps[grepl("2041", diversityMaps) & grepl("rcp85", diversityMaps)]
future70rcp26 <- diversityMaps[grepl("2061", diversityMaps) & grepl("rcp26", diversityMaps)]
future70rcp45 <- diversityMaps[grepl("2061", diversityMaps) & grepl("rcp45", diversityMaps)]
future70rcp85 <- diversityMaps[grepl("2061", diversityMaps) & grepl("rcp85", diversityMaps)]
f50r26Stack <- stack(paste0(relPath, diversityPath, future50rcp26))
f50r45Stack <- stack(paste0(relPath, diversityPath, future50rcp45))
f50r85Stack <- stack(paste0(relPath, diversityPath, future50rcp85))
f70r26Stack <- stack(paste0(relPath, diversityPath, future70rcp26))
f70r45Stack <- stack(paste0(relPath, diversityPath, future70rcp45))
f70r85Stack <- stack(paste0(relPath, diversityPath, future70rcp85))

# Calculate mean richness maps
beginCluster(detectCores() - 1)
f50r26Mean <- clusterR(f50r26Stack, calc, args = list(mean, na.rm = T))
writeRaster(f50r26Mean, paste0(relPath, outPath, "mean_2041.2060_rcp26.tif"))
f50r45Mean <- clusterR(f50r45Stack, calc, args = list(mean, na.rm = T))
writeRaster(f50r45Mean, paste0(relPath, outPath, "mean_2041.2060_rcp45.tif"))
f50r85Mean <- clusterR(f50r85Stack, calc, args = list(mean, na.rm = T))
writeRaster(f50r85Mean, paste0(relPath, outPath, "mean_2041.2060_rcp85.tif"))
f70r26Mean <- clusterR(f70r26Stack, calc, args = list(mean, na.rm = T))
writeRaster(f70r26Mean, paste0(relPath, outPath, "mean_2061.2080_rcp26.tif"))
f70r45Mean <- clusterR(f70r45Stack, calc, args = list(mean, na.rm = T))
writeRaster(f70r45Mean, paste0(relPath, outPath, "mean_2061.2080_rcp45.tif"))
f70r85Mean <- clusterR(f70r85Stack, calc, args = list(mean, na.rm = T))
writeRaster(f70r85Mean, paste0(relPath, outPath, "mean_2061.2080_rcp85.tif"))
endCluster()