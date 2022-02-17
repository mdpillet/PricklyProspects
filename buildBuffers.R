library(raster)
library(rgeos)
library(sp)
library(rgdal)
library(foreach)
library(doParallel)

# Set directories
relPath <- "F:/Chapter1/"
occPath <- "Data/Occurrences/Species"
outPath <- "Out/buildBuffers.txt"
bufferPath <- "Data/Buffers"
tempPath <- "Temp/"

rasterOptions(tmpdir = paste0(relPath, tempPath))

# List observation files
obs_files <- list.files(path = paste0(relPath, occPath), pattern = "shp$")

# Set up cluster
cl <- makeCluster(detectCores() - 1, outfile = paste0(relPath, outPath))
registerDoParallel(cl)

# Compute convex hull and construct buffers
foreach(i = 1:length(obs_files), .packages = c("sp", "raster", "rgdal", "rgeos")) %dopar% {
  print(obs_files[i])
  taxon <- strsplit(obs_files[i], ".", fixed = T)[[1]][1]
  tmp <- readOGR(dsn = paste0(relPath, occPath), layer = taxon, verbose = F)
  # Find convex hull
  convexhull <- chull(coordinates(tmp))
  convexhull <- c(convexhull, convexhull[1])
  samplingarea <- Polygons(list(Polygon(coordinates(tmp)[convexhull, ], hole = F)), ID = "range")
  samplingarea <- SpatialPolygons(list(samplingarea), proj4string = CRS("+init=epsg:6933"))
  # Construct buffers
  print(paste0(obs_files[i], ": constructing buffers"))
  for (j in c(0, 100000, 500000)) {
    bufferedarea <- buffer(samplingarea, width = j)
    bufferedarea <- SpatialPolygonsDataFrame(bufferedarea, data.frame(ID = "buffer", row.names = row.names(bufferedarea)))
    writeOGR(bufferedarea, dsn = paste0(relPath, bufferPath), paste0(taxon, "_", j / 1000, "km"), driver="ESRI Shapefile", overwrite_layer = T)
  }
  # Remove temporary files
  removeTmpFiles(h = 0)
}

# Close cluster connection
stopCluster(cl)