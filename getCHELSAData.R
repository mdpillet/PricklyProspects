library(RCurl)

#
# Get current data
#

# Set directory structure
ftpURL <- "ftp://envidatrepo.wsl.ch/uploads/chelsa/chelsa_V1/climatologies/bio/"
relPath <- "F:/Chapter1/"
outPath <- "Data/CHELSA/Raw/Current/"

# Get file names
fileNames <- getURL(ftpURL, verbose = TRUE, dirlistonly = TRUE, ftp.use.epsv = FALSE) 
fileNames <- unlist(strsplit(fileNames, "\r\n"))

# Download files
for (i in fileNames) {
  print(i)
  download.file(paste0(ftpURL, i), paste0(relPath, outPath, i), mode = "wb")
}

#
# Get future data
#

# Set directory structure
outPath <- "Data/CHELSA/Raw/Future/"

# Get file names
ftpURL <- c("ftp://envidatrepo.wsl.ch/uploads/chelsa/chelsa_V1/cmip5/2041-2060/bio/", "ftp://envidatrepo.wsl.ch/uploads/chelsa/chelsa_V1/cmip5/2061-2080/bio/")
fileNames <- unlist(strsplit(getURL(ftpURL[1], verbose = TRUE, dirlistonly = TRUE, ftp.use.epsv = FALSE), "\r\n"))
fileNames <- c(fileNames, unlist(strsplit(getURL(ftpURL[2], verbose = TRUE, dirlistonly = TRUE, ftp.use.epsv = FALSE), "\r\n")))
fileNamesSubset <- fileNames[grep("rcp26|rcp45|rcp85", fileNames)]
# These GCMs are chosen based on https://doi.org/10.1175/JCLI-D-14-00362.1
fileNamesSubset <- fileNamesSubset[grep("NorESM1-M|GFDL-ESM2M|HadGEM2-AO|CCSM4|BNU-ESM", fileNamesSubset)]
fileNamesFull <- paste0("ftp://envidatrepo.wsl.ch/uploads/chelsa/chelsa_V1/cmip5/", c(rep("2041-2060/bio/", length(fileNamesSubset) / 2), rep("2061-2080/bio/", length(fileNamesSubset) / 2)), fileNamesSubset)

# Download files
for (i in 1:length(fileNamesFull)) {
  if (!file.exists(paste0(relPath, outPath, fileNamesSubset[i]))) {
    print(i / length(fileNamesFull) * 100)
    download.file(fileNamesFull[i], paste0(relPath, outPath, fileNamesSubset[i]), mode = "wb")  
  }
}