# THIS ANALYSIS HAS BEEN REPEATED REPLACING OPUNTIADS WITH SHRUBS/TREES AS GROWTH FORM

library(sp)
library(raster)
library(rgdal)
library(dplyr)
library(lme4)
library(stringr)

# Set directory structure
relPath <- "F:/Chapter1/"
rangeChangePath <- "Data/RangeChanges/rangeChanges.csv"
rangeChangeSummPath <- "Data/RangeChanges/rangeChangesSummary.csv"
iucnPath <- "Data/RedList/"
iucnDivPath <- "Data/DiversityMaps/IUCN/richness.tif"
richnessPath <- "Data/DiversityMaps/Summary/mean_current.tif"
occPath <- "Data/Occurrences/Species"
habitPath <- "Data/Habits/habits_reanalysis.csv"
heightPath <- "Data/Habits/heights.csv"
exportPath <- "Data/EcologicalCovariates/"

# Read range changes
rangeChanges <- read.csv(paste0(relPath, rangeChangePath), header = T, stringsAsFactors = F)

# Group range changes by species
rangeChanges_grouped <- rangeChanges %>% group_by(Taxon) %>% summarize(medianRangeChange = median(RangeChange, na.rm = T),
                                                                         meanRangeChange = mean(RangeChange, na.rm = T),
                                                                         sdRangeChange = sd(RangeChange, na.rm = T),
                                                                         meanCurrentSize = mean(CurrentCount1, na.rm = T),
                                                                         medianCurrentSize = median(CurrentCount1, na.rm = T),
                                                                       minRangeChange = min(RangeChange, na.rm = T),
                                                                       maxRangeChange = max(RangeChange, na.rm = T))

# Export results by species for manuscript
write.csv(rangeChanges_grouped, paste0(relPath, rangeChangeSummPath), row.names = F)

# Read in IUCN data
iucnAss <- read.csv(paste0(relPath, iucnPath, "assessments.csv"), header = T, stringsAsFactors = F)
iucnAss <- iucnAss[, c("scientificName", "redlistCategory")]
iucnSyn <- read.csv(paste0(relPath, iucnPath, "synonyms.csv"), header = T, stringsAsFactors = F)
iucnSyn$Taxon <- paste(iucnSyn$genusName, iucnSyn$speciesName, sep = " ")

# Add IUCN data
rangeChanges_grouped$Species <- gsub("_", " ", rangeChanges_grouped$Taxon)
rangeChanges_grouped$Species <- gsub(".", "-", rangeChanges_grouped$Species, fixed = T)
for (i in 1:nrow(rangeChanges_grouped)) {
  tmp <- subset(iucnAss, scientificName == as.character(rangeChanges_grouped[i, "Species"]))
  if (nrow(tmp) == 1) rangeChanges_grouped[i, "Assessment"] <- tmp$redlistCategory
  # Check for synonyms
  else {
    tmp <- subset(iucnSyn, Taxon == as.character(rangeChanges_grouped[i, "Species"]))
    if (nrow(tmp) == 1) {
      tmp <- subset(iucnAss, scientificName == tmp$scientificName)
      rangeChanges_grouped[i, "Assessment"] <- tmp$redlistCategory
    }
  }
}

# Export assessment data
write.csv(rangeChanges_grouped[, c("Species", "Assessment")], paste0(relPath, iucnPath, "extrAssessment.csv"), row.names = F)

# Read in manual assessment data
assessments <- read.csv(paste0(relPath, iucnPath, "extrAssessmentManual.csv"), header = T, stringsAsFactors = F)
rangeChanges_grouped$Assessment <- assessments$Assessment

# Find current species richness for each species
# Read in richness raster
modelDiv <- raster(paste0(relPath, richnessPath))
iucnDiv <- raster(paste0(relPath, iucnDivPath))
for (i in 1:nrow(rangeChanges_grouped)) {
  # Get species name
  speciesName <- as.character(pull(rangeChanges_grouped[i, "Taxon"]))
  print(speciesName)
  # Read in occurrence data
  occ <- readOGR(dsn = paste0(relPath, occPath), layer = gsub(".", "-", speciesName, fixed = T), verbose = F)
  # Extract species richness
  modelExtract <- extract(modelDiv, occ)
  iucnExtract <- extract(iucnDiv, occ)
  rangeChanges_grouped[i, "richnessModeled_mean"] <- mean(modelExtract, na.rm = T)
  rangeChanges_grouped[i, "richnessIUCN_mean"] <- mean(iucnExtract, na.rm = T)
  rangeChanges_grouped[i, "richnessModeled_median"] <- median(modelExtract, na.rm = T)
  rangeChanges_grouped[i, "richnessIUCN_median"] <- median(iucnExtract, na.rm = T)
}

# Read in habit data
habits <- read.csv(paste0(relPath, habitPath), header = T, stringsAsFactors = F)
# habits[habits$Habit == "PricklyPear", "Habit"] <- "Shrub"
# rangeChanges_grouped$Habit <- habits$HabitReanalysis
habits$Species <- str_replace(habits$Species, " ", "_")
habits$Species <- str_replace(habits$Species, "-", ".")
rangeChanges_grouped$Habit <- character(length = 408)
for (i in 1:nrow(rangeChanges_grouped)) {
  print(i)
  tmp <- subset(habits, Species == as.character(rangeChanges_grouped[i, "Taxon"]))
  if (nrow(tmp) == 1) rangeChanges_grouped[i, "Habit"] <- tmp$HabitReanalysis
  else rangeChanges_grouped[i, "Habit"] <- NA
}

# Export covariates
write.csv(rangeChanges_grouped, paste0(relPath, exportPath, "covariates.csv"), row.names = F)

# Add height/diameter data
heights <- read.csv(paste0(relPath, heightPath), header = T, stringsAsFactors = F)
names(heights) <- c("Genus", "Species", "BodyHeight", "BodyDiameter", "BranchLength", "BranchDiameter")
heights$scientificName <- paste(heights$Genus, heights$Species, sep = "_")
rangeChanges_grouped$Height <- numeric(length = 408)
for (i in 1:nrow(rangeChanges_grouped)) {
  print(i)
  tmp <- subset(heights, scientificName == as.character(rangeChanges_grouped[i, "Taxon"]))
  if (nrow(tmp) == 1) rangeChanges_grouped[i, "Height"] <- tmp$BodyHeight
  else rangeChanges_grouped[i, "Height"] <- NA
}

# Perform regression
# summary(lm(meanRangeChange ~ meanCurrentSize + Habit + Assessment + richnessModeled_mean, data = rangeChanges_grouped))
summary(lm(meanRangeChange ~ meanCurrentSize + Habit + Assessment + richnessIUCN_mean, data = rangeChanges_grouped))
confint(lm(meanRangeChange ~ meanCurrentSize + Habit + Assessment + richnessIUCN_mean, data = rangeChanges_grouped))
# summary(lm(meanRangeChange ~ meanCurrentSize + Habit + Assessment + richnessIUCN_mean + Height + I(Height^2), data = rangeChanges_grouped))