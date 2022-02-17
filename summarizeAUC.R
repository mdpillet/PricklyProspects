library(ggplot2)

# Set directory structure
relPath <- "F:/Chapter1/"
AUCPath <- "Data/ModelStats/modelStats.csv"
AUCcvPath <- "Data/ModelStats/modelStatsCV.csv"
subsetPath <- "Data/Habits/habits.csv"
figPath <- "Figures/"

# Read in AUC
AUCs <- read.csv(paste0(relPath, AUCPath), header = T, stringsAsFactors = F)
AUCsCV <- read.csv(paste0(relPath, AUCcvPath), header = T, stringsAsFactors = F)

# Add columns for taxon, sampling distance, variable selection, and correlation filter
AUCs$Taxon <- sapply(lapply(strsplit(AUCs$Model, "_"), "[", 1:2), paste, collapse = "_")
AUCs$SamplingDistance <- sapply(strsplit(AUCs$Model, "_"), "[", 3)
AUCs$VariableSelection <- sapply(strsplit(AUCs$Model, "_"), "[", 4)
AUCs$CorrelationFilter <- sapply(sapply(lapply(strsplit(AUCs$Model, "_"), "[", 5), strsplit, ".", fixed = T), "[", 1)
AUCsCV$Taxon <- sapply(lapply(strsplit(AUCsCV$Model, "_"), "[", 1:2), paste, collapse = "_")
AUCsCV$SamplingDistance <- sapply(strsplit(AUCsCV$Model, "_"), "[", 3)
AUCsCV$VariableSelection <- sapply(strsplit(AUCsCV$Model, "_"), "[", 4)
AUCsCV$CorrelationFilter <- sapply(strsplit(AUCsCV$Model, "_"), "[", 5)
AUCsCV$FullModel <- paste0(paste(AUCsCV$Taxon, AUCsCV$SamplingDistance, AUCsCV$VariableSelection, AUCsCV$CorrelationFilter, sep = "_"), ".rda")

# Add cross-validated AUC to full models
for (i in 1:nrow(AUCs)) {
  print(i)
  tmp <- subset(AUCsCV, FullModel == AUCs[i, "Model"])
  tmp <- subset(tmp, Features == AUCs[i, "Features"])
  if (!is.na(nrow(tmp))) AUCs[i, "AUCcv"] <- mean(tmp$AUC, na.rm = T)
}

# Calculate difference between total and cross-validated AUC
AUCs$deltaAUC <- AUCs$AUC - AUCs$AUCcv

# List species to be cross-validated
cvSpecies <- read.csv(paste0(relPath, subsetPath), header = T, stringsAsFactors = F)
cvSpecies <- gsub(" ", "_", cvSpecies$Species)

# Check which species/models were not cross-validated (should be 408 species x 12 x 3 models = 14,688 models)
AUCs <- subset(AUCs, Taxon %in% cvSpecies)

# Get AUC statistics
summary(AUCs$AUC)
summary(AUCs$AUCcv)
summary(AUCs$deltaAUC)

# Compare model AUC
summary(lm(AUC ~ SamplingDistance + VariableSelection + CorrelationFilter, data = AUCs))

### ADD PUBLICATION THEME

theme_Publication <- function(base_size=14, base_family="") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

# Plot AUCs
ggplot(data = AUCs, aes(x = AUC)) + geom_histogram(binwidth = 0.01) +
  theme(axis.text = element_text(size = 26), axis.title = element_text(size = 28, face = "bold")) + 
  xlab("AUC") + ylab("Number of models") + theme_Publication()
ggsave(paste0(relPath, figPath, "AUCHist.png"))