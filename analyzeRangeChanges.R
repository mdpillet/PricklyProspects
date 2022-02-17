library(dplyr)
library(ggplot2)
library(lme4)
library(gridExtra)

# Set directory structure
relPath <- "F:/Chapter1/"
rangeChangePath <- "Data/RangeChanges/rangeChanges.csv"
figPath <- "Figures/"

# Read range changes
rangeChanges <- read.csv(paste0(relPath, rangeChangePath), header = T, stringsAsFactors = F)
rangeChanges <- subset(rangeChanges, Time != "current")

### REGRESSION

# Regress range change on factors
reg_noCorr <- lm(RangeChange ~ Time + RCP + GCM + Threshold + VariableSelection + Model + ProjectionDistance + SamplingDistance, data = rangeChanges)
summary(reg_noCorr)
reg_corr <- lm(RangeChange ~ Time + RCP + GCM + Threshold + VariableSelection + Model + ProjectionDistance + SamplingDistance + CorrelationFilter, data = rangeChanges)
summary(reg_corr)
confint(reg_corr)
reg_corr_mixed <- lmer(RangeChange ~ Time + RCP + GCM + Threshold + VariableSelection + Model + ProjectionDistance + SamplingDistance + CorrelationFilter + (1|Taxon), data = rangeChanges)
summary(reg_corr_mixed)

### OVERALL STATISTICS

# Group range change table by species and calculate summary statistics for each species
rangeChanges_S <- rangeChanges %>% group_by(Taxon) %>% summarize(medianRangeChange = median(RangeChange, na.rm = T),
                                                                   meanRangeChange = mean(RangeChange, na.rm = T),
                                                                   sdRangeChange = sd(RangeChange, na.rm = T))
rangeChanges_S$cvRangeChange <- rangeChanges_S$sdRangeChange / rangeChanges_S$meanRangeChange
rangeChanges_S$lowBound <- rangeChanges_S$meanRangeChange - rangeChanges_S$sdRangeChange
rangeChanges_S$upBound <- rangeChanges_S$meanRangeChange + rangeChanges_S$sdRangeChange

# Calculate mean of means and mean of medians
1 - mean(rangeChanges_S$meanRangeChange)
1 - median(rangeChanges_S$meanRangeChange)

# Calculate confidence interval
mean(rangeChanges_S$meanRangeChange) + 1.96 * sd(rangeChanges_S$meanRangeChange)
mean(rangeChanges_S$meanRangeChange) - 1.96 * sd(rangeChanges_S$meanRangeChange)

# Get counts of species with any loss and >25% loss
table(rangeChanges_S$meanRangeChange < 0.75) / nrow(rangeChanges_S) * 100
table(rangeChanges_S$medianRangeChange < 0.75) / nrow(rangeChanges_S) * 100
table(rangeChanges_S$meanRangeChange < 1) / nrow(rangeChanges_S) * 100
table(rangeChanges_S$medianRangeChange < 1) / nrow(rangeChanges_S) * 100

# Get counts of species with >25% gain
table(rangeChanges_S$meanRangeChange > 1.25) / nrow(rangeChanges_S) * 100
table(rangeChanges_S$medianRangeChange > 1.25) / nrow(rangeChanges_S) * 100

# Summarize robustness of projections
table(rangeChanges_S$meanRangeChange < 1)
speciesDecline <- subset(rangeChanges_S, meanRangeChange < 1)
table(speciesDecline$upBound < 1)
table(rangeChanges_S$meanRangeChange > 1)
speciesIncrease <- subset(rangeChanges_S, meanRangeChange > 1)
table(speciesIncrease$lowBound > 1)

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

### RANGE CHANGE PLOTS

# Plot range changes (histogram)
p0 <- ggplot(data = rangeChanges_S, aes(x = meanRangeChange)) + geom_histogram(binwidth = 0.05) +
  theme(axis.text = element_text(size = 26), axis.title = element_text(size = 28, face = "bold")) + 
  xlab("Mean proportional range size change") + ylab("Number of species") + xlim(0, 2.5) + ylim(0, 40) + theme_Publication()
ggsave(paste0(relPath, figPath, "RangeChangeHist_S.png"))

# Plot range changes (boxplot)
p9 <- ggplot(data = rangeChanges, aes(x = Model, y = RangeChange)) + geom_boxplot() + ylim(0, 2) + theme_Publication() +
  scale_x_discrete(breaks = c("modelL", "modelLQ", "modelLQH"),
                   labels = c("Linear", "Linear/quadratic", "Linear/quadratic/hinge")) + xlab("Model complexity") + ylab("") +
  scale_y_continuous(labels = NULL,
                     limits = c(0, 2))
p8 <- ggplot(data = rangeChanges, aes(x = Time, y = RangeChange)) + geom_boxplot() + ylim(0, 2) + theme_Publication() +
  scale_x_discrete(breaks = c("2041.2060", "2061.2080"),
                   labels = c("2041-2060", "2061-2080")) + ylab("") +
  scale_y_continuous(labels = NULL,
                     limits = c(0, 2))
p7 <- ggplot(data = rangeChanges, aes(x = RCP, y = RangeChange)) + geom_boxplot() + ylim(0, 2) + theme_Publication() +
  scale_x_discrete(breaks = c("rcp26", "rcp45", "rcp85"),
                   labels = c("2.6", "4.5", "8.5")) + ylab("") +
  scale_y_continuous(labels = c("0%", "50%", "100%", "150%", "200%"),
                     limits = c(0, 2))
p6 <- ggplot(data = rangeChanges, aes(x = Threshold, y = RangeChange)) + geom_boxplot() + ylim(0, 2) + theme_Publication() +
  scale_x_discrete(breaks = c("MaxTSS.tif", "Omiss5.tif"),
                   labels = c("Maximum TSS", "5% omission")) + xlab("Thresholding method") + ylab("") +
  scale_y_continuous(labels = NULL,
                     limits = c(0, 2))
p5 <- ggplot(data = rangeChanges, aes(x = GCM, y = RangeChange)) + geom_boxplot() + ylim(0, 2) + theme_Publication() + ylab("") +
  scale_y_continuous(labels = NULL,
                     limits = c(0, 2))
p4 <- ggplot(data = rangeChanges, aes(x = VariableSelection, y = RangeChange)) + geom_boxplot() + ylim(0, 2) + theme_Publication() +
  scale_x_discrete(breaks = c("varBIEN", "varPCAraw", "varPCAtrans", "varRandom"),
                   labels = c("One-size-fits-all", "PCA (raw)", "PCA (transformed)", "Random")) + xlab("Variable selection") + ylab("Proportional range change") +
  scale_y_continuous(labels = c("0%", "50%", "100%", "150%", "200%"),
                     limits = c(0, 2))
p3 <- ggplot(data = rangeChanges, aes(x = ProjectionDistance, y = RangeChange)) + geom_boxplot() + ylim(0, 2) + theme_Publication() +
  scale_x_discrete(breaks = c("d0", "d100", "d500"),
                   labels = c("0km", "100 km", "500 km")) + xlab("Projection distance") + ylab("") +
  scale_y_continuous(labels = NULL,
                     limits = c(0, 2))
p2 <- ggplot(data = rangeChanges, aes(x = SamplingDistance, y = RangeChange)) + geom_boxplot() + ylim(0, 2) + theme_Publication() + 
  scale_x_discrete(breaks = c("100km", "500km"),
                   labels = c("100 km", "500 km")) + xlab("Sampling distance") + ylab("") +
  scale_y_continuous(labels = NULL,
                     limits = c(0, 2))
p1 <- ggplot(data = rangeChanges, aes(x = CorrelationFilter, y = RangeChange)) + geom_boxplot() + ylim(0, 2) + theme_Publication() + 
  scale_x_discrete(breaks = c("corrF", "corrNA", "corrT"),
                   labels = c("False", "NA", "True")) + xlab("Correlation filter") + ylab("") +
  scale_y_continuous(labels = c("0%", "50%", "100%", "150%", "200%"),
                     limits = c(0, 2))
fig1 <- grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3)
ggsave(paste0(relPath, figPath, "Fig1.png"), fig1, height = 20, width = 18)

### BREAKDOWN BY RCP

# Group range change table by species and RCP and calculate summary statistics for each species
rangeChanges_SR <- rangeChanges %>% group_by(Taxon, RCP) %>% summarize(medianRangeChange = median(RangeChange, na.rm = T),
                                                                 meanRangeChange = mean(RangeChange, na.rm = T),
                                                                 sdRangeChange = sd(RangeChange, na.rm = T))
# Get mean change by RCP
RCPmeans <- rangeChanges_SR %>% group_by(RCP) %>% summarize(meanRangeChangeRCP = mean(meanRangeChange))

# Plot range changes
ggplot(data = rangeChanges_SR, aes(x = meanRangeChange)) + geom_histogram() +
  theme(axis.text = element_text(size = 26), axis.title = element_text(size = 28, face = "bold")) + 
  xlab("Mean proportional range size change") + ylab("Number of species") + 
  facet_wrap(~RCP) + geom_vline(data = RCPmeans, aes(xintercept = meanRangeChangeRCP), col = "red", lwd = 2) +
  geom_text(data = RCPmeans, aes(x = 1.4, y = 60, label = paste0("Mean: ", round(meanRangeChangeRCP, 2))))
ggsave(paste0(relPath, figPath, "RangeChangeHist_SR.png"))

### BREAKDOWN BY TIME

# Group range change table by species and time and calculate summary statistics for each species
rangeChanges_ST <- rangeChanges %>% group_by(Taxon, Time) %>% summarize(medianRangeChange = median(RangeChange, na.rm = T),
                                                                       meanRangeChange = mean(RangeChange, na.rm = T),
                                                                       sdRangeChange = sd(RangeChange, na.rm = T))
# Get mean change by time
timemeans <- rangeChanges_ST %>% group_by(Time) %>% summarize(meanRangeChangeTime = mean(meanRangeChange))

# Plot range changes
ggplot(data = rangeChanges_ST, aes(x = meanRangeChange)) + geom_histogram() +
  theme(axis.text = element_text(size = 26), axis.title = element_text(size = 28, face = "bold")) + 
  xlab("Mean proportional range size change") + ylab("Number of species") + 
  facet_wrap(~Time) + geom_vline(data = timemeans, aes(xintercept = meanRangeChangeTime), col = "red", lwd = 2) +
  geom_text(data = timemeans, aes(x = 1.4, y = 60, label = paste0("Mean: ", round(meanRangeChangeTime, 2))))
ggsave(paste0(relPath, figPath, "RangeChangeHist_ST.png"))

### BREAKDOWN BY SAMPLING DISTANCE

# Group range change table by species and sampling distance and calculate summary statistics for each species
rangeChanges_SS <- rangeChanges %>% group_by(Taxon, SamplingDistance) %>% summarize(medianRangeChange = median(RangeChange, na.rm = T),
                                                                        meanRangeChange = mean(RangeChange, na.rm = T),
                                                                        sdRangeChange = sd(RangeChange, na.rm = T))
# Get mean change by sampling distance
sdistmeans <- rangeChanges_SS %>% group_by(SamplingDistance) %>% summarize(meanRangeChangeSamplingDistance = mean(meanRangeChange))

# Plot range changes
ggplot(data = rangeChanges_SS, aes(x = meanRangeChange)) + geom_histogram() +
  theme(axis.text = element_text(size = 26), axis.title = element_text(size = 28, face = "bold")) + 
  xlab("Mean proportional range size change") + ylab("Number of species") + 
  facet_wrap(~SamplingDistance) + geom_vline(data = sdistmeans, aes(xintercept = meanRangeChangeSamplingDistance), col = "red", lwd = 2) +
  geom_text(data = sdistmeans, aes(x = 1.4, y = 60, label = paste0("Mean: ", round(meanRangeChangeSamplingDistance, 2))))
ggsave(paste0(relPath, figPath, "RangeChangeHist_SS.png"))

### BREAKDOWN BY VARIABLE SELECTION

# Group range change table by species and variable selection and calculate summary statistics for each species
rangeChanges_SV <- rangeChanges %>% group_by(Taxon, VariableSelection) %>% summarize(medianRangeChange = median(RangeChange, na.rm = T),
                                                                                    meanRangeChange = mean(RangeChange, na.rm = T),
                                                                                    sdRangeChange = sd(RangeChange, na.rm = T))
# Get mean change by variable selection
varmeans <- rangeChanges_SV %>% group_by(VariableSelection) %>% summarize(meanRangeChangeVar = mean(meanRangeChange))

# Plot range changes
ggplot(data = rangeChanges_SV, aes(x = meanRangeChange)) + geom_histogram() +
  theme(axis.text = element_text(size = 26), axis.title = element_text(size = 28, face = "bold")) + 
  xlab("Mean proportional range size change") + ylab("Number of species") + 
  facet_wrap(~VariableSelection) + geom_vline(data = varmeans, aes(xintercept = meanRangeChangeVar), col = "red", lwd = 2) +
  geom_text(data = varmeans, aes(x = 1.4, y = 60, label = paste0("Mean: ", round(meanRangeChangeVar, 2))))
ggsave(paste0(relPath, figPath, "RangeChangeHist_SV.png"))

### BREAKDOWN BY CORRELATION FILTER

# Group range change table by species and correlation filter and calculate summary statistics for each species
rangeChanges_SC <- rangeChanges %>% group_by(Taxon, CorrelationFilter) %>% summarize(medianRangeChange = median(RangeChange, na.rm = T),
                                                                                     meanRangeChange = mean(RangeChange, na.rm = T),
                                                                                     sdRangeChange = sd(RangeChange, na.rm = T))
# Get mean change by correlation filter
corrmeans <- rangeChanges_SC %>% group_by(CorrelationFilter) %>% summarize(meanRangeChangeCorr = mean(meanRangeChange))

# Plot range changes
ggplot(data = rangeChanges_SC, aes(x = meanRangeChange)) + geom_histogram() +
  theme(axis.text = element_text(size = 26), axis.title = element_text(size = 28, face = "bold")) + 
  xlab("Mean proportional range size change") + ylab("Number of species") + 
  facet_wrap(~CorrelationFilter) + geom_vline(data = corrmeans, aes(xintercept = meanRangeChangeCorr), col = "red", lwd = 2) +
  geom_text(data = corrmeans, aes(x = 1.4, y = 60, label = paste0("Mean: ", round(meanRangeChangeCorr, 2))))
ggsave(paste0(relPath, figPath, "RangeChangeHist_SC.png"))

### BREAKDOWN BY GCM

# Group range change table by species and GCM and calculate summary statistics for each species
rangeChanges_SG <- rangeChanges %>% group_by(Taxon, GCM) %>% summarize(medianRangeChange = median(RangeChange, na.rm = T),
                                                                                     meanRangeChange = mean(RangeChange, na.rm = T),
                                                                                     sdRangeChange = sd(RangeChange, na.rm = T))
# Get mean change by GCM
gcmmeans <- rangeChanges_SG %>% group_by(GCM) %>% summarize(meanRangeChangeGCM = mean(meanRangeChange))

# Plot range changes
ggplot(data = rangeChanges_SG, aes(x = meanRangeChange)) + geom_histogram() +
  theme(axis.text = element_text(size = 26), axis.title = element_text(size = 28, face = "bold")) + 
  xlab("Mean proportional range size change") + ylab("Number of species") + 
  facet_wrap(~GCM) + geom_vline(data = gcmmeans, aes(xintercept = meanRangeChangeGCM), col = "red", lwd = 2) +
  geom_text(data = gcmmeans, aes(x = 1.4, y = 60, label = paste0("Mean: ", round(meanRangeChangeGCM, 2))))
ggsave(paste0(relPath, figPath, "RangeChangeHist_SG.png"))

### BREAKDOWN BY PROJECTION DISTANCE

# Group range change table by species and projection distance and calculate summary statistics for each species
rangeChanges_SP <- rangeChanges %>% group_by(Taxon, ProjectionDistance) %>% summarize(medianRangeChange = median(RangeChange, na.rm = T),
                                                                       meanRangeChange = mean(RangeChange, na.rm = T),
                                                                       sdRangeChange = sd(RangeChange, na.rm = T))

# Get species counts
d0 <- subset(rangeChanges_SP, ProjectionDistance == "d0")
d100 <- subset(rangeChanges_SP, ProjectionDistance == "d100")
d500 <- subset(rangeChanges_SP, ProjectionDistance == "d500")
table(d0$meanRangeChange < 1) / nrow(d0) * 100
table(d0$medianRangeChange < 1) / nrow(d0) * 100
table(d100$meanRangeChange < 1) / nrow(d100) * 100
table(d100$medianRangeChange < 1) / nrow(d100) * 100
table(d500$meanRangeChange < 1) / nrow(d500) * 100
table(d500$medianRangeChange < 1) / nrow(d500) * 100

# Get mean change by projection distance
projmeans <- rangeChanges_SP %>% group_by(ProjectionDistance) %>% summarize(meanRangeChangeProjDist = mean(meanRangeChange))

# Plot range changes
ggplot(data = rangeChanges_SP, aes(x = meanRangeChange)) + geom_histogram() +
  theme(axis.text = element_text(size = 26), axis.title = element_text(size = 28, face = "bold")) + 
  xlab("Mean proportional range size change") + ylab("Number of species") + 
  facet_wrap(~ProjectionDistance) + geom_vline(data = projmeans, aes(xintercept = meanRangeChangeProjDist), col = "red", lwd = 2) +
  geom_text(data = projmeans, aes(x = 1.4, y = 60, label = paste0("Mean: ", round(meanRangeChangeProjDist, 2))))
ggsave(paste0(relPath, figPath, "RangeChangeHist_SP.png"))

### BREAKDOWN BY MODEL

# Group range change table by species and model and calculate summary statistics for each species
rangeChanges_SM <- rangeChanges %>% group_by(Taxon, Model) %>% summarize(medianRangeChange = median(RangeChange, na.rm = T),
                                                                                      meanRangeChange = mean(RangeChange, na.rm = T),
                                                                                      sdRangeChange = sd(RangeChange, na.rm = T))
# Get mean change by model
modelmeans <- rangeChanges_SM %>% group_by(Model) %>% summarize(meanRangeChangeModel = mean(meanRangeChange))

# Plot range changes
ggplot(data = rangeChanges_SM, aes(x = meanRangeChange)) + geom_histogram() +
  theme(axis.text = element_text(size = 26), axis.title = element_text(size = 28, face = "bold")) + 
  xlab("Mean proportional range size change") + ylab("Number of species") + 
  facet_wrap(~Model) + geom_vline(data = modelmeans, aes(xintercept = meanRangeChangeModel), col = "red", lwd = 2) +
  geom_text(data = modelmeans, aes(x = 1.4, y = 60, label = paste0("Mean: ", round(meanRangeChangeModel, 2))))
ggsave(paste0(relPath, figPath, "RangeChangeHist_SM.png"))

### BREAKDOWN BY THRESHOLD

# Group range change table by species and threshold and calculate summary statistics for each species
rangeChanges_STh <- rangeChanges %>% group_by(Taxon, Threshold) %>% summarize(medianRangeChange = median(RangeChange, na.rm = T),
                                                                         meanRangeChange = mean(RangeChange, na.rm = T),
                                                                         sdRangeChange = sd(RangeChange, na.rm = T))
# Get mean change by threshold
thresholdmeans <- rangeChanges_STh %>% group_by(Threshold) %>% summarize(meanRangeChangeThreshold = mean(meanRangeChange))

# Plot range changes
ggplot(data = rangeChanges_STh, aes(x = meanRangeChange)) + geom_histogram() +
  theme(axis.text = element_text(size = 26), axis.title = element_text(size = 28, face = "bold")) + 
  xlab("Mean proportional range size change") + ylab("Number of species") + 
  facet_wrap(~Threshold) + geom_vline(data = thresholdmeans, aes(xintercept = meanRangeChangeThreshold), col = "red", lwd = 2) +
  geom_text(data = thresholdmeans, aes(x = 1.4, y = 60, label = paste0("Mean: ", round(meanRangeChangeThreshold, 2))))
ggsave(paste0(relPath, figPath, "RangeChangeHist_STh.png"))