library(ggplot2)

# Set directory structure
relPath <- "F:/Chapter1/"
rangeChangePath <- "Data/RangeChanges/rangeChanges.csv"
figPath <- "Figures/"

# Read range changes
rangeChanges <- read.csv(paste0(relPath, rangeChangePath), header = T, stringsAsFactors = F)
rangeChanges <- subset(rangeChanges, Time != "current")

# Run ANOVA
anov_corr <- summary(aov(RangeChange ~ Time + RCP + GCM + Threshold + VariableSelection + Model + ProjectionDistance + SamplingDistance + CorrelationFilter, data = rangeChanges))
anov_noCorr <- summary(aov(RangeChange ~ Time + RCP + GCM + Threshold + VariableSelection + Model + ProjectionDistance + SamplingDistance, data = rangeChanges))
uncComponents <- data.frame(Source = c("Time", "RCP", "GCM", "Threshold", "Variable selection", "Model complexity", "Projection distance", "Sampling distance", "Correlation filter"),
                            Uncertainty = (anov_corr[[1]]$`Sum Sq` / sum(anov_corr[[1]]$`Sum Sq`[1:9]) * 100)[1:9])

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

# Plot uncertainty
uncComponents$Source <- factor(uncComponents$Source, levels = uncComponents$Source[order(uncComponents$Uncertainty, decreasing = T)])
uncComponents$Type <- c("Climate", "Climate", "Climate", "Model", "Model", "Model", "Model", "Model", "Model")
ggplot(data = uncComponents, aes(x = Source, y = Uncertainty)) + 
  geom_bar(stat = "identity", aes(fill = Type)) +
  ylab("% variance explained (excl. residual variance)") +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 28, face = "bold")) + theme_Publication()
ggsave(paste0(relPath, figPath, "UncertaintyRangeChanges.png"), width = 12.5)
