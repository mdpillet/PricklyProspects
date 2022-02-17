### THIS REANALYSIS IS BASED ON PROBABILISTIC ESTIMATION OF RICHNESS

library(raster)
library(RColorBrewer)
library(rasterVis)
library(gridExtra)
library(ggplot2)
library(cowplot)

# Set directory structure
relPath <- "F:/Chapter1/"
diversityPath <- "Data/DiversityMaps/SummaryReanalysis/"
iucnPath <- "Data/DiversityMaps/IUCN/"
outPath <- "Figures/"
footprintPath <- "Data/HumanFootprint/Dryadv3/Maps/HFP2009.tif"
figPath <- "Figures/"
uncPath <- "Data/UncertaintyRichness/uncertaintyMap_reanalysis.tif"

# Script for divergent color ramp (https://gist.github.com/johnbaums/306e4b7e69c87b1826db)
diverge0 <- function(p, ramp) {
  # p: a trellis object resulting from rasterVis::levelplot
  # ramp: the name of an RColorBrewer palette (as character), a character 
  #       vector of colour names to interpolate, or a colorRampPalette.
  require(RColorBrewer)
  require(rasterVis)
  if(length(ramp)==1 && is.character(ramp) && ramp %in% 
     row.names(brewer.pal.info)) {
    ramp <- suppressWarnings(colorRampPalette(brewer.pal(11, ramp)))
  } else if(length(ramp) > 1 && is.character(ramp) && all(ramp %in% colors())) {
    ramp <- colorRampPalette(ramp)
  } else if(!is.function(ramp)) 
    stop('ramp should be either the name of a RColorBrewer palette, ', 
         'a vector of colours to be interpolated, or a colorRampPalette.')
  rng <- range(p$legend[[1]]$args$key$at)
  s <- seq(-max(abs(rng)), max(abs(rng)), len=1001)
  i <- findInterval(rng[which.min(abs(rng))], s)
  zlim <- switch(which.min(abs(rng)), `1`=i:(1000+1), `2`=1:(i+1))
  p$legend[[1]]$args$key$at <- s[zlim]
  p[[grep('^legend', names(p))]][[1]]$args$key$col <- ramp(1000)[zlim[-length(zlim)]]
  p$panel.args.common$col.regions <- ramp(1000)[zlim[-length(zlim)]]
  p
}

# Read mean richness maps
current <- raster(paste0(relPath, diversityPath, "mean_current.tif"))
f50 <- raster(paste0(relPath, diversityPath, "mean_2041.2060.tif"))
f70 <- raster(paste0(relPath, diversityPath, "mean_2061.2080.tif"))

# Display maps
pal <- rasterTheme(brewer.pal(n = 11, "RdYlGn"))
plotA <- levelplot(current, par.settings = pal, at = seq(0, 70, 1), main = "Present")
plotB <- levelplot(f50, par.settings = pal, at = seq(0, 70, 1), main = "2050")
plotC <- levelplot(f70, par.settings = pal, at = seq(0, 70, 1), main = "2070")
grid.arrange(plotA, plotB, plotC, ncol = 3)

# Read IUCN richness map
iucn <- raster(paste0(relPath, iucnPath, "richness.tif"))

# Plot modeled current richness and IUCN richness maps
plotCurrent <- gplot(current) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "none") + 
  geom_tile(aes(fill = value)) + 
  scale_fill_gradient(low = "red", high = "green", na.value = NA, limits = c(0, 80)) + 
  coord_equal() +
  ggtitle("Modeled richness")
plotIUCN <- gplot(iucn) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  geom_tile(aes(fill = value)) + 
  scale_fill_gradient(low = "red", high = "green", na.value = NA, limits = c(0, 80)) + 
  coord_equal() +
  ggtitle("IUCN richness") +
  labs(fill = "Number of species")
fig3 <- plot_grid(plotCurrent, plotIUCN, labels = c("A", "B"), label_size = 12, align = "v", axis = "r")
save_plot(paste0(relPath, figPath, "ComparisonModeledIUCN_reanalysis.png"), fig3, ncol = 2)

# Plot correlation map
corRaster <- raster(paste0(relPath, iucnPath, "correlation_windowSize9.tif"))
plotCorr <- gplot(corRaster) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  geom_tile(aes(fill = value)) + 
  scale_fill_gradient(low = "red", high = "green", na.value = NA, limits = c(-1, 1)) + 
  coord_equal() +
  labs(fill = "Correlation")
ggsave(paste0(relPath, figPath, "CorrelationModeledIUCN_reanalysis.png"), plotCorr)

# Plot difference maps
diff50 <- f50 - current
writeRaster(diff50, paste0(relPath, diversityPath, "changeAbs_2050-current_reanalysis.tif"))
diff70 <- f70 - current
writeRaster(diff70, paste0(relPath, diversityPath, "changeAbs_2070-current_reanalysis.tif"))
plotD <- levelplot(diff50, main = "2050 - present")
plotD_div <- diverge0(plotD, "PiYG")
plotE <- levelplot(diff70, main = "2070 - present")
plotE_div <- diverge0(plotE, "PiYG")
grid.arrange(plotD_div, plotE_div, ncol = 2)

# Plot relative difference maps
diffrel50 <- diff50 / current
writeRaster(diffrel50, paste0(relPath, diversityPath, "changeRel_2050-current_reanalysis.tif"))
diffrel70 <- diff70 / current
writeRaster(diffrel70, paste0(relPath, diversityPath, "changeRel_2070-current_reanalysis.tif"))
plotF <- levelplot(diffrel50, par.settings = pal, at = seq(-1, 1, 0.2), main = "2050 - present (rel.)")
plotF_div <- diverge0(plotF, "PiYG")
plotG <- levelplot(diffrel70, par.settings = pal, at = seq(-1, 1, 0.2), main = "2070 - present (rel.)")
plotG_div <- diverge0(plotG, "PiYG")
grid.arrange(plotF_div, plotG_div, ncol = 2)

# Plot change maps
plot50A <- gplot(diff50) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5, size = 16),
        panel.background = element_rect(fill = NA), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "none") + 
  geom_tile(aes(fill = value)) + 
  scale_fill_gradient2(low = "red", mid = "yellow", high = "green", na.value = NA, limits = c(-15, 6)) + 
  coord_equal() +
  ggtitle("Absolute richness change (2041-2060)")
plot70A <- gplot(diff70) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5, size = 16),
        panel.background = element_rect(fill = NA), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  geom_tile(aes(fill = value)) + 
  scale_fill_gradient2(low = "red", mid = "yellow", high = "green", na.value = NA, limits = c(-15, 6)) + 
  coord_equal() +
  ggtitle("Absolute richness change (2061-2080)") +
  labs(fill = "Number of species")
plot50R <- gplot(diffrel50) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5, size = 16),
        panel.background = element_rect(fill = NA), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "none") + 
  geom_tile(aes(fill = value)) + 
  scale_fill_gradient2(low = "red", mid = "yellow", high = "green", na.value = NA, limits = c(-1, 1)) + 
  coord_equal() +
  ggtitle("Relative richness change (2041-2060)")
plot70R <- gplot(diffrel70) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5, size = 16),
        panel.background = element_rect(fill = NA), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  geom_tile(aes(fill = value)) + 
  scale_fill_gradient2(low = "red", mid = "yellow", high = "green", na.value = NA, limits = c(-1, 1)) + 
  coord_equal() +
  ggtitle("Relative richness change (2061-2080)") +
  labs(fill = "Proportional change\n(number of species)")
fig4 <- plot_grid(plot50A, plot70A, plot50R, plot70R, labels = c("A", "B", "C", "D"), label_size = 12, align = "v", axis = "r")
save_plot(paste0(relPath, figPath, "RichnessChanges_reanalysis.png"), fig4, ncol = 2, nrow = 2)

# Read SD richness maps
currentsd <- raster(paste0(relPath, diversityPath, "sd_current.tif"))
f50sd <- raster(paste0(relPath, diversityPath, "sd_2041.2060.tif"))
f70sd <- raster(paste0(relPath, diversityPath, "sd_2061.2080.tif"))

# Plot standard deviation
plotSD <- gplot(currentsd) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  geom_tile(aes(fill = value)) + 
  scale_fill_gradient(low = "green", high = "red", na.value = NA, limits = c(0, 30)) + 
  coord_equal() +
  labs(fill = "Standard deviation\n(number of species)")
ggsave(paste0(relPath, figPath, "SDCurrent_reanalysis.png"), plotSD)

# Display maps
palRev <- rasterTheme(rev(brewer.pal(n = 11, "RdYlGn")))
plotH <- levelplot(currentsd, main = "Present")
plotI <- levelplot(f50sd, main = "2050")
plotJ <- levelplot(f70sd, main = "2070")
grid.arrange(plotH, plotI, plotJ, ncol = 3)

# Read CV richness maps
currentcv <- raster(paste0(relPath, diversityPath, "cv_current.tif"))
f50cv <- raster(paste0(relPath, diversityPath, "cv_2041.2060.tif"))
f70cv <- raster(paste0(relPath, diversityPath, "cv_2061.2080.tif"))

# Plot CV
plotCV <- gplot(currentcv) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  geom_tile(aes(fill = value)) + 
  scale_fill_gradient(low = "green", high = "red", na.value = NA, trans = "log1p", labels = scales::trans_format("identity", function(x) round(x, 0))) + 
  coord_equal() +
  labs(fill = "Coefficient of variation")
ggsave(paste0(relPath, figPath, "CVCurrent_reanalysis.png"), plotCV)

plotCV50 <- gplot(f50cv) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "right") + 
  scale_fill_gradient(low = "green", high = "red", na.value = NA, trans = "log", labels = scales::trans_format("identity", function(x) round(x, 0))) + 
  coord_equal() +
  geom_tile(aes(fill = value)) + 
  ggtitle("2041-2060") +
  labs(fill = "Coefficient of variation")
plotCV70 <- gplot(f70cv) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  geom_tile(aes(fill = value)) + 
  scale_fill_gradient(low = "green", high = "red", na.value = NA, trans = "log", labels = scales::trans_format("identity", function(x) round(x, 0))) + 
  coord_equal() +
  geom_tile(aes(fill = value)) + 
  ggtitle("2061-2080") +
  labs(fill = "Coefficient of variation")
fig5 <- plot_grid(plotCV50, plotCV70, labels = c("A", "B"), label_size = 12, align = "v", axis = "r", ncol = 2)
save_plot(paste0(relPath, figPath, "FutureCV_reanalysis.png"), fig5, ncol = 2, nrow = 1)

# Display maps
plotK <- levelplot(currentcv, main = "Present")
plotL <- levelplot(f50cv, main = "2050")
plotM <- levelplot(f70cv, main = "2070")
grid.arrange(plotK, plotL, plotM, ncol = 3)

# Read in human footprint map
footprint <- raster(paste0(relPath, footprintPath))
footprint <- crop(footprint, current)
footprint <- projectRaster(footprint, current)
footprint_resampled <- resample(footprint, current)

# Calculate footprint and richness correlation
cor.test(values(footprint_resampled), values(current))
cor.test(values(footprint_resampled), values(diff70))
cor.test(values(footprint_resampled), values(diffrel70))

# Plot uncertainty map
uncMap <- raster(paste0(relPath, uncPath))
uncMap <- ratify(uncMap)
lc_levels <- c("Correlation filter", "Model complexity", "Projection distance", "Threshold", "Variable selection")
rat <- levels(uncMap)[[1]]
rat <- cbind(rat, lc_levels)
levels(uncMap) <- rat
png(file = paste0(relPath, figPath, "UncertaintyRichness.png"))
levelplot(uncMap, xlab = NULL, ylab = NULL, scales = list(draw = FALSE),
          par.settings = rasterTheme(region = c("blue", "green", "orange", "lightblue", "red"), axis.line = list(col = "transparent")),
          colorkey = list(legend.text = list(fontfamily = "mono")))
dev.off()