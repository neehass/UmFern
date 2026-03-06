# notes: methods
# -----------------------------------------
# El Khalil Cherif, 2019 
# Sentinel-3 reflectances at bands 6 (560 nm) and 7 (620 nm). 
# Since algorithms for Chlorophyll-a and Total Suspended Matter Concentrations 
# rely on these wavelengths, they too were found to be closely related to high As, Fe and P concentrations
# -----------------------------------------

# formulas:
# Suspenden partical matter (SPM), R = reflectance
# SPM = R(665 nm)/ R(490)
# SPM = a + b (R(red)/ R(blue))
# Wu 2025, Where a,b are empirical calibration coefficients fit using in-situ SPM data.

# Chlorophyl-a (Ch-a)
# Ch-a = a+b (R(710 nm)/ R(670))
# R(710) = reflectance peak from vegetation fluorescence.
# R(670) = absorption trough from Chl-a.
# a,b are calibrated with field data.

# landsat 4/5 (micrometer, resolustion in meters)
# Band 1 - Blue	0.45-0.52	30 <<<<<
# Band 2 - Green	0.52-0.60	30 
# Band 3 - Red	0.63-0.69	30 <<<<<<<<<<<<<<<<<<<<<
# Band 4 - Near Infrared (NIR)	0.76-0.90	30
# Band 5 - Shortwave Infrared (SWIR) 1	1.55-1.75	30
# Band 6 - Thermal	10.40-12.50	120 (resampled to 30)*
# Band 7 - Shortwave Infrared (SWIR) 2	2.08-2.35	30

# ----------------------------------------

library(terra)
library(ggplot2)
library(dplyr)
library(sf)
library(mapview)

getwd()
setwd("./Phosphate")
output <- "./R-scripts/outputs"
dir.create(output)

data_dir <- "./landsat_daten"

files_band <- list.files(data_dir, recursive = TRUE, pattern="T1_SR_B[0-9]{1,2}\\.TIF$")

# rastlist1984 <- 

parts <- strsplit(files_band, "_")

# Extract acquisition date (4th element)
dates <- sapply(parts, `[`, 4)

unique_dates <- unique(dates)

files_1984 <- files_band[grep("1984", files_band)]
dates_1984 <- unique_dates[grep("1984", unique_dates)]
scenes <- length(dates_1984)

parts <- strsplit(files_1984, "_")
scene_id <- sapply(parts, function(x) paste(x[1:4], collapse = "_"))

unique_scenes <- unique(scene_id)

rasters1989 <- lapply(unique_scenes, function(sc) {
  scene_files <- files_1984[scene_id == sc]
  rast(file.path(data_dir, scene_files))
})

windows()
plot(rasters1989[[1]], main = "Landsat 4/5 - Scene 1") # <<  relevante scene 
plot(rasters1989[[2]], main = "Landsat 4/5 - Scene 2")
# für 1986 = scene 6, für 1989 = scene 2

# relevante scene
coast <- rasters1989[[1]]

parts <- strsplit(names(coast), "_")
bands <- sapply(parts, `[`, length(parts[[1]]))
names(coast) <- bands

windows()
plot(coast$B2, main = "Band 2 - Green")
plot(coast$B3, main = "Band 3 - Red")

# Bounding Box definieren (xmin, xmax, ymin, ymax)
ext(coast)
gabes_ext <- ext(590000, 703515, 3724485, 3900000)

coast_crop <- crop(coast, gabes_ext)

plot(coast_crop$B2)

# phosphat industry locations
ind_loc <- t(data.frame(  # lon, lat
  "Gabes" = c(10.095439, 33.916444), #since 1979, expaned 1985
                "Sfax" = c(10.723889, 34.702778), # Triple Super Phosphate (TSP), not DAP, since 1952
                "Skhira" = c(10.148889, 34.3475))) # seit 2013 in betireb
colnames(ind_loc) <- c("lon", "lat")

# transforming to sf object
pt <- vect(matrix(ind_loc, ncol=2), crs="EPSG:4326")
ind_loc_sf <- project(pt, crs(coast))

crs(coast)
crs(ind_loc_sf)

# plot the cropped raster and overlay the industry locations
plot(coast_crop$B2)
points(ind_loc_sf, col="red", pch=4, cex=1)
text(ind_loc_sf, labels=rownames(ind_loc), pos=3, col="red", cex=0.8)

# plot per industry
crop_ind <- lapply(1:nrow(ind_loc_sf), function(i) {
  
  if(i == 3){
    ind_point <- c(10.148889, 34.3375)
  } else {
    ind_point <- ind_loc_sf[i]
  }
  
  # create buffer (CRS must be projected in meters!)
  ind_buffer <- buffer(ind_point, width = 8000) 
  
  # crop raster using buffer
  crop(coast_crop, ind_buffer)
})

windows()
png(filename = file.path(output, "overviewIndustry.png"), width = 1200, height = 400)
par(mfrow=c(1,3))
# Gabes
# plot(crop_ind[[1]]$B2, main = "Gabes - Band 2")
plotRGB(crop_ind[[1]], r = 3, g = 2, b = 1, main = "Gabes 11.07.1984 - RGB", stretch = "lin") # 19840711
points(ind_loc_sf[1], col="red", pch=4, cex=1)
text(ind_loc_sf[1], labels=rownames(ind_loc)[1], pos=3  , col="red", cex=2)
#Sfax
# plot(crop_ind[[2]]$B2, main = "Sfax - Band 2")
plotRGB(crop_ind[[2]], r = 3, g = 2, b = 1, main = "Sfax 11.07.1984 - RGB", stretch = "lin")
points(ind_loc_sf[2], col="red", pch=4, cex=1)
text(ind_loc_sf[2], labels=rownames(ind_loc)[2], pos=3  , col="red", cex=2)
# Skhira
# plot(crop_ind[[3]]$B2, main = "Skhira - Band 2")
plotRGB(crop_ind[[3]], r = 3, g = 2, b = 1, main = "Skhira 11.07.1984 - RGB", stretch = "lin")
points(ind_loc_sf[3], col="red", pch=4, cex=1)
text(ind_loc_sf[3], labels=rownames(ind_loc)[3], pos=3  , col="red", cex=2)
dev.off()

# ----------------------------------------------------------------
# ----------------------------------------------------------------
# ----------------------------------------------------------------

# DAP Produktion
dap_prod <- read.csv("C:/Users/maiim/Documents/25-26WS/Umweltfernerkundung/Phosphate/Production_DAP_tunisia.csv", sep = ";", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
dap_prod <- as.data.frame(t(dap_prod))
dap_prod$Year <- as.numeric(gsub("X", "", rownames(dap_prod)))
rownames(dap_prod) <- NULL 
colnames(dap_prod) <- c("DAP", "Year")
dap_prod$DAP <- as.numeric(dap_prod$DAP)
View(dap_prod)
mean98_2010 <- mean(dap_prod$DAP[dap_prod$Year >= 1998 & dap_prod$Year <= 2010])
year98_2010 <- dap_prod$Year[dap_prod$Year >= 1998 & dap_prod$Year <= 2010]

mean2011_2017 <- mean(dap_prod$DAP[dap_prod$Year > 2010 & dap_prod$Year <= 2017])
year2011_2017 <- dap_prod$Year[dap_prod$Year > 2010 & dap_prod$Year <= 2017]

dev.off()   # reset device (run until "null device")
png(filename = file.path(output, "DAP_production.png"), width = 800, height = 600)
plot(dap_prod$Year, dap_prod$DAP, type = "l", col = "black",
 main = "DAP Production in Gabes\nin 1000t", ylab = "Production of DAP (1000t)", xlab = "Year",  xaxt = "n",)
grid() 
axis(1, at = dap_prod$Year, labels = dap_prod$Year, las = 2)
lines(year98_2010, rep(mean98_2010, length(year98_2010)), col = "blue")
lines(year2011_2017, rep(mean2011_2017, length(year2011_2017)), col = "red")

usr <- par("usr")  # get plot limits
rect(xleft = 1998, ybottom = usr[3],   xright = 2010, ytop = usr[4],  col = rgb(173/255, 216/255, 230/255, 0.4), border = NA)
abline(v = c(1998, 2010), col="blue", lty=2)

rect(xleft = 2011, ybottom = usr[3],   xright = 2017, ytop = usr[4],  col = rgb(255/255, 0/255, 0/255, 0.2), border = NA)
abline(v = c(2011, 2017), col="red", lty=2)

legend(2012, max(dap_prod$DAP), legend = c("DAP", "Mean 1998-2010", "Mean 2011-2017"), col = c("black", "blue", "red"), 
  lty = 1, bty = "n", bg = "white")
dev.off()






