library(ggplot2)
library(ggpubr)
library(patchwork)
library(terra)
library(dplyr)
library(sf)
library(mapview)
library(tidyr)
library(tidyverse)
library(ggtext)
library(ncdf4)

getwd()
setwd("./Phosphate")
source('./R-scripts/hel-func.R', chdir = TRUE)

output <- "./R-scripts/RFpred-MODIS_TERRA"
dir.create(output)

data_elZrelli_shp <- "./elZrelli2018_shp"
data_elZrelli <- "./elZrelli2018"
data_elKateb_shap <- "./elKateb2018_shp"
# load interpolate data
data_dir_valid_masekd <- "./landsat2013-masked"

output_RF <- "./RF_TERRA"
output_MODLE <- file.path(output_RF, "RF_models")
raster_outRF_pred <- file.path(output_RF, "RF_pred") 

# Modis 
# data_modis <- "./ModisTerra-Chla/MODIST_L3m_CHL_NRT_2022.0-20260310_110439"
data_modis_Chla <- "./MODIS_SEPTEMBER_DATA_Chla"

# MODIS TERRA Chla
files_modis_Chla <-  list.files(data_modis_Chla , pattern = "_a\\.4km\\.nc$")
modisTERRA_chla <- rast(file.path(data_modis_Chla, files_modis_Chla))
dates_modis <- as.Date(substr(basename(files_modis_Chla), 13, 20), format="%Y%m%d")
time(modisTERRA_chla) <- dates_modis
crs(modisTERRA_chla) <- "EPSG:4326" 

# ----------------------------------------------------------------
# phosphat industry locations ------
# ------------------------------------------------
ind_loc <- t(data.frame(  # lon, lat
  "Gabes" = c(10.095439, 33.916444), # since 1979, expaned 1985
                "Sfax" = c(10.723889, 34.702778), # Triple Super Phosphate (TSP), not DAP, since 1952
                "Skhira" = c(10.148889, 34.3475))) # seit 2013 in betireb
colnames(ind_loc) <- c("lon", "lat")

# GABES and Sfax shp
shp_dir <- "C:/Users/maiim/Documents/25-26WS/Umweltfernerkundung/Phosphate/tunesia-shapefiles"
shp_files <- list.files(shp_dir, pattern = "\\.shp$", full.names = TRUE)
world_shp_list <- lapply(shp_files, vect)
tunisia_shp <- world_shp_list[[207]]
gulf_shp <- tunisia_shp[(tunisia_shp$name == "Gabès") | (tunisia_shp$name == "Sfax") | (tunisia_shp$name == "Médenine"), ]

# ---------------------------------------------------------
# load samplelocation data ------------------
# ---------------------------------------------------------
sampleloc_extent <- vect(file.path(data_elKateb_shap, "sampleLocation_extent.shp"))
sampleLoc_extent <- vect(file.path(data_elZrelli_shp, "elZrelli2018_sampleloc_EXTENT-Skhira.shp"))

sampleLoc_points <- vect(file.path(data_elZrelli_shp, "elZrelli2018_sampleloc2.shp"))
plot(sampleLoc_points)
plot(gulf_shp)
plot(sampleLoc_extent, add = TRUE)

# ---------------------------------------
# ---------------------------------------------------------
# crop CHla data  Modis Terra  ------------------
# ---------------------------------------------------------
sampleloc_extent_pj <- project(sampleLoc_extent, crs(modisTERRA_chla))
sampleloc_points_pj <- project(sampleLoc_points, crs(modisTERRA_chla))
indloc_pj <- project(vect(ind_loc, crs = "EPSG:4326"), crs(modisTERRA_chla))
gulf_shp_pj <- project(gulf_shp, crs(modisTERRA_chla))

chla_crop <- crop(modisTERRA_chla, sampleloc_extent_pj)
plot(chla_crop[[1]])
dates <- time(chla_crop)

september_dates <- dates[format(dates, "%m") == "09"]

chla_crop_sep <- chla_crop[[which(dates %in% september_dates)]]
time(chla_crop_sep)
plot(chla_crop_sep[[2]])
nan_prop <- sapply(1:nlyr(chla_crop_sep), function(i) {
  vals <- values(chla_crop_sep[[i]])         # extract raster values
  mean(is.na(vals))                          # proportion of NA
})

chla_crop_sepCLEAN <- chla_crop_sep[[ which(nan_prop <= 0.35) ]]
plot(chla_crop_sepCLEAN)
september_dates_clean <- time(chla_crop_sepCLEAN)

september2013 <- september_dates_clean[format(september_dates_clean, "%y") == "13"]

chla_2013 <- chla_crop_sepCLEAN[[which(september_dates_clean %in% september2013)]]
chla_2013_max10 <- chla_2013

plot(chla_2013[[1]])
plot(gulf_shp_pj, add = TRUE)
plot(sampleloc_extent_pj, add = TRUE)

buffer_points <- buffer(sampleloc_points_pj, width = 1000) # meter
chla_2013_sep_points <- mask(chla_2013, buffer_points)

WPI_interpolate <- rast(file.path(data_dir_valid_masekd, "WPi_interpolate.tif")) 
WPI_interpolate_pj <- project(WPI_interpolate, crs(chla_2013))
chla_2013_masked <- mask(chla_2013_sep_points, resample(WPI_interpolate_pj$var1.pred, chla_2013_sep_points))
plot(chla_2013_masked[[1]])
plot(sampleloc_extent_pj, add = TRUE)
# ----------------------------------------------------------------------
# Predictions mit Klassifkation Model ----------------------
# ----------------------------------------------------------------------
# load model 
RF_WPI_KLASS <- readRDS(file.path(output_MODLE, paste0("rf_model_WPI_R_KLASS1.rds")))
RF_WPI_KLASS <- readRDS(file.path(output_MODLE, paste0("rf_model_WPI_R_KLASS2.rds")))

chla_WPI_masked <- mask(chla_crop_sepCLEAN, resample(WPI_interpolate_pj$var1.pred, chla_crop_sepCLEAN))
plot(chla_WPI_masked[[1]])
# other time steps test -----
years <- format(as.Date(paste0(2000:2025, "-01-01")), "%y")

png(file.path(output, paste0("RF_Klass3_ALL_pred.png")), height = 3000, width = 5000)
par(mfrow=c(7,5), mar=c(1,1,2,1)) 

for(y in 1:length(years)){
    pred_date <- september_dates_clean[format(september_dates_clean, "%y") == years[y]]

    data_pred <- chla_crop_sepCLEAN[[which(september_dates_clean %in% pred_date)]]

    stacked_pred <- func_median_mean(data_pred)
    # plot(stacked_pred)

    pred_ALL <- predict(stacked_pred, RF_WPI_KLASS)
    pred_ALL <- mask(pred_ALL, stacked_pred$median)
    levels(pred_ALL) <- data.frame(ID=1:5, label=labs)

    #png(file.path(output, paste0("RF_Klass2_pred_20",years[y],".png")), height = 400, width = 800)
    
    plot(pred_ALL, main = paste0("20", years[y], " WPI (RF-KLASS1)"))
    plot(gulf_shp_pj, add = TRUE)
    points(sampleloc_points_pj, col = "cyan")
    points(indloc_pj, col = "red", pch = 15, cex = 1)
    # dev.off()

      rm(pred_ALL, stacked_pred)
  gc()


}
dev.off()
