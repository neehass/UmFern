library(ggplot2)
library(ggpubr)
library(patchwork)
library(terra)
library(dplyr)
library(sf)
library(mapview)

getwd()
setwd("./Phosphate")
source('./R-scripts/hel-func.R', chdir = TRUE)

output <- "./R-scripts/2014Valid-outputs"
dir.create(output)

data_elKateb_shap <- "./elKateb2018_shp"

data_CHla <- "./NOAAMSL12_chla/chlora/monthly/WW00"
data_nlw <- "./NOAAMSL12_chla/nlw/monthly/WW00"

files <- list.files(data_CHla)
chla <- rast(file.path(data_CHla, files))
plot(chla[[1]])

# nlw 
files_nlw <- list.files(data_nlw)
nlw <- rast(file.path(data_nlw, files_nlw))
plot(nlw[[1]])
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
sampleLoc_points <- vect(file.path(data_elKateb_shap, "elKateb2018_sampleLocations.shp"))
names(sampleLoc_points) <- c("Station", "Lat", "Lon")

# ---------------------------------------------------------
# crop CHla data  ------------------
# ---------------------------------------------------------
sampleloc_extent_pj <- project(sampleloc_extent, crs(chla))

chla_crop <- crop(chla, sampleloc_extent_pj)
plot(chla_crop)
time(chla_crop)

chla_crop_201408 <- chla_crop[[32]]
plot(chla_crop_201408)
time(chla_crop_201408)
