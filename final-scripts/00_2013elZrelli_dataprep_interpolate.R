# ------------------------------------------------------------------------
# Interpolation der Messungen 
# elZrelli is vaildation data, Messungen September 2013
# >> räumliche Verteilungen
# ------------------------------------------------------------------------

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

getwd()
setwd("./Phosphate")
source('./R-scripts/hel-func.R', chdir = TRUE)

output <- "./R-scripts/output/00_2013elZrelli_dataprep_interpolate"
dir.create(output, recursive = TRUE)

data_elZrelli_shp <- "./elZrelli2018_shp"
data_elZrelli <- "./elZrelli2018"

data_pif <- "./landsat-SEPTEMBER_PIF_LC08-09"

data_dir_valid_masekd <- "./landsat-interpolation"
dir.create(data_dir_valid_masekd)
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
sampleloc_extent_interpol <- vect(file.path(data_elZrelli_shp, "elZrelli2018_sampleloc_SMALL-Interpol.shp"))
sampleloc_extent_interpol <- erase(sampleloc_extent_interpol, gulf_shp)
sampleLoc_points <- vect(file.path(data_elZrelli_shp, "elZrelli2018_sampleloc2.shp"))

# ---------------------------------------------------------
# load elZrelli 2018, sep 2013 ------------------
# ---------------------------------------------------------
# load el Kateb 2018 data
# list.files(data_elZrelli)
# hevymetal_sl <- read.csv(file.path(data_elZrelli, "elZrelli2018_heavymetal.csv"), sep = ";", header = TRUE)

# split_F <- do.call(rbind, strsplit(hevymetal_sl$F..mg.L..., "< | ± "))
# hevymetal$F_mean <- as.numeric(split_F[,1])
# hevymetal$F_sd   <- as.numeric(split_F[,2])

# split_F <- do.call(rbind, strsplit(hevymetal$Cd..µg.L..., "< | ± "))
# hevymetal$Cd_sd_thres <- as.numeric(split_F[,2])
# hevymetal$Cd_mean <- as.numeric(split_F[,2])

# colnames(hevymetal_sel)
# hevymetal_sel <- hevymetal[,!names(hevymetal) %in% c("P_sd") ]
# hevymetal <- hevymetal_sel[, c(1:3, 11:23)]

# write.csv(hevymetal, file.path(data_elZrelli, "elZrelli2018_heavymetal_SORTED.csv"))

# heavymetal <-  read.csv(file.path(data_elZrelli, "elZrelli2018_heavymetal_SORTED.csv"), sep = ",", header = TRUE, row.names = 1)
# waterPolIndex <-  read.csv(file.path(data_elZrelli, "elZrelli2018_waterPollutionIndex.csv"), sep = ";", header = TRUE)

# sample_pints_values <- values(sampleLoc_points)
# sample_pints_values <- sample_pints_values[,names(sample_pints_values) %in% c("Sector", "Station", "Latitude", "Longitude") ]

# # el Zrelli 2018, sep 2013 Messungen
# heavymetal <- merge(heavymetal, sample_pints_values, by = "Station")
# colnames(waterPolIndex) <- c("Station" ,  "WPi", "Seawater.Quality.Classification")
# waterPolIndex <- merge(waterPolIndex, sample_pints_values, by = "Station")

# write.csv(heavymetal, file.path(data_elZrelli, "elZrelli2018_heavymetal_POINTS.csv"))
# write.csv(waterPolIndex, file.path(data_elZrelli, "elZrelli2018_waterPollutionIndex_POINTS.csv"))

heavymetal <-  read.csv(file.path(data_elZrelli, "elZrelli2018_heavymetal_POINTS.csv"), sep = ",", header = TRUE, row.names = 1)
waterPolIndex <-  read.csv(file.path(data_elZrelli, "elZrelli2018_waterPollutionIndex_POINTS.csv"), sep = ",", header = TRUE)

# ----------------------------------------------------------------------
# Interpolate Messungen >> räumlcihe Verteilung
# ----------------------------------------------------------------------

# # auf Ausdehnung der Maske zuschneiden
# plot(sampleloc_extent_interpol)

# WPI
wpi_interpolate_forward <- func_interpolation_forward(heavymetal = waterPolIndex, VAR = "WPi", maske = sampleloc_extent_interpol, 
        name = "WPi_interpolate_forwarde", data_dir_valid_masekd) # geringere Auflösung
gulf_shp_pj <- project(gulf_shp, wpi_interpolate_forward)
sampleLoc_points_pj <- project(sampleLoc_points, wpi_interpolate_forward)
indloc_pj <- project(vect(matrix(ind_loc, ncol=2), crs="EPSG:4326"), wpi_interpolate_forward)

brks <- c( 0,1,2,3,4, 6, 10, 15, 30)#seq(0.3, 20,  by = 4)
png(file.path(output,  paste0("wpi_interpolate_forward.png")), height = 1000, width = 1000)
plot(wpi_interpolate_forward$var1.pred, plg = list(title = "WPI"), , main = "WPI", breaks = brks)
plot(gulf_shp_pj, add = TRUE)
points(sampleLoc_points_pj, col = "red")
points(indloc_pj, col = "black", pch = 15, cex = 1)
dev.off()

wpi_interpolate <- func_interpoltate(waterPolIndex, VAR = "WPi", maske = sampleloc_extent_interpol, name = "WPi_interpolate", data_dir_valid_masekd)
brks <- c( 0,1,2,3,4, 6, 10, 15, 30)#seq(0.3, 20,  by = 4)

png(file.path(output,  paste0("wpi_interpolate.png")), height = 1000, width = 1000)
plot(wpi_interpolate$var1.pred, plg = list(title = "WPI"), , main = "WPI", breaks = brks)
plot(gulf_shp_pj, add = TRUE)
points(sampleLoc_points_pj, col = "red")
points(indloc_pj, col = "black", pch = 15, cex = 1)
dev.off()

# F
F_interpolate <- func_interpoltate(heavymetal, VAR = "F_mean",  maske = sampleloc_extent_interpol , name = "F_interpolate", data_dir_valid_masekd)
brks <- c( 0,0.3,5,10,15, 20)#seq(0.3, 20,  by = 4)
brks <- c( 0,1,2,3,4,5,10,15,20)
png(file.path(output,  paste0("F_interpolate",".png")), height = 1000, width = 1000)
plot(F_interpolate$var1.pred, plg = list(title = "[mg/l]"), 
  breaks = brks, main = "F")
plot(gulf_shp_pj, add = TRUE)
points(sampleLoc_points_pj, col = "red")
points(indloc_pj, col = "black", pch = 15, cex = 1)
dev.off()

# P
lm_data_p <- heavymetal
lm_data_p$P_mean[is.na(lm_data_p$P_mean)] <- 130
p_interpolate <- func_interpoltate(lm_data_p, VAR = "P_mean", maske = sampleloc_extent_interpol, name = "P_interpolate", data_dir_valid_masekd)
plot(p_interpolate)
brks <- c( 0, 0.5,1,2,3,4, 5, 6)#seq(0.3, 20,  by = 4)

png(file.path(output,  paste0("P_interpolate",".png")), height = 1000, width = 1000)
plot(p_interpolate$var1.pred / 1000, plg = list(title = "[mg/l]"), 
  breaks = brks, main = "P")
plot(gulf_shp_pj, add = TRUE)
points(sampleloc_points_pj, col = "red")
points(indloc_pj, col = "black", pch = 15, cex = 1)
dev.off()

# Pb
Pb_interpolate <- func_interpoltate(heavymetal, VAR = "Pb_mean", maske = sampleloc_extent_interpol, name = "Pb_interpolate", data_dir_valid_masekd)
brks <- c(0,0.2,0.4,0.6,0.8,1)
png(file.path(output,  paste0("Pb_interpolate",".png")), height = 1000, width = 1000)
plot(Pb_interpolate$var1.pred, plg = list(title = "[μg/l]"), 
  breaks = brks, main = "Zn")
plot(gulf_shp_pj, add = TRUE)
points(sampleLoc_points_pj, col = "red")
points(indloc_pj, col = "black", pch = 15, cex = 1)
dev.off()

# Zn
Zn_interpolate <- func_interpoltate(heavymetal, VAR = "Zn_mean", maske = sampleloc_extent_interpol, name = "Zn_interpolate", data_dir_valid_masekd)
brks <- c(5,8,10,15,20, 25)
png(file.path(output,  paste0("Zn_interpolate",".png")), height = 1000, width = 1000)
plot(Zn_interpolate$var1.pred, plg = list(title = "[μg/l]"), 
  breaks = brks, main = "Zn")
plot(gulf_shp_pj, add = TRUE)
points(sampleLoc_points_pj, col = "red")
points(indloc_pj, col = "black", pch = 15, cex = 1)
dev.off()

# Cu
Cu_interpolate <- func_interpoltate(heavymetal, VAR = "Cu_mean", maske = sampleloc_extent_interpol, name = "Cu_interpolate", data_dir_valid_masekd)
brks <- c(1,2,3,4,5,6,7,8)
png(file.path(output,  paste0("Cu_interpolate",".png")), height = 1000, width = 1000)
plot(Cu_interpolate$var1.pred, plg = list(title = "[μg/l]"), 
  breaks = brks, main = "Zn")
plot(gulf_shp_pj, add = TRUE)
points(sampleLoc_points_pj, col = "red")
points(indloc_pj, col = "black", pch = 15, cex = 1)
dev.off()

### END ########################################