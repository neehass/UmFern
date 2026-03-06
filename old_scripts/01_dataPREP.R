# DataPREP
# Produktion info of DAP (1998-2017)

# drei Abschnitte 
# vor 1998: Produktion steigt an, aber es gibt keine Datenpunkte
# 1998-2010: Produktion steigt weiter an
# 2011-2017: Produktion sinkt, 2011 einsturz durch Revolution

# erstmal konzentration auf Januar und Juli to reduce seasonal effect 
# in July microbal acitivity / breakdown higher

# packages
library(terra)
library(ggplot2)
library(dplyr)
library(sf)
library(mapview)

getwd()
setwd("./Phosphate")
source('./R-scripts/hel-func.R', chdir = TRUE)

output <- "./R-scripts/outputs"
dir.create(output)

data_dir <- "./landsat_daten"
data_crop_dir <- "./landsat_daten_crop"
dir.create(data_crop_dir)

# ----------------------------------------------------------------
# phosphat industry locations ------
# ------------------------------------------------
ind_loc <- t(data.frame(  # lon, lat
  "Gabes" = c(10.095439, 33.916444), # since 1979, expaned 1985
                "Sfax" = c(10.723889, 34.702778), # Triple Super Phosphate (TSP), not DAP, since 1952
                "Skhira" = c(10.148889, 34.3475))) # seit 2013 in betireb
colnames(ind_loc) <- c("lon", "lat")

# ------------------------------------------------
# load data  ------
files_band <- list.files(data_dir, recursive = TRUE, pattern="T1_SR_B[0-9]{1,2}\\.TIF$")
parts <- strsplit(files_band, "_")
satNR <- unique(sapply(parts, `[`, 1))
dates <- sapply(parts, `[`, 4)
unique_dates <- unique(dates)
# band overview:
satNR_bands <- list(
   "LT05" = list("B1" = "Blue", "B2" = "Green", "B3" = "Red", "B4" = "NIR", "B5" = "SWIR1", "B6" = "Thermal", "B7" = "SWIR2"),
   "LE07" =  list("B1" = "Blue", "B2" = "Green", "B3" = "Red", "B4" = "NIR", "B5" = "SWIR1", "B6" = "Thermal", "B7" = "SWIR2", "B8" = "Panchromatic"),
   "LC08" =  list("B1" = "CoastalAerosol", "B2" = "Blue", "B3" = "Green", "B4" = "Red", "B5" = "NIR", "B6" = "SWIR1", "B7" = "SWIR2", "B8" = "Panchromatic", "B9" = "TIRS1", "B10" = "TIRS2"),
   "LC09" =  list("B1" = "CoastalAerosol", "B2" = "Blue", "B3" = "Green", "B4" = "Red", "B5" = "NIR", "B6" = "SWIR1", "B7" = "SWIR2", "B8" = "Panchromatic", "B9" = "TIRS1", "B10" = "TIRS2")
)
# -------------------------------------------------
# 1. time step < 1998 -----------
# -------------------------------------------------
timestep1 <- 1984:1997

rasters_timestep1_july <- func_timestep_sel(timestep1, files_band, satNR_bands, unique_dates, data_dir, time = "JULY") # help-func # 20 scenens

# inspect JULY scenes timestep1
plot(rasters_timestep1_july[[8]]) # 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18
rasters_timestep1_july <- rasters_timestep1_july[8:18]

# select GABES location
pt <- vect(matrix(ind_loc, ncol=2), crs="EPSG:4326")
ind_loc_sf <- project(pt, crs(rasters_timestep1_july[[1]]))
# crs(rasters_timestep1_july[[1]])
# crs(ind_loc_sf)

# crop rasters to ind_loc
crop_timestep1_JULY <- crop_func(rasters_timestep1_july, ind_loc_sf) # help-func 
# strukur = [[raster scene]][[industry location]][[bands]]
# dev.off()
# plot(crop_timestep1_JULY[[1]][[1]]$B1)
# points(ind_loc_sf[1], col="red", pch=4, cex=1)
# text(ind_loc_sf[1], labels=rownames(ind_loc)[1], pos=3  , col="red", cex=2)

# mean für GABES über alle Szenen
mean_timestep1_JULY <- func_mean_loc(crop_timestep1_JULY) 
plot(mean_timestep1_JULY[[1]])

# std 
std_timestep1_JULY <- func_std_loc(crop_timestep1_JULY) 
col_fun <- colorRampPalette(c("white", "red"))
plot(std_timestep1_JULY[[1]], col = col_fun(100)) 
# B3 (Green) Hohe SD → Vegetation oder Wasseränderung
# Weiß → stabil
# Hellrot → moderate Veränderung
# Dunkelrot → starke zeitliche Dynamik

# coefficent of Variation (CV) = std/mean
cv_timestep1_JULY <- func_CV_loc(std_timestep1_JULY, mean_timestep1_JULY)
plot(cv_timestep1_JULY[[1]],  col = col_fun(100), range = c(0, 100))

# filter files for timestep1 - WINTER ---------------------------------------
rasters_timestep1_winter <- func_timestep_sel(timestep1, files_band, satNR_bands, unique_dates, data_dir, time = "WINTER") # help-func # 6 scenens

# inspect WINTER scenes timestep1
plot(rasters_timestep1_winter[[5]]) # 2, 3, 4, 5
rasters_timestep1_winter <- rasters_timestep1_winter[2:5]

# crop rasters to ind_loc
crop_timestep1_WINTER <- crop_func(rasters_timestep1_winter, ind_loc_sf) # help-func 

# mean für GABES über alle Szenen
mean_timestep1_WINTER <- func_mean_loc(crop_timestep1_WINTER) 
plot(mean_timestep1_WINTER[[1]])

# std 
std_timestep1_WINTER <- func_std_loc(crop_timestep1_WINTER) 
col_fun <- colorRampPalette(c("white", "red"))
plot(std_timestep1_WINTER[[1]], col = col_fun(100)) 

# coefficent of Variation (CV) = std/mean
cv_timestep1_WINTER <- func_CV_loc(std_timestep1_WINTER, mean_timestep1_WINTER)
plot(cv_timestep1_WINTER[[1]],  col = col_fun(100), range = c(0, 100))

# save mean croped raster 
timestep1_dir <- file.path(data_crop_dir, "timestep1")
dir.create(timestep1_dir)
func_save(mean_timestep1_JULY, std_timestep1_JULY, cv_timestep1_JULY, 
    mean_timestep1_WINTER, std_timestep1_WINTER, cv_timestep1_WINTER, 
    ind_loc, timestep1_dir)

# -------------------------------------------------
# 2. time step 1998 - 2010 ----------
# -------------------------------------------------
timestep2 <- 1998:2010

rasters_timestep2_july <- func_timestep_sel(timestep2, files_band, satNR_bands, unique_dates, data_dir, time = "JULY") # help-func # 22 scenen
rasters_timestep2_winter <- func_timestep_sel(timestep2, files_band, satNR_bands, unique_dates, data_dir, time = "WINTER") # help-func # 22 scenens

# summer
plot(rasters_timestep2_july[[9]]) # 1, 9:22
rasters_timestep2_july <- rasters_timestep2_july[c(1, 9:22)]
crop_timestep2_JULY <- crop_func(rasters_timestep2_july, ind_loc_sf) # help-func 

# winter 
plot(rasters_timestep2_winter[[1]]) # 1:7, 16, 17, 18, 19, 20, 21, 22
rasters_timestep2_winter <- rasters_timestep2_winter[c(1:7, 16:22)]
crop_timestep2_WINTER <- crop_func(rasters_timestep2_winter, ind_loc_sf) # help-func 

# mean 
mean_timestep2_JULY <- func_mean_loc(crop_timestep2_JULY)
mean_timestep2_WINTER <- func_mean_loc(crop_timestep2_WINTER)
plot(mean_timestep2_JULY[[1]])
plot(mean_timestep2_WINTER[[1]])

# std 
std_timestep2_JULY <- func_std_loc(crop_timestep2_JULY) 
std_timestep2_WINTER <- func_std_loc(crop_timestep2_WINTER) 

# coefficent of Variation (CV) = std/mean in %
cv_timestep2_JULY <- func_CV_loc(std_timestep2_JULY, mean_timestep2_JULY)
cv_timestep2_WINTER <- func_CV_loc(std_timestep2_WINTER, mean_timestep2_WINTER)
col_fun <- colorRampPalette(c("white", "red"))
plot(cv_timestep2_JULY[[1]],  col = col_fun(100), range = c(0, 100))
plot(cv_timestep2_WINTER[[1]],  col = col_fun(100), range = c(0, 100))

# speicher gecropted rasters timestep2
timestep2_dir <- file.path(data_crop_dir, "timestep2")
dir.create(timestep2_dir)
func_save(mean_timestep2_JULY, std_timestep2_JULY, cv_timestep2_JULY, 
    mean_timestep2_WINTER, std_timestep2_WINTER, cv_timestep2_WINTER, 
    ind_loc, timestep2_dir)

# -------------------------------------------------
# 3. time step 2011 - 2017  --------
# -------------------------------------------------
# juli and august added
timestep3 <- 2011:2017
rasters_timestep3_july <- func_timestep_sel(timestep3, files_band, satNR_bands, unique_dates, data_dir, time = "SUMMER") # help-func # 10 scenen
rasters_timestep3_winter <- func_timestep_sel(timestep3, files_band, satNR_bands, unique_dates, data_dir, time = "WINTER") # help-func # 13

# SUMMER
plot(rasters_timestep3_july[[1]]) # 1:7
rasters_timestep3_july <- rasters_timestep3_july[1:7]
crop_timestep3_JULY <- crop_func(rasters_timestep3_july, ind_loc_sf) # help-func 
# strukur = [[raster scene]][[industry location]][[bands]]

# WINTER
plot(rasters_timestep3_winter[[13]]) # 1:7, 13
rasters_timestep3_winter <- rasters_timestep3_winter[c(1:7, 13)]
crop_timestep3_WINTER <- crop_func(rasters_timestep3_winter, ind_loc_sf) # help-func 

# mean 
mean_timestep3_JULY <- func_mean_loc(crop_timestep3_JULY)
mean_timestep3_WINTER <- func_mean_loc(crop_timestep3_WINTER)
plot(mean_timestep3_JULY[[1]])
plot(mean_timestep3_WINTER[[1]])

# std 
std_timestep3_JULY <- func_std_loc(crop_timestep3_JULY) 
std_timestep3_WINTER <- func_std_loc(crop_timestep3_WINTER) 

# coefficent of Variation (CV) = std/mean in %
cv_timestep3_JULY <- func_CV_loc(std_timestep3_JULY, mean_timestep3_JULY)
cv_timestep3_WINTER <- func_CV_loc(std_timestep3_WINTER, mean_timestep3_WINTER)
col_fun <- colorRampPalette(c("white", "red"))
plot(cv_timestep3_JULY[[1]],  col = col_fun(100), range = c(0, 100))
plot(cv_timestep3_WINTER[[1]],  col = col_fun(100), range = c(0, 100))

# speicher gecropted rasters timestep3
timestep3_dir <- file.path(data_crop_dir, "timestep3")
dir.create(timestep3_dir)
func_save(mean_timestep3_JULY, std_timestep3_JULY, cv_timestep3_JULY, 
    mean_timestep3_WINTER, std_timestep3_WINTER, cv_timestep3_WINTER, 
    ind_loc, timestep3_dir)

# -------------------------------------------------
# 4. time step 2018 - 2026 --------
# -------------------------------------------------
timestep4 <- 2018:2026
rasters_timestep4_july <- func_timestep_sel(timestep4, files_band, satNR_bands, unique_dates, data_dir, time = "JULY") # help-func # 9 scenen
rasters_timestep4_winter <- func_timestep_sel(timestep4, files_band, satNR_bands, unique_dates, data_dir, time = "WINTER") # help-func # 11 

plot(rasters_timestep4_july[[9]]) # 1:6, 8:9
rasters_timestep4_july <- rasters_timestep4_july[c(1:6, 8:9)]
crop_timestep4_JULY <- crop_func(rasters_timestep4_july, ind_loc_sf) # help-func 

plot(rasters_timestep4_winter[[11]]) # 1:11
rasters_timestep4_winter <- rasters_timestep4_winter[1:11]
crop_timestep4_WINTER <- crop_func(rasters_timestep4_winter, ind_loc_sf) # help-func 

# mean 
mean_timestep4_JULY <- func_mean_loc(crop_timestep4_JULY)
mean_timestep4_WINTER <- func_mean_loc(crop_timestep4_WINTER)
plot(mean_timestep4_JULY[[1]])
plot(mean_timestep4_WINTER[[1]])

# std 
std_timestep4_JULY <- func_std_loc(crop_timestep4_JULY) 
std_timestep4_WINTER <- func_std_loc(crop_timestep4_WINTER) 

# coefficent of Variation (CV) = std/mean in %
cv_timestep4_JULY <- func_CV_loc(std_timestep4_JULY, mean_timestep4_JULY)
cv_timestep4_WINTER <- func_CV_loc(std_timestep4_WINTER, mean_timestep4_WINTER)
col_fun <- colorRampPalette(c("white", "red"))
plot(cv_timestep4_JULY[[1]],  col = col_fun(100), range = c(0, 100))
plot(cv_timestep4_WINTER[[1]],  col = col_fun(100), range = c(0, 100))

# speicher gecropted rasters timestep4
timestep4_dir <- file.path(data_crop_dir, "timestep4")
dir.create(timestep4_dir)
func_save(mean_timestep4_JULY, std_timestep4_JULY, cv_timestep4_JULY, 
    mean_timestep4_WINTER, std_timestep4_WINTER, cv_timestep4_WINTER, 
    ind_loc, timestep4_dir)
# -------------------------------------------------------------------------------------------------------
# mean per year (1984-2026)
# mean im zwei Jahres rythmus
# -------------------------------------------------------------------------------------------------------
data_dir_year <- "./landsat_daten_year"
dir.create(data_dir_year)
years <- 1984:2026 

y <- 43 ## 43 Jahre time periode
year <- years[y]
# SUMMER
raster_yearJULY <- func_timestep_sel(year, files_band, satNR_bands, unique_dates, data_dir, time = "JULY")
raster_yearWINTER <- func_timestep_sel(year, files_band, satNR_bands, unique_dates, data_dir, time = "WINTER")
# SUMMER
plot(raster_yearJULY[[2]]) 
raster_yearJULY <- raster_yearJULY # [c(2)]
crop_yearJULY <- crop_func(raster_yearJULY, ind_loc_sf) # help-func 

# WINTER
plot(raster_yearWINTER[[3]]) 
raster_yearWINTER <- raster_yearWINTER # [c(1, 2)] #[c(1,2,5, 6)]
crop_yearWINTER <- crop_func(raster_yearWINTER, ind_loc_sf) # help-func 

dir_year <- file.path(data_dir_year, year)
dir.create(dir_year)
if(length(crop_yearJULY) > 1 & length(crop_yearWINTER) > 1){
   # SUMMER 
   mean_year_JULY <- func_mean_loc(crop_yearJULY)
   std_years_JULY <- func_std_loc(crop_yearJULY) 
   cv_year_JULY <- func_CV_loc(std_years_JULY, mean_year_JULY)

    # Winter 
   mean_year_WINTER <- func_mean_loc(crop_yearWINTER)
   std_years_WINTER <- func_std_loc(crop_yearWINTER) 
   cv_year_WINTER <- func_CV_loc(std_years_WINTER, mean_year_WINTER)

   # save
   func_save(mean_year_JULY, std_years_JULY, cv_year_JULY, 
    mean_year_WINTER, std_years_WINTER, cv_year_WINTER, 
    ind_loc, dir_year)

} else {print("save manuel")}

# multiple scenes just for one month
# SUMMER 
mean_year_JULY <- func_mean_loc(crop_yearJULY)
std_years_JULY <- func_std_loc(crop_yearJULY) 
cv_year_JULY <- func_CV_loc(std_years_JULY, mean_year_JULY)
for(loc in 1:3){
    print(paste("Saving rasters for location:", rownames(ind_loc)[loc]))
        # SUMMER 
        r <- mean_year_JULY[[loc]]
        raster_name <- paste0("mean_",year,"_JULY_", rownames(ind_loc)[loc], ".tif")
        writeRaster(r, file.path(dir_year, raster_name), overwrite=TRUE)

         r <- std_years_JULY[[loc]]
        raster_name <- paste0("std_",year,"_JULY_", rownames(ind_loc)[loc], ".tif")
        writeRaster(r, file.path(dir_year, raster_name), overwrite=TRUE)

         r <- cv_year_JULY[[loc]]
        raster_name <- paste0("cv_",year,"_JULY_", rownames(ind_loc)[loc], ".tif")
        writeRaster(r, file.path(dir_year, raster_name), overwrite=TRUE)
}
# Winter 
mean_year_WINTER <- func_mean_loc(crop_yearWINTER)
std_years_WINTER <- func_std_loc(crop_yearWINTER) 
cv_year_WINTER <- func_CV_loc(std_years_WINTER, mean_year_WINTER)
for(loc in 1:3){
    print(paste("Saving rasters for location:", rownames(ind_loc)[loc]))
        # WINTER 
        r <- mean_year_WINTER[[loc]]
        raster_name <- paste0("mean_",year,"_WINTER_", rownames(ind_loc)[loc], ".tif")
        writeRaster(r, file.path(dir_year, raster_name), overwrite=TRUE)

         r <- std_years_WINTER[[loc]]
        raster_name <- paste0("std_",year,"_WINTER_", rownames(ind_loc)[loc], ".tif")
        writeRaster(r, file.path(dir_year, raster_name), overwrite=TRUE)

         r <- cv_year_WINTER[[loc]]
        raster_name <- paste0("cv_",year,"_WINTER_", rownames(ind_loc)[loc], ".tif")
        writeRaster(r, file.path(dir_year, raster_name), overwrite=TRUE)
}

# just 1 scene
for(loc in 1:3){
    print(paste("Saving rasters for location:", rownames(ind_loc)[loc]))
        # SUMMER 
        r <- crop_yearJULY[[1]][[loc]]
        raster_name <- paste0("mean_",year,"_JULY_", rownames(ind_loc)[loc], ".tif")
        writeRaster(r, file.path(dir_year, raster_name), overwrite=TRUE)
}
for(loc in 1:3){
    print(paste("Saving rasters for location:", rownames(ind_loc)[loc]))
        # WINTER 
        r <- crop_yearWINTER[[1]][[loc]]
        raster_name <- paste0("mean_",year,"_WINTER_", rownames(ind_loc)[loc], ".tif")
        writeRaster(r, file.path(dir_year, raster_name), overwrite=TRUE)
}


# -------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------
# water selection mit modified Normalized Difference Water Index (MNDWI)
# MNDWI= (Green-SWIR1)/(Green+SWIR1)

MNDWI <- func_MNDWI(GREEN = crop_timestep1_JULY[[1]][[3]]$B2, SWIR1 = crop_timestep1_JULY[[1]][[3]]$B5, threshold = -0.01) # help-func

# plot(MNDWI$mndwi, main = "MNDWI - July 2018-2026 - Scene 1")
water_only <- mask(crop_timestep1_JULY[[1]][[3]], MNDWI$mask_water, maskvalue=0)
plot(water_only)

SPM <- func_SPM(R_665 = water_only$B3, R_490 = water_only$B1) # help-func
plot(SPM, main = "SPM")

SPM <- func_SPM(R_665 = crop_timestep1_JULY[[1]][[3]]$B3, R_490 = crop_timestep1_JULY[[1]][[3]]$B1) # help-func
plot(SPM, main = "SPM")


MNDWI <- func_MNDWI(GREEN = crop_timestep4_JULY[[1]][[3]]$B3, SWIR1 = crop_timestep4_JULY[[1]][[3]]$B6, threshold = -0.01) # help-func
# plot(MNDWI$mndwi, main = "MNDWI - July 2018-2026 - Scene 1")
water_only <- mask(crop_timestep4_JULY[[1]][[3]], MNDWI$mask_water, maskvalue=0)
plot(water_only)

SPM <- func_SPM(R_665 = water_only$B3, R_490 = water_only$B1) # help-func
plot(SPM, main = "SPM")
