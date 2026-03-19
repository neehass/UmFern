# use models from validation sep 2013 für predictions 

# packages
library(terra)
library(ggplot2)
library(dplyr)
library(sf)
library(mapview)
library(ranger)

getwd()
setwd("./Phosphate")
source('./R-scripts/final-scripts/hel-func.R', chdir = TRUE)

output <- "./R-scripts/output/03_timeseries_RF_pred_LC08-09"
dir.create(output, recursive = TRUE)

data_dir_cloudmask <- "./landsat-SEPTEMBER_cloudmask"
data_crop_dir <- "./final_landsat-SEPTEMBER_PIF_LC08-09"
# data_crop_dir <- "./2final_landsat-SEPTEMBER_PIF_LC08-09"

folder_MODELS <- file.path("./final_RF_NORM_models_LC08-09")

folder_PRED <- "./final_RF_NORM_pred_LC08-09"

data_elZrelli_shp <- "./elZrelli2018_shp"
data_elZrelli <- "./elZrelli2018"

# windows()
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
sampleloc_extent3 <- vect(file.path(data_elZrelli_shp, "elZrelli2018_sampleloc_EXTENT-Skhira.shp"))
sampleloc_extent3 <- erase(sampleloc_extent3, gulf_shp)
sampleloc_extent4_land <- vect(file.path(data_elZrelli_shp, "elZrelli2018_sampleloc_EXTENT4_land.shp"))
sampleloc_extent4_land <- erase(sampleloc_extent4_land, gulf_shp)

# plot industrie loc and sample locs
# plot(sampleloc_extent3)
# points(ind_loc, col = "red", pch = 15)

# ---------------------------------------------------------
# load RF modles ------------------
# ---------------------------------------------------------
rf_model_WPI_RF <- readRDS(file.path(folder_MODELS, paste0("rf_model_WPI_RF.rds")))
rf_model_WPI_RF_CLASS <- readRDS(file.path(folder_MODELS, paste0("rf_model_WPI_RF_class.rds")))

rf_model_WPI_RF_CLASS_lessVAR <- readRDS(file.path(folder_MODELS, paste0("rf_model_WPI_RF_class_lessVAR.rds")))
rf_model_WPI_RF_CLASS_onlyTSMChla <- readRDS(file.path(folder_MODELS, paste0("rf_model_WPI_RF_class_TSMChla.rds")))

rf_model_F_RF <- readRDS(file.path(folder_MODELS, paste0("rf_model_F_RF.rds")))
rf_model_F_RF_CLASS <- readRDS(file.path(folder_MODELS, paste0("rf_model_F_RF_class.rds")))

# ---------------------------------------------------------
# load landsat data ------------------
# ---------------------------------------------------------
files_raster <- list.files(data_crop_dir, recursive = TRUE, pattern="^WATER_NORMIERT")
parts <- strsplit(files_raster, "_|\\.")
years <- sapply(parts, `[`, 4) 
unique_years <- unique(years)
print(unique_years) # von 2013 bis 2025 >> 12Jahre
length(unique_years) # 10 Scenen 

raster_years <- lapply(files_raster, function(x){rast(file.path(data_crop_dir, x))})
names(raster_years) <- paste0("year_",unique_years)
raster_years <- lapply(raster_years, function(x) {x[[ !names(x) %in% c("CoastalAerosol", "SWIR2") ]]})
length(raster_years)

png(file.path(output, paste0("RAW-Blue", ".png")), height = 800, width = 800)
par(mfrow=c(5,2))
for(y in 1:length(unique_years)){
    # max_abs <- max(abs(global(baseline_diff[[y]], "max", na.rm=TRUE)),
    #             abs(global(baseline_diff[[y]], "min", na.rm=TRUE)))

  plot(raster_years[[y]]$Blue,
       main = paste("Blue", unique_years[y]))
}
dev.off()

# ---------------------------------------------------------
# make predictions with WPI model ------------------
# ---------------------------------------------------------
name_RF <- "WPI_NORM"
folder_PRED_WPI <- file.path(folder_PRED, name_RF)
dir.create(folder_PRED_WPI)

rf_model_WPI_RF$forest$independent.variable.name

all_pred_WPI <- func_pred_RF(raster_years, rf_model_WPI_RF, sampleloc_extent4_land, outpur_dir = folder_PRED_WPI)
all_pred_WPI <- lapply(list.files(folder_PRED_WPI), function(x){rast(file.path(folder_PRED_WPI, x))})

brks <- c( 0,1,2,3,4, 6, 10, 15, 30)#seq(0.3, 20,  by = 4)
png(file.path(output, paste0(name_RF, "_interpolate_RF_pred.png")), height = 800, width = 800)
par(mfrow=c(3,4))
for(y in 1:length(unique_years)){
    # max_abs <- max(abs(global(baseline_diff[[y]], "max", na.rm=TRUE)),
    #             abs(global(baseline_diff[[y]], "min", na.rm=TRUE)))

  plot(all_pred_WPI[[y]],
       main = paste("WPI", unique_years[y]),  breaks = brks)
}
dev.off()

# Klassifikation and NORMIERT ------------------------
name_RF <- "WPI_class_NORM"
folder_PRED_WPI_class <- file.path(folder_PRED, name_RF)
dir.create(folder_PRED_WPI_class)

rf_model_WPI_RF_CLASS$forest$independent.variable.names # variable names

all_pred_WPI_Class <- func_pred_RF(raster_years, rf_model_WPI_RF_CLASS, 
                    sampleloc_extent4_land, outpur_dir = folder_PRED_WPI_class)
# all_pred_WPI_Class_lessVAR <- lapply(list.files(folder_PRED_WPI_class, pattern = ".tif$"), function(x){rast(file.path(folder_PRED_WPI_class, x))})

labs <- c("Not affected", "Slightly affected","Moderately affected","Strongly affected","Seriously affected")
png(file.path(output, paste0(name_RF, "_interpolate_RF_pred.png")), height = 800, width = 800)
par(mfrow=c(3,4))
for(y in 1:length(unique_years)){
    # max_abs <- max(abs(global(baseline_diff[[y]], "max", na.rm=TRUE)),
    #             abs(global(baseline_diff[[y]], "min", na.rm=TRUE)))
       levels(all_pred_WPI_Class[[y]]) <- data.frame(ID=1:5, label=labs)
  plot(all_pred_WPI_Class[[y]],
       main = paste("WPI", unique_years[y]), )
}
dev.off()

# less Variables Klassifikation and NORMIERT ------------------------
name_RF <- "WPI_class_NORM_lessVAR"
folder_PRED_WPI_class_lessVAR <- file.path(folder_PRED, name_RF)
dir.create(folder_PRED_WPI_class_lessVAR)

rf_model_WPI_RF_CLASS_lessVAR$forest$independent.variable.names # variable names
names(raster_years)

# all_pred_WPI_Class_lessVAR <- func_pred_RF(raster_years, rf_model_WPI_RF_CLASS_lessVAR, sampleloc_extent4_land, outpur_dir = folder_PRED_WPI_class_lessVAR)
all_pred_WPI_Class_lessVAR <- lapply(list.files(folder_PRED_WPI_class_lessVAR, pattern = ".tif$"), function(x){rast(file.path(folder_PRED_WPI_class_lessVAR, x))})

labs <- c("Not affected", "Slightly affected","Moderately affected","Strongly affected","Seriously affected")
png(file.path(output, paste0(name_RF, "_interpolate_RF_pred.png")), height = 800, width = 800)
par(mfrow=c(3,4))
for(y in 1:length(unique_years)){
    # max_abs <- max(abs(global(baseline_diff[[y]], "max", na.rm=TRUE)),
    #             abs(global(baseline_diff[[y]], "min", na.rm=TRUE)))
       levels(all_pred_WPI_Class_lessVAR[[y]]) <- data.frame(ID=1:5, label=labs)
  plot(all_pred_WPI_Class_lessVAR[[y]],
       main = paste("WPI", unique_years[y]), )
}
dev.off()

# only TSM und Chla Klassifikation and NORMIERT ------------------------
name_RF <- "WPI_class_NORM_TSMChla"
folder_PRED_WPI_class_TSMChla <- file.path(folder_PRED, name_RF)
dir.create(folder_PRED_WPI_class_TSMChla)

rf_model_WPI_RF_CLASS_onlyTSMChla$forest$independent.variable.names # variable names
names(raster_years)

all_pred_WPI_Class_TSMChla <- func_pred_RF(raster_years, rf_model_WPI_RF_CLASS_onlyTSMChla, sampleloc_extent4_land,
                     outpur_dir = folder_PRED_WPI_class_TSMChla)
# all_pred_WPI_Class_TSMChla <- lapply(list.files(folder_PRED_WPI_class_TSMChla, pattern = ".tif$"), 
          # function(x){rast(file.path(folder_PRED_WPI_class_TSMChla, x))})

labs <- c("Not affected", "Slightly affected","Moderately affected","Strongly affected","Seriously affected")
png(file.path(output, paste0(name_RF, "_interpolate_RF_pred.png")), height = 800, width = 800)
par(mfrow=c(3,4))
for(y in 1:length(unique_years)){
    # max_abs <- max(abs(global(baseline_diff[[y]], "max", na.rm=TRUE)),
    #             abs(global(baseline_diff[[y]], "min", na.rm=TRUE)))
       levels(all_pred_WPI_Class_TSMChla[[y]]) <- data.frame(ID=1:5, label=labs)
  plot(all_pred_WPI_Class_TSMChla[[y]],
       main = paste("WPI", unique_years[y]), )
}
dev.off()
# ---------------------------------------------------------
# make predictions with F model ------------------
# ---------------------------------------------------------
folder_PRED_F <- file.path(folder_PRED, "2F_NORM")
dir.create(folder_PRED_F)

all_pred_F <- func_pred_RF(raster_years, rf_model_F_RF, sampleloc_extent4_land,  outpur_dir = folder_PRED_F)
# all_pred_F <- lapply(list.files(folder_PRED_F), function(x){rast(file.path(folder_PRED_F, x))})

brks <- c(0,1,2,3,4,5,10,15,20)
png(file.path(output, "F_interpolate_RF_pred.png"), height = 800, width = 800)
par(mfrow=c(7,4))
for(y in 1:length(unique_years)){
    # max_abs <- max(abs(global(baseline_diff[[y]], "max", na.rm=TRUE)),
    #             abs(global(baseline_diff[[y]], "min", na.rm=TRUE)))

  plot(all_pred_F[[y]],
       main = paste("F", unique_years[y]),plg = list(title = "[mg/l]"),   breaks = brks)
}
dev.off()

# Klassifikation and NORMIERT ------------------------
name_RF <- "2F_class_NORM"
folder_PRED_WPI_class <- file.path(folder_PRED, name_RF)
dir.create(folder_PRED_F_class)

# rf_model_WPI_RF_CLASS$forest$independent.variable.names # variable names

all_pred_F_Class <- func_pred_RF(raster_years, rf_model_F_RF_CLASS, sampleloc_extent4_land, outpur_dir = folder_PRED_F_class)

labs <- c("F 0.0-0.3", "F 0.3-5", "F 5-10", "F 10-15", "F >= 15")
png(file.path(output, paste0(name_RF, "_interpolate_RF_pred.png")), height = 800, width = 800)
par(mfrow=c(6,4))
for(y in 1:length(unique_years)){
    # max_abs <- max(abs(global(baseline_diff[[y]], "max", na.rm=TRUE)),
    #             abs(global(baseline_diff[[y]], "min", na.rm=TRUE)))
       levels(all_pred_F_Class[[y]]) <- data.frame(ID=1:5, label=labs)
  plot(all_pred_F_Class[[y]],
       main = paste("WPI", unique_years[y]), )
}
dev.off()


# ---------------------------------------------------------
# make predictions with P model ------------------
# ---------------------------------------------------------
folder_PRED_P <- file.path(folder_PRED, "P3")
dir.create(folder_PRED_P)

all_pred_P <- func_pred_RF(raster_years, rf_model_P_RF, sampleloc_extent4_land, 
              outpur_dir = folder_PRED_P)
# all_pred_P <- lapply(list.files(folder_PRED_P), function(x){rast(file.path(folder_PRED_P, x))})

brks <- c( 0, 0.5,1,2,3,4, 5, 6) 
png(file.path(output, "P3_interpolate_RF_pred.png"), height = 800, width = 800)
par(mfrow=c(7,4))
for(y in 1:length(unique_years)){
    # max_abs <- max(abs(global(baseline_diff[[y]], "max", na.rm=TRUE)),
    #             abs(global(baseline_diff[[y]], "min", na.rm=TRUE)))

  plot(all_pred_P[[y]] / 1000,
       main = paste("P", unique_years[y]), plg = list(title = "[mg/l]"), breaks = brks)
}
dev.off()
