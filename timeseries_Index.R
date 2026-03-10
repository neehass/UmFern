# Veränderungen der Indices 
# packages
library(terra)
library(ggplot2)
library(dplyr)
library(sf)
library(mapview)
library(ranger)
library(purrr)

getwd()
setwd("./Phosphate")
source('./R-scripts/hel-func.R', chdir = TRUE)

output <- "./R-scripts/timeseries-outputs_NORM_RF"
dir.create(output)

data_crop_dir <- "./landsat-SEPTEMBER_PIF"
data_dir_cloudmask <- "./landsat-SEPTEMBER_cloudmask"
folder_out <- "timeseries_Index_NORM"
dir.create(folder_out)

data_elZrelli_shp <- "./elZrelli2018_shp"
data_elZrelli <- "./elZrelli2018"

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
# load landsat data ------------------
# ---------------------------------------------------------
files_raster <- list.files(data_crop_dir, recursive = TRUE, pattern="^WATER_NORMIERT")

parts <- strsplit(files_raster, "_|\\.")
years <- sapply(parts, `[`, 4)
unique_years <- unique(years)
print(unique_years) # von 1985 bis 2025 >> 40Jahre
length(unique_years) # 24 Scenen 

raster_years <- lapply(files_raster, function(x){rast(file.path(data_crop_dir, x))})
names(raster_years) <- paste0("year_",unique_years)
raster_years <- lapply(raster_years, function(x) {x[[ !names(x) %in% c("CoastalAerosol", "SWIR2") ]]})

png(file.path(output, paste0("RAW", ".png")), height = 800, width = 800)
par(mfrow=c(6,4))
for(y in 1:length(unique_years)){
    # max_abs <- max(abs(global(baseline_diff[[y]], "max", na.rm=TRUE)),
    #             abs(global(baseline_diff[[y]], "min", na.rm=TRUE)))

  plot(all_pred_WPI[[y]],
       main = paste("WPI", unique_years[y]),  breaks = brks)
}
dev.off()

# Indices ----------------------------------
# TSM_laili, logChla, NDWI, TI, sediment, NDSSI, SWIR1/NIR
raster_years_index <- lapply(raster_years, function(x){func_index(x)})
names(raster_years_index) <- paste0("year_",unique_years)

idx <- c("TSM_laili", "logChla", "NDWI", "TI", "sediment", "NDSSI", "SWIR1_NIR")

for(i in 1:lenght(idx)){
    print(idx[i])
    all_max_abs <- Inf
    for(y in 1:length(unique_years)){
    r <- raster_years_index[[y]][[idx[i]]]
    max_abs <- max(abs(global(r, "max", na.rm=TRUE)),
                    abs(global(r, "min", na.rm=TRUE)))
        all_max_abs <- max(max_abs)
    
    }

    png(file.path(folder_out, paste0(idx[i], "NORM2013.png")), height = 1000, width = 1000)
    par(mfrow=c(6,4))
    for(y in 1:length(unique_years)){
    plot(raster_years_index[[y]][[idx[i]]],
         zlim = c(-all_max_abs, all_max_abs),
        main = paste(idx[i], "\n", unique_years[y], "- 1985 (baseline)"))
        #plot(gabes_shp, add = TRUE, border="black", lwd=1)
    }
    dev.off()
}

# baseline 1985
baseline1985 <- raster_years_index[[1]]
namesTRUE <- lapply(raster_years_index, function(x){names(x) != names(baseline1985)})
any(unlist(namesTRUE))

basline_diff <- func_baseline_diff(raster_years_index[2:length(raster_years_index)], baseline1985)

idx <- c("TSM_laili", "logChla", "NDWI", "TI", "sediment", "NDSSI", "SWIR1_NIR")

for(i in 1:lenght(idx)){
    print(idx[i])
    all_max_abs <- Inf
    for(y in 1:length(unique_years)){
    r <- basline_diff[[y]][[idx[i]]]
    max_abs <- max(abs(global(r, "max", na.rm=TRUE)),
                    abs(global(r, "min", na.rm=TRUE)))
        all_max_abs <- max(all_max_abs, max_abs)
    
    }
    col_fun_diff <- colorRampPalette(c("red","white", "blue")) 

    par(mfrow=c(7,5))
    for(y in 1:length(unique_years)){
    plot(basline_diff[[y]][[idx[i]]], col = col_fun_diff(100), zlim = c(-all_max_abs, all_max_abs),
        main = paste(idx[i], "\n", unique_years[y], "- 1985 (baseline)"))
        #plot(gabes_shp, add = TRUE, border="black", lwd=1)
    }
    dev.off()
}
