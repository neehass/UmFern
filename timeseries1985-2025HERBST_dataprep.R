# september timeseries vorbereiten
# elZrelli sep 2013 is vaildationdata

# packages
library(terra)
library(ggplot2)
library(dplyr)
library(sf)
library(mapview)

getwd()
setwd("./Phosphate")
source('./R-scripts/hel-func.R', chdir = TRUE)

output <- "./R-scripts/timeserie1985-2025HERBST-outputs"
dir.create(output)

data_dir <- "./landsat-1995_2025"
data_crop_dir <- "./landsat-1995_2025_masked"
dir.create(data_crop_dir)
data_dir_cloudmask <- "./llandsat-1995_2025_cloudmask"
dir.create(data_dir_cloudmask)

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
# load samplelocation data ------------------
# ---------------------------------------------------------
list.files(data_elZrelli_shp)
# sampleloc_extent <- vect(file.path(data_elZrelli_shp, "elZrelli2018_sampleloc_EXTENT.shp"))
# sampleloc_extent2 <- vect(file.path(data_elZrelli_shp, "elZrelli2018_sampleloc_EXTENT2.shp"))
# sampleloc_extent2 <- erase(sampleloc_extent2, gulf_shp)
# sampleloc_extent3 <- vect(file.path(data_elZrelli_shp, "elZrelli2018_sampleloc_EXTENT-Skhira.shp"))
# sampleloc_extent3 <- erase(sampleloc_extent3, gulf_shp)
sampleloc_extent4_land <- vect(file.path(data_elZrelli_shp, "elZrelli2018_sampleloc_EXTENT4_land.shp"))
gabes_hafen <- vect(file.path(data_elZrelli_shp, "Gabes_hafen.shp"))

# plot industrie loc and sample locs
# windows()
plot(sampleloc_extent4_land)
# lines(gulf_shp, add = TRUE, col = "darkgrey", alpha = 0.2)
# points(ind_loc, col = "red", pch = 15)

# ---------------------------------------------------------
# load landsat data ------------------
# ---------------------------------------------------------
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

# ---------------------------------------------------------
# mean per year in September ------------------
# ---------------------------------------------------------
# 1986 should be skipped ! too many clouds 
years <- 1985:2026 
# windows()
y <- 1
year <- years[y] 
raster_name <- paste0("mean_raster_yearSEP_masked", "_", year, ".tif")
print(year)

# SEPTEMBER
raster_yearSEP <- func_timestep_sel(year, files_band, satNR_bands, unique_dates, data_dir, time = "SEPTEMBER")
dates <- datesyear[substr(datesyear, 5, 6) %in% c("09")] # lasse "10" weg wegen wolken maske 

# mask cloud and extent + mit Land für PIF
raster_yearSEP_masked <- func_mask_cloud_EXTENT(raster_yearSEP, data_dir, sampleloc_extent4_land, gabes_hafen)

# mean wenn mehere scenen verfügabr
if(length(raster_yearSEP_masked) > 1){
    print(paste(length(raster_yearSEP_masked), "scenes"))
    # mean_raster_yearSEP_masked <- func_mean_raster(raster_yearSEP_masked)

    idx0 <- sapply(raster_yearSEP_masked, function(x) {
    class(x) == "numeric" # numeric und damit 0, dh zu hohen Wolken Anteil in HAfen region 
    })

    validALL_masked_list_clean <- raster_yearSEP_masked[!idx0]
    dates <- dates[!idx0]

    # plot(validALL_masked_list_clean[[3]], main = dates[3])

    # make raster for all months (name months)
    validALL_masked_rast <- rast(validALL_masked_list_clean)
    names(validALL_masked_rast)

    # meadian über August, Sep, Okt
    all_bands <- names(validALL_masked_rast)
    bands <- unique(all_bands)

    median_list <- lapply(bands, function(b) {
    sel <- which(all_bands == b)
    app(validALL_masked_rast[[sel]], median, na.rm = TRUE)
    })

    median_raster <- rast(median_list)
    names(median_raster) <- bands

    plot(mean_raster_yearSEP_masked)
} else { 
    print("just one scene")
    mean_raster_yearSEP_masked <- raster_yearSEP_masked[[1]]
    plot(mean_raster_yearSEP_masked)
}

# save
writeRaster(mean_raster_yearSEP_masked, file.path(data_crop_dir, raster_name), overwrite=TRUE)


# save all wolken masken als tiff 
files_qa <- list.files(data_dir, recursive = TRUE, pattern="QA_PIXEL.TIF$")

for(y in 1985:2026){
    print(y)
    parts <- strsplit(files_qa, "_")
    dates <- sapply(parts, `[`, 4)

    year_file <- files_qa[grep(as.character(y), dates)]

    if(length(year_file) == 0 ){
        print("skip")
        next
    } else {
        print("maks exists")
        if(length(year_file) > 1){
              qa_list <- lapply(year_file, function(x){rast(file.path(data_dir, x))})
              qa_aligned <- lapply(qa_list, function(r) {resample(r, qa_list[[1]], method = "near")})
              qa <- rast(qa_aligned)
        } else {
              qa <- rast(file.path(data_dir, year_file))
        }

         # QA Band laden
      
        sampleloc_extent3_pj <- project(sampleloc_extent3, crs(qa))

        # plot(qa)

        cloud  <- get_bit(qa, 3)
        shadow <- get_bit(qa, 4)
        cloud_mask <- cloud | shadow 
        # plot(cloud_mask)

        # CLoud und Extent maskieren
        masked <- crop(cloud_mask, sampleloc_extent3_pj)
        masked <- mask(masked, sampleloc_extent3_pj)

        raster_name <- paste0("cloudExtent_mask_", as.character(y), ".tif")
        writeRaster(masked, file.path(data_dir_cloudmask, raster_name), overwrite=TRUE)

    }
}

# ---------------------------------------------------------
# end ------------------
# ---------------------------------------------------------
lapply(qa_list, ext)
lapply(qa_list, res)
lapply(qa_list, crs)
