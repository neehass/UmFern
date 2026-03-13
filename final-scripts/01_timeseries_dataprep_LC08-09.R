# ------------------------------------------------------------------------
# timeseries vorbereiten Landsatdaten LC08-09
# elZrelli is vaildation data, Messungen September 2013
# ------------------------------------------------------------------------

# packages
library(terra)
library(ggplot2)
library(dplyr)
library(sf)
library(mapview)

getwd()
setwd("./Phosphate")
source('./R-scripts/final-scripts/hel-func.R', chdir = TRUE)

output <- "./R-scripts/output/01_timeseries2013_2025_LC08"
dir.create(output)

data_dir <- "./landsat-SEPTEMBER"
data_crop_dir <- "./final_landsat-SEPTEMBER_land_LC08-09"
dir.create(data_crop_dir)
data_pif <- "./final_landsat-SEPTEMBER_PIF_LC08-09"
dir.create(data_pif)

# Messdaten SEp 2013
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
# plot(sampleloc_extent3)
# lines(gulf_shp, add = TRUE, col = "darkgrey", alpha = 0.2)
# points(ind_loc, col = "red", pch = 15)

# ---------------------------------------------------------
# load landsat data ------------------
# ---------------------------------------------------------
files_band <- list.files(data_dir, recursive = TRUE, pattern="T1_SR_B[0-9]{1,2}\\.TIF$")

# nur "LC08", "LC09
parts <- strsplit(files_band, "_")
idx <- which(sapply(parts, `[`, 1) %in% c("LC08", "LC09"))
parts <- parts[idx]
satNR <- unique(sapply(parts, `[`, 1))
dates <- sapply(parts, `[`, 4)
unique_dates <- unique(dates)
print(sort(unique_dates))
length(unique_dates)

# band overview:
satNR_bands <- list(
   "LT05" = list("B1" = "Blue", "B2" = "Green", "B3" = "Red", "B4" = "NIR", "B5" = "SWIR1", "B6" = "Thermal", "B7" = "SWIR2"),
   "LE07" =  list("B1" = "Blue", "B2" = "Green", "B3" = "Red", "B4" = "NIR", "B5" = "SWIR1", "B6" = "Thermal", "B7" = "SWIR2", "B8" = "Panchromatic"),
   "LC08" =  list("B1" = "CoastalAerosol", "B2" = "Blue", "B3" = "Green", "B4" = "Red", "B5" = "NIR", "B6" = "SWIR1", "B7" = "SWIR2", "B8" = "Panchromatic", "B9" = "TIRS1", "B10" = "TIRS2"),
   "LC09" =  list("B1" = "CoastalAerosol", "B2" = "Blue", "B3" = "Green", "B4" = "Red", "B5" = "NIR", "B6" = "SWIR1", "B7" = "SWIR2", "B8" = "Panchromatic", "B9" = "TIRS1", "B10" = "TIRS2")
)

# ---------------------------------------------------------
# medain per year in autum (Aug-Nov) ------------------
# --------------------------------------------------------- 
years <- 2013:2025 
# windows()
y <- 1
year <- years[y] 
print(year)
raster_year <- func_timestep_sel(year, files_band, satNR_bands, unique_dates, data_dir, time = "ALL") #, staNR_exclude = staNR_exclude)

years <- 2013
# Autum
start <- Sys.time()
for(y in 1:length(years)){
    year <- years[y] 
    # raster_name <- paste0("median_raster_year_land", "_", year, ".tif")
    raster_name <- paste0("mean_raster_year_land", "_", year, ".tif")
    print(year)

    raster_year <- func_timestep_sel(year, files_band, satNR_bands, unique_dates, data_dir, time = "ALL")

    if(class(raster_year[[1]]) == "SpatRaster") { # nur wenn szene existiert 
        # mask cloud and extent 
        # if(year == 2013){ # nur für 2013 ohne Hafen Region, sonst zu wenig pixel
        #    raster_year_masked <- func_mask_cloud_EXTENT(raster_year, data_dir, 
        #         sampleloc_extent4_land)
        # } else {
          raster_year_masked <- func_mask_cloud_EXTENT(raster_year, data_dir, 
                sampleloc_extent4_land, gabes_hafen)
        # }
         
        idx0 <- sapply(raster_year_masked, function(x) {
        class(x) == "numeric" # numeric und damit 0, dh zu hohen Wolken Anteil in HAfen region 
        })

        raster_year_masked_clean <- raster_year_masked[!idx0]

        # mean wenn mehere scenen verfügabr
        if(length(raster_year_masked_clean) != 0) {
            
            if(length(raster_year_masked) > 1){
                print(paste(length(raster_year_masked), "scenes"))
                median_raster_year_masked <- func_mean_raster(raster_year_masked_clean)
                # median_raster_year_masked <- func_median_raster(raster_year_masked_clean)
                plot(median_raster_year_masked)

                # save
                print(paste("save", year))
                writeRaster(median_raster_year_masked, file.path(data_crop_dir, raster_name), overwrite=TRUE)
            } else { 
                print("just one scene")
                median_raster_year_masked <- raster_year_masked[[1]]
                plot(median_raster_year_masked)

                # save
                print(paste("save", year))
                writeRaster(median_raster_year_masked, file.path(data_crop_dir, raster_name), overwrite=TRUE)
            }
        } else {
            print("Szene ist NUlL - Woleken Anteil zu hoch") 
            next
            }
    } else {
        print(paste("no szene", year))
        next
     }
}
end <- Sys.time()
diff_sec <- as.numeric(difftime(end, start, units = "secs"))
hours   <- floor(diff_sec / 3600)
minutes <- floor((diff_sec %% 3600) / 60)
seconds <- round(diff_sec %% 60, 2)
message(sprintf("Dauer: %02d:%02d:%05.2f (hh:mm:ss)", hours, minutes, seconds))

# ---------------------------------------------------------
# permanent stabile Pixel für Normierung filtern (pP)------------------
# ---------------------------------------------------------
# Referenz: 
# Schott, J.R., C. Salvaggio & W.J. Volchok (1988). 
# Radiometric scene normalization using pseudoinvariant features. Remote Sensing of Environment 26:1–16.
# Hessel et al. (2020)
# Relative Radiometric Normalization Using Automatically Extracted PIFs in Multi‑Sensor, Multi‑Temporal Series:
# → Anwendung über lange Zeitserien und PIF‑Auswahl automatisch anhand radiometrischer Konsistenz zwischen Bilder

# look at:  # https://search.r-project.org/CRAN/refmans/landsat/html/PIF.html

# -------------------------------------------
# # load alle raster  ----------------------------------
# -------------------------------------------
files <- list.files(data_crop_dir, pattern = "^mean")
parts <- strsplit(files, "_|\\.")
years <- sapply(parts, `[`, 5)
unique_years <- unique(years)
print(unique_years) # von 2013 2025

list_ALLyears <- lapply(files, function(x){rast(file.path(data_crop_dir, x))})
list_ALLyears <- lapply(list_ALLyears, function(x) {x[[ !names(x) %in% c("CoastalAerosol") ]]})

plot(list_ALLyears[[1]]$Blue,
       main = paste(unique_years[1]))

# png(file.path(output, "exampleTimeseries2013-2025.png"), height = 800, width = 800)
# par(mfrow=c(5,2))
# for(y in 1:length(unique_years)){
#   plot(list_ALLyears[[y]]$Green,
#        main = paste(unique_years[y], "- Green"))
# }
# dev.off()

# extent anpassen
extents <- lapply(list_ALLyears, ext)
# smallest_extent <- Reduce(intersect, extents)
areas <- sapply(extents, function(e) (xmax(e) - xmin(e)) * (ymax(e) - ymin(e)))
id_smallest <- which.min(areas)

# ref <- list_ALLyears[[15]] # 2013
ref_smallest <- list_ALLyears[[id_smallest]] 

list_ALLyears_aligned <- lapply(list_ALLyears, function(r) {
  r <- crop(r, ext(ref_smallest))  # auf Ref-Extent zuschneiden
  r <- mask(r, ref_smallest) # auf kleinsten Auschnitt zuschneiden 
  return(r)
})

# png(file.path(output, "MASKIERT_smallestTimeseries2013-2025.png"), height = 800, width = 800)
# par(mfrow=c(5,2))
# for(y in 1:length(unique_years)){
#   plot(list_ALLyears_aligned[[y]]$Green,
#        main = paste(unique_years[y], "- Green"))
# }
# dev.off()
# names(list_ALLyears_aligned[[1]])

# # alle zusammen
# rasterALL <- rast(list_ALLyears_aligned)
# bands <- unique(names(rasterALL))
# n_bands <- length(bands)
# n_years <- length(unique_years)

# if (n_bands * n_years != nlyr(rasterALL)) {
#   stop("Layeranzahl stimmt nicht mit Bändern * Jahren überein!")
# }

# gulf_shp_pj <- project(gulf_shp, crs(rasterALL))
# rasterALL_masked <- mask(rasterALL, gulf_shp_pj)
# names(rasterALL_masked)
# writeRaster(rasterALL_masked, file.path(data_pif, "MASKED_LAND_all_raster2013-2025_masked_smallest.tif"), overwrite=TRUE)

# layer_names <- unlist(lapply(unique_years, function(y) {
#   paste0(y, "_", bands)
# }))
# names(rasterALL) <- layer_names
# names(rasterALL) 
# writeRaster(rasterALL, file.path(data_pif, "all_raster2013-2025_masked_smallest.tif"), overwrite=TRUE)

# -------------------------------------------
# Permanente Pixel ----------------------------------
# -------------------------------------------
rasterALL <- rast(file.path(data_pif, "all_raster2013-2025_masked_smallest.tif"))
rasterALL_masked <- rast(file.path(data_pif, "MASKED_LAND_all_raster2013-2025_masked_smallest.tif"))

# Standardabweichung über alle bänder 
sd_raster <- app(rasterALL_masked, sd, na.rm = TRUE)
col_fun <- colorRampPalette(c("white", "red"))

png(file.path(output, "SD_smallestTimeseries2013-2025.png"), height = 800, width = 800)
plot(sd_raster, col = col_fun(100000), plg = list(title = "SD 2013-2025"))
dev.off()

cv <- app(rasterALL_masked, function(x){
  sd(x, na.rm=TRUE) / mean(x, na.rm=TRUE)
})

png(file.path(output, "CV_smallestTimeseries2013-2025.png"), height = 800, width = 800)
plot(cv, col = col_fun(10), plg = list(title = "CV 2013-2025"))
dev.off()

# 10% der stabilsten pixel 
vals_sd <- values(sd_raster)
threshold_sd <- quantile(vals_sd, 0.10, na.rm = TRUE)

vals_cv <- values(cv)
threshold_cv <- quantile(vals_cv, 0.10, na.rm = TRUE)

pP_mask_sd <- sd_raster <= threshold_sd
plot(pP_mask_sd)

pP_mask_cv <- cv <= threshold_cv
plot(pP_mask_cv)

pP_mask_sdCV <- pP_mask_sd & pP_mask_cv
png(file.path(output, "PIF_CV_SD_Timeseries2013-2025.png"), height = 800, width = 800)
plot(pP_mask_sdCV, main = "Permanent Pixel Mask (10% quantile CV & SD)")
dev.off()

writeRaster(pP_mask_sd, file.path(data_pif, "permanentPixel_sd.tif"), overwrite=TRUE)
writeRaster(pP_mask_cv, file.path(data_pif, "permanentPixel_CV.tif"), overwrite=TRUE)
writeRaster(pP_mask_sdCV, file.path(data_pif, "permanentPixel_sd_CV.tif"), overwrite=TRUE)
pP_mask_sdCV <- rast(file.path(data_pif, "permanentPixel_sd_CV.tif"))

# sd und cv pro band ------------------------------------------
bands_unique <- unique(names(rasterALL_masked))
sd_per_band <- list()
cv_per_band <- list()
pP_mask_sdCV_perband <- list()

for(b in bands_unique){
  print(paste(b, "in", c(bands_unique)))
  # alle Bänder mit diesem Namen auswählen
  sel <- grep(paste0("^", b, "$"), names(rasterALL_masked))
  band_stack <- rasterALL_masked[[sel]]
  
  # SD pro Pixel über die Jahre
  sd_per_band[[b]] <- app(band_stack, fun = sd, na.rm = TRUE)
 cv_per_band[[b]] <- app(band_stack, fun = function(x){
    sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)})
}

for(b in bands_unique){
  print(paste(b, "in", c(bands_unique)))
  # SD und CV für dieses Band
  sd_r <- sd_per_band[[b]]
  cv_r <- cv_per_band[[b]]
  
  # Schwellenwerte: 10% der stabilsten Pixel
  vals_sd <- values(sd_r)
  threshold_sd <- quantile(vals_sd, 0.10, na.rm = TRUE)
  
  vals_cv <- values(cv_r)
  threshold_cv <- quantile(vals_cv, 0.10, na.rm = TRUE)
  
  # Masken: SD und CV kleiner als Schwellenwert
  mask_sd <- sd_r <= threshold_sd
  mask_cv <- cv_r <= threshold_cv
  
  # kombinierte PIF-Maske
  pP_mask_sdCV_perband[[b]] <- mask_sd & mask_cv
  
  # Optional: plotten
  png(file.path(output, paste0(b, "_permanentPixel_sd_CV.png")), height = 800, width = 800)
  plot(pP_mask_sdCV_perband[[b]], main = paste("10% stabilste Pixel:", b))
  dev.off()
  
  writeRaster(pP_mask_sdCV_perband[[b]], file.path(data_pif, paste0(b, "_permanentPixel_sd_CV.tif")), overwrite=TRUE)
}

# -------------------------------------------
# normiere Alle szenen anhand 2013 -------------------------
# -------------------------------------------
ref_scene <- list_ALLyears_aligned[[which(unique_years == "2013")]] # 2013
# windows()
plot(ref_scene)
gulf_shp_pj <- project(gulf_shp, crs(ref_scene))

# normieren mit library(lmodel2) ---------
norm_ALLyears_list_band <- vector("list", length(list_ALLyears_aligned))
bandnames <- names(list_ALLyears_aligned[[1]])

for(i in seq_along(list_ALLyears_aligned)) {
  print(paste(i, "von", length(list_ALLyears_aligned)))
  scene <- list_ALLyears_aligned[[i]]

 # Maskierete Szene mit PIF
 scene_norm_bands <- list()
  for(b in seq_along(bandnames)) {

    # Einzelnes Band extrahieren
    ref_mask <- mask(ref_scene[[b]], pP_mask_sdCV_perband[[b]], maskvalues = 0)
    ref_mask <- mask(ref_mask, gulf_shp_pj)

    tar_band <- scene[[b]]

    # Gemeinsame PIF-Pixel für Regression
    common_mask <- !is.na(values(ref_mask)) & !is.na(values(tar_band))
    ref_vals <- values(ref_mask)[common_mask]
    tar_vals <- values(tar_band)[common_mask]

      #   # lineare Regression (nur PIF-Pixel)
    #   model <- lm(ref_vals ~ tar_vals)
    #   a <- coef(model)[1]
    #   b <- coef(model)[2]

    # Major-Axis Regression wie bei https://search.r-project.org/CRAN/refmans/landsat/html/PIF.html
    model <- lmodel2(ref_vals ~ tar_vals) # zeigt bloß hinweise
    a <- model$regression.results[2, "Intercept"]
    b_coef <- model$regression.results[2, "Slope"]

    # Normalisierung des Bands
    band_norm <- tar_band * b_coef + a

    # speichern
    scene_norm_bands[[b]] <- band_norm
  }

  # Normalisierte Bänder wieder zu einem Rasterstack zusammenfügen
  scene_norm <- rast(scene_norm_bands)
  names(scene_norm) <- bandnames

  # speichern
  writeRaster(scene_norm, file.path(data_pif, paste0("NORMIERT2013_proBand_",unique_years[i],"_permanentPixel_sd_CV.tif")), overwrite=TRUE)
  norm_ALLyears_list_band[[i]] <- scene_norm

  # Speicher freigeben
  rm(scene_masked, scene_norm_bands)
  gc()
}

# norm_files <- list.files(data_pif, pattern = "^NORMIERT2013_proBand")
norm_files <- list.files(data_pif, pattern = "^NORMIERT2013_proBand")
norm_ALLyears_list_band <- lapply(norm_files, function(x){rast(file.path(data_pif, x))})

# windows()
# plot(norm_ALLyears_list_band[[1]])
norm_dir <- file.path(output, "normiert_proBAND_sd_CV")
dir.create(norm_dir)
for(b in seq_along(bandnames)) {
    png(file.path(norm_dir, 
      paste0(bandnames[b],"_NORMIERT2013_proBand_pP_smallestTimeseries2013-2025_2.png")), 
      height = 800, width = 800)
    par(mfrow=c(3,4))
    for(y in 1:length(unique_years)){
    plot(norm_ALLyears_list_band[[y]][[b]],
        main = paste(unique_years[y] , "-", bandnames[b]))
    }
    dev.off()
}

# maskiere water --------------
gulf_shp_pj <- project(gulf_shp, crs(norm_ALLyears_list_band[[1]]))
norm_ALLyears_list_band_water <- lapply(norm_ALLyears_list_band, function(x){
    masked <- mask(x, gulf_shp_pj, inverse = TRUE)
    names(masked) <- names(x)
    return(masked)})
plot(norm_ALLyears_list_band_water[[1]])

# save 
for(i in seq_along(norm_ALLyears_list_band_water)){
  x <- norm_ALLyears_list_band_water[[i]]
  writeRaster(
    x,
    file.path(data_pif, paste0("WATER_NORMIERT2013_proBand_", unique_years[i], "_permanentPixel_sd_CV.tif")),
    overwrite = TRUE
  )
}

### END #####################################################

# # -------------------------------------------
# # normiere Alle szenen anhand 2013 -------------------------
# # -------------------------------------------
# ref_scene <- list_ALLyears_aligned[[which(unique_years == "2013")]] # 2013

# # permanent pixel maskieren 
# ref_mask <- mask(ref_scene, pP_mask_sdCV, maskvalues = 0)
# ref_mask <- mask(ref_mask, gulf_shp_pj)
# plot(ref_mask)

# tar_mask <- vector("list", length(list_ALLyears_aligned))

# for(i in seq_along(list_ALLyears_aligned)){

#   x <- list_ALLyears_aligned[[i]]

#   x_masked <- mask(x, pP_mask_sdCV, maskvalues = 0)

#   x_masked <- mask(x_masked, gulf_shp_pj)


#   tar_mask[[i]] <- x_masked

#   # Speicher freigeben
#   rm(mask_proj, gulf_pj, x_masked)
#   gc()
# }

# plot(tar_mask[[1]])

# # normieren ---------
# library(lmodel2)
# norm_ALLyears_list <- vector("list", length(list_ALLyears_aligned))
# bandnames <- names(list_ALLyears_aligned[[1]])

# for(i in seq_along(list_ALLyears_aligned)) {

#   scene <- list_ALLyears_aligned[[i]]

#  # Maskierete Szene mit PIF
#   scene_masked <- tar_mask[[i]]
 
#  scene_norm_bands <- list()
#   for(b in seq_along(bandnames)) {

#     # Einzelnes Band extrahieren
#     ref_band <- ref_mask[[b]]
#     tar_band <- scene[[b]]

#     # Gemeinsame PIF-Pixel für Regression
#     common_mask <- !is.na(values(ref_band)) & !is.na(values(tar_band))
#     ref_vals <- values(ref_band)[common_mask]
#     tar_vals <- values(tar_band)[common_mask]

#       #   # lineare Regression (nur PIF-Pixel)
#     #   model <- lm(ref_vals ~ tar_vals)
#     #   a <- coef(model)[1]
#     #   b <- coef(model)[2]

#     # Major-Axis Regression wie bei https://search.r-project.org/CRAN/refmans/landsat/html/PIF.html
#     model <- lmodel2(ref_vals ~ tar_vals) # zeigt bloß hinweise
#     a <- model$regression.results[2, "Intercept"]
#     b_coef <- model$regression.results[2, "Slope"]

#     # Normalisierung des Bands
#     band_norm <- tar_band * b_coef + a

#     # speichern
#     scene_norm_bands[[b]] <- band_norm
#   }

#   # Normalisierte Bänder wieder zu einem Rasterstack zusammenfügen
#   scene_norm <- rast(scene_norm_bands)
#   names(scene_norm) <- bandnames

#   # speichern
#   writeRaster(scene_norm, file.path(data_pif, paste0("NORMIERT2013_",unique_years[i],"_permanentPixel_sd_CV.tif")), overwrite=TRUE)
#   norm_ALLyears_list[[i]] <- scene_norm

#   # Speicher freigeben
#   rm(scene_masked, scene_norm_bands)
#   gc()
# }

# plot(norm_ALLyears_list[[1]])
# norm_dir <- file.path(output, "normiert_sd_CV")
# dir.create(norm_dir)
# for(b in seq_along(bandnames)) {
#     png(file.path(norm_dir, paste0(bandnames[b],"_NORMIERT2013_pP_smallestTimeseries1985-2025.png")), height = 800, width = 800)
#     par(mfrow=c(6,4))
#     for(y in 1:length(unique_years)){
#     plot(norm_ALLyears_list[[y]][[b]],
#         main = paste(unique_years[y] , "-", bandnames[b]))
#     }
#     dev.off()
# }

# # PIF Pseudo-Invariant Features # https://search.r-project.org/CRAN/refmans/landsat/html/PIF.html
# library(landsat)
# pif_list <- lapply(list_ALLyears_aligned, function(x){PIF(x[["Red"]], x[["NIR"]], x[["SWIR2"]], level = 0.99)}) # (band3 = Rot, band4 = NIR, band7 = SWIR2
# pif_list[[1]]
# # use major axis regression: error in both x and y
# nov.correction <- lmodel2:::lmodel2(july3@data[july.pif@data[,1] == 1, 1] ~ 
# nov3@data[july.pif@data[,1] == 1, 1])$regression.results[2, 2:3]
# nov3.corrected <- nov3
# # nov3.corrected@data[,1] <- nov3@data[,1] * nov.correction[2] + nov.correction[1]
# ref_scene <- list_ALLyears_aligned[[which(unique_years == "2013")]] # 2013

# pif_ref <- PIF(ref_scene[["Red"]], ref_scene[["NIR"]], ref_scene[["SWIR2"]], level = 0.1)
# sum(pif_ref == 1, na.rm = TRUE)

# tar <- list_ALLyears_aligned[[1]] 

# ref_vals <- values(ref_scene)   # gibt Vektor aller Pixel
# tar_vals <- values(tar)

# # nur die Pixel, die PIF sind
# idx <- pif_ref == 1

# var(values(ref_scene)[pif_ref == 1], na.rm = TRUE)
# var(values(tar)[pif_ref == 1], na.rm = TRUE)

# ref_vals <- ref_vals[idx]
# tar_vals <- tar_vals[idx]

# # NA entfernen
# valid_idx <- !is.na(ref_vals) & !is.na(tar_vals)  & is.finite(ref_vals) & is.finite(tar_vals)
# ref_vals <- ref_vals[valid_idx]
# tar_vals <- tar_vals[valid_idx]
# if(length(ref_vals) < 10 || var(ref_vals) == 0 || var(tar_vals) == 0){
#   stop("Nicht genügend valide PIF-Pixel mit Varianz > 0")
# }

# # Regression
# model <- lmodel2(ref_vals ~ tar_vals)
# tar.correction <- model$regression.results[2, 2:3] 

# ---------------------------------------------------------
# # save all wolken masken als tiff ------------------
# ---------------------------------------------------------

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

