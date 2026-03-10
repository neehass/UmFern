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

output <- "./R-scripts/2013Valid-outputs-MODIS_TERRA"
dir.create(output)

data_elZrelli_shp <- "./elZrelli2018_shp"
data_elZrelli <- "./elZrelli2018"
data_elKateb_shap <- "./elKateb2018_shp"
# load interpolate data
data_dir_valid_masekd <- "./landsat2013-masked"

# Modis 
# data_modis <- "./ModisTerra-Chla/MODIST_L3m_CHL_NRT_2022.0-20260310_110439"
data_modis <- "./MODIS_SEPTEMBER_DATA"

# data_CHla <- "./NOAAMSL12_chla/chlora/monthly/WW00"
# data_nlw <- "./NOAAMSL12_chla/nlw/monthly/WW00"

# # NOAA Chlorophyl and # # nlw 
# files <- list.files(data_CHla)
# chla <- rast(file.path(data_CHla, files))
# plot(chla[[1]])
# files_nlw <- list.files(data_nlw)
# nlw <- rast(file.path(data_nlw, files_nlw))
# plot(nlw[[1]])

# MODIS TERRA 
files_modis <-  list.files(data_modis , pattern = "_a\\.4km\\.nc$")
modisTERRA_chla <- rast(file.path(data_modis, files_modis))
dates_modis <- as.Date(substr(basename(files_modis), 13, 20), format="%Y%m%d")
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

sampleLoc_points <- vect(file.path(data_elZrelli_shp, "elZrelli2018_sampleloc2.shp"))
plot(sampleLoc_points)

# ---------------------------------------
# ---------------------------------------------------------
# crop CHla data  Modis Terra  ------------------
# ---------------------------------------------------------
sampleloc_extent_pj <- project(sampleloc_extent, crs(modisTERRA_chla))
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
chla_2013_max10[chla_2013_max10 >= 10] <- NA

plot(chla_2013_max10, main = as.character(september2013))

buffer_points <- buffer(sampleloc_points_pj, width = 1000) # meter
chla_2013_sep_points <- mask(chla_2013, buffer_points)

plot(chla_2013_sep_points, main = as.character(time(chla_crop_sep)))

# ---------------------------------------------------------
# load elZrelli sep 2013 data --------
# ---------------------------------------------------------
heavymetal <-  read.csv(file.path(data_elZrelli, "elZrelli2018_heavymetal_SORTED.csv"), sep = ",", header = TRUE, row.names = 1)
waterPolIndex <-  read.csv(file.path(data_elZrelli, "elZrelli2018_waterPollutionIndex.csv"), sep = ";", header = TRUE)

WPI_interpolate <- rast(file.path(data_dir_valid_masekd, "WPi_interpolate.tif")) 
F_interpolate <- rast(file.path(data_dir_valid_masekd, "F_interpolate.tif")) 
P_interpolate <- rast(file.path(data_dir_valid_masekd, "P_interpolate.tif")) 
Zn_interpolate <- rast(file.path(data_dir_valid_masekd, "ZN_interpolate.tif")) 
Cu_interpolate <- rast(file.path(data_dir_valid_masekd, "Cu_interpolate.tif")) 

# CORRELATION Chla with wpi -------------------
WPI_interpolate_pj <- project(WPI_interpolate, crs(chla_2013))
plot(WPI_interpolate)
# chla_2013_masked <- mask(resample(chla_2013, WPI_interpolate_pj$var1.pred), WPI_interpolate_pj$var1.pred)
chla_2013_masked <- mask(resample(chla_2013_sep_points, WPI_interpolate_pj$var1.pred), WPI_interpolate_pj$var1.pred)
png(file.path(output, "MODIS_chla_2013SEP.png"), height = 1000, width = 1000)
plot(chla_2013_masked)
dev.off()

valuesChla_2013masekd <- values(chla_2013_masked)
valuesWPI <- values(WPI_interpolate_pj$var1.pred)

# correlation 
df_corr <- data.frame(
  Chla = valuesChla_2013masekd,
  WPI  = valuesWPI
)
names(df_corr) <- c("Chla", "WPI")

perCOR <- cor(valuesChla_2013masekd, valuesWPI, use = "complete.obs")
R2_WPI <- mean(perCOR^2)

p_cor <- ggplot(df_corr, aes(Chla, WPI))
p_cor <- p_cor +
  geom_bin2d() +
  scale_fill_viridis_c() + 
  geom_smooth(method = "lm", color = "red") +
  labs(
    x = "Chl-a [mg/m3]",
    y = "Water Pollution Index (WPI)",
    title = "Correlation between Chl-a and WPI"
  ) +
  annotate(
    "text", 
    x = max(df_corr$Chla, na.rm = TRUE)*0.7,   # Position x
    y = max(df_corr$WPI, na.rm = TRUE)*0.9,    # Position y
    label = paste0("median R² = ", round(median(perCOR^2), 3), 
                  "\nmax R² = ", round(max(perCOR^2), 3),
                  "\nmin R² = ", round(min(perCOR^2), 3)),     # Text
    color = "red",
    size = 5
  ) + theme_minimal()
p_cor
ggsave(file.path(output, "MODIScoor_WPi_Chl-a.png"), p_cor, height = 6, width = 6, scale = 1.2)

# ------------------------
# Flourine 
F_interpolate_pj <- project(F_interpolate, crs(chla_2013))
valuesF <- values(F_interpolate_pj$var1.pred)

# correlation 
df_corrF <- data.frame(
  Chla = valuesChla_2013masekd,
  F  = valuesF
)
names(df_corrF) <- c("Chla", "F")

perCOR_F <- cor(valuesChla_2013masekd, valuesF, use = "complete.obs")
R2_F<- mean(perCOR_F^2)

p_cor <- ggplot(df_corrF, aes(Chla, F))
p_cor <- p_cor +
  geom_bin2d() +
  scale_fill_viridis_c() + 
  geom_smooth(method = "lm", color = "red") +
  labs(
    x = "Chl-a [mg/m3]",
    y = "Flourine (F) [mg/l]",
    title = "Correlation between Chl-a and F"
  ) +
  annotate(
    "text", 
    x = max(df_corr$Chla, na.rm = TRUE)*0.7,   # Position x
    y = max(df_corr$WPI, na.rm = TRUE)*0.9,    # Position y
    label = paste0("median R² = ", round(median(perCOR_F^2), 3), 
                  "\nmax R² = ", round(max(perCOR_F^2), 3),
                  "\nmin R² = ", round(min(perCOR_F^2), 3)),     # Text
    color = "red",
    size = 5
  ) + theme_minimal()
p_cor
ggsave(file.path(output, "MODIScoor_F_Chl-a.png"), p_cor, height = 6, width = 6, scale = 1.2)

# ------------------------
# Phosphore 
P_interpolate_pj <- project(P_interpolate, crs(chla_2013))
# plot(P_interpolate_pj)
# chla_2013_masked <- mask(resample(chla_2013, P_interpolate_pj$var1.pred), P_interpolate_pj$var1.pred)
chla_2013_masked <- mask(resample(chla_2013_sep_points, P_interpolate_pj$var1.pred), P_interpolate_pj$var1.pred)
plot(chla_2013_masked)

valuesChla <- values(chla_2013_masked)
valuesP <- values(P_interpolate_pj$var1.pred)

df_corrP <- data.frame(
  Chla = valuesChla,
  P    = valuesP
)

names(df_corrP) <- c("Chla", "P")
head(df_corrP)
df_corrP <- na.omit(df_corrP)

perCOR_P <- cor(valuesChla, valuesP, use = "complete.obs")
R2_P <- perCOR_P^2

p_cor <- ggplot(df_corrP, aes(Chla, P))
p_cor <- p_cor +
  geom_bin2d() +
  scale_fill_viridis_c() + 
  geom_smooth(method = "lm", color = "red") +
  labs(
    x = "Chl-a [mg/m3]",
    y = "Phosphore (P) [µg/l]",
    title = "Correlation between Chl-a and P"
  ) +
  annotate(
    "text", 
    x = max(df_corrP$Chla, na.rm = TRUE)*0.7,   # Position x
    y = max(df_corrP$P, na.rm = TRUE),    # Position y
    label =  paste0("median R² = ", round(median(perCOR_P^2), 3), 
                  "\nmax R² = ", round(max(perCOR_P^2), 3),
                  "\nmin R² = ", round(min(perCOR_P^2), 3)),     # Text
    color = "red",
    size = 5
  ) + theme_minimal()
p_cor
ggsave(file.path(output, "MODIScoor_P_Chl-a.png"), p_cor, height = 6, width = 6, scale = 1.2)

# ------------------------
# Zinc 
Cu_interpolate_pj <- project(Cu_interpolate, crs(chla_2013))
# plot(Cu_interpolate_pj)
# chla_2013_masked <- mask(resample(chla_2013, Cu_interpolate_pj$var1.pred), Cu_interpolate_pj$var1.pred)
chla_2013_masked <- mask(resample(chla_2013_sep_points, Cu_interpolate_pj$var1.pred), Cu_interpolate_pj$var1.pred)
plot(chla_2013_masked)

valuesChla <- values(chla_2013_masked)
valuesCu <- values(Cu_interpolate_pj$var1.pred)

df_corrCu <- data.frame(
  Chla = valuesChla,
  Zn   = valuesCu
)

names(df_corrCu) <- c("Chla", "Cu")
# head(df_corrCu)
df_corrCu <- na.omit(df_corrCu)

perCOR_Cu <- cor(valuesChla, valuesCu, use = "complete.obs")
R2_Cu <- perCOR_Cu^2

p_cor <- ggplot(df_corrCu, aes(Chla, Cu))
p_cor <- p_cor +
  geom_bin2d() +
  scale_fill_viridis_c() + 
  geom_smooth(method = "lm", color = "red") +
  labs(
    x = "Chl-a [mg/m3]",
    y = "Cupper (Cu) [µg/l]",
    title = "Correlation between Chl-a and Zn"
  ) +
  annotate(
    "text", 
    x = max(df_corrCu$Chla, na.rm = TRUE)*0.7,   # Position x
    y = max(df_corrCu$Cu, na.rm = TRUE),    # Position y
    label = paste0("median R² = ", round(median(perCOR_Cu), 3), 
                  "\nmax R² = ", round(max(perCOR_Cu^2), 3),
                  "\nmin R² = ", round(min(perCOR_Cu^2), 3)),     # Text
    color = "red",
    size = 5
  ) + theme_minimal()
p_cor
ggsave(file.path(output, "MODIScoor_Cu_Chl-a.png"), p_cor, height = 6, width = 6, scale = 1.2)

# ------------------------
# Zinc 
Zn_interpolate_pj <- project(Zn_interpolate, crs(chla_2013))
plot(Zn_interpolate_pj)
# chla_2013_masked <- mask(resample(chla_2013, Zn_interpolate_pj$var1.pred), Zn_interpolate_pj$var1.pred)
chla_2013_masked <- mask(resample(chla_2013_sep_points, Zn_interpolate_pj$var1.pred), Zn_interpolate_pj$var1.pred)
plot(chla_2013_masked)

valuesChla <- values(chla_2013_masked)
valuesZn <- values(Zn_interpolate_pj$var1.pred)

df_corrZn <- data.frame(
  Chla = valuesChla,
  Zn   = valuesZn
)

names(df_corrZn) <- c("Chla", "Zn")
head(df_corrZn)
df_corrZn <- na.omit(df_corrZn)

perCOR_Zn <- cor(valuesChla, valuesZn, use = "complete.obs")
R2_Zn <- perCOR_Zn^2

p_cor <- ggplot(df_corrZn, aes(Chla, Zn))
p_cor <- p_cor +
  geom_bin2d() +
  scale_fill_viridis_c() + 
  geom_smooth(method = "lm", color = "red") +
  labs(
    x = "Chl-a [mg/m3]",
    y = "zinc (Zn) [µg/l]",
    title = "Correlation between Chl-a and Zn"
  ) +
  annotate(
    "text", 
    x = max(df_corrZn$Chla, na.rm = TRUE)*0.7,   # Position x
    y = max(df_corrZn$Zn, na.rm = TRUE),    # Position y
    label =  paste0("median R² = ", round(median(perCOR_Zn), 3), 
                  "\nmax R² = ", round(max(perCOR_Zn^2), 3),
                  "\nmin R² = ", round(min(perCOR_Zn^2), 3)),     # Text
    color = "red",
    size = 5
  ) + theme_minimal()
p_cor
ggsave(file.path(output, "MODIScoor_Zn_Chl-a.png"), p_cor, height = 6, width = 6, scale = 1.2)
# # ---------------------------------------------------------------------------------
# # ---------------------------------------------------------
# # RF trainig ------------------
# # ---------------------------------------------------------
output_RF <- "./RF_TERRA"
dir.create(output_RF)
output_MODLE <- file.path(output_RF, "./RF_models")
dir.create(output_MODLE)
raster_outRF_pred <- file.path(output_RF, "./RF_pred") 
dir.create(raster_outRF_pred)

set.seed(42)
# windows()
# -------------------------------------------
# Random Forest für Water Pollution Index ----
WPI_interpolate_pj <- project(WPI_interpolate, crs(chla_2013_sep_points))

chla_2013_masked <- mask(chla_2013_sep_points, resample(WPI_interpolate_pj$var1.pred, chla_2013_sep_points))
plot(chla_2013_masekd_stacked_RF)

WPI_interpolate_masked <- mask(WPI_interpolate_pj, terra::resample(chla_2013_sep_points[[1]], WPI_interpolate_pj))
plot(WPI_interpolate_masked)

WPI_interpolate_masked2 <- terra::resample(WPI_interpolate_pj, chla_2013[[1]])
plot(WPI_interpolate_masked2)
chla_2013stacked_RF <- mask(chla_2013, WPI_interpolate_masked2)
chla_2013stacked_RF <- func_median_mean(chla_2013stacked_RF)

brks <- c(0,1,2,3,4, 6, 10, 15, 30)#seq(0.3, 20,  by = 4)
RF_WPI <- func_RF_ranger_MODIS(dat_interpolate = WPI_interpolate_masked2, RSdata_valid = chla_2013stacked_RF, 
        model_name = "WPI_RF2", output_MODLE, output_RF, raster_outRF_pred,
                brks = brks, unite = "WPI")  


rcl <- matrix(c(-Inf, 1, 1,
  1, 2, 2,
  2, 3, 3,
  3, 5, 4,
  5, Inf, 5), ncol = 3, byrow = TRUE)


WPI_interpolate_masked2 <- terra::resample(WPI_interpolate_pj, chla_2013[[1]])
WPI_interpolate_masked2$WPI_class <- classify(WPI_interpolate_masked2$var1.pred, rcl)
WPI_interpolate_masked2$WPI_class <- as.factor(WPI_interpolate_masked2$WPI_class)
levels(WPI_interpolate_masked2$WPI_class) <- data.frame(ID=1:5, label=labs)
plot(WPI_interpolate_masked2$WPI_class)

chla_2013stacked_RF <- mask(chla_2013, WPI_interpolate_masked2)
chla_2013stacked_RF <- func_median_mean(chla_2013stacked_RF)
chla_2013stacked_RF <- chla_2013stacked_RF[[c("median", "sd")]]
plot(chla_2013stacked_RF)
ncol(chla_2013stacked_RF)

labs <- c("Not affected", "Slightly affected","Moderately affected","Strongly affected","Seriously affected")  
RF_WPI_KLASS <- func_RF_ranger_classMODIS(dat_interpolate = WPI_interpolate_masked2, 
                RSdata_valid = chla_2013stacked_RF, model_name = "WPI_R_KLASS", output_MODLE, output_RF, raster_outRF_pred,
                rcl = brks, unite = "WPI", labs)

chla_2013stacked_pred <- func_median_mean(chla_2013)
plot(chla_2013stacked_pred)

pred_ALL2013 <- predict(chla_2013stacked_pred, RF_WPI_KLASS$RF)
levels(pred_ALL2013) <- data.frame(ID=1:5, label=labs)
plot(pred_ALL2013, main = model_name)
plot(gulf_shp_pj, add = TRUE)
points(sampleloc_points_pj, col = "cyan")
points(indloc_pj, col = "red", pch = 15, cex = 1)

# other time steps test


pred_date <- september_dates_clean[format(september_dates_clean, "%y") == "25"]

data_pred <- chla_crop_sepCLEAN[[which(september_dates_clean %in% pred_date)]]

stacked_pred <- func_median_mean(data_pred)
plot(stacked_pred)

pred_ALL <- predict(stacked_pred, RF_WPI_KLASS$RF)
levels(pred_ALL) <- data.frame(ID=1:5, label=labs)
plot(pred_ALL, main = paste(model_name, "2000"))
plot(gulf_shp_pj, add = TRUE)
points(sampleloc_points_pj, col = "cyan")
points(indloc_pj, col = "red", pch = 15, cex = 1)


# # ---------------------------------------------------------------------------------
# # ---------------------------------------------------------
# # crop CHla data  NOAA ------------------
# # ---------------------------------------------------------
# sampleloc_extent_pj <- project(sampleloc_extent, crs(chla))
# sampleloc_points_pj <- project(sampleLoc_points, crs(chla))
# gulf_shp_pj <- project(gulf_shp, crs(chla))

# chla_crop <- crop(chla, sampleloc_extent_pj)
# plot(chla_crop[[1]])
# dates <- time(chla_crop)

# september_dates <- dates[format(dates, "%m") == "09"]

# chla_crop_sep <- chla_crop[[which(dates %in% september_dates)]]
# plot(chla_crop_sep, main = as.character(time(chla_crop_sep)))

# max_abs <- max(global(chla_crop_sep, "max", na.rm=TRUE))
# min_abs <- min(global(chla_crop_sep, "min", na.rm=TRUE))
# cols <- hcl.colors(50, "RdBu", rev = TRUE)
# breaks <- seq(min_abs, max_abs, length=20)
# par(mfrow=c(4,4))
# for(y in 1:length(september_dates)){
#   plot(chla_crop_sep[[y]], main = as.character(september_dates[y] ), 
#   zlim = c(min_abs, max_abs), col = cols)
# }
# dev.off()

# september2013 <- september_dates[format(september_dates, "%y") == "13"]

# chla_2013 <- chla_crop_sep[[which(september_dates %in% september2013)]]
# chla_2013_max10 <- chla_2013
# chla_2013_max10[chla_2013_max10 >= 10] <- NA

# buffer_points <- buffer(sampleloc_points_pj, width = 1000) # meter
# chla_2013_sep_points <- mask(chla_2013, buffer_points)

# plot(chla_2013_sep_points, main = as.character(time(chla_crop_sep)))
# plot(chla_crop_sep, main = as.character(time(chla_crop_sep)))
# plot(chla_2013_max10, main = as.character(september2013))
# plot(gulf_shp_pj, add = TRUE)

# # ---------------------------------------------------------
# # load elZrelli sep 2013 data --------
# # ---------------------------------------------------------
# heavymetal <-  read.csv(file.path(data_elZrelli, "elZrelli2018_heavymetal_SORTED.csv"), sep = ",", header = TRUE, row.names = 1)
# waterPolIndex <-  read.csv(file.path(data_elZrelli, "elZrelli2018_waterPollutionIndex.csv"), sep = ";", header = TRUE)

# WPI_interpolate <- rast(file.path(data_dir_valid_masekd, "WPi_interpolate.tif")) 
# F_interpolate <- rast(file.path(data_dir_valid_masekd, "F_interpolate.tif")) 
# P_interpolate <- rast(file.path(data_dir_valid_masekd, "P_interpolate.tif")) 
# Zn_interpolate <- rast(file.path(data_dir_valid_masekd, "ZN_interpolate.tif")) 
# Cu_interpolate <- rast(file.path(data_dir_valid_masekd, "Cu_interpolate.tif")) 

# # CORRELATION Chla with wpi -------------------
# WPI_interpolate_pj <- project(WPI_interpolate, crs(chla_2013))
# plot(WPI_interpolate)
# # chla_2013_masked <- mask(resample(chla_2013, WPI_interpolate_pj$var1.pred), WPI_interpolate_pj$var1.pred)
# chla_2013_masked <- mask(resample(chla_2013_sep_points, WPI_interpolate_pj$var1.pred), WPI_interpolate_pj$var1.pred)
# plot(chla_2013_masked)

# valuesChla_2013masekd <- values(chla_2013_masked)
# valuesWPI <- values(WPI_interpolate_pj$var1.pred)

# # correlation 
# df_corr <- data.frame(
#   Chla = valuesChla_2013masekd,
#   WPI  = valuesWPI
# )
# names(df_corr) <- c("Chla", "WPI")

# perCOR <- cor(valuesChla_2013masekd, valuesWPI, use = "complete.obs")
# R2_WPI <- perCOR^2

# p_cor <- ggplot(df_corr, aes(Chla, WPI))
# p_cor <- p_cor +
#   geom_bin2d() +
#   scale_fill_viridis_c() + 
#   geom_smooth(method = "lm", color = "red") +
#   labs(
#     x = "Chl-a [mg/m3]",
#     y = "Water Pollution Index (WPI)",
#     title = "Correlation between Chl-a and WPI"
#   ) +
#   annotate(
#     "text", 
#     x = max(df_corr$Chla, na.rm = TRUE)*0.7,   # Position x
#     y = max(df_corr$WPI, na.rm = TRUE)*0.9,    # Position y
#     label = paste0("R² = ", round(R2_WPI, 3)),     # Text
#     color = "red",
#     size = 5
#   ) + theme_minimal()
# p_cor
# ggsave(file.path(output, "NOAAcoor_WPi_Chl-a.png"), p_cor, height = 6, width = 6, scale = 1.2)

# # entlang longitude --------------------------------
# # calculate scaling factor ------------------
# scale_factor <- max(valuesChla_2013masekd, na.rm = TRUE) /
#                 max(valuesWPI, na.rm = TRUE)

# valuesWPI_scaled <- valuesWPI * scale_factor

# coords_chla <- xyFromCell(chla_2013_masked, 1:ncell(chla_2013_masked))
# values_long_INDEX <- data.frame(
#   Longitude = coords_chla[,1],
#   Latitude  = coords_chla[,2],
#   chla = valuesChla_2013masekd
# )
# names(values_long_INDEX) <- c("Longitude", "Latitude", "Chla")

# coords_WPI <- xyFromCell(WPI_interpolate_pj, 1:ncell(WPI_interpolate_pj))
# waterPolIndex <- data.frame(
#   Longitude = coords_WPI[,1],
#   Latitude = coords_WPI[,2],
#   WPi_scaled = valuesWPI * scale_factor
# )
# names(waterPolIndex)<- c("Longitude", "Latitude", "WPI")

# gabesLON <- ind_loc[[1]]
# p_WPI_INDEX <- ggplot(values_long_INDEX, aes(Longitude, Chla)) +
  
#   # Chl-a
#   geom_point(color = "coral") +

#   # WPI
#   geom_point(data = waterPolIndex, aes(x = Longitude, y = WPI),
#              color = "blue", size = 3) +

#   # Gabes vertical line
#   geom_vline(xintercept = gabesLON, colour = "black", linetype = "dashed") +

#   # Achsen
#   scale_y_continuous(
#     name = "Chl-a [mg/m3]",
#     sec.axis = sec_axis(~ . / scale_factor, name = "Water Pollution Index (WPi)"), 
#   ) +

#   # Themes
#   theme(
#     axis.title.y.right = element_text(color = "blue"),
#     axis.text.y.right  = element_text(color = "blue")
#   ) +

#   # Labels
#   labs(
#     x = "Longitude",
#     title = "Spectral profile across Gulf of Gabès - Water Pollution Index (WPi)"
#   ) + theme_minimal()

# p_WPI_INDEX

# label_text <- paste0(
#   "<span style='color:coral;'>Chl-a: R = ", round(perCOR,2),
#   " | R² = ", round(R2_WPI,2)
# )

# p_WPI_INDEX2 <- p_WPI_INDEX  +
#   annotate("richtext",     x = max(values_long_INDEX$Longitude, na.rm = TRUE)-0.1,
#          y = max(values_long_INDEX$Chla, na.rm = TRUE)/2,
#            label = label_text,
#            hjust = 0,
#            size = 3)
# p_WPI_INDEX2
# ggsave(file.path(output, "WPi_Chla.png"), p_WPI_INDEX2, height = 6, width = 6, scale = 1.2)

# # ------------------------
# # Flourine 
# F_interpolate_pj <- project(F_interpolate, crs(chla_2013))
# valuesF <- values(F_interpolate_pj$var1.pred)

# # correlation 
# df_corrF <- data.frame(
#   Chla = valuesChla_2013masekd,
#   F  = valuesF
# )
# names(df_corrF) <- c("Chla", "F")

# perCOR_F <- cor(valuesChla_2013masekd, valuesF, use = "complete.obs")
# R2_F<- perCOR_F^2

# p_cor <- ggplot(df_corrF, aes(Chla, F))
# p_cor <- p_cor +
#   geom_bin2d() +
#   scale_fill_viridis_c() + 
#   geom_smooth(method = "lm", color = "red") +
#   labs(
#     x = "Chl-a [mg/m3]",
#     y = "Flourine (F) [mg/l]",
#     title = "Correlation between Chl-a and WPI"
#   ) +
#   annotate(
#     "text", 
#     x = max(df_corr$Chla, na.rm = TRUE)*0.7,   # Position x
#     y = max(df_corr$WPI, na.rm = TRUE)*0.9,    # Position y
#     label = paste0("R² = ", round(R2_F, 3)),     # Text
#     color = "red",
#     size = 5
#   ) + theme_minimal()
# p_cor
# ggsave(file.path(output, "NOAAcoor_F_Chl-a.png"), p_cor, height = 6, width = 6, scale = 1.2)

# # ------------------------
# # Phosphore 
# P_interpolate_pj <- project(P_interpolate, crs(chla_2013))
# plot(P_interpolate_pj)
# # chla_2013_masked <- mask(resample(chla_2013, P_interpolate_pj$var1.pred), P_interpolate_pj$var1.pred)
# chla_2013_masked <- mask(resample(chla_2013_sep_points, P_interpolate_pj$var1.pred), P_interpolate_pj$var1.pred)
# plot(chla_2013_masked)

# valuesChla <- values(chla_2013_masked)
# valuesP <- values(P_interpolate_pj$var1.pred)

# df_corrP <- data.frame(
#   Chla = valuesChla,
#   P    = valuesP
# )

# names(df_corrP) <- c("Chla", "P")
# head(df_corrP)
# df_corrP <- na.omit(df_corrP)

# perCOR_P <- cor(valuesChla, valuesP, use = "complete.obs")
# R2_P <- perCOR_P^2

# p_cor <- ggplot(df_corrP, aes(Chla, P))
# p_cor <- p_cor +
#   geom_bin2d() +
#   scale_fill_viridis_c() + 
#   geom_smooth(method = "lm", color = "red") +
#   labs(
#     x = "Chl-a [mg/m3]",
#     y = "Phosphore (P) [µg/l]",
#     title = "Correlation between Chl-a and P"
#   ) +
#   annotate(
#     "text", 
#     x = max(df_corrP$Chla, na.rm = TRUE)*0.7,   # Position x
#     y = max(df_corrP$P, na.rm = TRUE),    # Position y
#     label = paste0("R² = ", round(R2_P, 3)),     # Text
#     color = "red",
#     size = 5
#   ) + theme_minimal()
# p_cor
# ggsave(file.path(output, "NOAAcoor_P_Chl-a.png"), p_cor, height = 6, width = 6, scale = 1.2)

# # ------------------------
# # Zinc 
# Cu_interpolate_pj <- project(Cu_interpolate, crs(chla_2013))
# plot(Cu_interpolate_pj)
# # chla_2013_masked <- mask(resample(chla_2013, Cu_interpolate_pj$var1.pred), Cu_interpolate_pj$var1.pred)
# chla_2013_masked <- mask(resample(chla_2013_sep_points, Cu_interpolate_pj$var1.pred), Cu_interpolate_pj$var1.pred)
# plot(chla_2013_masked)

# valuesChla <- values(chla_2013_masked)
# valuesCu <- values(Cu_interpolate_pj$var1.pred)

# df_corrCu <- data.frame(
#   Chla = valuesChla,
#   Zn   = valuesCu
# )

# names(df_corrCu) <- c("Chla", "Cu")
# head(df_corrCu)
# df_corrCu <- na.omit(df_corrCu)

# perCOR_Cu <- cor(valuesChla, valuesCu, use = "complete.obs")
# R2_Cu <- perCOR_Cu^2

# p_cor <- ggplot(df_corrCu, aes(Chla, Cu))
# p_cor <- p_cor +
#   geom_bin2d() +
#   scale_fill_viridis_c() + 
#   geom_smooth(method = "lm", color = "red") +
#   labs(
#     x = "Chl-a [mg/m3]",
#     y = "Cupper (Cu) [µg/l]",
#     title = "Correlation between Chl-a and Zn"
#   ) +
#   annotate(
#     "text", 
#     x = max(df_corrCu$Chla, na.rm = TRUE),   # Position x
#     y = max(df_corrCu$Cu, na.rm = TRUE),    # Position y
#     label = paste0("R² = ", round(R2_Cu, 3)),     # Text
#     color = "red",
#     size = 5
#   ) + theme_minimal()
# p_cor
# ggsave(file.path(output, "NOAAcoor_Cu_Chl-a.png"), p_cor, height = 6, width = 6, scale = 1.2)

# # ------------------------
# # Zinc 
# Zn_interpolate_pj <- project(Zn_interpolate, crs(chla_2013))
# plot(Zn_interpolate_pj)
# # chla_2013_masked <- mask(resample(chla_2013, Zn_interpolate_pj$var1.pred), Zn_interpolate_pj$var1.pred)
# chla_2013_masked <- mask(resample(chla_2013_sep_points, Zn_interpolate_pj$var1.pred), Zn_interpolate_pj$var1.pred)
# plot(chla_2013_masked)

# valuesChla <- values(chla_2013_masked)
# valuesZn <- values(Zn_interpolate_pj$var1.pred)

# df_corrZn <- data.frame(
#   Chla = valuesChla,
#   Zn   = valuesZn
# )

# names(df_corrZn) <- c("Chla", "Zn")
# head(df_corrZn)
# df_corrZn <- na.omit(df_corrZn)

# perCOR_Zn <- cor(valuesChla, valuesZn, use = "complete.obs")
# R2_Zn <- perCOR_Zn^2

# p_cor <- ggplot(df_corrZn, aes(Chla, Zn))
# p_cor <- p_cor +
#   geom_bin2d() +
#   scale_fill_viridis_c() + 
#   geom_smooth(method = "lm", color = "red") +
#   labs(
#     x = "Chl-a [mg/m3]",
#     y = "zinc (Zn) [µg/l]",
#     title = "Correlation between Chl-a and Zn"
#   ) +
#   annotate(
#     "text", 
#     x = max(df_corrZn$Chla, na.rm = TRUE),   # Position x
#     y = max(df_corrZn$Zn, na.rm = TRUE),    # Position y
#     label = paste0("R² = ", round(R2_Zn, 5)),     # Text
#     color = "red",
#     size = 5
#   ) + theme_minimal()
# p_cor
# ggsave(file.path(output, "NOAAcoor_Zn_Chl-a.png"), p_cor, height = 6, width = 6, scale = 1.2)

