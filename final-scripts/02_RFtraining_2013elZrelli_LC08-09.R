# In September 2013, 
# surface seawater samples were collected from 16 stations located along the central coastal area of GG

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

output <- "./R-scripts/02_2013NORM_Valid_elZrelli-outputs_LC08-09"
dir.create(output)

data_elZrelli_shp <- "./elZrelli2018_shp"
data_elZrelli <- "./elZrelli2018"

data_pif <- "./landsat-SEPTEMBER_PIF_LC08-09"

data_dir_valid_masekd <- "./landsat-interpolation"
dir.create(data_dir_valid_masekd)

# VAlidation Chla Reana
data_CHla <- "./NOAAMSL12_chla/chlora/monthly/WW00"

# files <- list.files(data_CHla)
# REANA_chla <- rast(file.path(data_CHla, files))

# Validation nwl Reana
data_nlw <- "./NOAAMSL12_chla/nlw/monthly/WW00"

# files_nlw <- list.files(data_nlw)
# REANA_nlw <- rast(file.path(data_nlw, files_nlw))

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

sampleloc_extent2 <- vect(file.path(data_elZrelli_shp, "elZrelli2018_sampleloc_EXTENT2.shp"))
sampleloc_extent2 <- erase(sampleloc_extent2, gulf_shp)
sampleLoc_points <- vect(file.path(data_elZrelli_shp, "elZrelli2018_sampleloc2.shp"))

sampleloc_extent4_land <- vect(file.path(data_elZrelli_shp, "elZrelli2018_sampleloc_EXTENT4_land.shp"))
sampleloc_extent4_land <- erase(sampleloc_extent4_land, gulf_shp)

# plot industrie loc and sample locs
# plot(sampleloc_extent_interpol)
# lines(gulf_shp, add = TRUE, col = "darkgrey", alpha = 0.2)
# points(sampleLoc_points, col = "cyan")
# points(ind_loc, col = "red", pch = 15)

# ---------------------------------------------------------
# load elZrelli 2018, sep 2013 ------------------
# ---------------------------------------------------------
# load el Kateb 2018 data >> vorbereitet in 00_elZrelli_dataprep_interpolate.R
# list.files(data_elZrelli)
heavymetal <-  read.csv(file.path(data_elZrelli, "elZrelli2018_heavymetal_POINTS.csv"), sep = ",", header = TRUE, row.names = 1)
waterPolIndex <-  read.csv(file.path(data_elZrelli, "elZrelli2018_waterPollutionIndex_POINTS.csv"), sep = ";", header = TRUE)

# ---------------------------------------------------------
# load landsat data, sep 2013 ------------------
# ---------------------------------------------------------
# crop sample locations ------------------
# # load normierte scene 
# norm_files <- list.files(data_pif, pattern = "^WATER_NORMIERT2013_proBand_")
norm_files <- list.files(data_pif, pattern = "^2_WATER_NORMIERT2013_proBand_")
parts <- strsplit(norm_files, "_")
years <- sapply(parts, `[`, 5) #4 
idxyear2013 <- which(years == "2013")

valid_SEP_Extent <- rast(file.path(data_pif, norm_files[idxyear2013]))
plot(valid_SEP_Extent)

# # crop sample loc extent -----
# sampleloc_points_pj <- project(sampleLoc_points, crs(valid_SEP_Extent))

buffer <- buffer(sampleloc_points_pj, width = 200) # meter
# valid_SEP_points <- mask(valid_SEP_Extent, buffer)
# valid_SEP_points <- crop(valid_SEP_points, valid_SEP_Extent)
# plot(valid_SEP_points$Blue)

# # writeRaster(valid_SEP_points, file.path(data_pif,"NORM13_perBand_POINTS_valid_sampleLoc_SEP2013.tif") , overwrite = TRUE)
# writeRaster(valid_SEP_points, file.path(data_pif,"2_NORM13_perBand_POINTS_valid_sampleLoc_SEP2013.tif") , overwrite = TRUE)

# load 
# valid_SEP_points <- rast(file.path(data_pif,"NORM13_perBand_POINTS_valid_sampleLoc_SEP2013.tif"))
valid_SEP_points <- rast(file.path(data_pif,"2_NORM13_perBand_POINTS_valid_sampleLoc_SEP2013.tif"))

sampleloc_points_pj <- project(sampleLoc_points, crs(valid_SEP_Extent))
indloc_pj <- project(vect(matrix(ind_loc, ncol=2), crs="EPSG:4326"), crs(valid_SEP_Extent))
gulf_shp_pj <- project(gulf_shp, crs(valid_SEP_Extent))

# plot 
# png(file.path(output, "Extent.png"), height = 1000, width = 1000)
# plot(valid_SEP_Extent$SWIR1)
# points(indloc_pj, col = "red")
# points(sampleloc_points_pj, col = "cyan")
# dev.off()

# png(file.path(output, "points.png"), height = 1000, width = 1000)
# plot(valid_SEP_points$SWIR1)
# points(indloc_pj, col = "red")
# lines(gulf_shp_pj, col = "black", alpha = 0.2)
# dev.off()


# ---------------------------------------------------------
# ANALYSE alle variables  --------------------------------------------
# ---------------------------------------------------------
sample_pints_values <- values(sampleloc_points_pj)
sample_pints_values <- sample_pints_values[,names(sample_pints_values) %in% c("Sector", "Station", "Latitude", "Longitude") ]

# landsat values
valid_SEP_points_VALUES <- terra::extract(valid_SEP_Extent, buffer, fun = mean)
valid_SEP_points_VALUES$Station <- buffer$Station
valid_SEP_points_VALUES <- na.omit(valid_SEP_points_VALUES)
valid_SEP_points_VALUES <- merge(valid_SEP_points_VALUES, sample_pints_values, by = "Station")
colnames(valid_SEP_points_VALUES)
# chl-a TSM
valid_SEP_points_VALUES$TMS_tassan <- func_TSM_tassan(valid_SEP_points_VALUES$Blue, valid_SEP_points_VALUES$Green, valid_SEP_points_VALUES$Red)
valid_SEP_points_VALUES$TMS_laili <- func_TMS_Laili(valid_SEP_points_VALUES$Blue, valid_SEP_points_VALUES$Red)
valid_SEP_points_VALUES$log_Chl  <- func_log_Chl(valid_SEP_points_VALUES$Red, valid_SEP_points_VALUES$Green)
# water indices
valid_SEP_points_VALUES$NDWI  <- func_NDWI(valid_SEP_points_VALUES$Green, valid_SEP_points_VALUES$NIR)
valid_SEP_points_VALUES$TI  <- func_TI(valid_SEP_points_VALUES$Red, valid_SEP_points_VALUES$Green)
valid_SEP_points_VALUES$sediment  <- func_sedi(valid_SEP_points_VALUES$Red, valid_SEP_points_VALUES$Blue)
valid_SEP_points_VALUES$NDSSI  <- func_NDSSI(valid_SEP_points_VALUES$Red, valid_SEP_points_VALUES$NIR)

# plot
gabesLON <- ind_loc[[1]]

#####################################################
# watrer index  Red , Green
colnames(valid_SEP_points_VALUES)
values_long <- valid_SEP_points_VALUES %>%
    #   select(Longitude, Green, Red, SWIR1) %>%
  select(Longitude, Blue, Green, Red, SWIR1, SWIR2, NIR) %>%
  pivot_longer(# cols = c(Green, Red, SWIR1),
    cols = c( Blue, Green, Red, SWIR1, SWIR2, NIR),
               names_to = "Band",
               values_to = "Value")

# calculate scaling factor
scale_factor <- max(values_long$Value, na.rm = TRUE) /
                max(waterPolIndex$WPi, na.rm = TRUE)

waterPolIndex <- waterPolIndex %>%
  mutate(WPi_scaled = WPi * scale_factor)
# windows()
p_WPI_GR <- ggplot(values_long, aes(Longitude, Value, colour = Band)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = gabesLON,
             colour = "black",
             linetype = "dashed") +
      scale_colour_manual(values = c(Green = "darkgreen",
                                     Red = "red", 
                                     SWIR1 = "darkorange", 
                                     SWIR2 = "orange", NIR = "magenta",
                                     Blue = "blue")) + 
    # Water Pollution Index
  geom_point(data = waterPolIndex,
             aes(x = Longitude, y = WPi_scaled),
             colour = "blue",
             size = 3) +

  geom_line(data = waterPolIndex,
            aes(x = Longitude, y = WPi_scaled),
            colour = "blue",
            linewidth = 1) +

  geom_vline(xintercept = gabesLON,
             colour = "black",
             linetype = "dashed") +

  scale_y_continuous(
    name = "Reflectance",
    sec.axis = sec_axis(~ . / scale_factor,
                        name = "Water Pollution Index (WPi)")) +
theme(
  axis.title.y.right = element_text(color = "blue"),
  axis.text.y.right  = element_text(color = "blue")
) +
  labs(x = "Longitude",
       y = "Reflectance",
       title = "Spectral profile across Gulf of Gabès - Water Pollution Index (WPi)")
p_WPI_GR

df_corr <- values_long %>%
  pivot_wider(names_from = Band, values_from = Value) %>%
  left_join(waterPolIndex[, c("Longitude", "WPi")],
            by = "Longitude")
# WPI vs TSM
cor_Green  <- cor(df_corr$WPi, df_corr$Green, use = "complete.obs")
R2_Green   <- cor_Green^2
cor_Red  <- cor(df_corr$WPi, df_corr$Red, use = "complete.obs")
R2_Red   <- cor_Red^2
cor_nir  <- cor(df_corr$WPi, df_corr$NIR, use = "complete.obs")
R2_nir  <- cor_nir^2
cor_Sw1  <- cor(df_corr$WPi, df_corr$SWIR1, use = "complete.obs")
R2_SW1  <- cor_Sw1^2
cor_Sw2  <- cor(df_corr$WPi, df_corr$SWIR2, use = "complete.obs")
R2_SW2  <- cor_Sw2^2

label_text <- paste0(
  "<span style='color:darkgreen;'>Green: R = ", round(cor_Green,2),
  " | R² = ", round(R2_Green,2), "</span><br>",
  "<span style='color:red;'>Red: R = ", round(cor_Red,2),
  " | R² = ", round(R2_Red,2),"</span><br>",
     "<span style='color:magenta;'>NIR: R = ", round(cor_nir,2),
  " | R² = ", round(R2_nir,2),
  "</span><br>",
   "<span style='color:darkorange;'>SWIR1: R = ", round(cor_Sw1,2),
  " | R² = ", round(R2_SW1,2),
  "</span><br>",
   "<span style='color:orange;'>SWIR2: R = ", round(cor_Sw2,2),
  " | R² = ", round(R2_SW2,2)
)

p_WPI_GR2 <- p_WPI_GR  +
  annotate("richtext",     x = max(df_corr$Longitude, na.rm = TRUE)-0.1,
         y = max(values_long$Value, na.rm = TRUE)/2,
           label = label_text,
           hjust = 0,
           size = 3)

p_WPI_GR2
ggsave(file.path(output, "WPi_Green_Red.png"), p_WPI_GR2, height = 6, width = 6, scale = 1.2)

#####################################################
# watrer index  TSM, chl-a
values_long_INDEX <- valid_SEP_points_VALUES %>%
  select(Longitude, TMS_tassan, TMS_laili, log_Chl, NDWI, TI, sediment, NDSSI) %>%
  pivot_longer(cols = c(TMS_tassan, TMS_laili, log_Chl,  NDWI, TI, sediment, NDSSI),
               names_to = "Index",
               values_to = "Value")

# calculate scaling factor
scale_factor <- max(values_long_INDEX$Value, na.rm = TRUE) /
                max(waterPolIndex$WPi, na.rm = TRUE)

waterPolIndex <- waterPolIndex %>%
  mutate(WPi_scaled = WPi * scale_factor)

p_WPI_INDEX <- ggplot(values_long_INDEX, aes(Longitude, Value, colour = Index)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = gabesLON,
             colour = "black",
             linetype = "dashed") +
   scale_colour_manual(values = c(TMS_tassan = "darkgreen",
                                    TMS_laili = "darkolivegreen",
                                 log_Chl = "coral",
                                 NDWI = "magenta",
                                 TI = "red",
                                 sediemnt = "pink",
                                 NDSSI = "purple")) + 
 # Water Pollution Index
  geom_point(data = waterPolIndex,
             aes(x = Longitude, y = WPi_scaled),
             colour = "blue",
             size = 3) +

  geom_line(data = waterPolIndex,
            aes(x = Longitude, y = WPi_scaled),
            colour = "blue",
            linewidth = 1) +

  geom_vline(xintercept = gabesLON,
             colour = "black",
             linetype = "dashed") +

  scale_y_continuous(
    name = "Reflectance",
    sec.axis = sec_axis(~ . / scale_factor,
                        name = "Water Pollution Index (WPi)")) +
theme(
  axis.title.y.right = element_text(color = "blue"),
  axis.text.y.right  = element_text(color = "blue")
) +
  labs(x = "Longitude",
       y = "Index TSM [mg/l], Chl-a [mg/m3]",
       title = "Spectral profile across Gulf of Gabès - Water Pollution Index (WPi)")


df_corr <- values_long_INDEX %>%
  pivot_wider(names_from = Index, values_from = Value) %>%
  left_join(waterPolIndex[, c("Longitude", "WPi")],
            by = "Longitude")
# WPI vs TSM
cor_TSM_tassan  <- cor(df_corr$WPi, df_corr$TMS_tassan, use = "complete.obs")
R2_TSM_tassan   <- cor_TSM_tassan^2

cor_TSM_laili  <- cor(df_corr$WPi, df_corr$TMS_laili, use = "complete.obs")
R2_TSM_laili   <- cor_TSM_laili^2

# WPI vs Chl-a
cor_CHLA <- cor(df_corr$WPi, df_corr$log_Chl, use = "complete.obs")
R2_CHLA  <- cor_CHLA^2

# WPI vs CNDWI
cor_NDWI <- cor(df_corr$WPi, df_corr$NDWI, use = "complete.obs")
R2_NDWI <- cor_NDWI^2

# WPI vs TI
cor_TI <- cor(df_corr$WPi, df_corr$TI, use = "complete.obs")
R2_TI <- cor_TI^2
# WPI vs CNDWI
cor_sedi <- cor(df_corr$WPi, df_corr$sediment, use = "complete.obs")
R2_sedi <- cor_sedi^2
# WPI vs CNDWI
cor_NDSSI <- cor(df_corr$WPi, df_corr$NDSSI, use = "complete.obs")
R2_NDSSI <- cor_NDSSI^2

label_text <- paste0(
  "<span style='color:darkgreen;'>TSM (tassan): R = ", round(cor_TSM_tassan,2),
  " | R² = ", round(R2_TSM_tassan,2), "</span><br>",
  "<span style='color:darkolivegreen;'>TSM (laili): R = ", round(cor_TSM_laili,2),
  " | R² = ", round(R2_TSM_laili,2), "</span><br>",
  "<span style='color:coral;'>Chl-a: R = ", round(cor_CHLA,2),
  " | R² = ", round(R2_CHLA,2), "</span><br><br>",

 "<span style='color:magenta;'>NDWI: R = ", round(cor_NDWI,2),
  " | R² = ", round(cor_NDWI,2), "</span><br>",
   "<span style='color:red;'>TI: R = ", round(cor_TI,2),
  " | R² = ", round(cor_TI,2), "</span><br>",
  "<span style='color:pink;'>sediment: R = ", round(cor_sedi,2),
  " | R² = ", round(R2_sedi,2), "</span><br>",
   "<span style='color:purple;'>NDSSI: R = ", round(cor_NDSSI,2),
  " | R² = ", round(R2_NDSSI,2), "</span>"
)

p_WPI_INDEX2 <- p_WPI_INDEX  +
  annotate("richtext",     x = max(df_corr$Longitude, na.rm = TRUE)-0.15,
         y = max(values_long_INDEX$Value, na.rm = TRUE)/2,
           label = label_text,
           hjust = 0,
           size = 3)
p_WPI_INDEX2
ggsave(file.path(output, "WPi_INDEX.png"), p_WPI_INDEX2, height = 6, width = 6, scale = 1.2)

#####################################################
# Florine F und phosphorus P Red , Green
#####################################################
values_long <- valid_SEP_points_VALUES %>%
  select(Longitude, Green, Red, SWIR1, NIR) %>%
  pivot_longer(cols = c(Green, Red, SWIR1, NIR),
               names_to = "Band",
               values_to = "Value")

# calculate scaling factor
scale_factorF <- max(values_long$Value, na.rm = TRUE) /
                max(heavymetal$F_mean , na.rm = TRUE) # da in mg/l die ander in ym/l
scale_factorP <- max(values_long$Value, na.rm = TRUE) /
                max(heavymetal$P_mean, na.rm = TRUE)
scale_factorCU <- max(values_long$Value, na.rm = TRUE) /
        max(heavymetal$Cu_mean, na.rm = TRUE)           

scale_factorZN <- max(values_long$Value, na.rm = TRUE) /
        max(heavymetal$Zn_mean, na.rm = TRUE)  
heavymetal <- heavymetal %>%
  mutate(F = F_mean * scale_factorF,
        P = P_mean * scale_factorP,
        Cu = Cu_mean*scale_factorCU ,
         Zn = Zn_mean*scale_factorZN)

heavy_long <- heavymetal %>%
  select(Longitude, F, P, Cu, Zn) %>%
  pivot_longer(cols = c(F, P, Cu, Zn),
               names_to = "Band",
               values_to = "Value")

library(ggnewscale)
p_FP_GR <- ggplot(values_long,
                  aes(Longitude, Value, colour = Band)) +

  # ---------- BANDS ----------
  geom_point() +
  geom_line() +

  scale_colour_manual(
    name = "Spectral bands",
    values = c(
      Green = "darkgreen",
      Red = "red",
      SWIR1 = "darkorange",
      NIR = "magenta"
    )
  ) +

  geom_vline(xintercept = gabesLON,
             colour = "black",
             linetype = "dashed") +

  # 👉 NEUE FARBSKALA STARTEN
  ggnewscale::new_scale_colour() +

  # ---------- HEAVY METALS ----------
  geom_point(data = heavy_long,
             aes(Longitude, Value, colour = Band),
             size = 3) +

  geom_line(data = heavy_long,
            aes(Longitude, Value, colour = Band),
            linewidth = 1) +

  scale_colour_manual(
    name = "Heavy metals",
    values = c(
      F  = "#08519C",
      P  = "#3182BD",
      Cu = "#6BAED6",
      Zn = "#BDD7E7"
    )
  ) +

     scale_y_continuous(
    name = "Reflectance",
    sec.axis = sec_axis(~ . / scale_factorP,
                        name = "scaled Heavy metal values ")) + # [μg/l] F in [mg/l]
theme(
  axis.title.y.right = element_text(color = "blue"),
  axis.text.y.right  = element_text(color = "blue")
) +
  labs(
    x = "Longitude",
    y = "Reflectance",
    title = "Spectral profile across Gulf of Gabès"
  ) 

p_FP_GR

df_corr <- values_long %>%
  pivot_wider(names_from = Band, values_from = Value) %>%
  left_join(heavymetal[, c("Longitude", "F_mean", "P_mean", "Zn_mean", "Cu_mean")],
            by = "Longitude")
# WPI vs TSM
cor_Green  <- cor(df_corr$F_mean, df_corr$Green, use = "complete.obs")
R2_Green   <- cor_Green^2
cor_Red  <- cor(df_corr$F_mean, df_corr$Red, use = "complete.obs")
R2_Red   <- cor_Red^2
cor_Sw  <- cor(df_corr$F_mean, df_corr$SWIR1, use = "complete.obs")
R2_SW  <- cor_Sw^2
cor_Nir <- cor(df_corr$F_mean, df_corr$NIR, use = "complete.obs")
R2_nir  <- cor_Nir^2
cor_Green2  <- cor(df_corr$P_mean, df_corr$Green, use = "complete.obs")
R2_Green2   <- cor_Green2^2
cor_Red2  <- cor(df_corr$P_mean, df_corr$Red, use = "complete.obs")
R2_Red2  <- cor_Red2^2
cor_Sw2  <- cor(df_corr$P_mean, df_corr$SWIR1, use = "complete.obs")
R2_SW2  <- cor_Sw2^2
cor_Nir2 <- cor(df_corr$P_mean, df_corr$NIR, use = "complete.obs")
R2_nir2  <- cor_Nir2^2

cor_Green3  <- cor(df_corr$Zn_mean, df_corr$Green, use = "complete.obs")
R2_Green3  <- cor_Green3^2
cor_Red3  <- cor(df_corr$Zn_mean, df_corr$Red, use = "complete.obs")
R2_Red3  <- cor_Red3^2
cor_Sw3  <- cor(df_corr$Zn_mean, df_corr$SWIR1, use = "complete.obs")
R2_SW3  <- cor_Sw3^2
cor_Nir3 <- cor(df_corr$Zn_mean, df_corr$NIR, use = "complete.obs")
R2_nir3  <- cor_Nir3^2

cor_Green4  <- cor(df_corr$Cu_mean, df_corr$Green, use = "complete.obs")
R2_Green4  <- cor_Green4^2
cor_Red4  <- cor(df_corr$Cu_mean, df_corr$Red, use = "complete.obs")
R2_Red4  <- cor_Red4^2
cor_Sw4  <- cor(df_corr$Cu_mean, df_corr$SWIR1, use = "complete.obs")
R2_SW4  <- cor_Sw4^2
cor_Nir4 <- cor(df_corr$Cu_mean, df_corr$NIR, use = "complete.obs")
R2_nir4  <- cor_Nir4^2

label_text <- paste0(
     "<span style='color:deepskyblue3;'><b>Copper (Cu):</b></span><br>",
  "<span style='color:darkgreen;'>Cu Green: R = ", round(cor_Green4,2),
  " | R² = ", round(R2_Green4,2), "</span><br>",
  "<span style='color:red;'>Cu Red: R = ", round(cor_Red4,2),
  " | R² = ", round(R2_Red4,2),"</span><br>",
   "<span style='color:darkorange;'>Cu SWIR1: R = ", round(cor_Sw4,2),
  " | R² = ", round(R2_SW4,2), "</span><br>",   
   "<span style='color:magenta;'>Cu NIR: R = ", round(cor_Nir4,2),
  " | R² = ", round(R2_nir4,2), "</span><br><br>",

  "<span style='color:navy;'><b>Fluorine (F):</b></span><br>",
  "<span style='color:darkgreen;'>F Green: R = ", round(cor_Green,2),
  " | R² = ", round(R2_Green,2), "</span><br>",
  "<span style='color:red;'>F Red: R = ", round(cor_Red,2),
  " | R² = ", round(R2_Red,2),"</span><br>",
   "<span style='color:darkorange;'>F SWIR1: R = ", round(cor_Sw,2),
  " | R² = ", round(R2_SW,2), "</span><br>",
     "<span style='color:magenta;'>F NIR: R = ", round(cor_Nir,2),
  " | R² = ", round(R2_nir,2), "</span><br><br>",

  "<span style='color:royalblue;'><b>Phosphorus (P):</b></span><br>",
  "<span style='color:darkgreen;'>P Green: R = ", round(cor_Green2,2),
  " | R² = ", round(R2_Green2,2), "</span><br>",
  "<span style='color:red;'>P Red: R = ", round(cor_Red2,2),
  " | R² = ", round(R2_Red2,2),"</span><br>",
   "<span style='color:darkorange;'>P SWIR1: R = ", round(cor_Sw2,2),
  " | R² = ", round(R2_SW2,2),"</span><br>",
   "<span style='color:magenta;'>P NIR: R = ", round(cor_Nir2,2),
  " | R² = ", round(R2_nir2,2), "</span><br><br>",

 "<span style='color:lightblue;'><b>Zinc (Zn):</b></span><br>",
  "<span style='color:darkgreen;'>Zn Green: R = ", round(cor_Green3,2),
  " | R² = ", round(R2_Green3,2), "</span><br>",
  "<span style='color:red;'>Zn Red: R = ", round(cor_Red3,2),
  " | R² = ", round(R2_Red3,2),"</span><br>",
   "<span style='color:darkorange;'>Zn SWIR1: R = ", round(cor_Sw3,2),
  " | R² = ", round(R2_SW3,2),"</span><br>",
   "<span style='color:magenta;'>Zn NIR: R = ", round(cor_Nir3,2),
  " | R² = ", round(R2_nir3,2)

)

p_FP_GR2 <- p_FP_GR  +
  annotate("richtext",     x = max(df_corr$Longitude, na.rm = TRUE)*0.985,
         y = (max(values_long$Value, na.rm = TRUE)/2) - 1500,
           label = label_text,
           hjust = 0,
           size = 3)

p_FP_GR2
ggsave(file.path(output, "heavymetal_Green_Red.png"), p_FP_GR2, height = 6, width = 6, scale = 1.2)

#####################################################
# Florine F und phosphorus P INDEX
#####################################################
values_long_INDEX <- valid_SEP_points_VALUES %>%
  select(Longitude, TMS_tassan, TMS_laili, log_Chl, NDWI, TI, sediment, NDSSI) %>%
  pivot_longer(cols = c(TMS_tassan, TMS_laili, log_Chl, NDWI, TI, sediment, NDSSI),
               names_to = "Band",
               values_to = "Value")

# calculate scaling factor
scale_factorF <- max(values_long_INDEX$Value, na.rm = TRUE) /
                max(heavymetal$F_mean , na.rm = TRUE) # da in mg/l die ander in ym/l
scale_factorP <- max(values_long_INDEX$Value, na.rm = TRUE) /
                max(heavymetal$P_mean, na.rm = TRUE)
scale_factorCU <- max(values_long_INDEX$Value, na.rm = TRUE) /
        max(heavymetal$Cu_mean, na.rm = TRUE)           

scale_factorZN <- max(values_long_INDEX$Value, na.rm = TRUE) /
        max(heavymetal$Zn_mean, na.rm = TRUE)  
heavymetal_IDX <- heavymetal %>%
  mutate(F = F_mean * scale_factorF,
        P = P_mean * scale_factorP,
        Cu = Cu_mean*scale_factorCU ,
         Zn = Zn_mean*scale_factorZN)

heavy_long <- heavymetal_IDX %>%
  select(Longitude, F, P, Cu, Zn) %>%
  pivot_longer(cols = c(F, P, Cu, Zn),
               names_to = "Band",
               values_to = "Value")

library(ggnewscale)
p_FP_IDX <- ggplot(values_long_INDEX,
                  aes(Longitude, Value, colour = Band)) +

  # ---------- BANDS ----------
  geom_point() +
  geom_line() +

  scale_colour_manual(
    name = "Index",
    values = c(TMS_tassan = "darkgreen",
                                    TMS_laili = "darkolivegreen",
                                 log_Chl = "coral",
                                 NDWI = "magenta",
                                 TI = "red",
                                 sediemnt = "pink",
                                 NDSSI = "purple")) + 
  geom_vline(xintercept = gabesLON,
             colour = "black",
             linetype = "dashed") +

  # 👉 NEUE FARBSKALA STARTEN
  ggnewscale::new_scale_colour() +

  # ---------- HEAVY METALS ----------
  geom_point(data = heavy_long,
             aes(Longitude, Value, colour = Band),
             size = 3) +

  geom_line(data = heavy_long,
            aes(Longitude, Value, colour = Band),
            linewidth = 1) +

  scale_colour_manual(
    name = "Pollutants",
    values = c(
      F  = "#08519C",
      P  = "#3182BD",
      Cu = "#6BAED6",
      Zn = "#BDD7E7"
    )
  ) +

     scale_y_continuous(
    name = "Reflectance",
    sec.axis = sec_axis(~ . / scale_factorP,
                        name = "scaled Heavy metals values ")) + # [μg/l] F in [mg/l]
theme(
  axis.title.y.right = element_text(color = "blue"),
  axis.text.y.right  = element_text(color = "blue")
) +
  labs(
    x = "Longitude",
    y = "Reflectance",
    title = "Spectral profile across Gulf of Gabès"
  ) 

p_FP_IDX

df_corr <- values_long_INDEX %>%
  pivot_wider(names_from = Band, values_from = Value) %>%
  left_join(heavymetal[, c("Longitude", "F_mean", "P_mean", "Zn_mean", "Cu_mean")],
            by = "Longitude")
# WPI vs TSM
# WPI vs TSM
element_cols <- c("F_mean", "P_mean", "Zn_mean", "Cu_mean")

cor_results_tassan <- sapply(element_cols, function(col) {
  r <- cor(df_corr[[col]], df_corr$TMS_tassan, use = "complete.obs")
  c(R = r, R2 = r^2)
})

t(cor_results_tassan)

cor_results_laili <- sapply(element_cols, function(col) {
  r <- cor(df_corr[[col]], df_corr$TMS_laili, use = "complete.obs")
  c(R = r, R2 = r^2)
})

t(cor_results_laili)

cor_results_chla <- sapply(element_cols, function(col) {
  r <- cor(df_corr[[col]], df_corr$log_Chl, use = "complete.obs")
  c(R = r, R2 = r^2)
})

t(cor_results_chla)

cor_results_NDWI <- sapply(element_cols, function(col) {
  r <- cor(df_corr[[col]], df_corr$NDWI, use = "complete.obs")
  c(R = r, R2 = r^2)
})

t(cor_results_NDWI)

cor_results_TI <- sapply(element_cols, function(col) {
  r <- cor(df_corr[[col]], df_corr$TI, use = "complete.obs")
  c(R = r, R2 = r^2)
})

t(cor_results_TI)

cor_results_sedi <- sapply(element_cols, function(col) {
  r <- cor(df_corr[[col]], df_corr$sediment, use = "complete.obs")
  c(R = r, R2 = r^2)
})

t(cor_results_sedi)

cor_results_NDSSI <- sapply(element_cols, function(col) {
  r <- cor(df_corr[[col]], df_corr$NDSSI, use = "complete.obs")
  c(R = r, R2 = r^2)
})

t(cor_results_NDSSI)

ggsave(file.path(output, "heavymetal_INDEX.png"), p_FP_IDX, height = 6, width = 6, scale = 1.2)

# --------------------------------------------------------
# RandomForest Modelle mit Interpolated WPI, F und P 
# --------------------------------------------------------
# RF F
library(caret)
library(randomForest)
library(ranger)

output_RF <- file.path(output, "2_RF_NORM_interpolate_LC08-09")
dir.create(output_RF)
output_MODEL <- file.path(output_RF, "RF_NORM_models")
dir.create(output_MODEL)
raster_outRF_pred <- "./2_RF_NORM_pred_LC08-09"
dir.create(raster_outRF_pred)

# load interpolate data << 00_2013elZrelli_dataprep_interpolate.R 
WPI_interpolate <- rast(file.path(data_dir_valid_masekd, "WPi_interpolate.tif")) 
F_interpolate <- rast(file.path(data_dir_valid_masekd, "F_interpolate.tif")) 
P_interpolate <- rast(file.path(data_dir_valid_masekd, "P_interpolate.tif")) 
Cu_interpolate <- rast(file.path(data_dir_valid_masekd, "Cu_interpolate.tif")) 

# prep landsat remot sensing data
sampleloc_extent2_pj <- project(sampleloc_extent2, crs(valid_SEP_Extent))

set.seed(42)
# windows()
# -------------------------------------------
# Random Forest für Water Pollution Index ----
WPI_interpolate_pj <- project(WPI_interpolate,  crs(valid_SEP_Extent))
WPI_interpolate_pj_res <- resample(WPI_interpolate_pj$var1.pred, valid_SEP_Extent)

RSdata_valid <- mask(valid_SEP_Extent, WPI_interpolate_pj_res)
names(RSdata_valid) <- names(valid_SEP_Extent) 
plot(RSdata_valid)

brks <- c( 0,1,2,3,4, 6, 10, 15, 30)#seq(0.3, 20,  by = 4)
RF_WPI <- func_RF_ranger(WPI_interpolate, RSdata_valid, model_name = "WPI_RF", output_MODLE, output_RF, raster_outRF_pred,
                brks = brks, unite = "WPI")       

# -------------------------------------------
# RandomForest für Florine F ----
F_interpolate_pj <- project(F_interpolate,  crs(valid_SEP_Extent))
F_interpolate_pj_res <- resample(F_interpolate_pj$var1.pred, valid_SEP_Extent)

RSdata_valid <- mask(valid_SEP_Extent, F_interpolate_pj_res)
names(RSdata_valid) <- names(valid_SEP_Extent) 
plot(RSdata_valid)

unite <- "[mg/l]"
brks <- c( 0,1,2,3,4,5,10,15,20)
RF_F <- func_RF_ranger(F_interpolate, RSdata_valid, model_name = "F_RF", 
                output_MODLE, output_RF, raster_outRF_pred,
                brks = brks, unite = unite)

# -------------------------------------------
# RandomForest für Phsophor P ----
P_interpolate_pj <- project(P_interpolate,  crs(valid_SEP_Extent))
P_interpolate_pj_res <- resample(P_interpolate_pj$var1.pred, valid_SEP_Extent)

RSdata_valid <- mask(valid_SEP_Extent, P_interpolate_pj_res)
names(RSdata_valid) <- names(valid_SEP_Extent) 
plot(RSdata_valid)

unite <- "[mg/l]"
brks <- c( 0, 0.5,1,2,3,4, 5, 6)#
RF_P <- func_RF_ranger(P_interpolate, RSdata_valid, model_name = "P_RF", 
                output_MODLE, output_RF, raster_outRF_pred,
                brks = brks, unite = unite)

 
#---------------------------------------------------------------
#---------------------------------------------------------------
#--------------------------------------------------------------- 
# KLASSIFIKATION ------------------
# RandomForest Klassifikation Model mit Interpolated F und P 
set.seed(42)
# windows()
# -------------------------------------------
# Water Pollution Index ----
# Random Forest KLASSIFIKATION für Water Pollution Index ----
WPI_interpolate_pj <- project(WPI_interpolate,  crs(valid_SEP_Extent))
WPI_interpolate_pj_res <- resample(WPI_interpolate_pj$var1.pred, valid_SEP_Extent)

RSdata_valid <- mask(valid_SEP_Extent, WPI_interpolate_pj_res)
names(RSdata_valid) <- names(valid_SEP_Extent) 
plot(RSdata_valid)

rcl <- matrix(c(-Inf, 1, 1,
  1, 2, 2,
  2, 3, 3,
  3, 5, 4,
  5, Inf, 5), ncol = 3, byrow = TRUE)

labs <- c("Not affected", "Slightly affected","Moderately affected","Strongly affected","Seriously affected")
RF_class_WPI <- func_RF_ranger_class(WPI_interpolate, RSdata_valid, "WP_RF_class", output_MODEL,
                                  output_RF, raster_outRF_pred,
                                  rcl, unite = "WPI", labs)   



# -------------------------------------------
# Random Forest KLASSIFIKATION für Flourine ----

RF_class_F <- func_RF_ranger_class(F_interpolate, RSdata_valid_norm, "WP_RF_class_NORM", output_MODEL,
                                  output_RF, raster_outRF_pred,
                                  rcl, unite = "WPI")   


# -------------------------------------------
# Random Forest KLASSIFIKATION für Water Pollution Index ----

RF_class_WPI <- func_RF_ranger_class(P_interpolate, RSdata_valid_norm, "WP_RF_class_NORM",
                                  sampleloc_Rvalid_pj, output_MODEL,
                                  output_RF, raster_outRF_pred,
                                  rcl, unite = "WPI")   


# #---------------------------------------------------------------
# #---------------------------------------------------------------
# #--------------------------------------------------------------- 

# # Ridge Regression 
# # windows()
# library(glmnet)

# output_RR <- file.path(output, "RidgeR_interpolate")
# dir.create(output_RR)
# output_MODLE_RR <- "./RidgeR_models"
# dir.create(output_MODLE_RR)
# raster_outRidgeR_pred <- file.path(data_dir_valid_masekd, "RidgeR_pred")
# dir.create(raster_outRidgeR_pred)

# # load interpolate data
# WPI_interpolate <- rast(file.path(data_dir_valid_masekd, "WPi_interpolate.tif")) 
# F_interpolate <- rast(file.path(data_dir_valid_masekd, "F_interpolate.tif")) 
# P_interpolate <- rast(file.path(data_dir_valid_masekd, "P_interpolate.tif")) 

# # prep landsat remot sensing data
# sampleloc_extent2_pj <- project(sampleloc_extent2, crs(valid_SEP_Extent_masked))
# mask_interpolate <- project(F_interpolate,  crs(valid_SEP_Extent_masked))

# RSdata_valid <- resample(valid_SEP_Extent_masked, mask_interpolate)
# RSdata_valid <- mask(RSdata_valid, sampleloc_extent2_pj)
# names(RSdata_valid) <- names(valid_SEP_Extent_masked) 
# #plot(RSdata_valid)

# WPI_interpolate_pj <- project(WPI_interpolate,  crs(RSdata_valid))
# RSdata_valid <- resample(RSdata_valid, WPI_interpolate_pj)
# all_RF_data <- c(RSdata_valid, WPI_interpolate_pj$var1.pred)

# rf_df_sample <- spatSample(all_RF_data, size = 50000, # sample 
#         method = "random", as.df = TRUE, na.rm = TRUE)
# rf_df_sample <- rf_df_sample[, !colnames(rf_df_sample) %in% c("CoastalAerosol", "SWIR2")]

# RSdata_valid <- RSdata_valid[[ !names(RSdata_valid) %in%  c("CoastalAerosol", "SWIR2") ]]

# # Data
# x <- as.matrix(rf_df_sample[, !colnames(rf_df_sample) %in% c("var1.pred")])   # Prädiktoren
# colnames(x)
# y <- rf_df_sample[, "var1.pred"]               # Zielvariable

# foldid <- sample(1:5, size = nrow(x), replace = TRUE) # räumliche CV
# cv_ridge <- cv.glmnet( x,  y, alpha = 0.3, standardize = TRUE, foldid = foldid, nfolds = 10) # cv für optimale regulariseirung

# # R² berechnen
# pred <- predict(cv_ridge, newx = x, s = "lambda.min")
# R2 <- 1 - sum((y - pred)^2) / sum((y - mean(y))^2)
# R2
# rmse <- sqrt(min(cv_ridge$cvm))
# r2   <- 1 - min(cv_ridge$cvm) / var(y)
# cat("CV RMSE:", rmse, "\n")
# cat("CV R2:", r2, "R2:", R2)

# cv_ridge$lambda.min
# plot(cv_ridge) # CV-Fehler vs. Lambda.

# # terra perdiction 
# ridge_fun <- function(model, data) {predict(model, newx = as.matrix(data), s = "lambda.min")}
# pred_raster <- terra::predict(RSdata_valid, cv_ridge, fun = ridge_fun)

# brks <- c( 0,1,2,3,4, 6, 10, 15, 30)
# plot(pred_raster, breaks = brks)

# coef(cv_ridge, s = "lambda.min") # koeffizenten / sensitvity

# ------------------------------
# ------------------------------------------
# ########################################################
# # first try 
# # linear Regression Model für
# # F anhand Green, Red, Swir, NIR, ,TSMlaili
# # P Red, Swir1 TSM laili
# # WPI Chl-a 

# lm_data <- merge(valid_SEP_points_VALUES, heavymetal, by = "Station")
# lm_data <- merge(lm_data, waterPolIndex, by = "Station")
# titel <- c(letters, LETTERS)
# # F  
# var_names <- c("Green", "Red", "SWIR1", "NIR" ,"TMS_laili", "log_Chl")
# col <- c("darkgreen", "red", "orange", "magenta","darkolivegreen", "coral")
# p_F <- list()
# for(i in 1:length(var_names)){
#   p_TMS <- ggplot(lm_data, aes(x = .data[[var_names[i]]], y = F_mean))+
#     geom_point(color = col[i]) +
#       geom_smooth(method = "lm", color = "black", se = FALSE) +
#     stat_cor(
#         aes(label = after_stat( paste0("italic(R)==", round(r, 2), "*','~~italic(R^2)==", round(rr, 2)))),# after_stat(paste0("italic(R)==", round(r, 2)))), # R [-1,1], zeigt richtung + stärke, jenäher 1/-1 desto stärker positiver/negativer Zusammenhang
#         parse = TRUE,  label.x.npc = "left", label.y.npc = "top") +

#     theme_bw() +
#     labs(
#       y = "F [mg/l]",
#       x = var_names[i] ,
#       title = paste0("(",titel[i], ") ", var_names[i])
#     ) 
#    p_TMS 
#   p_F[[i]] <-  p_TMS

# }
# p_F_all <- wrap_plots(p_F) +
#   plot_layout(guides = "collect") &
#   theme(legend.position = "bottom")
# p_F_all
# ggsave(file.path(output, "F_regression.png"), p_F_all, height = 10, width = 10, scale = 1.2)

# # P 
# var_names <- c("Green", "Red", "SWIR1", "NIR" ,"TMS_laili", "log_Chl")
# col <- c("darkgreen", "red", "orange", "magenta","darkolivegreen", "coral")
# p_P <- list()
# for(i in 1:length(var_names)){
#   p_TMS <- ggplot(lm_data, aes(x = .data[[var_names[i]]], y = P_mean))+
#     geom_point(color = col[i]) +
#       geom_smooth(method = "lm", color = "black", se = FALSE) +
#     stat_cor(
#         aes(label = after_stat( paste0("italic(R)==", round(r, 2), "*','~~italic(R^2)==", round(rr, 2)))),# after_stat(paste0("italic(R)==", round(r, 2)))), # R [-1,1], zeigt richtung + stärke, jenäher 1/-1 desto stärker positiver/negativer Zusammenhang
#         parse = TRUE,  label.x.npc = "left", label.y.npc = "top") +

#     theme_bw() +
#     labs(
#       y = "P [μg/l]",
#       x = var_names[i] ,
#       title = paste0("(",titel[i], ") ", var_names[i])
#     ) 
#    p_TMS 
#   p_P[[i]] <-  p_TMS

# }
# p_P_all <- wrap_plots(p_P) +
#   plot_layout(guides = "collect") &
#   theme(legend.position = "bottom")
# p_P_all
# ggsave(file.path(output, "P_regression.png"), p_P_all, height = 10, width = 10, scale = 1.2)


# # WPI
# var_names <- c("Green", "Red", "SWIR1", "NIR" ,"TMS_laili", "log_Chl")
# col <- c("darkgreen", "red", "orange", "magenta","darkolivegreen", "coral")
# p_Wpi <- list()
# for(i in 1:length(var_names)){
#   p_TMS <- ggplot(lm_data, aes(x = .data[[var_names[i]]], y = WPi))+
#     geom_point(color = col[i]) +
#       geom_smooth(method = "lm", color = "black", se = FALSE) +
#     stat_cor(
#         aes(label = after_stat( paste0("italic(R)==", round(r, 2), "*','~~italic(R^2)==", round(rr, 2)))),# after_stat(paste0("italic(R)==", round(r, 2)))), # R [-1,1], zeigt richtung + stärke, jenäher 1/-1 desto stärker positiver/negativer Zusammenhang
#         parse = TRUE,  label.x.npc = "left", label.y.npc = "top") +

#     theme_bw() +
#     labs(
#       y = "WPi",
#       x = var_names[i] ,
#       title = paste0("(",titel[i], ") ", var_names[i])
#     ) 
#    p_TMS 
#   p_Wpi[[i]] <-  p_TMS

# }
# p_wpi_all <- wrap_plots(p_Wpi) +
#   plot_layout(guides = "collect") &
#   theme(legend.position = "bottom")
# p_wpi_all
# ggsave(file.path(output, "Wpi_regression.png"), p_wpi_all, height = 10, width = 10, scale = 1.2)


# ######
# # linear Regression Model für:
# # F und P

# # F model 
# model_F_band <- lm(F_mean ~ Green + Red + SWIR1 + NIR + Blue +CoastalAerosol ,
#               data = lm_data)
# summary(model_F_band)

# model_F_TSM_laili <- lm(F_mean ~ TMS_laili,
#               data = lm_data)
# summary(model_F_TSM_laili) # beinhaltet blau

# model_F_band_blue <- lm(F_mean ~ Green + Red + SWIR1 + Blue,
#               data = lm_data)
# summary(model_F_band_blue) # blau verbessert das Model deutlich

# model_F_band_mini <- lm(F_mean ~ SWIR1 + NIR,
#               data = lm_data)
# summary(model_F_band_mini)
# AIC(model_F_band_blue, model_F_band_mini)

# # predict F

# F_map <- predict(valid_SEP_Extent, model_F_band)

# brks <- c( 0.3,5,10,15, 20)#seq(0.3, 20,  by = 4)
# plot(F_map, plg = list(title = "[mg/l]"))#,  breaks = brks)

# # tms laili 
# TSM <- lapp(valid_SEP_Extent[[c("Blue","Green","Red")]], fun = function(b, g, r){func_TSM_tassan(b, g, r)})

# names(TSM) <- "TMS_laili"
# F_map <- predict(TSM, model_F_TSM_laili)
# #brks <- c( 0.3,5,10,15, 20)#seq(0.3, 20,  by = 4)
# plot(F_map, plg = list(title = "[mg/l]"))#,  breaks = brks)

# # P model 
# lm_data_p <- lm_data
# lm_data_p$P_mean[is.na(lm_data_p$P_mean)] <- 130

# model_P_band <- lm(P_mean ~ Red + SWIR1 + Green + Blue + CoastalAerosol +NIR,
#               data = lm_data_p)
# summary(model_P_band)

# P_map <- predict(valid_SEP_Extent, model_P_band)
# brks <- c(0, 0.1,2,3,4, 6,20,40)#seq(0.3, 20,  by = 4)
# plot(P_map/1000, plg = list(title = "[mg/l]"),  breaks = brks)

# plot(valid_SEP_Extent_masked)

# plot(merge_valid_SUMMER08_EXTENT)
