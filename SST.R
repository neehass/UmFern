library(ggplot2)
library(ggpubr)
library(patchwork)
library(terra)
library(dplyr)
library(sf)
library(mapview)
library(ncdf4)

getwd()
setwd("./Phosphate")
source('./R-scripts/hel-func.R', chdir = TRUE)

output <- "./R-scripts/SST-outputs"
dir.create(output)

data_elKateb_shap <- "./elKateb2018_shp"

SST_data <- "./NOAA_OI_SST_V2_High_Resolution"
SST_monthly_hr <- rast(file.path(SST_data, "sst.mon.mean.nc"))
windows()
plot(SST_monthly_hr[[1]])

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
# crop SST data  ------------------
# ---------------------------------------------------------

SST_crop_monthly <- crop(SST_monthly_hr, sampleloc_extent)

plot(SST_crop_monthly[[1]])
plot(gulf_shp, add = TRUE, col = "darkgrey", alpha = 0.2, type = "l")
points(sampleLoc_points, col = "cyan")
points(ind_loc, col = "red", pch = 15)

writeRaster(SST_crop_monthly, file.path(SST_data,"sst.mon.mean_cropped_GABES.tif") , overwrite = TRUE)

time(SST_crop_monthly)
class(time(SST_crop_monthly))

dates <- time(SST_crop_monthly)

idx <- which(
  dates >= as.Date("1984-01-01") &
  dates <= as.Date("2025-12-01")
)
SST_crop_1984_2025 <- SST_crop_monthly[[idx]]
plot(SST_crop_1984_2025)

# ---------------------------------------------------------
# trendanlyse SST data  ------------------
# ---------------------------------------------------------
months <- as.integer(format(time(SST_crop_1984_2025), "%m"))

SST_clim <- tapp(
  SST_crop_1984_2025,
  index = months,
  fun = mean,
  na.rm = TRUE
)

# anomalie:
SST_anom <- SST_crop_1984_2025 - SST_clim[[months]]

# trend
t <- 1:nlyr(SST_anom)

trend_fun <- function(x) {
  ok <- !is.na(x)
  
  # Mindestanzahl gültiger Werte (z.B. > 50 %)
  if (sum(ok) < length(x) * 0.5) return(NA)
  
  coef(lm(x[ok] ~ t[ok]))[2]
}

SST_trend <- app(SST_anom, trend_fun)

# sign 
trend_fun_sig <- function(x) {
  ok <- !is.na(x)
  if (sum(ok) < 24) return(c(NA, NA))
  
  fit <- lm(x[ok] ~ t[ok])
  c(coef(fit)[2], summary(fit)$coefficients[2,4])
}

SST_trend_sig <- app(SST_anom, trend_fun_sig)
names(SST_trend_sig) <- c("slope", "p_value")

# sign trend pro Jahr: 

SST_sig <- SST_trend_sig[[2]] < 0.05
SST_trend_year <- SST_trend * 12
SST_trend_sig_masked <- mask(SST_trend_year, SST_sig, maskvalues = FALSE)

cols <- colorRampPalette(c("blue", "white", "red"))(100)
max_val <- max(abs(values(SST_trend_year)), na.rm = TRUE)

png(filename = file.path(output,"SST_trend_significant.png"), width = 2000, height = 2500, res = 300)
plot(SST_trend_sig_masked, col = cols, main = paste("sign. SST Trend (°C pro Jahr)\n von", time(SST_crop_1984_2025)[1], "bis" , time(SST_crop_1984_2025)[length(time(SST_crop_1984_2025))]),
 zlim = c(-max_val, max_val))
lines(gulf_shp, col = "black", alpha = 0.2)
# points(sampleLoc_points, col = "cyan")
points(ind_loc, col = "black", pch = 15)
dev.off()
