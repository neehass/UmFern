# ---------------------------------------------------------
# Differences over time series ------------------
# ---------------------------------------------------------

# packages
library(terra)
library(ggplot2)
library(dplyr)
library(sf)
library(mapview)
library(ranger)
library(tidyr)

getwd()
setwd("./Phosphate")
source('./R-scripts/final-scripts/hel-func.R', chdir = TRUE)

output <- "./R-scripts/output/04_timeseries_results_diff_LC08-09"
dir.create(output, recursive = TRUE)

folder_PRED <- "./final_RF_NORM_pred_LC08-09"
# ----------------------------------------------------------------
# phosphat industry locations ------
# ------------------------------------------------
ind_loc <- t(data.frame(  # lon, lat
  "Gabes" = c(10.095439, 33.916444), # since 1979, expaned 1985
                "Sfax" = c(10.723889, 34.702778), # Triple Super Phosphate (TSP), not DAP, since 1952
                "Skhira" = c(10.148889, 34.3475))) # seit 2013 in betireb
colnames(ind_loc) <- c("lon", "lat")

# ---------------------------------------------------------
# load predictions with WPI model ------------------
# ---------------------------------------------------------
# REGRESSION -------------------
name_RF <- "WPI_NORM"
folder_PRED_WPI <- file.path(folder_PRED, name_RF)

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

all_pred_WPI_Class <- lapply(list.files(folder_PRED_WPI_class, pattern = ".tif$"), function(x){rast(file.path(folder_PRED_WPI_class, x))})

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

# -----------------------------------------------------------
# CV REGRESSION ------------------
# ---------------------------------------------------------
rast_all_pred_WPI <- rast(all_pred_WPI)

cv_pred_REG <- app(rast_all_pred_WPI, function(x){
  sd(x, na.rm=TRUE) / mean(x, na.rm=TRUE)
})
col_fun <- colorRampPalette(c("white", "red"))
png(file.path(output, "CV_PredictionTimeseries2013-2025.png"), height = 350, width = 300)
plot(cv_pred_REG, col = col_fun(10), plg = list(title = "CV Predictions\n 2013-2025"))
dev.off()

# -----------------------------------------------------------
# baseline 2013 difference ------------------
# ---------------------------------------------------------
baseline2013 <- all_pred_WPI[[1]]
diff <- func_baseline_diff(all_pred_WPI, baseline2013)
years <- unique_years[2:length(unique_years)]

df_list <- lapply(1:(length(unique_years)-1), function(y) {
  r <- diff[[y+1]]
  df <- as.data.frame(r, xy = TRUE)
  df <- df %>% 
    pivot_longer(cols = names(r), names_to = "band", values_to = "value")   
     
df$year <- unique_years[y+1]
  return(df)
})

plot_df <- bind_rows(df_list)
lim <- round(max(abs(plot_df$value)))

p_diff_reg <- ggplot(plot_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-lim, lim)) +  # gleiche Farbskala
     coord_equal() +
  
  facet_wrap(~ year, ncol = 5) +
  theme_minimal() +
  labs(fill = "Difference\ to 2013", x = "", y= "") +
  theme(legend.position = "bottom", 
  axis.text.x = element_text(angle = 45, hjust = 1)) 


p_diff_reg
ggsave(file.path(output, "diffPredciton_baseline2013.png"), p_diff_reg, 
          width = 7, height = 6, scale = 1.2)


# -----------------------------------------------------------
# raw Index unterschiede ------------------
# ---------------------------------------------------------
data_crop_dir <- "./final_landsat-SEPTEMBER_PIF_LC08-09"
files_raster <- list.files(data_crop_dir, recursive = TRUE, pattern="^WATER_NORMIERT")
parts <- strsplit(files_raster, "_|\\.")
years <- sapply(parts, `[`, 4) 
unique_years <- unique(years)
print(unique_years) # von 2013 bis 2025 >> 12Jahre
length(unique_years) # 10 Scenen 

raster_years <- lapply(files_raster, function(x){rast(file.path(data_crop_dir, x))})
names(raster_years) <- paste0("year_",unique_years)
raster_years <- lapply(raster_years, function(x) {x[[ !names(x) %in% c("CoastalAerosol", "SWIR2") ]]})

raster_yearsINDEX <- lapply(raster_years, function(x){
        x$TSM_laili <- func_TMS_Laili(x$Blue, x$Red)
     x$logChla   <- func_log_Chl(x$Red, x$Green)

     x$NDWI <- func_NDWI(x$Green, x$NIR) 
     x$TI   <- func_TI(x$NIR, x$Green)

     x$sediment <- func_sedi(x$Red, x$Blue)
     x$NDSSI    <- func_NDSSI(x$Red, x$NIR)

     x$SWIR1_NIR <- x$SWIR1 / x$NIR

     x
})

# Baseline 2013 difference ------------------
baseline2013_INDEX <- raster_yearsINDEX[[1]]
diff_INDEX <- func_baseline_diff(raster_yearsINDEX, baseline2013_INDEX)
years <- unique_years[2:length(unique_years)]

dfINDEx_list <- lapply(1:(length(unique_years)-1), function(y) {
  r <- diff_INDEX[[y+1]]
  df <- as.data.frame(r, xy = TRUE)
  df <- df %>% 
    pivot_longer(cols = names(r), names_to = "band", values_to = "value")   
     
df$year <- unique_years[y+1]
  return(df)
})

plot_dfINDEx <- bind_rows(dfINDEx_list)

index <- c("TSM_laili", "logChla", "NDWI", "TI", "sediment", "NDSSI", "SWIR1_NIR")
units <- c("[mg/l]", "[mg/l]", "", "", "", "", "")
for(i in 1:length(index)){
     print(index[i])
     TSM_df <- plot_dfINDEx[plot_dfINDEx$band == index[i], ]
     if(max(abs(TSM_df$value)) <= 1){
               lim <- round(max(abs(TSM_df$value)), 2)
     } else {
          lim <- round(max(abs(TSM_df$value))) 
     }

     p_diff_TSM <- ggplot(TSM_df, aes(x = x, y = y, fill = value)) +
     geom_raster() +
     scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-lim, lim)) +  # gleiche Farbskala
          coord_equal() +
     
     facet_wrap(~ year, ncol = 5) +
     theme_minimal() +
     labs(fill = paste("Difference",units[i], "\n to 2013"), x = "", y = "") +
     theme(legend.position = "bottom", 
     axis.text.x = element_text(angle = 45, hjust = 1)) 

     p_diff_TSM
     ggsave(file.path(output, paste0(index[i], "_diffINDEX_baseline2013.png")), p_diff_TSM, 
               width = 7, height = 6, scale = 1.2)
}


