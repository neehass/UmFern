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
library(cowplot)

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
# 2015 02.08. exluded !
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
png(file.path(output, paste0(name_RF, "_interpolate_RF_pred.png")), height = 500, width = 800)
par(mfrow=c(2,5))
for(y in 1:length(unique_years)){
    # max_abs <- max(abs(global(baseline_diff[[y]], "max", na.rm=TRUE)),
    #             abs(global(baseline_diff[[y]], "min", na.rm=TRUE)))
       levels(all_pred_WPI_Class[[y]]) <- data.frame(ID=1:5, label=labs)
  plot(all_pred_WPI_Class[[y]],
       main = paste(unique_years[y]), legend = FALSE)
}
dev.off()

png(file.path(output, paste0(name_RF, "_interpolate_RF_pred_LEGEND.png")), height = 200, width = 300)
levels(all_pred_WPI_Class[[1]]) <- data.frame(ID=1:5, label=labs)
plot(all_pred_WPI_Class[[1]],
       main = paste("WPI", unique_years[y]), legend = TRUE)

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
# CV Classification ------------------
# ---------------------------------------------------------
rast_all_pred_WPI_class <- rast(all_pred_WPI_Class)

cv_pred_Class<- app(rast_all_pred_WPI_class, function(x){
  sd(x, na.rm=TRUE) / mean(x, na.rm=TRUE)
})
col_fun <- colorRampPalette(c("white", "red"))
png(file.path(output, "CV_Class_PredictionTimeseries2013-2025.png"), height = 350, width = 300)
plot(cv_pred_Class, col = col_fun(10), plg = list(title = "CV Predictions\n 2013-2025"))
dev.off()

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
     
     facet_wrap(~ year, ncol = 4) +
     theme_minimal() +
     labs(fill = paste("Difference",units[i], "\n to 2013"), x = "", y = "") +
     theme(legend.position = "bottom", 
     axis.text.x = element_text(angle = 45, hjust = 1)) 

     # p_diff_TSM
     ggsave(file.path(output, paste0(index[i], "_diffINDEX_baseline2013.png")), p_diff_TSM, 
               width = 4, height = 6, scale = 1.2)
}

# -----------------------------------------------------------
# baseline 2013 difference PREDICTIONS Regression ------------------
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
  
  facet_wrap(~ year, ncol = 4) +
  theme_minimal() +
  labs(fill = "Difference\ to 2013", x = "", y= "") +
  theme(legend.position = "bottom", 
  axis.text.x = element_text(angle = 45, hjust = 1)) 


# p_diff_reg
ggsave(file.path(output, "diffPredciton_baseline2013.png"), p_diff_reg, 
          width = 4, height = 6, scale = 1.2)

# -----------------------------------------------------------
# baseline 2013 difference PREDICTIONS Klassifikation  ------------------
# ---------------------------------------------------------
baseline2013_class <- all_pred_WPI_Class[[1]]

ind_loc <- data.frame(
   name = c("Gabes", "Skhira","Sfax"),
  lon = c(10.09544, 10.14889, 10.72389),
  lat = c(33.91644, 34.34750, 34.70278)
)
ind_sf <- st_as_sf(ind_loc, coords = c("lon","lat"), crs = 4326)  # WGS84
ind_sf_proj <- st_transform(ind_sf, crs = st_crs(baseline2013_class)$wkt)
df_indloc <- cbind(ind_sf_proj, st_coordinates(ind_sf_proj))

# comaprison 
p_comp_class <- list()
fixed_colors <- c( # chat gpt
  "1 → 2" = "#fee5d9",
  "1 → 3" = "#fcae91",
  "1 → 4" = "#fb6a4a",
  "1 → 5" = "#de2d26",
  "2 → 1" = "#edf8fb",
  "2 → 3" = "#b2e2e2",
  "2 → 4" = "#66c2a4",
  "2 → 5" = "#238b45",
  "3 → 1" = "#f2f0f7",
  "3 → 2" = "#cbc9e2",
  "3 → 4" = "#9e9ac8",
  "3 → 5" = "#6a51a3",
  "4 → 1" = "#ffffcc",
  "4 → 2" = "#c2e699",
  "4 → 3" = "#78c679",
  "4 → 5" = "#238443",
  "5 → 1" = "#fef0d9",
  "5 → 2" = "#fdcc8a",
  "5 → 3" = "#fc8d59",
  "5 → 4" = "#d7301f"
)
dir.create(file.path(output, "diffplots_propTab"))
for(i in 1:(length(all_pred_WPI_Class)-1)){
  i <- 1
  year1 <- all_pred_WPI_Class[[i+1]]
  compareGeom(baseline2013_class, year1)

  df <- as.data.frame(year1, xy = TRUE)  # year1 Raster
  df <- df[!is.na(values(baseline2013_class)) & !is.na(values(year1)), ]
  names(df)[3] <- "pred"
  df$baseline <- values(baseline2013_class)[!is.na(values(baseline2013_class))]
  df <- df[complete.cases(df), ]

  df$change <- df$baseline != df$pred
  df_change <- df[df$change, ]
  # df_change <- na.omit(df_change)
  df_change$transition <- paste(df_change$baseline, "→", df_change$pred)
  df_change$transition <- as.factor(df_change$transition)
  df_change$transition_num <- as.numeric(df_change$transition)
  unique(df_change$transition)
  names(df_change)

  # rast_change <- rast(df_change[, c("x","y","transition_num")], type = "xyz")
  # cats_df <- data.frame(
  # ID = 1:length(levels(df_change$transition)),   # numeric codes
  # category = levels(df_change$transition))       # factor labels
  # levels(rast_change) <- levels(df_change$transition)
  # plot(rast_change, col = rainbow(length(levels(df_change$transition))))

  # library(RColorBrewer)
  # par(mar=c(3,4,2,2))
  # display.brewer.all()

  # pal1 <- RColorBrewer::brewer.pal(8, "Blues")
  # pal2 <- RColorBrewer::brewer.pal(8, "YlGn")
  # pal3 <- rev(RColorBrewer::brewer.pal(4, "PuRd"))

  # combined <- colorRampPalette(c(pal1, pal2, pal3))(
  #   length(unique(df_change$transition))
  # )
  p <- ggplot(df_change, aes(x = x, y = y)) +
     geom_tile(aes(fill = transition, color = transition), linewidth = 1) +  # geom_raster() +
     scale_color_manual(values = fixed_colors) +
    scale_fill_manual(values = fixed_colors) +
      geom_point(data = df_indloc[1,], aes(x = X, y = Y), color = "red", size = 2) +
    coord_equal() + 
    theme_minimal() + 
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) + # wie terra)
    labs(title = unique_years[i+1], 
        color = "Class\ntransition", x = "", y= "", fill = "Class\ntransition")
#  p
  ggsave(file.path(output, "diffplots", paste0(unique_years[i+1], "_diffPredciton_CLASS_baseline2013.png")), p, 
          width = 5, height = 6, scale = 1.2) #3 , 4

  p_comp_class[[i]] <- p

  prop_tab <- prop.table(table(df$baseline, df$pred), 1) # Kontingenztabelle (Konfusionsmatrix)
  write.csv(prop_tab, file.path(output, "diffplots_propTab", 
      paste0(unique_years[i+1], "_diffPredciton_CLASS_baseline2013.csv")))

  if(i != (length(unique_years) -1)){
    rm(df, df_change) # behalte letzte
  }
  
}


prop.table(table(df$baseline, df$pred), 1) # Kontingenztabelle (Konfusionsmatrix)



# zusammenführen -----
xlim <- range(df_change$x)
ylim <- range(df_change$y)
legend <- get_legend(p_comp_class[[2]])
legend
plots_no_legend <- lapply(p_comp_class, function(p) 
  p + theme(legend.position = "none") +
      coord_cartesian(xlim = xlim, ylim = ylim)  # keine coord_equal()
)
combined <- plot_grid(plotlist = plots_no_legend, ncol = 4, align = "hv")
final <- plot_grid(combined, legend, ncol = 1, rel_heights = c(1, 0.1))
# final
ggsave(file.path(output, "diffPredciton_CLASS_baseline2013.png"), final, 
          width = 5, height = 6, scale = 1.2)
