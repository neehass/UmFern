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
indloc_pj <- project(vect(matrix(ind_loc, ncol=2), crs="EPSG:4326"), crs(all_pred_WPI[[1]]))

brks <- c( 0,1,2,3,4, 6, 10, 15, 30)#seq(0.3, 20,  by = 4)
png(file.path(output, paste0(name_RF, "_interpolate_RF_pred.png")), height = 800, width = 800)
par(mfrow=c(3,4))
for(y in 1:length(unique_years)){
    # max_abs <- max(abs(global(baseline_diff[[y]], "max", na.rm=TRUE)),
    #             abs(global(baseline_diff[[y]], "min", na.rm=TRUE)))

  plot(all_pred_WPI[[y]],
       main = paste("WPI", unique_years[y]),  breaks = brks, legend = FALSE
       #, plg=list(x="bottomleft", cex=1.3)
       )
  points(indloc_pj, col = "black", pch = 15, cex = 1)
}
dev.off()


# Klassifikation and NORMIERT ------------------------
name_RF <- "WPI_class_NORM"
folder_PRED_WPI_class <- file.path(folder_PRED, name_RF)

all_pred_WPI_Class <- lapply(list.files(folder_PRED_WPI_class, pattern = ".tif$"), function(x){rast(file.path(folder_PRED_WPI_class, x))})
indloc_pj <- project(vect(matrix(ind_loc, ncol=2), crs="EPSG:4326"), crs(all_pred_WPI[[1]]))

labs <- c("Not affected", "Slightly affected","Moderately affected","Strongly affected","Seriously affected")
png(file.path(output, paste0(name_RF, "_interpolate_RF_pred.png")), height = 500, width = 800)
par(mfrow=c(2,5))
for(y in 1:length(unique_years)){
    # max_abs <- max(abs(global(baseline_diff[[y]], "max", na.rm=TRUE)),
    #             abs(global(baseline_diff[[y]], "min", na.rm=TRUE)))
       levels(all_pred_WPI_Class[[y]]) <- data.frame(ID=1:5, label=labs)
  plot(all_pred_WPI_Class[[y]],
       main = paste(unique_years[y]), legend = FALSE)
    points(indloc_pj, col = "black", pch = 15, cex = 1)
}
dev.off()

png(file.path(output, paste0(name_RF, "_interpolate_RF_pred_LEGEND.png")), height = 200, width = 300)
levels(all_pred_WPI_Class[[1]]) <- data.frame(ID=1:5, label=labs)
plot(all_pred_WPI_Class[[1]],
       main = paste("WPI", unique_years[y]), legend = TRUE)

dev.off()

# less Variables Klassifikation and NORMIERT ------------------------
name_RF <- "WPI_class_NORM_lessVAR"
folder_PRED_WPI_class_lessVARS <- file.path(folder_PRED, name_RF)

all_pred_WPI_Class_lessVARS <- lapply(list.files(folder_PRED_WPI_class_lessVARS, pattern = ".tif$"), 
  function(x){rast(file.path(folder_PRED_WPI_class_lessVARS, x))})
indloc_pj <- project(vect(matrix(ind_loc, ncol=2), crs="EPSG:4326"), crs(all_pred_WPI_Class_lessVARS[[1]]))

labs <- c("Not affected", "Slightly affected","Moderately affected","Strongly affected","Seriously affected")
png(file.path(output, paste0(name_RF, "_interpolate_RF_pred.png")), height = 800, width = 800)
par(mfrow=c(3,4))
for(y in 1:length(unique_years)){
    # max_abs <- max(abs(global(baseline_diff[[y]], "max", na.rm=TRUE)),
    #             abs(global(baseline_diff[[y]], "min", na.rm=TRUE)))
       levels(all_pred_WPI_Class_lessVARS[[y]]) <- data.frame(ID=1:5, label=labs)
  plot(all_pred_WPI_Class_lessVARS[[y]],
       main = paste(unique_years[y]), legend = FALSE #,
       # plg=list(x="bottomleft", cex=1.3)
       )
      points(indloc_pj, col = "black", pch = 15, cex = 1)
}
dev.off()


# only TSM und Chla Klassifikation and NORMIERT ------------------------
name_RF <- "WPI_class_NORM_TSMChla"
folder_PRED_WPI_class_TSMChla <- file.path(folder_PRED, name_RF)

all_pred_WPI_Class_TSMChla <- lapply(list.files(folder_PRED_WPI_class_TSMChla, pattern = ".tif$"), 
          function(x){rast(file.path(folder_PRED_WPI_class_TSMChla, x))})
indloc_pj <- project(vect(matrix(ind_loc, ncol=2), crs="EPSG:4326"), 
  crs(all_pred_WPI_Class_TSMChla[[1]]))

labs <- c("Not affected", "Slightly affected","Moderately affected","Strongly affected","Seriously affected")
png(file.path(output, paste0(name_RF, "_interpolate_RF_pred.png")), height = 500, width = 800)
par(mfrow=c(2,5))
for(y in 1:length(all_pred_WPI_Class_TSMChla)){
    # max_abs <- max(abs(global(baseline_diff[[y]], "max", na.rm=TRUE)),
    #             abs(global(baseline_diff[[y]], "min", na.rm=TRUE)))
       levels(all_pred_WPI_Class_TSMChla[[y]]) <- data.frame(ID=1:5, label=labs)
  plot(all_pred_WPI_Class_TSMChla[[y]],
       main = paste(unique_years[y]), legend = FALSE)
  points(indloc_pj, col = "black", pch = 15, cex = 1)
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
baseline2013_class <- as.factor(baseline2013_class)

# comaprison 
comp_class <- list()
dir.create(file.path(output, "diffplots_propTab_Class"))
for(i in 1:(length(all_pred_WPI_Class)-1)){
  year1_org <- all_pred_WPI_Class[[i+1]]
  year1_org <- as.factor(year1_org)
  year1 <- year1_org
  names(year1) <- "pred"
  year1$baseline <- values(baseline2013_class)
  year1$change <- year1$baseline != year1$pred
  year1_change <- mask(year1, year1$change, maskvalues = FALSE)
  trans_vec <- paste(values(year1_change$baseline), "→", values(year1_change$pred))
  valid <- (!is.na(trans_vec)) & (trans_vec != "NaN → NaN") & (trans_vec != "NA → NA")
  year1_change$transition <- NA
  year1_change$transition[valid] <- as.factor(trans_vec[valid])
  cats_df <- data.frame(
  ID = 1:length(levels(factor(trans_vec[valid]))),
  category = levels(factor(trans_vec[valid]))
  )

  levels(year1_change$transition) <- cats_df

  finaltransition <- year1_change$transition

  comp_class[[i]] <- finaltransition

  # Prob tabelle pro Klasse
   df <- as.data.frame(year1_org, xy = TRUE)  # year1 Raster
  df <- df[!is.na(values(baseline2013_class)) & !is.na(values(year1_org)), ]
  names(df)[3] <- "pred"
  df$baseline <- values(baseline2013_class)[!is.na(values(baseline2013_class))]
  df <- df[complete.cases(df), ]

  prop_tab <- prop.table(table(df$baseline, df$pred), 1) # Kontingenztabelle (Konfusionsmatrix)
  write.csv(prop_tab, file.path(output, "diffplots_propTab_Class", 
      paste0(unique_years[i+1], "_diffPredciton_CLASS_baseline2013.csv")))

  if(i != (length(unique_years) -1)){
    rm(year1, year1_change, trans_vec, valid, finaltransition, cats_df) # behalte letzte
  }
  
}
fixed_colors2 <- c( # chat gpt
# to class 1 (oranges)
  "5 → 1" = "#993404",
  "4 → 1" = "#d95f0e",
  "3 → 1" = "#fe9929",
  "2 → 1" = "#fec44f",
  # to class 2 (greens)
  "5 → 2" = "#238b45",
  "4 → 2" = "#41ab5d",
  "3 → 2" = "#74c476",
  "1 → 2" = "#c7e9c0",
# to class 3 (purples)
  "5 → 3" = "#6a51a3",
  "4 → 3" = "#807dba",
  "2 → 3" = "#9e9ac8",
  "1 → 3" = "#dadaeb",
 # to class 4 (blues)
  "5 → 4" = "#2171b5",
  "3 → 4" = "#4292c6",
  "2 → 4" = "#6baed6",
  "1 → 4" = "#c6dbef",
  # to class 5 (magenta/pink)
"1 → 5" = "#8e0152",
"2 → 5" = "#c51b8a",
"3 → 5" = "#f768a1",
"4 → 5" = "#fde0ef"
)
cats <- levels(comp_class[[1]])[[1]]
cols <- fixed_colors2[cats$category]

png(file.path(output, paste0("diffWPI_CLASS", "_interpolate_RF_pred.png")), height = 500, width = 800)
par(mfrow=c(2,5))
for(y in 1:length(comp_class)){
  # Row 1, first slot empty
  if(y == 1){
    plot.new()  # leave first slot blank
  }
  
  # Row 2, first plot is class 5
  if(y == 2){
    plot(comp_class[[y]],
         main = paste(unique_years[y+1]),
         legend = FALSE,
         col = cols)
     points(indloc_pj, col = "black", pch = 15, cex = 1)
  } else {
    plot(comp_class[[y]],
         main = paste(unique_years[y+1]),
         legend = FALSE,
         col = cols)
     points(indloc_pj, col = "black", pch = 15, cex = 1)
  }
}
dev.off()

png(file.path(output, paste0("diffWPI_CLASS", "LEGEND.png")), height = 500, width = 800)
terra::plot(comp_class[[1]], col = cols,  plg=list( cex=2))
dev.off()

# -----------------------------------------------------------
# only TSM and Chla baseline 2013 difference PREDICTIONS Klassifikation  ------------------
# ---------------------------------------------------------
baseline2013_classTSMChla <- all_pred_WPI_Class_TSMChla[[1]]
baseline2013_classTSMChla <- as.factor(baseline2013_classTSMChla)

# comaprison 
comp_class_TSMChla <- list()
dir.create(file.path(output, "diffplots_propTab_TSMChla"))
for(i in 1:(length(all_pred_WPI_Class_TSMChla)-1)){
  year1_org <- all_pred_WPI_Class_TSMChla[[i+1]]
  year1_org <- as.factor(year1_org)
  year1 <- year1_org
  names(year1) <- "pred"
  year1$baseline <- values(baseline2013_classTSMChla)
  year1$change <- year1$baseline != year1$pred
  year1_change <- mask(year1, year1$change, maskvalues = FALSE)
  trans_vec <- paste(values(year1_change$baseline), "→", values(year1_change$pred))
  valid <- (!is.na(trans_vec)) & (trans_vec != "NaN → NaN") & (trans_vec != "NA → NA")
  year1_change$transition <- NA
  year1_change$transition[valid] <- as.factor(trans_vec[valid])
  cats_df <- data.frame(
  ID = 1:length(levels(factor(trans_vec[valid]))),
  category = levels(factor(trans_vec[valid]))
  )

  levels(year1_change$transition) <- cats_df

  finaltransition <- year1_change$transition

  comp_class_TSMChla[[i]] <- finaltransition

  # Prob tabelle pro Klasse
   df <- as.data.frame(year1_org, xy = TRUE)  # year1 Raster
  df <- df[!is.na(values(baseline2013_classTSMChla)) & !is.na(values(year1_org)), ]
  names(df)[3] <- "pred"
  df$baseline <- values(baseline2013_classTSMChla)[!is.na(values(baseline2013_classTSMChla))]
  df <- df[complete.cases(df), ]

  prop_tab <- prop.table(table(df$baseline, df$pred), 1) # Kontingenztabelle (Konfusionsmatrix)
  write.csv(prop_tab, file.path(output, "diffplots_propTab_TSMChla", 
      paste0(unique_years[i+1], "_diffPredciton_CLASS_baseline2013.csv")))

  if(i != (length(unique_years) -1)){
    rm(year1, year1_change, trans_vec, valid, finaltransition, cats_df) # behalte letzte
  }
  
}
fixed_colors <- c( # chat gpt
# to class 1 (oranges)
  "5 → 1" = "#993404",
  "4 → 1" = "#d95f0e",
  "3 → 1" = "#fe9929",
  "2 → 1" = "#fec44f",
  # to class 2 (greens)
  "5 → 2" = "#238b45",
  "4 → 2" = "#41ab5d",
  "3 → 2" = "#74c476",
  "1 → 2" = "#c7e9c0",
# to class 3 (purples)
  "5 → 3" = "#6a51a3",
  "4 → 3" = "#807dba",
  "2 → 3" = "#9e9ac8",
  "1 → 3" = "#dadaeb",
 # to class 4 (blues)
  "5 → 4" = "#2171b5",
  "3 → 4" = "#4292c6",
  "2 → 4" = "#6baed6",
  "1 → 4" = "#c6dbef",
  # to class 5 (magenta/pink)
"1 → 5" = "#8e0152",
"2 → 5" = "#c51b8a",
"3 → 5" = "#f768a1",
"4 → 5" = "#fde0ef"
)

cats <- levels(comp_class_TSMChla[[1]])[[1]]
cols <- fixed_colors_sorted[cats$category]
png(file.path(output, paste0("diffWPI_CLASS_TSMChla", "_interpolate_RF_pred.png")), height = 500, width = 800)
par(mfrow=c(2,5))
for(y in 1:length(comp_class_TSMChla)){
  # Row 1, first slot empty
  if(y == 1){
    plot.new()  # leave first slot blank
  }
  
  # Row 2, first plot is class 5
  if(y == 2){
    plot(comp_class_TSMChla[[y]],
         main = paste(unique_years[y+1]),
         legend = FALSE,
         col = cols)
     points(indloc_pj, col = "black", pch = 15, cex = 1)
  } else {
    plot(comp_class_TSMChla[[y]],
         main = paste(unique_years[y+1]),
         legend = FALSE,
         col = cols)
     points(indloc_pj, col = "black", pch = 15, cex = 1)
  }
}
dev.off()

png(file.path(output, paste0("diffWPI_CLASS_TSMChla", "LEGEND.png")), height = 500, width = 800)
terra::plot(comp_class_TSMChla[[1]], col = cols,  plg=list( cex=2))
dev.off()



# fixed_colors <- c( # chat gpt
# # From class 1 (magenta/pink)
# "1 → 2" = "#8e0152",
# "1 → 3" = "#c51b8a",
# "1 → 4" = "#f768a1",
# "1 → 5" = "#fde0ef", 
  
#   # From class 2 (blues)
#   "2 → 1" = "#2171b5",
#   "2 → 3" = "#4292c6",
#   "2 → 4" = "#6baed6",
#   "2 → 5" = "#c6dbef",
  
#   # From class 3 (purples)
#   "3 → 1" = "#6a51a3",
#   "3 → 2" = "#807dba",
#   "3 → 4" = "#9e9ac8",
#   "3 → 5" = "#dadaeb",
  
#   # From class 4 (greens)
#   "4 → 1" = "#238b45",
#   "4 → 2" = "#41ab5d",
#   "4 → 3" = "#74c476",
#   "4 → 5" = "#c7e9c0",
  
#   # From class 5 (oranges)
#   "5 → 1" = "#993404",
#   "5 → 2" = "#d95f0e",
#   "5 → 3" = "#fe9929",
#   "5 → 4" = "#fec44f"
# )