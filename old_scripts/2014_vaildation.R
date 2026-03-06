# El Kateb 2018
# Phosphorus and organic matter anlysis around Gabes
# Seven coastal stations were sampled at >1 m water depth


# packages
# update.packages(ask = FALSE, checkBuilt = TRUE)
# install.packages("rlang", dependencies = TRUE)
# install.packages("ggplot2", dependencies = TRUE)
# .libPaths()
# q()
# packageVersion("rlang")
library(ggplot2)
library(ggpubr)
library(patchwork)
library(terra)
library(dplyr)
library(sf)
library(mapview)
library(ggnewscale)

getwd()
setwd("./Phosphate")
source('./R-scripts/hel-func.R', chdir = TRUE)

output <- "./R-scripts/2014Valid-outputs"
dir.create(output)

data_elKateb_shap <- "./elKateb2018_shp"
data_elKateb <- "./elKateb2018_tabels"

data_dir <- "./landsat_daten"
data_dir_valid <- "./landsat2014"
data_dir_valid_masekd <- "./landsat2014-masked"

dir.create(data_dir_valid_masekd)
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

# plot industrie loc and sample locs
plot(sampleloc_extent)
lines(gulf_shp, add = TRUE, col = "darkgrey", alpha = 0.2)
points(sampleLoc_points, col = "cyan")
points(ind_loc, col = "red", pch = 15)
# ---------------------------------------------------------
# load landsat data ------------------
# ---------------------------------------------------------
files_band <- list.files(data_dir_valid, recursive = TRUE, pattern="T1_SR_B[0-9]{1,2}\\.TIF$")
parts <- strsplit(files_band, "_")
satNR <- unique(sapply(parts, `[`, 1))
dates <- sapply(parts, `[`, 4)
unique_dates <- unique(dates)

files_band2 <- list.files(data_dir, recursive = TRUE, pattern="T1_SR_B[0-9]{1,2}\\.TIF$")
parts2 <- strsplit(files_band2, "_")
satNR2 <- unique(sapply(parts2, `[`, 1))
dates2 <- sapply(parts, `[`, 4)
unique_dates2 <- unique(dates2)


# band overview:
satNR_bands <- list(
   "LT05" = list("B1" = "Blue", "B2" = "Green", "B3" = "Red", "B4" = "NIR", "B5" = "SWIR1", "B6" = "Thermal", "B7" = "SWIR2"),
   "LE07" =  list("B1" = "Blue", "B2" = "Green", "B3" = "Red", "B4" = "NIR", "B5" = "SWIR1", "B6" = "Thermal", "B7" = "SWIR2", "B8" = "Panchromatic"),
   "LC08" =  list("B1" = "CoastalAerosol", "B2" = "Blue", "B3" = "Green", "B4" = "Red", "B5" = "NIR", "B6" = "SWIR1", "B7" = "SWIR2", "B8" = "Panchromatic", "B9" = "TIRS1", "B10" = "TIRS2"),
   "LC09" =  list("B1" = "CoastalAerosol", "B2" = "Blue", "B3" = "Green", "B4" = "Red", "B5" = "NIR", "B6" = "SWIR1", "B7" = "SWIR2", "B8" = "Panchromatic", "B9" = "TIRS1", "B10" = "TIRS2")
)

# load SUMMER -----------------
timestep <- 2014

vaild_SUMMER <- func_timestep_sel(timestep, files_band, satNR_bands, unique_dates, data_dir_valid, "SUMMER")

vaild_SUMMER2 <- func_timestep_sel(timestep, files_band2, satNR_bands, unique_dates2, data_dir, "SUMMER")

plot(vaild_SUMMER[[1]])
plot(vaild_SUMMER[[2]])
plot(vaild_SUMMER[[3]])

merge_valid_SUMMER08 <- mosaic(vaild_SUMMER[[2]],vaild_SUMMER[[3]])
plot(merge_valid_SUMMER08)

# ---------------------------------------------------------
# crop sample locations ------------------
# ---------------------------------------------------------
# crop sample loc extent -----
sampleloc_extent_pj <- project(sampleloc_extent, crs(merge_valid_SUMMER08))

merge_valid_SUMMER08_EXTENT <- mask(merge_valid_SUMMER08, sampleloc_extent_pj)
plot(merge_valid_SUMMER08_EXTENT$CoastalAerosol)
points(sampleloc_pj, col = "cyan")
writeRaster(merge_valid_SUMMER08_EXTENT, file.path(data_dir_valid_masekd,"EXTENT_valid_sampleLoc_SUMMER08.tif") , overwrite = TRUE)

# crop sample loctions points with small buffer  -----
gulf_shp_pj <- project(gulf_shp, crs(merge_valid_SUMMER08_EXTENT))
sampleloc_pj <- project(sampleLoc_points, crs(merge_valid_SUMMER08_EXTENT))
buffer <- buffer(sampleloc_pj, width = 550) # meter

merge_valid_SUMMER08_SAMPELPOINTS <- mask(merge_valid_SUMMER08_EXTENT, buffer)
plot(merge_valid_SUMMER08_SAMPELPOINTS$Blue)
points(sampleloc_pj, col = "cyan", pch = 2)
lines(gulf_shp_pj, col = "darkgrey", alpha = 0.2)

writeRaster(merge_valid_SUMMER08_SAMPELPOINTS, file.path(data_dir_valid_masekd,"POINTS_valid_sampleLoc_SUMMER08.tif") , overwrite = TRUE)

# just values at sample station
vals <- extract(merge_valid_SUMMER08_EXTENT, buffer, fun = mean)
vals$Sample <- buffer$Station
vals <- na.omit(vals)
plot(vals$CoastalAerosol)

vals$TMS_tassan <- func_TSM_tassan(vals$Blue, vals$Green, vals$Red)
vals$TMS_laili <- func_TMS_Laili(vals$Blue, vals$Red)
vals$TMS_budhiman <- func_TMS_Budhiman(vals$Red)
vals$log_Chl  <- func_log_Chl(vals$Red, vals$Green)
vals$salin <- func_salin(vals$Blue, vals$Green, vals$Red)

write.csv(vals, file.path(data_dir_valid_masekd,"POINTS_values_valid_sampleLoc_SUMMER08.csv"))

# ---------------------------------------------------------
# Regression analyse samples and landsat refelctions ------------------
# PG-02 = July 2014
# Bi = July 2014
# ---------------------------------------------------------
#load data 
vals <- read.csv(file.path(data_dir_valid_masekd,"POINTS_values_valid_sampleLoc_SUMMER08.csv"), row.names = 1)
merge_valid_SUMMER08_SAMPELPOINTS <- rast(file.path(data_dir_valid_masekd,"POINTS_valid_sampleLoc_SUMMER08.tif"))
merge_valid_SUMMER08_EXTENT <- rast(file.path(data_dir_valid_masekd,"EXTENT_valid_sampleLoc_SUMMER08.tif"))

# load el Kateb 2018 data
list.files(data_elKateb)
mineralComp <- read.csv(file.path(data_elKateb, "elKateb2018_mineralComp.csv"), sep = ";", header = TRUE)
phosRES <- read.csv(file.path(data_elKateb, "elKateb2018_phosRes.csv"), sep = ";", skip = 1) # Total phosphorus (Ptotal) represents the sum of all phosphorus reservoirs (Ptotal = Pex + PFe + Pauthi + Pdetr + Porg)
geoCHEM <- read.csv(file.path(data_elKateb, "elKateb2018_geochem.csv"), sep = ";") 

# allKateb <- merge(mineralComp, phosRES, by = "Sample")
allKateb_chem <-  merge(phosRES, geoCHEM, by = "Sample")
allKateb_chem <-  merge(allKateb_chem, vals, by = "Sample")
allKateb_chem$Loc <- sapply(strsplit(allKateb_chem$Sample, "-"), `[`, 1)
allKateb_chem$LocNR <- sapply(strsplit(allKateb_chem$Sample, "-"), `[`, 2)
# TOC = Total organic carbon 
# delta 13 C -25% stärkeren Anteil terrestrischer organischer Stoffe >> oft europhierung anthorpogen >> abnehmend mit entfernung >> evtl korreliert mit total suspenden Matter oder Chlorphyl

# für analyse TSM
# Ptotal (µmol P/g) → Hohe Gesamt-P-Gehalte deuten oft auf hohen Partikeltransport / Eutrophierung / Resuspension hin
# >> Pdet (detrital P) → Repräsentiert mineralisches, terrigenes Material → Sehr guter Proxy für Sedimentzufuhr = hohe TSM
# Fe (µmol/g) → Eisen ist häufig an mineralische Partikel gebunden → Steigt bei terrigenem Eintrag und Resuspension
# Pauthi (authigenic P) → Kann bei intensiver Partikel-Wasser-Interaktion steigen
# >> Distance from coast (km) → TSM nimmt typischerweise offshore ab → Sehr guter Vergleich mit Landsat-Transekten
# Pex (exchangeable P) → Kann bei frischen Sedimenten / Resuspension erhöht sein
# >> Mineral Carbon (MINC)
# >> OI (Oxygen Index)
# >> HI (Hydrogen Index)

# für anlyse chl-a
# >> Porg (organic P) → Direkt mit biologischer Produktion verbunden → Bester P-basierter Proxy für Algen >> gut für Chl-a-Vergleich
# >> Ptotal → Gesamt-Nährstoffverfügbarkeit
# Pex (labiles P) → Sofort bioverfügbar → kann Algenwachstum antreiben
# >> TOC (Total Organic Carbon) → Korreliert oft mit produktiven Gewässern >> korreliert mit Chl-a indirekt
# >> C/N-Verhältnis >> Hinweise auf terrestrische vs. marine OM → Niedrig = eher phytoplanktonisch >> hoch → terrestrisch / Abwasser  

# ---------------------------------------------------------
# ANALYSE Für GABES und P/BI --------------------------------------------
# ---------------------------------------------------------

# wichtigste Faktoren um Industrie zu erkennen:
# **1Phosphor (Gesamt und Fraktionen)** — stärkster Indikator 
# **2️Organische Substanz / TOC**  
# **3️C/P-Verhältnis**  → niedrig = viel anorganischer P (Industrie)
# **4️δ¹³C der organischen Materie**  
# (# **5️Sedimentzusammensetzung**)
var_names <- c("Ptotal", "TOC", "C.P", "δ13C", "Fe", "OI")
ort <- "DJB" # "GBS"
col <- "blue"
allKateb_chem_GBS <- allKateb_chem[(allKateb_chem$Loc == "GBS")| (allKateb_chem$Loc == "DJB"), ] # | (allKateb_chem$Loc == "DJB")

# Plota of Variables per station 
p_vars <- list()
titel <- letters[seq_along(var_names)]
for(i in 1:length(var_names)){
  p_TMS <- ggplot(allKateb_chem_GBS, aes(y = .data[[var_names[i]]], x = Sample, col = Loc))+
    geom_point() +
      geom_smooth(method = "lm", color = "black", se = FALSE) +
    stat_cor(
        aes(label = after_stat( paste0("italic(R)==", round(r, 2), "*','~~italic(R^2)==", round(rr, 2)))),# after_stat(paste0("italic(R)==", round(r, 2)))), # R [-1,1], zeigt richtung + stärke, jenäher 1/-1 desto stärker positiver/negativer Zusammenhang
        parse = TRUE,  label.x.npc = "left", label.y.npc = "bottom") +
  scale_color_manual(values = c("GBS" = "blue", "DJB" = "red")) +
    theme_bw() +
    labs(
      y = var_names[i],
      x =  "Sample Locations",
      title = paste0("(",titel[i], ")")
    ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
   p_TMS 
  p_vars[[i]] <-  p_TMS

}

p_vars_all <- wrap_plots(p_vars) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p_vars_all
ggsave(file.path(output, "imPFacotr_variables.png"), p_vars_all, height = 10, width = 10, scale = 1.2)

# TSM Tassan
p_tassan <- list()
titel <- letters[seq_along(var_names)]
for(i in 1:length(var_names)){
  p_TMS <- ggplot(allKateb_chem_GBS, aes(x = .data[[var_names[i]]], y = TMS_tassan, col = Loc))+
    geom_point() +
      geom_smooth(method = "lm", color = "black", se = FALSE) +
    stat_cor(
        aes(label = after_stat( paste0("italic(R)==", round(r, 2), "*','~~italic(R^2)==", round(rr, 2)))),# after_stat(paste0("italic(R)==", round(r, 2)))), # R [-1,1], zeigt richtung + stärke, jenäher 1/-1 desto stärker positiver/negativer Zusammenhang
        parse = TRUE,  label.x.npc = "left", label.y.npc = "bottom") +

    theme_bw() +
    labs(
      y = "TSM (Tassan)",
      x = var_names[i] ,
      title = paste0("(",titel[i], ")")
    ) 
   p_TMS 
  p_tassan[[i]] <-  p_TMS

}
p_tassan_all <- wrap_plots(p_tassan) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p_tassan_all
ggsave(file.path(output, paste0(ort,"impFactor_TSM_tassan.png")), p_tassan_all, height = 10, width = 10, scale = 1.2)

# TMS Laili
p_laili <- list()
titel <- letters[seq_along(var_names)]
for(i in 1:length(var_names)){
  p_TMS <- ggplot(allKateb_chem_GBS, aes(x = .data[[var_names[i]]], y = TMS_laili, col = Loc))+
    geom_point() +
      geom_smooth(method = "lm", color = "black", se = FALSE) +
    stat_cor(
        aes(label = after_stat( paste0("italic(R)==", round(r, 2), "*','~~italic(R^2)==", round(rr, 2)))),# after_stat(paste0("italic(R)==", round(r, 2)))), # R [-1,1], zeigt richtung + stärke, jenäher 1/-1 desto stärker positiver/negativer Zusammenhang
        parse = TRUE,  label.x.npc = "left", label.y.npc = "bottom") +

    theme_bw() +
    labs(
      y = "TSM (Laili)",
      x = var_names[i] ,
      title = paste0("(",titel[i], ")")
    ) 
   p_TMS 
  p_laili[[i]] <-  p_TMS

}

p_laili_all <- wrap_plots(p_laili) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p_laili_all
ggsave(file.path(output, paste0(ort,"impFactor_TSM_laili.png")), p_laili_all, height = 10, width = 10, scale = 1.2)

# Chl-a
p_chl <- list()
titel <- letters[seq_along(var_names)]
for(i in 1:length(var_names)){
  p_TMS <- ggplot(allKateb_chem_GBS, aes(x = .data[[var_names[i]]], y = log_Chl, col = Loc))+
    geom_point() +
      geom_smooth(method = "lm", color = "black", se = FALSE) +
    stat_cor(
        aes(label = after_stat( paste0("italic(R)==", round(r, 2), "*','~~italic(R^2)==", round(rr, 2)))),# after_stat(paste0("italic(R)==", round(r, 2)))), # R [-1,1], zeigt richtung + stärke, jenäher 1/-1 desto stärker positiver/negativer Zusammenhang
        parse = TRUE,  label.x.npc = "left", label.y.npc = "bottom") +

    theme_bw() +
    labs(
      x = var_names[i],
      y =  "Chl-a",
      title = paste0("(",titel[i], ")")
    ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
   p_TMS 
  p_chl[[i]] <-  p_TMS

}

p_chl_all <- wrap_plots(p_chl) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p_chl_all 
ggsave(file.path(output, paste0(ort,"impFactor_Chl.png")), p_chl_all, height = 10, width = 10, scale = 1.2)

# ---------------------------------------------------------
# ANALYSE TMS --------------------------------------------
# ---------------------------------------------------------
# Pdetr, Distance from coast, , Ptotal <- phosRES
# Mineral Carbon (MINC), OI (Oxygen Index), HI (Hydrogen Index) <- geochem

TMS_phosRES <- data.frame(Sample = phosRES$Sample, DistanceCoast = phosRES$Distance..km., 
                Ptotal = phosRES$Ptotal, Pdetr = phosRES$Pdet, Fe = phosRES$Fe, PFe = phosRES$Pfe)
TMS_geochem <- data.frame(Sample = geoCHEM$Sample, MINC = geoCHEM$MINC, delta13C = geoCHEM$δ13C, OI = geoCHEM$OI, HI = geoCHEM$HI)
TMS_vals <- data.frame(Sample = vals$Sample, TMS_tassan = vals$TMS_tassan, TMS_laili = vals$TMS_laili, TMS_budhiman = vals$TMS_budhiman)


TMS_ANA <- merge(TMS_phosRES, TMS_geochem, by = "Sample")
TMS_ANA <- merge(TMS_ANA, TMS_vals, by = "Sample")
TMS_ANA$Loc <- sapply(strsplit(TMS_ANA$Sample, "-"), `[`, 1)
TMS_ANA$LocNR <- sapply(strsplit(TMS_ANA$Sample, "-"), `[`, 2)
TMS_ANA_clean <- na.omit(TMS_ANA)

## Coraltion table: 
var_names <- c("DistanceCoast", "Ptotal", "Pdetr", "Fe",  "PFe", "MINC", "delta13C", "OI" , "HI")
cor_table_tassan <- func_cor_tab_loc(var_names, TMS_ANA, y_name = "TMS_tassan", loc_name = "Loc")
cor_table_laili <- func_cor_tab_loc(var_names, TMS_ANA, y_name = "TMS_laili", loc_name = "Loc")

# TSM Tassan
p_tassan <- list()
var_names <- c("DistanceCoast", "Ptotal", "Pdetr", "Fe",  "PFe", "MINC", "delta13C", "OI" , "HI")
titel <- letters[seq_along(var_names)]
for(i in 1:length(var_names)){
  p_TMS <- ggplot(TMS_ANA, aes(x = .data[[var_names[i]]], y = TMS_tassan, col = Loc))+
    geom_point() +
      geom_smooth(method = "lm", color = "black", se = FALSE) +
    stat_cor(
        aes(label = after_stat( paste0("italic(R)==", round(r, 2), "*','~~italic(R^2)==", round(rr, 2)))),# after_stat(paste0("italic(R)==", round(r, 2)))), # R [-1,1], zeigt richtung + stärke, jenäher 1/-1 desto stärker positiver/negativer Zusammenhang
        parse = TRUE,  label.x.npc = "left", label.y.npc = "bottom") +

    theme_bw() +
    labs(
      y = "TSM (Tassan)",
      x = var_names[i] ,
      title = paste0("(",titel[i], ")")
    ) 
   p_TMS 
  p_tassan[[i]] <-  p_TMS

}
p_tassan_all <- wrap_plots(p_tassan) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p_tassan_all
ggsave(file.path(output, "TSM_tassan.png"), p_tassan_all, height = 10, width = 10, scale = 1.2)

# TMS Laili
p_laili <- list()
var_names <- c("DistanceCoast", "Ptotal", "Pdetr", "Fe",  "PFe", "MINC", "delta13C", "OI" , "HI")
titel <- letters[seq_along(var_names)]
for(i in 1:length(var_names)){
  p_TMS <- ggplot(TMS_ANA, aes(x = .data[[var_names[i]]], y = TMS_laili, col = Loc))+
    geom_point() +
      geom_smooth(method = "lm", color = "black", se = FALSE) +
    stat_cor(
        aes(label = after_stat( paste0("italic(R)==", round(r, 2), "*','~~italic(R^2)==", round(rr, 2)))),# after_stat(paste0("italic(R)==", round(r, 2)))), # R [-1,1], zeigt richtung + stärke, jenäher 1/-1 desto stärker positiver/negativer Zusammenhang
        parse = TRUE,  label.x.npc = "left", label.y.npc = "bottom") +

    theme_bw() +
    labs(
      y = "TSM (Laili)",
      x = var_names[i] ,
      title = paste0("(",titel[i], ")")
    ) 
   p_TMS 
  p_laili[[i]] <-  p_TMS

}

p_laili_all <- wrap_plots(p_laili) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p_laili_all
ggsave(file.path(output, "TSM_laili.png"), p_laili_all, height = 10, width = 10, scale = 1.2)

# Plota of Variables per station 
p_vars <- list()
var_names <- c("DistanceCoast", "Ptotal", "Pdetr", "Fe",  "PFe", "MINC", "delta13C", "OI" , "HI")
titel <- letters[seq_along(var_names)]
for(i in 1:length(var_names)){
  p_TMS <- ggplot(TMS_ANA, aes(y = .data[[var_names[i]]], x = Sample, col = Loc))+
    geom_point() +
      geom_smooth(method = "lm", color = "black", se = FALSE) +
    stat_cor(
        aes(label = after_stat( paste0("italic(R)==", round(r, 2), "*','~~italic(R^2)==", round(rr, 2)))),# after_stat(paste0("italic(R)==", round(r, 2)))), # R [-1,1], zeigt richtung + stärke, jenäher 1/-1 desto stärker positiver/negativer Zusammenhang
        parse = TRUE,  label.x.npc = "left", label.y.npc = "bottom") +

    theme_bw() +
    labs(
      y = var_names[i],
      x =  "Sample Locations",
      title = paste0("(",titel[i], ")")
    ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
   p_TMS 
  p_vars[[i]] <-  p_TMS

}

p_vars_all <- wrap_plots(p_vars) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p_vars_all
ggsave(file.path(output, "variables_TMS.png"), p_vars_all, height = 10, width = 10, scale = 1.2)


# ---------------------------------------------------------
# ANALYSE Chl-a --------------------------------------------
# ---------------------------------------------------------
# Porg, Distance from coast, , Ptotal, Pex <- phosRES
# C/N, TOC <- geochem

Chl_phosRES <- data.frame(Sample = phosRES$Sample, DistanceCoast = phosRES$Distance..km., 
                Ptotal = phosRES$Ptotal, Pex = phosRES$Pex, Porg = phosRES$Porg)
Chl_geochem <- data.frame(Sample = geoCHEM$Sample, TOC = geoCHEM$TOC, C.N = geoCHEM$C.N)
Chl_vals <- data.frame(Sample = vals$Sample, log_Chl = vals$log_Chl)

Chl_ANA <- merge(Chl_phosRES, Chl_geochem, by = "Sample")
Chl_ANA <- merge(Chl_ANA, Chl_vals, by = "Sample")
Chl_ANA$Loc <- sapply(strsplit(TMS_ANA$Sample, "-"), `[`, 1)
Chl_ANA$LocNR <- sapply(strsplit(TMS_ANA$Sample, "-"), `[`, 2)

## Coraltion table: 
var_names <- c("DistanceCoast",  "Ptotal", "Pex", "Porg", "TOC", "C.N")
cor_table_Chl <- func_cor_tab_loc(var_names, Chl_ANA, y_name = "log_Chl", loc_name = "Loc")

# Chl-a
p_chl <- list()
titel <- letters[seq_along(var_names)]
for(i in 1:length(var_names)){
  p_TMS <- ggplot(Chl_ANA, aes(x = .data[[var_names[i]]], y = log_Chl, col = Loc))+
    geom_point() +
      geom_smooth(method = "lm", color = "black", se = FALSE) +
    stat_cor(
        aes(label = after_stat( paste0("italic(R)==", round(r, 2), "*','~~italic(R^2)==", round(rr, 2)))),# after_stat(paste0("italic(R)==", round(r, 2)))), # R [-1,1], zeigt richtung + stärke, jenäher 1/-1 desto stärker positiver/negativer Zusammenhang
        parse = TRUE,  label.x.npc = "left", label.y.npc = "bottom") +

    theme_bw() +
    labs(
      x = var_names[i],
      y =  "Chl-a",
      title = paste0("(",titel[i], ")")
    ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
   p_TMS 
  p_chl[[i]] <-  p_TMS

}

p_chl_all <- wrap_plots(p_chl) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p_chl_all 
ggsave(file.path(output, "Chl.png"), p_chl_all, height = 10, width = 10, scale = 1.2)

# Plota of Variables per station 
p_vars <- list()
var_names <- c("DistanceCoast",  "Ptotal", "Pex", "Porg", "TOC", "C.N")
titel <- letters[seq_along(var_names)]
for(i in 1:length(var_names)){
  p_TMS <- ggplot(Chl_ANA, aes(y = .data[[var_names[i]]], x = Sample, col = Loc))+
    geom_point() +
      geom_smooth(method = "lm", color = "black", se = FALSE) +
    stat_cor(
        aes(label = after_stat( paste0("italic(R)==", round(r, 2), "*','~~italic(R^2)==", round(rr, 2)))),# after_stat(paste0("italic(R)==", round(r, 2)))), # R [-1,1], zeigt richtung + stärke, jenäher 1/-1 desto stärker positiver/negativer Zusammenhang
        parse = TRUE,  label.x.npc = "left", label.y.npc = "bottom") +

    theme_bw() +
    labs(
      y = var_names[i],
      x =  "Sample Locations",
      title = paste0("(",titel[i], ")")
    ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
   p_TMS 
  p_vars[[i]] <-  p_TMS

}

p_vars_all <- wrap_plots(p_vars) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p_vars_all
ggsave(file.path(output, "variables_Chl.png"), p_vars_all, height = 10, width = 10, scale = 1.2)

# ---------------------------------------------------------
# ANALYSE minerals  --------------------------------------------
# ---------------------------------------------------------

mineral <- merge(mineralComp, TMS_vals, by = "Sample")
mineral <- merge(mineral, Chl_vals, by = "Sample")
mineral$Loc <- sapply(strsplit(mineral$Sample, "-"), `[`, 1)
mineral$LocNR <- sapply(strsplit(mineral$Sample, "-"), `[`, 2)

var_names <- colnames(mineral)[2:(length(colnames(mineral))-6)]

# Chl-a
p_chl <- list()
titel <- letters[seq_along(var_names)]
for(i in 1:length(var_names)){
  p_TMS <- ggplot(mineral, aes(x = .data[[var_names[i]]], y = log_Chl, col = Loc))+
    geom_point() +
      geom_smooth(method = "lm", color = "black", se = FALSE) +
    stat_cor(
        aes(label = after_stat( paste0("italic(R)==", round(r, 2), "*','~~italic(R^2)==", round(rr, 2)))),# after_stat(paste0("italic(R)==", round(r, 2)))), # R [-1,1], zeigt richtung + stärke, jenäher 1/-1 desto stärker positiver/negativer Zusammenhang
        parse = TRUE,  label.x.npc = "left", label.y.npc = "bottom") +

    theme_bw() +
    labs(
      x = var_names[i],
      y =  "Chl-a",
      title = paste0("(",titel[i], ")")
    ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
   p_TMS 
  p_chl[[i]] <-  p_TMS

}

p_chl_all <- wrap_plots(p_chl) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p_chl_all 
ggsave(file.path(output, "minerals_Chl.png"), p_chl_all, height = 10, width = 10, scale = 1.2)

# TSM Tassan
p_tassan <- list()
titel <- letters[seq_along(var_names)]
for(i in 1:length(var_names)){
  p_TMS <- ggplot(mineral, aes(x = .data[[var_names[i]]], y = TMS_tassan, col = Loc))+
    geom_point() +
      geom_smooth(method = "lm", color = "black", se = FALSE) +
    stat_cor(
        aes(label = after_stat( paste0("italic(R)==", round(r, 2), "*','~~italic(R^2)==", round(rr, 2)))),# after_stat(paste0("italic(R)==", round(r, 2)))), # R [-1,1], zeigt richtung + stärke, jenäher 1/-1 desto stärker positiver/negativer Zusammenhang
        parse = TRUE,  label.x.npc = "left", label.y.npc = "bottom") +

    theme_bw() +
    labs(
      y = "TSM (Tassan)",
      x = var_names[i] ,
      title = paste0("(",titel[i], ")")
    ) 
   p_TMS 
  p_tassan[[i]] <-  p_TMS

}
p_tassan_all <- wrap_plots(p_tassan) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p_tassan_all
ggsave(file.path(output, "minerals_TSM_tassan.png"), p_tassan_all, height = 10, width = 10, scale = 1.2)

# ---------------------------------------------------------
# ANALYSE alle variables  --------------------------------------------
# ---------------------------------------------------------

allKateb_TMS <- merge(allKateb, TMS_vals, by = "Sample")
allKateb_CHL <- merge(allKateb, Chl_vals, by = "Sample")

titel <- c(letters, LETTERS)

## Coraltion table: 
var_names <- colnames(allKateb_TMS)[2:(length(colnames(allKateb_TMS))-5)]
all_cor_table_tassan <- func_cor_tab_loc(var_names, allKateb_TMS, y_name = "TMS_tassan", loc_name = "Loc")
all_cor_table_laili <- func_cor_tab_loc(var_names, allKateb_TMS, y_name = "TMS_laili", loc_name = "Loc")

var_names <- colnames(allKateb_CHL)[2:(length(colnames(allKateb_CHL))-3)]
all_cor_table_Chl <- func_cor_tab_loc(var_names, allKateb_CHL, y_name = "log_Chl", loc_name = "Loc")

var_names <- colnames(mineral)[2:(length(colnames(mineral))-6)]
all_cor_table_mineral_tassan <- func_cor_tab_loc(var_names, mineral, y_name = "TMS_tassan", loc_name = "Loc")
all_cor_table_mineral_laili <- func_cor_tab_loc(var_names, mineral, y_name = "TMS_laili", loc_name = "Loc")
all_cor_table_mineral_Chl <- func_cor_tab_loc(var_names, mineral, y_name = "log_Chl", loc_name = "Loc")

# r2 < 0.2
# Phosphatindustrie Gabès: Typische Einträge:
# Phosphogips-Schlämme Phosphatpartikel
# Eisen-PhosphatFeinstsedimente
# gelöstes Phosphat → Ausfällung
# hohe Trübung nahe der Quelle
all_cor_table_tassan[all_cor_table_tassan$R2 > 0.5, ] # starke negative correltation >> Viele deiner Variablen nehmen mit steigendem TSM ab >> Je größer die Distanz, desto geringer das TSM.
# Pdetr und Pauthi zeigt starke neg correlation bei GBs >> jehöher Var desto geringer TMS 
# Daher steigt MINC offshore, während TSM sinkt.
all_cor_table_laili[all_cor_table_laili$R2 > 0.6, ]
all_cor_table_Chl[all_cor_table_Chl$R2 > 0.7, ]

all_cor_table_mineral_tassan[all_cor_table_mineral_tassan$R2 > 0.7, ] # calcit magnesia aragonit in DJB hohe negative correlation >> Je höher die Vars, desto geringer das TSM.
all_cor_table_mineral_laili[all_cor_table_mineral_laili$R2 > 0.6, ]
all_cor_table_mineral_Chl[all_cor_table_mineral_Chl$R2 > 0.7, ]

# TMS Laili
var_names <- colnames(allKateb_TMS)[2:(length(colnames(allKateb_TMS))-5)]
p_laili <- list()
for(i in 1:length(var_names)){
  p_TMS <- ggplot(allKateb_TMS, aes(x = .data[[var_names[i]]], y = TMS_laili, col = Loc))+
    geom_point() +
      geom_smooth(method = "lm", color = "black", se = FALSE) +
    stat_cor(
        aes(label = after_stat( paste0("italic(R)==", round(r, 2), "*','~~italic(R^2)==", round(rr, 2)))),# after_stat(paste0("italic(R)==", round(r, 2)))), # R [-1,1], zeigt richtung + stärke, jenäher 1/-1 desto stärker positiver/negativer Zusammenhang
        parse = TRUE,  label.x.npc = "left", label.y.npc = "bottom") +

    theme_bw() +
    labs(
      y = "TSM (Laili)",
      x = var_names[i] ,
      title = paste0("(",titel[i], ")")
    ) 
   p_TMS 
  p_laili[[i]] <-  p_TMS

}

p_laili_all <- wrap_plots(p_laili) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p_laili_all
ggsave(file.path(output, "all_vars_TSM_laili.png"), p_laili_all, height = 20, width = 20, scale = 1.2)

# TMS tassan
p_tassan <- list()
for(i in 1:length(var_names)){
  p_TMS <- ggplot(allKateb_TMS, aes(x = .data[[var_names[i]]], y = TMS_tassan, col = Loc))+
    geom_point() +
      geom_smooth(method = "lm", color = "black", se = FALSE) +
    stat_cor(
        aes(label = after_stat( paste0("italic(R)==", round(r, 2), "*','~~italic(R^2)==", round(rr, 2)))),# after_stat(paste0("italic(R)==", round(r, 2)))), # R [-1,1], zeigt richtung + stärke, jenäher 1/-1 desto stärker positiver/negativer Zusammenhang
        parse = TRUE,  label.x.npc = "left", label.y.npc = "bottom") +

    theme_bw() +
    labs(
      y = "TSM (Tassan)",
      x = var_names[i] ,
      title = paste0("(",titel[i], ")")
    ) 
   p_TMS 
  p_tassan[[i]] <-  p_TMS

}

p_tassan_all <- wrap_plots(p_tassan) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p_tassan_all
ggsave(file.path(output, "all_vars_TSM_tassan.png"), p_tassan_all, height = 20, width = 20, scale = 1.2)

# Chl-a
var_names <- colnames(allKateb_CHL)[2:(length(colnames(allKateb_CHL))-3)]
p_chl <- list()
for(i in 1:length(var_names)){
  p_TMS <- ggplot(allKateb_CHL, aes(x = .data[[var_names[i]]], y = log_Chl, col = Loc))+
    geom_point() +
      geom_smooth(method = "lm", color = "black", se = FALSE) +
    stat_cor(
        aes(label = after_stat( paste0("italic(R)==", round(r, 2), "*','~~italic(R^2)==", round(rr, 2)))),# after_stat(paste0("italic(R)==", round(r, 2)))), # R [-1,1], zeigt richtung + stärke, jenäher 1/-1 desto stärker positiver/negativer Zusammenhang
        parse = TRUE,  label.x.npc = "left", label.y.npc = "bottom") +

    theme_bw() +
    labs(
      y = "TSM (Tassan)",
      x = var_names[i] ,
      title = paste0("(",titel[i], ")")
    ) 
   p_TMS 
  p_chl[[i]] <-  p_TMS

}

p_chl_all <- wrap_plots(p_chl) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p_chl_all
ggsave(file.path(output, "all_vars_Chl-a.png"), p_chl_all, height = 20, width = 20, scale = 1.2)

# Plota of Variables per station 
p_vars <- list()
var_names <- colnames(allKateb_CHL)[2:(length(colnames(allKateb_CHL))-3)]
for(i in 1:length(var_names)){
  p_TMS <- ggplot(allKateb_CHL, aes(y = .data[[var_names[i]]], x = Sample, col = Loc))+
    geom_point() +
      geom_smooth(method = "lm", color = "black", se = FALSE) +
    stat_cor(
        aes(label = after_stat( paste0("italic(R)==", round(r, 2), "*','~~italic(R^2)==", round(rr, 2)))),# after_stat(paste0("italic(R)==", round(r, 2)))), # R [-1,1], zeigt richtung + stärke, jenäher 1/-1 desto stärker positiver/negativer Zusammenhang
        parse = TRUE,  label.x.npc = "left", label.y.npc = "bottom") +

    theme_bw() +
    labs(
      y = var_names[i],
      x =  "Sample Locations",
      title = paste0("(",titel[i], ")")
    ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
   p_TMS 
  p_vars[[i]] <-  p_TMS

}

p_vars_all <- wrap_plots(p_vars) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p_vars_all
ggsave(file.path(output, "variables_ALL.png"), p_vars_all, height = 20, width = 20, scale = 1.2)


# ---------------------------------------------------------
# PLOT TMS  --------------------------------------------
# ---------------------------------------------------------
gulf_shp_pj <- project(gulf_shp, crs(merge_valid_SUMMER08_EXTENT))
ind_loc_pj <- project(vect(matrix(ind_loc, ncol=2), crs="EPSG:4326"), crs(merge_valid_SUMMER08_EXTENT))
sampleloc_pj <- project(sampleLoc_points, crs(merge_valid_SUMMER08_EXTENT))
sampleloc_extent_pj <- project(sampleloc_extent, crs(merge_valid_SUMMER08_EXTENT))

# merge_valid_SUMMER08_SAMPELPOINTS <- crop(merge_valid_SUMMER08_SAMPELPOINTS, sampleloc_extent_pj)
png(filename = file.path(output,"RAW_SamplePoints.png"), width = 2000, height = 2500, res = 300)
plot(merge_valid_SUMMER08_SAMPELPOINTS)
plot(ind_loc_pj, add = TRUE, col = "red")
plot(gulf_shp_pj, add = TRUE)
dev.off()

# TSM_points <- lapp(merge_valid_SUMMER08_SAMPELPOINTS[[c("Blue","Green","Red")]], fun = function(b, g, r){func_TSM_tassan(b, g, r)})
# writeRaster(TSM_points, file.path(data_dir_valid_masekd,"TSM_POINTS_valid_sampleLoc_SUMMER08.tif"), overwrite=TRUE)
TSM_points <- rast(file.path(data_dir_valid_masekd,"TSM_POINTS_valid_sampleLoc_SUMMER08.tif"))
png(filename = file.path(output,"TMS_SamplePoints.png"), width = 2000, height = 2500, res = 300)
plot(TSM_points, main = "2014-08 TMS [mg/l] at SamplePoints", plg = list(title = "[mg/l]"))
plot(ind_loc_pj, add = TRUE, col = "red")
plot(gulf_shp_pj, add = TRUE)
dev.off()

# Chla_points <- lapp(merge_valid_SUMMER08_SAMPELPOINTS[[c("Red","Green")]], fun = function(r, g){func_log_Chl(r, g)})
# writeRaster(Chla_points, file.path(data_dir_valid_masekd,"Chla_POINTS_valid_sampleLoc_SUMMER08.tif"),overwrite=TRUE)
Chla_points <- rast(file.path(data_dir_valid_masekd,"Chla_POINTS_valid_sampleLoc_SUMMER08.tif"))
png(filename = file.path(output,"Chla_SamplePoints.png"), width = 2000, height = 2500, res = 300)
plot(Chla_points, main = "2014-08 Chl-a [mg/m3] at SamplePoints", plg = list(title = "[mg/m3]"))
plot(ind_loc_pj, add = TRUE, col = "red")
plot(gulf_shp_pj, add = TRUE)
dev.off()

# --------------------------------------
all_cor_table_tassan[all_cor_table_tassan$R2 > 0.5, ] # starke negative correltation >> Viele deiner Variablen nehmen mit steigendem TSM ab >> Je größer die Distanz, desto geringer das TSM.
# Pdetr und Pauthi zeigt starke neg correlation bei GBs >> jehöher Var desto geringer TMS 
# Daher steigt MINC offshore, während TSM sinkt.
all_cor_table_laili[all_cor_table_laili$R2 > 0.5, ]
all_cor_table_Chl[all_cor_table_Chl$R2 > 0.7, ]

all_cor_table_mineral_tassan[all_cor_table_mineral_tassan$R2 > 0.7, ] # calcit magnesia aragonit in DJB hohe negative correlation >> Je höher die Vars, desto geringer das TSM.
all_cor_table_mineral_laili[all_cor_table_mineral_laili$R2 > 0.6, ]
all_cor_table_mineral_Chl[all_cor_table_mineral_Chl$R2 > 0.7, ]
# --------------------------------------------


# merge_valid_SUMMER08_EXTENT_water <- mask(merge_valid_SUMMER08_EXTENT, gulf_shp_pj, inverse = TRUE)
# names(merge_valid_SUMMER08_EXTENT_water) <- names(merge_valid_SUMMER08_EXTENT)
# merge_valid_SUMMER08_EXTENT_water <- crop(merge_valid_SUMMER08_EXTENT_water, sampleloc_extent_pj)
# writeRaster(merge_valid_SUMMER08_EXTENT_water, file.path(data_dir_valid_masekd,"Water_EXTENT_valid_sampleLoc_SUMMER08.tif") , overwrite = TRUE)
merge_valid_SUMMER08_EXTENT_water <- rast(file.path(data_dir_valid_masekd,"Water_EXTENT_valid_sampleLoc_SUMMER08.tif"))
plot(merge_valid_SUMMER08_EXTENT_water)

B <- merge_valid_SUMMER08_EXTENT_water$Blue
G <- merge_valid_SUMMER08_EXTENT_water$Green
R <- merge_valid_SUMMER08_EXTENT_water$Red

# TSM <- lapp(merge_valid_SUMMER08_EXTENT_water[[c("Blue","Green","Red")]], fun = function(b, g, r){func_TSM_tassan(b, g, r)})
# writeRaster(TSM, file.path(data_dir_valid_masekd,"TSM_EXTENT_valid_sampleLoc_SUMMER08.tif"), overwrite=TRUE)
TSM <- rast(file.path(data_dir_valid_masekd,"TSM_EXTENT_valid_sampleLoc_SUMMER08.tif"))
max_val <- round(max(values(TSM), na.rm = TRUE))
min_val <- round(min(values(TSM), na.rm = TRUE))

z <- c(min_val, max_val)
brks <- seq(min_val, max_val, by = 2)
# plot TSM
png(filename = file.path(output,"TMS_EXTENT.png"), width = 2000, height = 2500, res = 300)
plot(TSM, main = "Total Suspended Matter (Tassan), August 2014", zlim =  z,  plg = list(title = "[mg/m3]"), breaks = brks,  col =rev(hcl.colors(99, "YlGnBu")))
plot(ind_loc_pj, add = TRUE, col = "red")
# points(sampleloc_pj, col = "cyan", alpha = 0.3, pch = 1)
dev.off()

# Chl_a <- lapp(merge_valid_SUMMER08_EXTENT_water[[c("Red","Green")]], fun = function(r, g){func_log_Chl(r, g)})
# writeRaster(TSM, file.path(data_dir_valid_masekd,"Chl-a_EXTENT_valid_sampleLoc_SUMMER08.tif"), overwrite=TRUE)
Chl_a <- rast(file.path(data_dir_valid_masekd,"Chl-a_EXTENT_valid_sampleLoc_SUMMER08.tif"))

#####
# VAlidation Chla Reana
data_CHla <- "./NOAAMSL12_chla/chlora/monthly/WW00"

files <- list.files(data_CHla)
chla <- rast(file.path(data_CHla, files))
plot(chla[[1]])

# ---------------------------------------------------------
# crop CHla data  ------------------
# ---------------------------------------------------------
sampleloc_extent_pjchla <- project(sampleloc_extent, crs(chla))
ind_loc_pj_reana <- project(ind_loc_pj, crs(chla))
chla_crop <- crop(chla, sampleloc_extent_pjchla)
plot(chla_crop)
time(chla_crop)

chla_crop_201408 <- chla_crop[[32]]

time(chla_crop_201408)

# plot REANA
vals1 <- values(chla_crop_201408)
vals2 <- values(Chl_a)

max_val <- round(max(c(vals1, vals2), na.rm = TRUE))
min_val <- round(min(c(vals1, vals2), na.rm = TRUE))

z <- c(min_val, max_val)
brks <- seq(min_val, max_val, length.out = 19)

png(filename = file.path(output,"REANA_Chla_EXTENT.png"), width = 2000, height = 2500, res = 300)
plot(chla_crop_201408, main = "Reanalyse Chl-a 2014-08 [mg/m3]", zlim =  z,  plg = list(title = "[mg/m3]"), breaks = brks, col =  rev(hcl.colors(99, "YlGnBu")))
plot(ind_loc_pj_reana, add = TRUE, col = "red")
# points(sampleloc_pj, col = "cyan", alpha = 0.3, pch = 1)
dev.off()

## 

# plot Chl-a
png(filename = file.path(output,"Chla_EXTENT.png"), width = 2000, height = 2500, res = 300)
plot(Chl_a, main = "Chl-a [mg/m3] , August 2014", zlim = z, plg = list(title = "[mg/m3]"),  breaks = brks, col = rev(hcl.colors(99, "YlGnBu")))
plot(ind_loc_pj, add = TRUE, col = "red")
# points(sampleloc_pj, col = "cyan", alpha = 0.3, pch = 1)
dev.off()


summary(values(chla_crop_201408))
summary(values(Chl_a))


######
# nlw 
data_nlw <- "./NOAAMSL12_chla/nlw/monthly/WW00"

# nlw 
files_nlw <- list.files(data_nlw)
nlw <- rast(file.path(data_nlw, files_nlw))
time(nlw)
nwl_201408 <- nlw[[187]]
plot(nwl_201408)
time(nwl_201408)

# crop 
sampleloc_extent_pj <- project(sampleloc_extent, crs(nlw))
ind_loc_pj_nlw <- project(ind_loc_pj, sampleloc_extent)

nlw_crop201408 <- crop(nwl_201408, sampleloc_extent_pj)
time(nlw_crop201408)
plot(nlw_crop201408)

# plot nlw
max_val <- round(max(values(nlw_crop201408), na.rm = TRUE))
min_val <- round(min(values(nlw_crop201408), na.rm = TRUE))

z <- c(min_val, max_val)

png(filename = file.path(output,"REANA_nlw_EXTENT.png"), width = 2000, height = 2500, res = 300)
plot(nlw_crop201408, main = "nlw [W/ m^-2 um^-1 sr^-1] , August 2014", zlim = z, plg = list(title = "[W/ m^-2 um^-1 sr^-1]"),  
   col = rev(hcl.colors(99, "YlGnBu")))
plot(ind_loc_pj_nlw, add = TRUE, col = "red")
# points(sampleloc_pj, col = "cyan", alpha = 0.3, pch = 1)
dev.off()


##############################################################################
# Longituded nur für GBS
allKateb_chem_GBS <- allKateb_chem[(allKateb_chem$Loc == "GBS"), ]

values_long <- vals %>%
  select(Longitude, Green, Red, SWIR1) %>%
  pivot_longer(cols = c(Green, Red, SWIR1),
               names_to = "Band",
               values_to = "Value")

# calculate scaling factor
scale_factorF <- max(values_long$Value, na.rm = TRUE) /
                max(heavymetal$F_mean * 1000 , na.rm = TRUE) # da in mg/l die ander in ym/l
scale_factorP <- max(values_long$Value, na.rm = TRUE) /
                max(heavymetal$P_mean, na.rm = TRUE)
scale_factorCU <- max(values_long$Value, na.rm = TRUE) /
        max(heavymetal$Cu_mean, na.rm = TRUE)           

scale_factorZN <- max(values_long$Value, na.rm = TRUE) /
        max(heavymetal$Zn_mean, na.rm = TRUE)  
heavymetal <- heavymetal %>%
  mutate(F = F_mean * scale_factorF,
        P = P_mean * scale_factorP,
        Cu = Cu_mean*scale_factorCU,
         Zn = Zn_mean*scale_factorZN)

heavy_long <- heavymetal %>%
  select(Longitude, F, P, Cu, Zn) %>%
  pivot_longer(cols = c(F, P, Cu, Zn),
               names_to = "Band",
               values_to = "Value")

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
      SWIR1 = "darkorange"
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
                        name = "Heavy metals [μg/l]")) +
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
  left_join(heavymetal[, c("Longitude", "F_mean", "P_mean", "Zn_mean")],
            by = "Longitude")
# WPI vs TSM
cor_Green  <- cor(df_corr$F_mean, df_corr$Green, use = "complete.obs")
R2_Green   <- cor_Green^2
cor_Red  <- cor(df_corr$F_mean, df_corr$Red, use = "complete.obs")
R2_Red   <- cor_Red^2
cor_Sw  <- cor(df_corr$F_mean, df_corr$SWIR1, use = "complete.obs")
R2_SW  <- cor_Sw^2

cor_Green2  <- cor(df_corr$P_mean, df_corr$Green, use = "complete.obs")
R2_Green2   <- cor_Green2^2
cor_Red2  <- cor(df_corr$P_mean, df_corr$Red, use = "complete.obs")
R2_Red2  <- cor_Red2^2
cor_Sw2  <- cor(df_corr$P_mean, df_corr$SWIR1, use = "complete.obs")
R2_SW2  <- cor_Sw2^2

cor_Green3  <- cor(df_corr$Zn_mean, df_corr$Green, use = "complete.obs")
R2_Green3  <- cor_Green3^2
cor_Red3  <- cor(df_corr$Zn_mean, df_corr$Red, use = "complete.obs")
R2_Red3  <- cor_Red3^2
cor_Sw3  <- cor(df_corr$Zn_mean, df_corr$SWIR1, use = "complete.obs")
R2_SW3  <- cor_Sw3^2

label_text <- paste0(
    "Florine (F):", "</span><br>",
  "<span style='color:darkgreen;'>F Green: R = ", round(cor_Green,2),
  " | R² = ", round(R2_Green,2), "</span><br>",
  "<span style='color:red;'>F Red: R = ", round(cor_Red,2),
  " | R² = ", round(R2_Red,2),"</span><br>",
   "<span style='color:darkorange;'>F SWIR1: R = ", round(cor_Sw,2),
  " | R² = ", round(R2_SW,2),"</span><br>",
"</span><br>",
"Phosphorus (P):", "</span><br>",
  "<span style='color:darkgreen;'>P Green: R = ", round(cor_Green2,2),
  " | R² = ", round(R2_Green2,2), "</span><br>",
  "<span style='color:red;'>P Red: R = ", round(cor_Red2,2),
  " | R² = ", round(R2_Red2,2),"</span><br>",
   "<span style='color:darkorange;'>P SWIR1: R = ", round(cor_Sw2,2),
  " | R² = ", round(R2_SW2,2),"</span><br>",
 "</span><br>",
"Zinc (Zn):", "</span><br>",
  "<span style='color:darkgreen;'>Zn Green: R = ", round(cor_Green3,2),
  " | R² = ", round(R2_Green3,2), "</span><br>",
  "<span style='color:red;'>Zn Red: R = ", round(cor_Red3,2),
  " | R² = ", round(R2_Red3,2),"</span><br>",
   "<span style='color:darkorange;'>Zn SWIR1: R = ", round(cor_Sw3,2),
  " | R² = ", round(R2_SW3,2)

)

p_FP_GR2 <- p_FP_GR  +
  annotate("richtext",     x = max(df_corr$Longitude, na.rm = TRUE)-0.12,
         y = max(values_long$Value, na.rm = TRUE)-950,
           label = label_text,
           hjust = 0,
           size = 5)

p_FP_GR2
ggsave(file.path(output, "geochem_Green_Red.png"), p_FP_GR2, height = 10, width = 10, scale = 1.2)