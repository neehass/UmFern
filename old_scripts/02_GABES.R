# in GABES 
# lon, lat: "Gabes" = c(10.095439, 33.916444)
# production since 1979, expaned 1985

# El Khalil Cherif, 2019 
# Sentinel-3 reflectances at bands 6 (560 nm) and 7 (620 nm). 
# >> band 2 und 3 von landsat 4/5
# >> band 

# packages
library(terra)
library(ggplot2)
library(dplyr)
library(sf)
library(mapview)

getwd()
setwd("./Phosphate")
source('./R-scripts/hel-func.R', chdir = TRUE)

output <- "./R-scripts/GABES-outputs"
dir.create(output)

data_dir <- "./landsat_daten"
data_crop_dir <- "./landsat_daten_crop"
data_masked_dir <- "./landsat_daten_masked"
data_dir_year <- "./landsat_daten_year"

dir.create(data_crop_dir)
dir.create(data_masked_dir)

col_fun <- colorRampPalette(c("white", "red")) # color funcion for CV maps
# ----------------------------------------------------------------
# DAP Produktion ------
# ------------------------------------------------
dap_prod <- read.csv("C:/Users/maiim/Documents/25-26WS/Umweltfernerkundung/Phosphate/Production_DAP_tunisia.csv", sep = ";", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
dap_prod <- as.data.frame(t(dap_prod))
dap_prod$Year <- as.numeric(gsub("X", "", rownames(dap_prod)))
rownames(dap_prod) <- NULL 
colnames(dap_prod) <- c("DAP", "Year")
dap_prod$DAP <- as.numeric(dap_prod$DAP)
View(dap_prod)
mean98_2010 <- mean(dap_prod$DAP[dap_prod$Year >= 1998 & dap_prod$Year <= 2010])
year98_2010 <- dap_prod$Year[dap_prod$Year >= 1998 & dap_prod$Year <= 2010]

mean2011_2017 <- mean(dap_prod$DAP[dap_prod$Year > 2010 & dap_prod$Year <= 2017])
year2011_2017 <- dap_prod$Year[dap_prod$Year > 2010 & dap_prod$Year <= 2017]

dev.off()   # reset device (run until "null device")
png(filename = file.path(output, "DAP_production.png"), width = 800, height = 600)
plot(dap_prod$Year, dap_prod$DAP, type = "l", col = "black",
 main = "DAP Production in Gabes\nin 1000t", ylab = "Production of DAP (1000t)", xlab = "Year",  xaxt = "n",)
grid() 
axis(1, at = dap_prod$Year, labels = dap_prod$Year, las = 2)
lines(year98_2010, rep(mean98_2010, length(year98_2010)), col = "blue")
lines(year2011_2017, rep(mean2011_2017, length(year2011_2017)), col = "red")

usr <- par("usr")  # get plot limits
rect(xleft = 1998, ybottom = usr[3],   xright = 2010, ytop = usr[4],  col = rgb(173/255, 216/255, 230/255, 0.4), border = NA)
abline(v = c(1998, 2010), col="blue", lty=2)

rect(xleft = 2011, ybottom = usr[3],   xright = 2017, ytop = usr[4],  col = rgb(255/255, 0/255, 0/255, 0.2), border = NA)
abline(v = c(2011, 2017), col="red", lty=2)

legend(2012, max(dap_prod$DAP), legend = c("DAP", "Mean 1998-2010", "Mean 2011-2017"), col = c("black", "blue", "red"), 
  lty = 1, bty = "n", bg = "white")
dev.off()
##############################################################################################################
##############################################################################################################
# load croped files and industry locations
files_band <- list.files(data_crop_dir, recursive = TRUE)

# load shape files for tunisia
shp_dir <- "C:/Users/maiim/Documents/25-26WS/Umweltfernerkundung/Phosphate/tunesia-shapefiles"
shp_files <- list.files(shp_dir, pattern = "\\.shp$", full.names = TRUE)
world_shp_list <- lapply(shp_files, vect)
tunisia_shp <- world_shp_list[[207]]
gabes_shp <- tunisia_shp[tunisia_shp$name == "Gabès", ]

###############################################
# SUMMER ------------------------
###################################################
files_SUMMER_Gabes_mean <- files_band[grep("^.*/mean.*_JULY_Gabes\\.tif$", files_band)]
files_SUMMER_Gabes_cv <- files_band[grep("^.*/cv.*_JULY_Gabes\\.tif$", files_band)]

mean_SUMMER_Gabes <- lapply(files_SUMMER_Gabes_mean, function(f) {rast(file.path(data_crop_dir, f))})
cv_SUMMER_Gabes <- lapply(files_SUMMER_Gabes_cv, function(f) {rast(file.path(data_crop_dir, f))})
plot(mean_SUMMER_Gabes[[1]]) # timestep 1

# projection transformieren
gabes_shp <- project(gabes_shp, crs(mean_SUMMER_Gabes[[1]])) # project to raster CRS

nbands <- nlyr(cv_SUMMER_Gabes[[1]]) # timestep 1
par(mfrow=c(2,3))
for(b in 1:nbands){
  plot(cv_SUMMER_Gabes[[1]][[b]], 
       col = col_fun(100), 
       zlim = c(0,100),
       main = paste("CV SUMMER Gabes - Band", b))
  plot(gabes_shp, add = TRUE, border="black", lwd=2)
}
dev.off()

# use inverse gabes-shap als maske
masked_mean_SUMMER_Gabes <- lapply(mean_SUMMER_Gabes, function(x) {
  masked <- mask(x, gabes_shp,  inverse = TRUE)
  names(masked) <- names(x)
  return(masked)
})

masked_cv_SUMMER_Gabes <- lapply(cv_SUMMER_Gabes, function(x) {
  masked <- mask(x, gabes_shp,  inverse = TRUE)
    names(masked) <- names(x)
  return(masked)
})

# save masked rasters
timestep <- c("timestep1", "timestep2", "timestep3", "timestep4")
for(i in 1:length(masked_cv_SUMMER_Gabes)){
  masked_raster <- masked_mean_SUMMER_Gabes[[i]]
  raster_name <- paste0("GABES_mean_masked_", timestep[i], "_JULY.tif")
  writeRaster(masked_raster, file.path(data_masked_dir, raster_name), overwrite=TRUE)

  masked_raster <- masked_cv_SUMMER_Gabes[[i]]
  raster_name <- paste0("GABES_cv_masked_", timestep[i], "_JULY.tif")
  writeRaster(masked_raster, file.path(data_masked_dir, raster_name), overwrite=TRUE)

  png(filename = file.path(output, paste0("CV_SUMMER_Gabes_masked_", timestep[i], ".png")), width = 800, height = 500)
  plot(masked_cv_SUMMER_Gabes[[i]], col = col_fun(100), range = c(0,100))
  dev.off()
}

##############################################################################################################
# ANAlYSE SUMMER --------------------
##############################################################################################################
# Veränderung über alle timesteps: CV
output_ana <- file.path(output, "analyse")
dir.create(output_ana)
dev.off()
mean_SUMMER <- func_mean_timestep(masked_mean_SUMMER_Gabes)
std_SUMMER <- func_std_timestep(masked_mean_SUMMER_Gabes)
cv_SUMMER <- func_cv_timestep(std_SUMMER, mean_SUMMER)

png(filename = file.path(output_ana, paste0("Vergleich_CV_SUMMER_Gabes_masked_", timestep[i], ".png")), width = 800, height = 500)
plot(cv_SUMMER, col = col_fun(100), range = c(0,100))    
dev.off()

# difference between timesteps
diff_meanSummer <- masked_mean_SUMMER_Gabes[[1]] - masked_mean_SUMMER_Gabes[[4]]

col_fun_diff <- colorRampPalette(c("red","white", "blue")) 
max_abs <- max(abs(global(diff_meanSummer, "max", na.rm=TRUE)),
               abs(global(diff_meanSummer, "min", na.rm=TRUE)))

plot(diff_meanSummer, col = col_fun_diff(100), zlim = c(-max_abs, max_abs),
    main = "Difference of mean SUMMER Gabes\nbetween timestep 1 and 4")

##############################################################################################################
# Analyse baseline 1984 -------
# >> veränderung pro jahr (y <- 43 ## 43 Jahre time periode)
##############################################################################################################
files_year <- list.files(data_dir_year, recursive = TRUE)
files_year_SUMMER_mean  <- files_year[grep("^.*/mean_.*_JULY_Gabes\\.tif$", files_year)] # files_year[grep("^.*/mean_.*_WINTER_Gabes\\.tif$", files_year)] 

years_gabes_SUMMER <- sapply(strsplit(files_year_SUMMER_mean, "/"), `[`, 1)
years_gabes_SUMMER <- as.numeric(years_gabes_SUMMER)
length(years_gabes_SUMMER)

mean_year_SUMMER_Gabes <- lapply(files_year_SUMMER_mean, function(f) {rast(file.path(data_dir_year, f))})
# projection transformieren
gabes_shp <- project(gabes_shp, crs(mean_year_SUMMER_Gabes[[1]])) # project to raster CRS

# crop to gabes
masked_mean_year_SUMMER_Gabes <- lapply(mean_year_SUMMER_Gabes, function(x) {
  masked <- mask(x, gabes_shp,  inverse = TRUE)
  names(masked) <- names(x)
  return(masked)
})

# baseline difference
baselien_year_1984 <- masked_mean_year_SUMMER_Gabes[[1]] # 1984
plot(baselien_year_1984) # main = "Baseline year 1984 - SUMMER Gabes"

baseline_diff <- func_baseline_diff(masked_mean_year_SUMMER_Gabes, baselien_year_1984)


band_nr <- 3
all_max_abs <- Inf
for(y in 1:length(years_gabes_SUMMER)){
  r <- baseline_diff[[y]][[2]]
  max_abs <- max(abs(global(r, "max", na.rm=TRUE)),
                abs(global(r, "min", na.rm=TRUE)))
    all_max_abs <- max(all_max_abs, max_abs)
 
}
col_fun_diff <- colorRampPalette(c("red","white", "blue")) 

par(mfrow=c(7,5))
for(y in 1:length(years_gabes_SUMMER)){
    # max_abs <- max(abs(global(baseline_diff[[y]], "max", na.rm=TRUE)),
    #             abs(global(baseline_diff[[y]], "min", na.rm=TRUE)))
    band <- names(baseline_diff[[y]]$Red)
  plot(baseline_diff[[y]]$Red, col = col_fun_diff(100), zlim = c(-all_max_abs, all_max_abs),
       main = paste(band, "\n", years_gabes_SUMMER[y], "- 1984 (baseline)"))
      plot(gabes_shp, add = TRUE, border="black", lwd=1)
}
dev.off()


###########################
# Emiyati et al 2017 
# According to Danbara 2014 [15], the water leaving reflectance 𝜌𝑤(𝜆)
# retrieved from atmospheric correction can be converted to remote
# sensing reflectance [15] using formula:
# Rrs = pw(lambda)/ pi
# Tassan (1994) algorithms: Total supended MAtter TSM in mg/l:
# TSM (Blue = 490 nm, Green = 550 nm and Red = 670 nm)
# 𝑻𝑺𝑴 = (𝟏𝟎^(𝟏.𝟔+𝟎.𝟐𝟑 𝒍𝒐𝒈((𝑹(𝟓𝟓𝟎)−𝑹(𝟔𝟕𝟎))(𝑹(𝟓𝟓𝟎)𝑹(𝟒𝟗𝟎))^−𝟎.𝟏))− 𝟏0
TSM_tassan <- function(Blue, Green, Red){
  RBlue <- Blue / pi
  RGreen <- Green/ pi
  RRed <- Red / pi
  log_term <- (RGreen -RRed)*(RGreen/RBlue)^(-0.1)

  TSM <- (10^(1.6+0.23 * log10(log_term)))-10
  TSM <- TSM/10
  return(TSM) 
}


#test 
test_TSM <- lapply(masked_mean_year_SUMMER_Gabes, function(x) {TSM_tassan(x$Blue, x$Green, x$Red)})
y <- 26
plot(test_TSM[[y]], main = years_gabes_SUMMER[y])

# Sawitri Subiyanto 2017
# Laili algorithm (Band 2 = Blue and Band 4 = Red): 𝑇𝑆𝑆 (𝑚𝑔⁄l ) = 31.42 ∗ (log(𝑅𝑅𝑆2)/log(𝑅𝑅𝑆4)) − 12.719
TMS_Laili <- function(Blue, Red){
  RBlue = Blue / pi
  RRed = Red / pi
  log_term <- log10(RBlue)/log10(RRed)
  TMS <- 31.32 * log_term - 12.719
  return(TMS)
}

test_TSM_laili <- lapply(masked_mean_year_SUMMER_Gabes, function(x) {TMS_Laili(x$Blue, x$Red)})
y <- 31
plot(test_TSM_laili[[y]], main = years_gabes_SUMMER[y])


###
# Gilang Buditama et al 2017
# B2 = Blue B3 = Green , B4 = Red
# Log Chl = (2,41*B4⁄B3)+0,187 [mg/m3]
# Sln = 29,983+165,047(B2)-260,227(B3)+2,609(B4) >> Salinity [ppt]
# TSS = 7,9038*exp(23,942*B4) mg/l

log_Chl <- function(Red, Green){
   (2.41* (Red / Green))+0.187
} 
test_Chl <- lapply(masked_mean_year_SUMMER_Gabes, function(x) {log_Chl(x$Red, x$Green)})
y <- 27
plot(test_Chl[[y]], main = years_gabes_SUMMER[y])

TMS_Budhiman <- function(Red) {
  # RRed <- Red / pi
  TMS <- 7.9038 *exp(23.942*Red)
  return(TMS)
}
test_TSM_bud <- lapply(masked_mean_year_SUMMER_Gabes, function(x) {TMS_Budhiman(x$Red)})
y <- 1
plot(test_TSM_bud[[y]], main = years_gabes_SUMMER[y])

salin <- function(Blue, Green, Red){
  29.983+165.047*(Blue)-260.227*(Green)+2.609*(Red)
}
test_salin <- lapply(masked_mean_year_SUMMER_Gabes, function(x) {salin(x$Blue, x$Green, x$Red)})
y <- 31
plot(test_salin[[y]], main = years_gabes_SUMMER[y])


# quelle CHatgpt
# Nechad
A <- 289.29
C <- 0.1686

TSM <- (A * red) / (1 - (red / C))

TSM <- 327.7 * red - 2.7