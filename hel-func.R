# help- functions

# cloud mask 
# Bit-Funktion
get_bit <- function(r, bit){
  bit_mask <- bitwShiftL(1L, bit)
  app(r, function(x) (bitwAnd(x, bit_mask)) > 0)
}

func_timestep_sel <- function(timestep, files_band, satNR_bands, unique_dates, data_dir, time, staNR_exclude = NULL){
    print(time)
    timestep <- as.character(timestep)
    timestep <- unique_dates[grep(paste(timestep, collapse = "|"), unique_dates)]

    if(time == "JULY" | time == "SUMMER"){
            dates_time <- timestep[substr(timestep, 5, 6)  %in%  c("06","07", "08")]

    }else if (time == "WINTER") {
       
         dates_time <- timestep[substr(timestep, 5, 6) %in% c("01", "02", "12")]
         
    } else if (time == "SEPTEMBER") {
       
         dates_time <- timestep[substr(timestep, 5, 6) %in% c("09")] 

    }  else if (time == "HERBST") {
       
         dates_time <- timestep[substr(timestep, 5, 6) %in% c("08","09", "10", "11")] 

    } else if (time == "ALL") {
       
         dates_time <- timestep
    } else {
        stop("Invalid time argument. Use 'SUMMER|JULY' or 'WINTER'.")
    }

    if(length(dates_time) == 0){ # stop function return 0
      print("no scences")
        rasters_timestep <- list()
        rasters_timestep[[1]] <- 0
        return(rasters_timestep)
      
    }

    # filter files for timestep2 - JULY --------------------------
    files_timestep <- files_band[grep(paste(dates_time, collapse = "|"), files_band)]
    parts <- strsplit(files_timestep, "_")
    satNR_files <- unique(sapply(parts, `[`, 1))
    
    if(!is.null(staNR_exclude)){
        idx <- which(sapply(parts, `[`, 1) == staNR_exclude)
        parts <- parts[!idx]
        
        if(length(parts) == 0){
              print(paste(staNR_exclude, "only - STOP"))
             rasters_timestep <- list()
            rasters_timestep[[1]] <- 0
            return(rasters_timestep)
        }

    }
    
    scene_id <- sapply(parts, function(x) paste(x[1:4], collapse = "_"))
    unique_scenes <- unique(scene_id)
    
    # erstelle raster
    rasters_timestep <- lapply(unique_scenes, function(sc) {
        scene_files <- files_timestep[scene_id == sc]
        
        if(!all(file.exists(file.path(data_dir, scene_files)))){
            warning(paste("Some files for scene", sc, "do not exist. Skipping this scene."))
            # delete files which dont exists

            scene_files  <- file.path(data_dir, scene_files)[file.exists(file.path(data_dir, scene_files))]
        }

        rast(file.path(data_dir, scene_files))
        
        
    })

    rasters_timestep <- lapply(rasters_timestep, function(r) {
        parts <- strsplit(names(r), "_")
         sat <- unique(sapply(parts, `[`, 1)) # satellite 
         satNRband_names <- unlist(lapply(parts, function(p){satNR_bands[[sat]][[p[[length(p)]]]]}))

        names(r) <- satNRband_names
        return(r)
    })

    print(unique_scenes)
    print(paste("Number of rasters:", length(rasters_timestep)))

    return(rasters_timestep)
}

func_mask_cloud_EXTENT <- function(raster, data_dir, sampleloc_extent, gabes_hafen = NULL){
    raster_masked <- list()
    for(i in 1:length(raster)){
        print(paste("scene", i, "von", length(raster)))
        # Region corpne/ maskieren 
        # wolken maskierung QA_PIXEL

        secen_raster <- raster[[i]]
        name_parts <- strsplit(sources(secen_raster), "/")
        name_parts <-  strsplit(unique(sapply(name_parts, `[`, length(name_parts[[1]]))), "_")
        name <- paste(unique(as.vector(sapply(name_parts, `[`, 1:7))), collapse = "_")

        sampleloc_extent3_pj <- project(sampleloc_extent, crs(secen_raster))

        # QA Band laden
        qa <- rast(file.path(data_dir, paste0(name,"_QA_PIXEL.TIF")))
        # plot(qa)

        cloud  <- get_bit(qa, 1)
        cloud  <- get_bit(qa, 2)
        cloud  <- get_bit(qa, 3)
        shadow <- get_bit(qa, 4)

        cloud_mask <- cloud | shadow 
        # plot(cloud_mask)

        if(!is.null(gabes_hafen)){
            # check wolken percentage in hafen region
            gabes_hafen_pj <- project(gabes_hafen, crs(secen_raster))
            cloud_mask_extent <- crop(cloud_mask, gabes_hafen_pj)
            cloud_mask_extent <- mask(cloud_mask_extent, gabes_hafen_pj)
            #  plot(cloud_mask_extent)
            total_pixels <- global(!is.na(cloud_mask_extent), "sum")
            cloud_pixels <- global(cloud_mask_extent, "sum", na.rm = TRUE)
            cloud_percent <- (cloud_pixels / total_pixels) 
        } else {
            cloud_percent <- 0
        }

        # if clouds innerhalb hafen gebiet kleiner 10 %%
        if(cloud_percent < 0.1){
            # CLoud und Extent maskieren
            masked <- mask(secen_raster, cloud_mask, maskvalues=1)
            masked <- crop(masked, sampleloc_extent3_pj)
            masked <- mask(masked, sampleloc_extent3_pj)
        } else {
            print("Wolken Anteil zuhoch im Hafen gebiet")
            masked <- 0
        }
     
        raster_masked[[i]] <- masked
        # plot(masked)

    }

    return(raster_masked)
}

library(terra)
library(landsat)
library(lmodel2)

# Funktion: PIF-basierte MA-Normalisierung
normalize_landsat_time_series <- function(raster_list, ref_index = 1, level = 0.99) {
  
  ref <- raster_list[[ref_index]]
  
  # PIF berechnen
  ref_pif <- PIF(ref[["Red"]], ref[["NIR"]], ref[["SWIR2"]], level = level)
  
  # PIF-Maske in SpatRaster umwandeln
  pif_raster <- ref[[1]]
  values(pif_raster) <- as.numeric(ref_pif)  # 1 = stabil, NA = instabil
  
  # Referenzmaskierte Raster
  ref_mask <- ref
  for(b in seq_along(names(ref))) {
    ref_mask[[b]] <- mask(ref[[b]], pif_raster, maskvalues = 0)
  }
  
  norm_list <- vector("list", length(raster_list))
  
  for(i in seq_along(raster_list)) {
    scene <- raster_list[[i]]
    
    # PIF-Maske auf Szene projizieren
    pif_scene <- project(pif_raster, scene)
    scene_masked <- mask(scene, pif_scene, maskvalues = 0)
    
    scene_norm_bands <- list()
    
    for(b in seq_along(names(scene))) {
      ref_band <- ref_mask[[b]]
      tar_band <- scene_masked[[b]]
      
      common_mask <- !is.na(values(ref_band)) & !is.na(values(tar_band))
      ref_vals <- values(ref_band)[common_mask]
      tar_vals <- values(tar_band)[common_mask]
      
      model <- lmodel2(ref_vals ~ tar_vals)
      a <- model$regression.results[2, "Intercept"]
      b_coef <- model$regression.results[2, "Slope"]
      
      band_norm <- tar_band * b_coef + a
      scene_norm_bands[[b]] <- band_norm
    }
    
    scene_norm <- rast(scene_norm_bands)
    names(scene_norm) <- names(scene)
    
    norm_list[[i]] <- scene_norm
    
    rm(scene_masked, scene_norm_bands, pif_scene)
    gc()
  }
  
  return(norm_list)
}

crop_func <- function(raster, ind_loc_sf){

    raster_list <- list()
    for(r in 1:length(raster)){
        data <- raster[[r]]
        
        crop_ind <- lapply(1:nrow(ind_loc_sf), function(i) {
        
            if(i == 3){
                ind_point <- c(10.148889, 34.3375)
                pt <- vect(matrix(ind_point, ncol=2), crs="EPSG:4326")
                ind_point <- project(pt, crs(data))
            } else {
                ind_point <- ind_loc_sf[i]
            }
            
            # create buffer (CRS must be projected in meters!)
            ind_buffer <- buffer(ind_point, width = 8000) 
            
            # crop raster using buffer
            crop(data, ind_buffer)
            })
        raster_list[[r]] <- crop_ind
        }
     
    return(raster_list)
    # strukur = [[raster scene]][[industry location]][[bands]]
}

# mean über ratser
func_mean_raster <- function(raster) {
  
    mean_band <- list()
    common_bands <- Reduce(intersect, lapply(raster, function(x) names(x)))
    print(paste("Common bands:", paste(common_bands, collapse = ", ")))

    # ref <- raster[[2]] # 2013

    # raster_aligned <- lapply(raster, function(r) {
    # r <- resample(r, ref)   # Pixelgröße angleichen
    # r <- crop(r, ext(ref))  # auf Ref-Extent zuschneiden
    # return(r)
    # })

    for(b in common_bands){

        # Extrahiere Band aus allen Raster
        band_list <- lapply(raster, function(x) x[[b]])

        # Stack nur gültige Raster
        band_list <- band_list[!sapply(band_list, is.null)]
        band_stack <- do.call(c, band_list)

        mean_band[[b]] <- app(band_stack, fun = mean, na.rm = TRUE)
    }

    mean_raster <- rast(mean_band)
    names(mean_raster) <- common_bands
   
    return(mean_raster)
}

# mean per industry location
func_mean_loc <- function(raster) {
    mean_perloc <-list()
    for(loc in 1:3){
    raster_Loc <- list()
    for(i in 1:length(raster)){
        print(paste("Scene", i))
        print(raster[[i]][[loc]])
        raster_Loc[[i]] <- raster[[i]][[loc]]

    }
  
    mean_band <- list()
    common_bands <- Reduce(intersect, lapply(raster, function(x) names(x[[loc]])))
    print(paste("Common bands for location", loc, ":", paste(common_bands, collapse = ", ")))
    for(b in common_bands){

        # Extrahiere Band aus allen Raster
        band_list <- lapply(raster_Loc, function(x) x[[b]])

        # Stack nur gültige Raster
        band_list <- band_list[!sapply(band_list, is.null)]
        band_stack <- do.call(c, band_list)

        mean_band[[b]] <- app(band_stack, fun = mean, na.rm = TRUE)
    }

    mean_raster <- rast(mean_band)
    names(mean_raster) <- common_bands
    # Mittelwert über alle Raster
    mean_perloc[[loc]] <- mean_raster
    }
    names(mean_perloc) <- rownames(ind_loc)
    return(mean_perloc)
}

func_std_loc <- function(raster) {
    sd_perloc <-list()
    for(loc in 1:3){
        raster_Loc <- list()
        for(i in 1:length(raster)){
            print(paste("Scene", i))
            print(raster[[i]][[loc]])
            raster_Loc[[i]] <- raster[[i]][[loc]]

        }

        sd_band <- list()
        common_bands <- Reduce(intersect, lapply(raster, function(x) names(x[[loc]])))
        print(paste("Common bands for location", loc, ":", paste(common_bands, collapse = ", ")))
        for(b in common_bands){
            band_list <- lapply(raster_Loc, function(x) x[[b]])
            band_stack <- do.call(c, band_list)
            sd_band[[b]] <- app(band_stack, fun = sd,  na.rm = TRUE)
        }

        sd_raster <- rast(sd_band)
        names(sd_raster) <- common_bands
        # Mittelwert über alle Raster
        sd_perloc[[loc]] <- sd_raster
    }
    names(sd_perloc) <- rownames(ind_loc)
    return(sd_perloc)
}

# coefficent of Variation (CV) = std/mean *100 in %
# 0-5 % Sehr stabil
# 5–20% Moderat stabil
# 20–50% Stark variabel
# >50% Extrem dynamisch
func_CV_loc <- function(sd_perloc, mean_perloc) {
    cv_perloc <- list()
    for(loc in 1:3){
        cv_band <- list()
        for(b in 1:nlyr(mean_perloc[[loc]])){
            mean_band <- mean_perloc[[loc]][[b]]
            sd_band <- sd_perloc[[loc]][[b]]
            cv_band[[b]] <- (sd_band / mean_band) * 100
        }
        cv_raster <- rast(cv_band)
        names(cv_raster) <- names(mean_perloc[[loc]])
        cv_perloc[[loc]] <- cv_raster
    }
    names(cv_perloc) <- rownames(ind_loc)
    return(cv_perloc)
}

# save mean, std, CV
func_save <- function(mean_timestep_JULY, std_timestep_JULY, cv_timestep_JULY, 
    mean_timestep_WINTER, std_timestep_WINTER, cv_timestep_WINTER, 
    ind_loc, timestep_dir){
    timestep <- basename(timestep_dir)
    for(loc in 1:3){
    print(paste("Saving rasters for location:", rownames(ind_loc)[loc]))
        # SUMMER 
        r <- mean_timestep_JULY[[loc]]
        raster_name <- paste0("mean_",timestep,"_JULY_", rownames(ind_loc)[loc], ".tif")
        writeRaster(r, file.path(timestep_dir, raster_name), overwrite=TRUE)

        r <- std_timestep_JULY[[loc]]
        raster_name <- paste0("std_",timestep,"_JULY_", rownames(ind_loc)[loc], ".tif")
        writeRaster(r, file.path(timestep_dir, raster_name), overwrite=TRUE)

        r <- cv_timestep_JULY[[loc]]
        raster_name <- paste0("cv_",timestep,"_JULY_", rownames(ind_loc)[loc], ".tif")
        writeRaster(r, file.path(timestep_dir, raster_name), overwrite=TRUE)


        # WINTER
        r <- mean_timestep_WINTER[[loc]]
        raster_name <- paste0("mean_",timestep,"_WINTER_", rownames(ind_loc)[loc], ".tif")
        writeRaster(r, file.path(timestep_dir, raster_name), overwrite=TRUE)

        r <- std_timestep_WINTER[[loc]]
        raster_name <- paste0("std_",timestep,"_WINTER_", rownames(ind_loc)[loc], ".tif")
        writeRaster(r, file.path(timestep_dir, raster_name), overwrite=TRUE)

        r <- cv_timestep_WINTER[[loc]]
        raster_name <- paste0("cv_",timestep,"_WINTER_", rownames(ind_loc)[loc], ".tif")
        writeRaster(r, file.path(timestep_dir, raster_name), overwrite=TRUE)

    }

    print("Rasters saved successfully.")
}

# mean über alle timesteps 
func_mean_timestep <- function(raster){
    mean_band <- list()
    common_bands <- Reduce(intersect, lapply(raster, function(x) names(x)))
    for(b in common_bands){

        # Extrahiere Band aus allen Raster
        band_list <- lapply(raster, function(x) x[[b]])

        # Stack nur gültige Raster
        band_list <- band_list[!sapply(band_list, is.null)]
        band_stack <- do.call(c, band_list)

        mean_band[[b]] <- app(band_stack, fun = mean, na.rm = TRUE)
    }

    mean_raster <- rast(mean_band)
    names(mean_raster) <- common_bands
    print(paste(common_bands, collapse = ", "))
    return(mean_raster)
}

# std über alle timesteps
func_std_timestep <- function(raster){
    mean_band <- list()
    common_bands <- Reduce(intersect, lapply(raster, function(x) names(x)))
    for(b in common_bands){

        # Extrahiere Band aus allen Raster
        band_list <- lapply(raster, function(x) x[[b]])

        # Stack nur gültige Raster
        band_list <- band_list[!sapply(band_list, is.null)]
        band_stack <- do.call(c, band_list)

        mean_band[[b]] <- app(band_stack, fun = sd, na.rm = TRUE)
    }

    mean_raster <- rast(mean_band)
    names(mean_raster) <- common_bands
    print(paste(common_bands, collapse = ", "))
    return(mean_raster)
}

# CV über alle timesteps
func_cv_timestep <- function(std_raster, mean_raster){
    cv_band <- list()
    for(b in 1:nlyr(mean_raster)){
        mean_band <- mean_raster[[b]]
        sd_band <- std_raster[[b]]
        cv_band[[b]] <- (sd_band / mean_band) * 100
    }
    cv_raster <- rast(cv_band)
    names(cv_raster) <- names(mean_raster)
    return(cv_raster)
}

# diff pro jahr zu baseline
func_baseline_diff <- function(raster_list, baseline){
    diff_list <- lapply(raster_list, function(x) {

        x <- resample(x, baseline)
        diff <- x - baseline
        names(diff) <- names(x)
        return(diff)
    })
    return(diff_list)
}


# water selection mit modified Normalized Difference Water Index (MNDWI)
# MNDWI= (Green-SWIR1)/(Green+SWIR1)
func_MNDWI <- function(GREEN, SWIR1, threshold = 0){
    mndwi <- (GREEN - SWIR1) / (GREEN + SWIR1)
    mask_water <- mndwi > threshold
    return(list(mndwi = mndwi, mask_water = mask_water))
}

func_NDWI <- function(GREEN, NIR){
    ndwi <- (GREEN - NIR) / (GREEN + NIR)
    return(ndwi)
}

# Normalized Difference Turbidity Index
func_NDTI <- function(RED, GREEN){
    ndti <- (RED - GREEN) / (RED + GREEN)
    return(ndti)
}

# Turbidity Index
func_TI <- function(RED, GREEN){
    ti <- (RED) / (GREEN)
    return(ti)
}

# Sediemtn
func_sedi <- function(RED, BLUE){
    sedi <- (RED) / (BLUE)
    return(sedi)
}

# NDSSI=−NDVI 0 | Red stärker → Sediment / Trübung |, < 0 | NIR stärker → Vegetation |
# klares wasser NDSSI nahe 0
# Noramlized disolved supended Sediemtn
func_NDSSI <- function(RED, NIR){
    ndssi <- (RED - NIR) / (RED + NIR)
    return(ndssi)
}

# Suspenden partical matter (SPM), R = reflectance
# SPM = R(665 nm)/ R(490)
func_SPM <- function(R_665, R_490){
    spm <- R_665 / R_490
    return(spm)
}

#####
# Emiyati et al 2017 
# According to Danbara 2014 [15], the water leaving reflectance 𝜌𝑤(𝜆)
# retrieved from atmospheric correction can be converted to remote
# sensing reflectance [15] using formula:
# Rrs = pw(lambda)/ pi
# Tassan (1994) algorithms: Total supended MAtter TSM in mg/l:
# TSM (Blue = 490 nm, Green = 550 nm and Red = 670 nm)
# 𝑻𝑺𝑴 = (𝟏𝟎^(𝟏.𝟔+𝟎.𝟐𝟑 𝒍𝒐𝒈((𝑹(𝟓𝟓𝟎)−𝑹(𝟔𝟕𝟎))(𝑹(𝟓𝟓𝟎)𝑹(𝟒𝟗𝟎))^−𝟎.𝟏))− 𝟏0
func_TSM_tassan <- function(Blue, Green, Red){

  RBlue  <- Blue  / pi
  RGreen <- Green / pi
  RRed   <- Red   / pi

  log_term <- (RGreen - RRed) * (RGreen / RBlue)^(-0.1)

  valid <- log_term > 0 & is.finite(log_term)

  TSM <- rep(NA_real_, length(log_term))
  TSM[valid] <- (10^(1.6 + 0.23 * log10(log_term[valid]))) - 10

  TSM / 10
}

# Sawitri Subiyanto 2017
# Laili algorithm (Band 2 = Blue and Band 4 = Red): 𝑇𝑆𝑆 (𝑚𝑔⁄l ) = 31.42 ∗ (log(𝑅𝑅𝑆2)/log(𝑅𝑅𝑆4)) − 12.719
func_TMS_Laili <- function(Blue, Red){
  RBlue = Blue / pi
  RRed = Red / pi
  log_term <- log10(RBlue)/log10(RRed)
  TMS <- 31.32 * log_term - 12.719
  return(TMS)
}



###
# Gilang Buditama et al 2017
# B2 = Blue B3 = Green , B4 = Red
# Log Chl = (2,41*B4⁄B3)+0,187 [mg/m3]
# Sln = 29,983+165,047(B2)-260,227(B3)+2,609(B4) >> Salinity [ppt]
# TSS = 7,9038*exp(23,942*B4) mg/l

func_log_Chl <- function(Red, Green){
   (2.41* (Red / Green))+0.187
} 

func_TMS_Budhiman <- function(Red) {
  # RRed <- Red / pi
  TMS <- 7.9038 *exp(23.942*Red)
  return(TMS)
}

func_salin <- function(Blue, Green, Red){
  29.983+165.047*(Blue)-260.227*(Green)+2.609*(Red)
}

####
func_cor_tab <- function(var_names, df, dfY){
    cor_table <- data.frame()
    for (v in var_names) {

    x <- df[[v]]
    y <- dfY

    test <- cor.test(x, y, use = "complete.obs", method = "pearson")

    cor_table <- rbind(cor_table, data.frame(   Variable = v, r = unname(test$estimate),R2 = unname(test$estimate)^2,
        p_value = test$p.value, n = sum(complete.cases(x, y))))
    }
    return(cor_table)
}

func_cor_tab_loc <- function(var_names, df, y_name, loc_name = "Loc") {

  cor_table <- data.frame()

  loc_levels <- unique(df[[loc_name]])

  for (loc in loc_levels) {

    df_sub <- df[df[[loc_name]] == loc, ]

    for (v in var_names) {

      x <- df_sub[[v]]
      y <- df_sub[[y_name]]

      if (sum(complete.cases(x, y)) > 2) {

        test <- cor.test(x, y, use = "complete.obs", method = "pearson")

        cor_table <- rbind(
          cor_table,
          data.frame(
            Loc = loc,
            Variable = v,
            r = round(unname(test$estimate), 2),
            R2 = round(unname(test$estimate)^2, 2),
            p_value = test$p.value,
            n = sum(complete.cases(x, y))
          )
        )
      }
    }
  }

  cor_table
}

# interpolate Messungen anhand messpunkten 
library(terra)
library(sf)
library(gstat)
library(sp)

func_interpoltate <- function(heavymetal, VAR, mask, name, data_dir_valid_masekd){

  F_data <- heavymetal[, c(VAR, "Latitude", "Longitude")]

  coordinates(F_data) <- ~Longitude + Latitude
  proj4string(F_data) <- CRS("+proj=longlat +datum=WGS84")
  F_data_utm <- spTransform(
    F_data,
    CRS("+proj=utm +zone=32 +datum=WGS84")
  )

  r <- rast(ext(vect(F_data_utm)), resolution = 5000)
  bb <- bbox(F_data_utm)

  grid <- spsample(F_data_utm, type = "regular", cellsize = 30) 
  gridded(grid) <- TRUE

  # -----------------------------
  # IDW-Modell (gstat)
  # -----------------------------
  gs <- gstat(
    formula =   as.formula(paste(VAR, "~ 1")), # F_mean
    locations = F_data_utm,
    nmax = 20,
    set = list(idp = 2.0)
  )

  idw_out <- predict(gs, grid) 

  idw_out_rast <- rast(idw_out)
  plot(idw_out_rast)

  # Maskieren
  idw_out_rast_mask <- mask(idw_out_rast, mask)
  writeRaster(idw_out_rast_mask, file.path(data_dir_valid_masekd, paste0(name,".tif")) , overwrite = TRUE)

  return(idw_out_rast_mask)
}
library(ggplot2)
library(dplyr)

func_plot_RF_varImp <- function(RF_model) { # chatgpt

  # Importance extrahieren
  imp_raw <- ranger::importance(RF_model)

  # Falls Vektor (ranger) → Data Frame bauen
  if (is.vector(imp_raw)) {
    imp <- data.frame(
      Variable = names(imp_raw),
      Importance = imp_raw
    )
  } else {
    # Falls Matrix (randomForest)
    imp <- as.data.frame(imp_raw)
    imp$Variable <- rownames(imp)

    # Erste Spalte als Importance nehmen
    imp <- imp %>%
      rename(Importance = 1)
  }

  # Sortieren
  imp <- imp %>%
    arrange(desc(Importance))

  # Plot
  p <- ggplot(imp, aes(x = reorder(Variable, Importance), y = Importance)) +
    geom_col(fill = "#2C7FB8") +
    coord_flip() +
    labs(
      title = "Variable Importance (Random Forest)",
      x = "",
      y = "Importance"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.major.y = element_blank()
    )

  return(p)
}

func_saveRF_metric <- function(rf_model, model_name, output_dir){
      mse_oob <- rf_model$prediction.error
    r2_oob  <- rf_model$r.squared
    rmse_oob <- sqrt(mse_oob)

      # Text erstellen
    txt <- paste0(
        "Model: ", model_name, "\n",
        "----------------------------------\n",
        "OOB MSE:  ", round(mse_oob, 6), "\n",
        "OOB RMSE: ", round(rmse_oob, 6), "\n",
        "OOB R2:   ", round(r2_oob, 6), "\n",
        "Number of trees: ", rf_model$num.trees, "\n",
        "Mtry: ", rf_model$mtry, "\n",
        "Min node size: ", rf_model$min.node.size, "\n"
    )

    # Schreiben
    writeLines(txt, con = file.path(output_dir, paste0(model_name, ".txt")))

    
}

func_saveRF_metric_class <- function(rf_model, model_name, output_dir){

  # OOB Error
  oob_error <- rf_model$prediction.error
  oob_accuracy <- 1 - oob_error

  # Confusion Matrix (falls classification = TRUE)
  conf_mat <- rf_model$confusion.matrix
  
  # Text erstellen
  txt <- paste0(
    "Model: ", model_name, "\n",
    "----------------------------------\n",
    "Model type: Classification\n",
    "OOB Error Rate: ", round(oob_error, 4), "\n",
    "OOB Accuracy:   ", round(oob_accuracy, 4), "\n",
    "Number of trees: ", rf_model$num.trees, "\n",
    "Mtry: ", rf_model$mtry, "\n",
    "Min node size: ", rf_model$min.node.size, "\n",
    "Splitrule: ", rf_model$splitrule, "\n",
    "\nConfusion Matrix (OOB):\n",
    capture.output(print(conf_mat)) |> paste(collapse = "\n")
  )

  # Schreiben
  writeLines(txt, con = file.path(output_dir, paste0(model_name, "_metrics.txt")))
}
set.seed(42)
library(terra)
library(ranger)

# Funktion zum pixelweisen Z-Score
zscore <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

anomal <- function(x) {
  (x - mean(x, na.rm = TRUE))
}


func_RF_ranger <- function(dat_interpolate, RSdata_valid, model_name, output_MODLE, output_RF, raster_outRF_pred,
                              brks, unite){
        WPI_interpolate_pj <- project(dat_interpolate,  crs(RSdata_valid))
        WPI_interpolate_pj <- resample(WPI_interpolate_pj, RSdata_valid)
        # RSdata_valid <- resample(RSdata_valid, WPI_interpolate_pj)

        RSdata_valid$TSM_laili <- app( RSdata_valid[[c("Blue","Red")]], fun = function(x) {   func_TMS_Laili(x[,1], x[,2]) })
        RSdata_valid$logChla <- app(RSdata_valid[[c("Red","Green")]],fun = function(x) {  func_log_Chl(x[,1], x[,2]) })
        
        RSdata_valid$NDWI <- app(RSdata_valid[[c("Green","NIR")]],fun = function(x) { func_NDWI(x[,1] , x[,2])})
        RSdata_valid$TI <- app(RSdata_valid[[c("Red","Green")]],fun = function(x) { func_TI(x[,1] , x[,2])})
        RSdata_valid$sediment <- app(RSdata_valid[[c("Red","Blue")]],fun = function(x) { func_sedi(x[,1] , x[,2])})
        RSdata_valid$NDSSI <- app(RSdata_valid[[c("Red","NIR")]],fun = function(x) { func_NDSSI(x[,1] , x[,2])})
       
        RSdata_valid$SWIR1_NIR <- app(RSdata_valid[[c("SWIR1","NIR")]],fun = function(x) {  x[,1] / x[,2] })
        # RSdata_valid$B_G <- app(RSdata_valid[[c("Blue","Green")]],fun = function(x) {  x[,1] / x[,2] })
        # RSdata_valid$R_NIR_SWIR1 <- app(RSdata_valid[[c("Red","NIR", "SWIR1")]],fun = function(x) {  (x[,1] + x[,2]) /  x[,3]})
       
        all_RF_data <- c(RSdata_valid, WPI_interpolate_pj$var1.pred)
        # plot(all_RF_data)

        rf_df_sample <- spatSample(all_RF_data,size = 50000, # sample 
        method = "random", as.df = TRUE, na.rm = TRUE)

         if(unique(sapply(strsplit(model_name, "_"), `[`, 1)) == "WPI"){ # einheit von ym/l zu mg/l
         # ohne NDWI
           rf_df_sample <- rf_df_sample[, !colnames(rf_df_sample) %in% c("NDWI","CoastalAerosol", "SWIR2", "Blue" , "Green" , "Red", "NIR" ,"SWIR1")]
    
        } else {
              rf_df_sample <- rf_df_sample[, !colnames(rf_df_sample) %in% c("NDWI","CoastalAerosol", "SWIR2", "Blue" , "Green" , "Red", "NIR" ,"SWIR1")]
        }
        rf_df_sample <- na.omit(rf_df_sample)
        # # zscore normieren
        # # Alle Spalten außer evtl. Zielvariable normalisieren
        # cols_to_scale <- setdiff(names(rf_df_sample), "var1.pred") # falls var1.pred Ziel ist

        # rf_df_sample[cols_to_scale] <- lapply(rf_df_sample[cols_to_scale], zscore)

        # colnames(rf_df_sample)
        
        rf_model_WPI_RF <- ranger( # schneller da viele pixel
            dependent.variable.name = "var1.pred",
            data = rf_df_sample,
            num.trees = 500,
            mtry = floor(sqrt(ncol(rf_df_sample) - 1)),
            min.node.size = 5,
            importance = "permutation",  
            num.threads = parallel::detectCores() - 1
        )
        print("save model")
        saveRDS(rf_model_WPI_RF, file = file.path(output_MODEL, paste0("rf_model_", model_name, ".rds")))
        func_saveRF_metric(rf_model_WPI_RF, model_name, output_RF)

        # Variable Importance
        # ranger::importance(rf_model_WPI_RF)
        p <- func_plot_RF_varImp(rf_model_WPI_RF)
        ggsave(file.path(output_RF, paste0(model_name,"_pred_VarImp.png")), p, width = 5, height = 5, scale = 1.2)

        # predict 
        pred <- predict(RSdata_valid, rf_model_WPI_RF)
        pred_WPI <- mask(pred, WPI_interpolate_pj$var1.pred)#   
        names(pred_WPI) <- names(pred)
        writeRaster(pred_WPI, file.path(raster_outRF_pred, paste0(model_name,"_pred.tif")), overwrite = TRUE)

        # breaks definded 
        if(unique(sapply(strsplit(model_name, "_"), `[`, 1)) == "P"){ # einheit von ym/l zu mg/l
            pred_WPI <- pred_WPI/1000
            WPI_interpolate_pj[[1]] <- WPI_interpolate_pj[[1]] /1000
        }
        png(file.path(output_RF,  paste0(model_name, "_interpolate_RF_pred.png")), height = 1000, width = 1000)
        plot(pred_WPI, plg = list(title = unite), main = model_name, breaks = brks)
        plot(gulf_shp_pj, add = TRUE)
        points(sampleloc_points_pj, col = "red")
        points(indloc_pj, col = "black", pch = 15, cex = 1)
        dev.off()

        # ERROR RMSE
        pred_WPI_ERROR <- abs(pred_WPI - WPI_interpolate_pj$var1.pred)
        rmse <- sqrt(global((pred_WPI - WPI_interpolate_pj$var1.pred)^2,
                    fun = "mean", # Ein Pixel liegt im Mittel etwa ±RMSE Einheiten vom Referenzwert entfernt.
                    na.rm = TRUE))
        writeRaster(pred_WPI_ERROR, file.path(raster_outRF_pred, paste0("RMSE_pred_", model_name,".tif")), overwrite = TRUE)
   
        # R2
        SS_res <- global((pred_WPI - WPI_interpolate_pj[[1]])^2, "sum", na.rm = TRUE)[1,1]

        # Gesamtvarianz
        mean_obs <- global(WPI_interpolate_pj[[1]], "mean", na.rm = TRUE)[1,1]
        SS_tot <- global((WPI_interpolate_pj[[1]] - mean_obs)^2, "sum", na.rm = TRUE)[1,1]

        R2 <- 1 - (SS_res / SS_tot)

        png(file.path(output_RF,  paste0("RMSE_",model_name,"_interpolate_RF_pred",".png")), height = 1000, width = 1000)
        rng <- max(abs(values(pred_WPI_ERROR)), na.rm = TRUE)
        plot(pred_WPI_ERROR,col = hcl.colors(100, "RdBu"),zlim = c(-rng, rng), 
                main = paste(model_name, "\nRMSE:", round(rmse[1,], 2), "R2:", round(R2, 2)))
        points(sampleloc_points_pj, col = "red")
        points(indloc_pj, col = "black", pch = 15, cex = 1)
        dev.off()
    

    return(list(RF = rf_model_WPI_RF, pred = pred_WPI))
}

func_RF_ranger_class <- function(dat_interpolate, RSdata_valid, model_name, output_MODEL,
                                  output_RF, raster_outRF_pred,
                                  rcl, unite, labs){


    # Projektion & Resampling
    WPI_interpolate_pj <- project(dat_interpolate, crs(RSdata_valid))
    WPI_interpolate_pj <- resample(WPI_interpolate_pj, RSdata_valid)
    # RSdata_valid <- resample(RSdata_valid, WPI_interpolate_pj)

    # Feature Engineering
    RSdata_valid$TSM_laili <- app(RSdata_valid[[c("Blue","Red")]],
                                fun = function(x) func_TMS_Laili(x[,1], x[,2]))

    RSdata_valid$logChla <- app(RSdata_valid[[c("Red","Green")]],
                                fun = function(x) func_log_Chl(x[,1], x[,2]))


    RSdata_valid$NDWI <- app(RSdata_valid[[c("Green","NIR")]],fun = function(x) { func_NDWI(x[,1] , x[,2])})
    RSdata_valid$TI <- app(RSdata_valid[[c("Red","Green")]],fun = function(x) { func_TI(x[,1] , x[,2])})
    RSdata_valid$sediment <- app(RSdata_valid[[c("Red","Blue")]],fun = function(x) { func_sedi(x[,1] , x[,2])})
    RSdata_valid$NDSSI <- app(RSdata_valid[[c("Red","NIR")]],fun = function(x) { func_NDSSI(x[,1] , x[,2])})

    RSdata_valid$SWIR1_NIR <- app(RSdata_valid[[c("SWIR1","NIR")]],fun = function(x) {  x[,1] / x[,2] })
    # RSdata_valid$B_G <- app(RSdata_valid[[c("Blue","Green")]],fun = function(x) {  x[,1] / x[,2] })
    # RSdata_valid$R_NIR_SWIR1 <- app(RSdata_valid[[c("Red","NIR", "SWIR1")]],fun = function(x) {  (x[,1] + x[,2]) /  x[,3]})


    # WPI klassifizieren
    WPI_interpolate_pj$WPI_class <- classify(WPI_interpolate_pj$var1.pred, rcl)
    WPI_interpolate_pj$WPI_class <- as.factor(WPI_interpolate_pj$WPI_class)

    all_RF_data <- c(RSdata_valid, WPI_interpolate_pj$WPI_class)

    # Sampling
    rf_df_sample <- spatSample(all_RF_data, size = 50000,   method = "random", as.df = TRUE, na.rm = TRUE)

    rf_df_sample$WPI_class <- as.factor(rf_df_sample$WPI_class)

    if(unique(sapply(strsplit(model_name, "_"), `[`, 1)) == "WPI"){ # einheit von ym/l zu mg/l
            # ohne NDWI
            rf_df_sample <- rf_df_sample[, !colnames(rf_df_sample) %in% c("NDWI","CoastalAerosol", "SWIR2", "Blue" , "Green" , "Red", "NIR" ,"SWIR1")]

        } else {
                rf_df_sample <- rf_df_sample[, !colnames(rf_df_sample) %in% c("NDWI","CoastalAerosol", "SWIR2", "Blue" , "Green" , "Red", "NIR" ,"SWIR1")]
        }
    rf_df_sample <- na.omit(rf_df_sample)
    # # normieren
    # # Alle Spalten außer evtl. Zielvariable normalisieren
    # cols_to_scale <- setdiff(names(rf_df_sample), "WPI_class") # falls var1.pred Ziel ist

    # rf_df_sample[cols_to_scale] <- lapply(rf_df_sample[cols_to_scale], anomal)

    # RF Classification
    rf_model <- ranger(
        dependent.variable.name = "WPI_class",
        data = rf_df_sample,
        num.trees = 500,
        mtry = floor(sqrt(ncol(rf_df_sample) - 1)),
        importance = "permutation",
        classification = TRUE,
        num.threads = parallel::detectCores() - 1
    )
    print("save model")
    saveRDS(rf_model, file = file.path(output_MODEL, paste0("rf_model_", model_name, ".rds")))
    func_saveRF_metric_class(rf_model, model_name, output_RF)

    # Variable Importance
    # ranger::importance(rf_model_WPI_RF)
    p <- func_plot_RF_varImp(rf_model)
    ggsave(file.path(output_RF, paste0(model_name,"_pred_VarImp.png")), p, width = 5, height = 5, scale = 1.2)

    # prediction 
    print("make prediction")
    pred <- predict(RSdata_valid, rf_model)
    pred_WPI <- mask(pred, WPI_interpolate_pj$var1.pred)# 
    writeRaster(pred_WPI, file.path(raster_outRF_pred, paste0(model_name,"_pred.tif")), overwrite = TRUE)

    # breaks definded 
    if(unique(sapply(strsplit(model_name, "_"), `[`, 1)) == "P"){ # einheit von ym/l zu mg/l
        pred_WPI <- pred_WPI/1000
        WPI_interpolate_pj[[1]] <- WPI_interpolate_pj[[1]] /1000
    }

    levels(pred_WPI) <- data.frame(ID=1:5, label=labs)

    png(file.path(output_RF,  paste0(model_name, "_interpolate_RF_pred.png")), height = 1000, width = 1000)
    plot(pred_WPI, main = model_name)
    plot(gulf_shp_pj, add = TRUE)
    points(sampleloc_points_pj, col = "red")
    points(indloc_pj, col = "black", pch = 15, cex = 1)
    dev.off()

    # ERROR RMSE
    pred_WPI_ERROR <- abs(pred_WPI - WPI_interpolate_pj)
    levels(pred_WPI_ERROR) <- data.frame(ID=1:5, label=labs)

    rmse <- sqrt(global((pred_WPI - WPI_interpolate_pj)^2,
                fun = "mean", # Ein Pixel liegt im Mittel etwa ±RMSE Einheiten vom Referenzwert entfernt.
                na.rm = TRUE))
    writeRaster(pred_WPI_ERROR, file.path(raster_outRF_pred, paste0("RMSE_pred_", model_name,".tif")), overwrite = TRUE)


    # R2
    SS_res <- global((pred_WPI - WPI_interpolate_pj[[1]])^2, "sum", na.rm = TRUE)[1,1]

    # Gesamtvarianz
    mean_obs <- global(WPI_interpolate_pj[[1]], "mean", na.rm = TRUE)[1,1]
    SS_tot <- global((WPI_interpolate_pj[[1]] - mean_obs)^2, "sum", na.rm = TRUE)[1,1]

    R2 <- 1 - (SS_res / SS_tot)

    png(file.path(output_RF,  paste0("RMSE_",model_name,"_interpolate_RF_pred",".png")), height = 1000, width = 1000)
    rng <- max(abs(values(pred_WPI_ERROR)), na.rm = TRUE)
    plot(pred_WPI_ERROR,col = hcl.colors(100, "RdBu"),zlim = c(-rng, rng), 
            main = paste(model_name, "\nRMSE:", round(rmse[1,], 2), "R2:", round(R2, 2)))
    points(sampleloc_points_pj, col = "red")
    points(indloc_pj, col = "black", pch = 15, cex = 1)
    dev.off()


    return(list(RF = rf_model, pred = pred_WPI))
}

func_pred_RF <- function(raster, model, sampleloc_extent3, outpur_dir){

  tot_start <- Sys.time()
  all_pred <- list()

  for(i in 1:length(raster)){
    print(paste("year", i, "von", length(raster)))

    start <- Sys.time()
    rast <- raster[[i]]
     # Feature Engineering
    rast$TSM_laili <- app(rast[[c("Blue","Red")]],
                                fun = function(x) func_TMS_Laili(x[,1], x[,2]))

    rast$logChla <- app(rast[[c("Red","Green")]],
                                fun = function(x) func_log_Chl(x[,1], x[,2]))


    rast$NDWI <- app(rast[[c("Green","NIR")]],fun = function(x) { func_NDWI(x[,1] , x[,2])})
    rast$TI <- app(rast[[c("Red","Green")]],fun = function(x) { func_TI(x[,1] , x[,2])})
    rast$sediment <- app(rast[[c("Red","Blue")]],fun = function(x) { func_sedi(x[,1] , x[,2])})
    rast$NDSSI <- app(rast[[c("Red","NIR")]],fun = function(x) { func_NDSSI(x[,1] , x[,2])})

    rast$SWIR1_NIR <- app(rast[[c("SWIR1","NIR")]],fun = function(x) {  x[,1] / x[,2] })
    # rast$B_G <- app(rast[[c("Blue","Green")]],fun = function(x) {  x[,1] / x[,2] })
    # rast$R_NIR_SWIR1 <- app(rast[[c("Red","NIR", "SWIR1")]],fun = function(x) {  (x[,1] + x[,2]) /  x[,3]})
    

    # prediction
    pred <- terra::predict(rast, model)
    pred_mask <- pred

    # mask 
    #sampleloc_extent3_pj <- project(sampleloc_extent3, crs(pred))
    # pred_mask <- mask(pred, sampleloc_extent3_pj)
    # plot(pred_mask)
    
    stop <- Sys.time()
    time <- stop - start
    message("model:",time)

    # save
    writeRaster(pred_mask, file.path(outpur_dir, paste0("pred_", names(raster)[i],".tif")), overwrite = TRUE)
    all_pred[[i]] <- pred_mask

  }
  tot_stop <- Sys.time()
  time <- tot_stop - tot_start
  message("total model:",time)

  return(all_pred)

}


func_index <- function(rast){
    rast$TSM_laili <- app(rast[[c("Blue","Red")]],
                                fun = function(x) func_TMS_Laili(x[,1], x[,2]))

    rast$logChla <- app(rast[[c("Red","Green")]],
                                fun = function(x) func_log_Chl(x[,1], x[,2]))


    rast$NDWI <- app(rast[[c("Green","NIR")]],fun = function(x) { func_NDWI(x[,1] , x[,2])})
    rast$TI <- app(rast[[c("Red","Green")]],fun = function(x) { func_TI(x[,1] , x[,2])})
    rast$sediment <- app(rast[[c("Red","Blue")]],fun = function(x) { func_sedi(x[,1] , x[,2])})
    rast$NDSSI <- app(rast[[c("Red","NIR")]],fun = function(x) { func_NDSSI(x[,1] , x[,2])})

    rast$SWIR1_NIR <- app(rast[[c("SWIR1","NIR")]],fun = function(x) {  x[,1] / x[,2] })

    return(rast)
}

func_median_mean <- function(raster){
    median_chla_2013_masked <- app(raster, median, na.rm = TRUE)
    mean_chla_2013_masked <- app(raster, mean, na.rm = TRUE)
    sd_chla_2013_masked <- app(raster, sd, na.rm=TRUE)
    max_chla_2013_masked <- app(raster, max, na.rm=TRUE)
    min_chla_2013_masked <- app(raster, min, na.rm=TRUE)

    stacked <- c(median_chla_2013_masked, mean_chla_2013_masked, sd_chla_2013_masked, max_chla_2013_masked, min_chla_2013_masked)
    return(stacked)
}

func_RF_ranger_MODIS <- function(dat_interpolate, RSdata_valid, model_name, output_MODLE, output_RF, raster_outRF_pred,
                              brks, unite){
        WPI_interpolate_pj <- project(dat_interpolate,  crs(RSdata_valid))

        RSdata_valid <- resample(RSdata_valid, WPI_interpolate_pj)
        names(RSdata_valid)
        
        all_RF_data <- c(RSdata_valid, WPI_interpolate_pj$var1.pred)
        plot(all_RF_data)

        if(ncol(all_RF_data >= 50000)){
             rf_df_sample <- spatSample(all_RF_data,size = 50000, # sample 
        method = "random", as.df = TRUE, na.rm = TRUE)
        } else {
            rf_df_sample <- as.dataframe(all_RF_data)
        }
       
        
        rf_model_WPI_RF <- ranger( # schneller da viele pixel
            dependent.variable.name = "var1.pred",
            data = rf_df_sample,
            num.trees = 500,
            mtry = floor(sqrt(ncol(rf_df_sample) - 1)),
            min.node.size = 5,
            importance = "permutation",  
            num.threads = parallel::detectCores() - 1
        )
        
        print("save model")
        saveRDS(rf_model_WPI_RF, file = file.path(output_MODLE, paste0("rf_model_",model_name,".rds")))
        func_saveRF_metric(rf_model_WPI_RF, model_name, output_RF)

        # Variable Importance
        # ranger::importance(rf_model_WPI_RF)
        p <- func_plot_RF_varImp(rf_model_WPI_RF)
        ggsave(file.path(output_RF, paste0(model_name,"_pred_VarImp.png")), p, width = 5, height = 5, scale = 1.2)

        print("prediction")
        pred_WPI <- predict(RSdata_valid, rf_model_WPI_RF)
        writeRaster(pred_WPI, file.path(raster_outRF_pred, paste0(model_name,"_pred.tif")), overwrite = TRUE)

        # breaks definded 
        if(unique(sapply(strsplit(model_name, "_"), `[`, 1)) == "P"){ # einheit von ym/l zu mg/l
            pred_WPI <- pred_WPI/1000
            WPI_interpolate_pj[[1]] <- WPI_interpolate_pj[[1]] /1000
        }
        png(file.path(output_RF,  paste0(model_name, "_interpolate_RF_pred.png")), height = 1000, width = 1000)
        plot(pred_WPI, plg = list(title = unite), main = model_name, breaks = brks)
        plot(gulf_shp_pj, add = TRUE)
        points(sampleloc_points_pj, col = "red")
        points(indloc_pj, col = "black", pch = 15, cex = 1)
        dev.off()

        # ERROR RMSE
        pred_WPI_ERROR <- abs(pred_WPI - WPI_interpolate_pj)
        rmse <- sqrt(global((pred_WPI - WPI_interpolate_pj)^2,
                    fun = "mean", # Ein Pixel liegt im Mittel etwa ±RMSE Einheiten vom Referenzwert entfernt.
                    na.rm = TRUE))
        writeRaster(pred_WPI_ERROR, file.path(raster_outRF_pred, paste0("RMSE_pred_", model_name,".tif")), overwrite = TRUE)
   

        # R2
        SS_res <- global((pred_WPI - WPI_interpolate_pj[[1]])^2, "sum", na.rm = TRUE)[1,1]

        # Gesamtvarianz
        mean_obs <- global(WPI_interpolate_pj[[1]], "mean", na.rm = TRUE)[1,1]
        SS_tot <- global((WPI_interpolate_pj[[1]] - mean_obs)^2, "sum", na.rm = TRUE)[1,1]

        R2 <- 1 - (SS_res / SS_tot)

        png(file.path(output_RF,  paste0("RMSE_",model_name,"_interpolate_RF_pred",".png")), height = 1000, width = 1000)
        rng <- max(abs(values(pred_WPI_ERROR)), na.rm = TRUE)
        plot(pred_WPI_ERROR,col = hcl.colors(100, "RdBu"),zlim = c(-rng, rng), 
                main = paste(model_name, "\nRMSE:", round(rmse[1,], 2), "R2:", round(R2, 2)))
        points(sampleloc_points_pj, col = "red")
        points(indloc_pj, col = "black", pch = 15, cex = 1)
        dev.off()
    

    return(list(RF = rf_model_WPI_RF, pred = pred_WPI))
}


func_RF_ranger_classMODIS <- function(dat_interpolate, RSdata_valid, model_name, output_MODEL,
                                  output_RF, raster_outRF_pred,
                                  rcl, unite, labs){


    # Projektion & Resampling
    WPI_interpolate_pj <- project(dat_interpolate, crs(RSdata_valid))
    RSdata_valid <- resample(RSdata_valid, WPI_interpolate_pj)

    # WPI klassifizieren
    WPI_interpolate_pj$WPI_class <- classify(WPI_interpolate_pj$var1.pred, rcl)
    WPI_interpolate_pj$WPI_class <- as.factor(WPI_interpolate_pj$WPI_class)

    all_RF_data <- c(RSdata_valid, WPI_interpolate_pj$WPI_class)

    valid_mask <- app(RSdata_valid, function(x) all(!is.na(x)))
    predictor_stack <- mask(all_RF_data, valid_mask)

    # Sampling
    if(ncol(all_RF_data) >= 50000){
            rf_df_sample <- spatSample(predictor_stack,size = 50000, # sample 
    method = "random", as.df = TRUE, na.rm = TRUE)
    } else {
        rf_df_sample <- as.data.frame(predictor_stack)
    }
    
    rf_df_sample$WPI_class <- as.factor(rf_df_sample$WPI_class)
    rf_df_sample <- na.omit(rf_df_sample)
    
    # # normieren
    # # Alle Spalten außer evtl. Zielvariable normalisieren
    # cols_to_scale <- setdiff(names(rf_df_sample), "WPI_class") # falls var1.pred Ziel ist

    # rf_df_sample[cols_to_scale] <- lapply(rf_df_sample[cols_to_scale], anomal)

    # RF Classification
    rf_model <- ranger(
        dependent.variable.name = "WPI_class",
        data = rf_df_sample,
        num.trees = 500,
        mtry = floor(sqrt(ncol(rf_df_sample) - 1)),
        importance = "permutation",
        classification = TRUE,
        num.threads = parallel::detectCores() - 1
    )

    saveRDS(rf_model, file = file.path(output_MODEL, paste0("rf_model_", model_name, ".rds")))
    func_saveRF_metric_class(rf_model, model_name, output_RF)

    # Variable Importance
    # ranger::importance(rf_model_WPI_RF)
    p <- func_plot_RF_varImp(rf_model)
    print(ranger::importance(rf_model))
    ggsave(file.path(output_RF, paste0(model_name,"_pred_VarImp.png")), p, width = 5, height = 5, scale = 1.2)

    # predict 
    pred_WPI <- predict(predictor_stack, rf_model)
    pred_WPI <- mask(pred_WPI, WPI_interpolate_pj$var1.pred)
    levels(pred_WPI) <- data.frame(ID=1:5, label=labs)

    writeRaster(pred_WPI, file.path(raster_outRF_pred, paste0(model_name,"_pred.tif")), overwrite = TRUE)

    # breaks definded 
    png(file.path(output_RF,  paste0(model_name, "_interpolate_RF_pred.png")), height = 1000, width = 1000)
    plot(pred_WPI, main = model_name)
    plot(gulf_shp_pj, add = TRUE)
    points(sampleloc_points_pj, col = "red")
    points(indloc_pj, col = "black", pch = 15, cex = 1)
    dev.off()

    # ERROR RMSE
    pred_WPI_ERROR <- abs(pred_WPI - WPI_interpolate_pj)
    levels(pred_WPI_ERROR) <- data.frame(ID=1:5, label=labs)

    rmse <- sqrt(global((pred_WPI - WPI_interpolate_pj)^2,
                fun = "mean", # Ein Pixel liegt im Mittel etwa ±RMSE Einheiten vom Referenzwert entfernt.
                na.rm = TRUE))
    writeRaster(pred_WPI_ERROR, file.path(raster_outRF_pred, paste0("RMSE_pred_", model_name,".tif")), overwrite = TRUE)


    # R2
    SS_res <- global((pred_WPI - WPI_interpolate_pj[[1]])^2, "sum", na.rm = TRUE)[1,1]

    # Gesamtvarianz
    mean_obs <- global(WPI_interpolate_pj[[1]], "mean", na.rm = TRUE)[1,1]
    SS_tot <- global((WPI_interpolate_pj[[1]] - mean_obs)^2, "sum", na.rm = TRUE)[1,1]

    R2 <- 1 - (SS_res / SS_tot)

    png(file.path(output_RF,  paste0("RMSE_",model_name,"_interpolate_RF_pred",".png")), height = 1000, width = 1000)
    rng <- max(abs(values(pred_WPI_ERROR)), na.rm = TRUE)
    plot(pred_WPI_ERROR,col = hcl.colors(100, "RdBu"),zlim = c(-rng, rng), 
            main = paste(model_name, "\nRMSE:", round(rmse[1,], 2), "R2:", round(R2, 2)))
    points(sampleloc_points_pj, col = "red")
    points(indloc_pj, col = "black", pch = 15, cex = 1)
    dev.off()


    return(list(RF = rf_model, pred = pred_WPI))
}
