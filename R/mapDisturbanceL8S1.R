

#' Map forest disturbance using Landsat 8 and Sentinel-1
#'
#' Detect and map forest disturbance using Landsat 8 and Sentinel-1 time series.
#'
#' This function can be mainly used for three ways.
#' 1. build RF models and then detect & map disturbance.
#' 2. build RF models (no mapping nor detection of disturbance)
#' 3. detect & map disturbance (mandatory with RF models)
#'
#' The use of both Landsat 8 and Sentinel-1 is full implementation of this algorithm.
#' However, you can use either Landsat 8 or Sentinle-1 if you want. The process takes very long time.
#' This use parallel process to save time.
#'
#' @note Only 50\% of sample ID is used for RF modeling.Subfolders are automatically genereted for saving RF models
#' and mapping temporally tiles of a disturbance map.
#'
#' @param ls_l8 a file list of Landsat 8 with 6 bands. All files should have the same extent.
#' @param ls_s1 a file list of Sentinel-1 with 2 bands (VV and VH, in this order). If this is used with Landsat,
#' all files should have the same extent as Landsat (i.e. same origin, crs and resolution)
#' @param l8_doys a list of julian day of Landsat 8 (numeric). This should be the same order and length as \code{ls_l8}.
#' @param s1_doys a list of julian day of Sentinel-1 (numeric). This should be the same order and length as \code{ls_s1}.
#' @param dt_ref a dataframe of reference data for RF. This should contain column 'x', 'y', and 'date'. x and y should
#' indicate locations in the same crs as Landsat or Sentinel. 'date' should indicate the timing of disturbance.
#' If there is no disturbance at that location, use NA.
#' @param ls_dem a file name of DEM with 2 bands (elevation and slope in this order). This is used for RF model.
#' If this is not to be used for model, set NULL.
#' @param dir_save a directory for save the results.
#' @param VI names of spectral index used for Landsat.
#' @param rf_model a list of two RF models (as list) to skip RF model. Set NULL to build RF model.
#' @param startDOY a julian day of the initial timing of disturbance detection (e.g. 2016-02-01 should be 2016.08767)
#' @param endDOY a julian day of the terminate day of disturbance detection.
#' @param mmu numeric. minimum mapping unit (pixel) for final maps. New tif image will be additionally generated.
#' @param only_rf logical. If TRUE, only build RF models. If FALSE, build RF models and map disturbance detection.
#' @param max_cores numeric. Maximum numbers of cores used for parallel processing.
#' @param threshold numeric from 0 to 1. A threshold to detect disturbance from time series disturbance probabilities.
#'
#' @importFrom dplyr rename %>% filter mutate
#' @importFrom raster raster extent crop writeRaster ncell
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom gdalUtils mosaic_rasters
#' @importFrom foreach foreach %dopar%
#'
#' @return a filename of mapped disturbance.
#'
#' @export
#' @examples
#'

mapDisturbanceL8S1 <- function(ls_l8, ls_s1, l8_doys, s1_doys, dt_ref, ls_dem = NULL, dir_save, VI, rf_model = NULL,
                                startDOY, endDOY, mmu = NULL, only_rf = F, max_cores = 20, threshold = 0.5){

  # convert list ---
  if(class(ls_l8) == "character") ls_l8 <- as.list(ls_l8)
  if(class(ls_s1) == "character") ls_s1 <- as.list(ls_s1)
  if(class(ls_dem) == "character") ls_dem <- as.list(ls_dem)

  # buid RF model if not specified ---
  if(is.null(rf_model)){
    cat(catTime(), "Building RF model...")

    # add DOYs and ID column in reference data ---
    if(is.null(dt_ref$date)) stop("Please prepare reference data with 'date' column.")
    doys_each  <- as.numeric(strftime(dt_ref$date, format = "%j"))
    dt_ref <- dt_ref %>%
      rename(date_dis = date) %>%
      mutate(DOY_dis = as.numeric(substring(dt_ref$date, 1, 4)) + doys_each / 365,
             ID = 1:nrow(dt_ref))

    # prepare dirs ---
    dir_rf <- file.path(dir_save, "rf_result")


    # L8 ---
    if(!is.null(ls_l8)){
      cat(catTime(), " extract L8 reference data...")
      ## extract values from each image ---
      col_names_l8 <- paste0("B", 2:7)
      dt_l8 <- extractParallel(ls_l8, l8_doys, dt_ref, col_names = col_names_l8, max_cores)

      # implement CCDC for each location (ID) ---
      cat(catTime(), " fit harmonic models to L8 reference data...")
      ## IDs ---
      id_uniq <- sort(unique(dt_l8$ID))
      ## parallel ---
      n_cl <- min(max_cores, detectCores())
      cluster <- makeCluster(n_cl)
      registerDoParallel(cluster)
      ## implement at each location ---
      dt_l8_c <-
        foreach(i = 1:length(id_uniq), .packages = c("raster", "dplyr"), .combine = rbind,
                .export = c("ccdcTimeSeries", "init_rirls", "rirls", "ccdc", "init_ccdc", "rmse")) %dopar% {

          ### extract data for ID ---
          dt_l8_i <- dt_l8[dt_l8$ID == id_uniq[i],]
          if(dim(dt_l8_i)[1] < 12) return(NULL)
          ### CCDC ---
          dt_l8_ccdc_i <- ccdcTimeSeries(dt_l8_i, VI = VI, startDOY = startDOY, fillNA = T, init_rirls,
                                         rirls, ccdc, init_ccdc, rmse)
          return(dt_l8_ccdc_i)
        }
      stopCluster(cluster)

      # assign elevation and slope ---
      if(!is.null(ls_dem)){
        dt_topo <- extractTopo(ls_dem, dt_ref, col_names = c("elevation", "slope"))
        dt_l8_c$elevation <- dt_topo$elevation[match(dt_l8_c$ID, dt_topo$ID)]
        dt_l8_c$slope     <- dt_topo$slope[match(dt_l8_c$ID, dt_topo$ID)]
      }

      # select cols used for RF model ---
      index <- sub("[[:digit:]]+", "", VI)
      colUse_l8 <- c("manual",
                     paste0(rep(index, each = 4), paste0("_coef", 1:4)),
                     paste0(index, "_RMSE"),
                     index)
      if(!is.null(ls_dem)) colUse_l8 <- c(colUse_l8, "slope", "elevation")

      # sample selection ---
      n_sample_id <- round(length(id_uniq)*0.5)
      sample_train_id <- sort(sample(id_uniq, n_sample_id))

      sample_train_l8 <- match(dt_l8_c$ID, sample_train_id)
      sample_train_l8 <- which(!is.na(sample_train_l8))

      # RF ---
      cat(catTime(), " tune a RF model of L8 reference data...")
      res_rf_l8 <- implementRF(dt_train = dt_l8_c[sample_train_l8,  colUse_l8],
                               dt_test  = dt_l8_c[-sample_train_l8, colUse_l8],
                               dir_rf = dir_rf, name_rdata = "rf_L8.Rdata", do_prallel = T, max_cores)
    }


    # S1 ---
    if(!is.null(ls_s1)){
      cat(catTime(), " extract S1 reference data...")
      ## extract values from each image ---
      col_names_s1 <- c("VV", "VH")
      if(raster::nlayers(stack(ls_s1[1])) == 1) col_names_s1 <- c("VV") # temporally. use VV if only 1 band available

      dt_s1 <- extractParallel(ls_s1, s1_doys, dt_ref, col_names = col_names_s1, max_cores)


      cat(catTime(), " fit harmonic models to S1 reference data...")
      # implement CCDC for each location (ID) ---
      id_uniq <- sort(unique(dt_s1$ID))

      ## extract values from each image ---
      n_cl <- min(max_cores, detectCores())
      cluster <- makeCluster(n_cl)
      registerDoParallel(cluster)
      ## implement at each location ---
      dt_s1_c <-
        foreach(i = 1:length(id_uniq), .packages = c("raster", "dplyr"), .combine = rbind,
                .export = c("ccdcTimeSeriesSAR", "init_rirls", "rirls", "ccdc", "init_ccdc", "rmse")) %dopar% {
          ### extract data for ID ---
          dt_s1_i <- dt_s1[dt_s1$ID == id_uniq[i],]
          if(dim(dt_s1_i)[1] < 12) return(NULL)
          ### CCDC ---
          dt_s1_ccdc_i <- ccdcTimeSeriesSAR(dt_s1_i, VI = col_names_s1, startDOY = startDOY, init_rirls,
                                            rirls, ccdc, init_ccdc, rmse)
          return(dt_s1_ccdc_i)
        }
      stopCluster(cluster)


      # assign elevation and slope ---
      if(!is.null(ls_dem)){
        dt_topo <- extractTopo(ls_dem, dt_ref, col_names = c("elevation", "slope"))
        dt_s1_c$elevation <- dt_topo$elevation[match(dt_s1_c$ID, dt_topo$ID)]
        dt_s1_c$slope     <- dt_topo$slope[match(dt_s1_c$ID, dt_topo$ID)]
      }


      # select cols used for RF model ---
      index <- col_names_s1
      colUse_s1 <- c("manual",
                     paste0(rep(index, each = 4), paste0("_coef", 1:4)),
                     paste0(index, "_RMSE"),
                     index)
      if(!is.null(ls_dem)) colUse_s1 <- c(colUse_s1, "slope", "elevation")

      # sample selection ---
      n_sample_id <- round(length(id_uniq)*0.5)
      sample_train_id <- sort(sample(id_uniq, n_sample_id))

      sample_train_s1 <- match(dt_s1_c$ID, sample_train_id)
      sample_train_s1 <- which(!is.na(sample_train_s1))

      # RF ---
      cat(catTime(), " tune a RF model of S1 reference data...")
      res_rf_s1 <- implementRF(dt_train = dt_s1_c[sample_train_s1,  colUse_s1],
                               dt_test  = dt_s1_c[-sample_train_s1, colUse_s1],
                               dir_rf = dir_rf, name_rdata = "rf_S1.Rdata", do_prallel = T, max_cores)
    }


    if(only_rf) return(NULL) # if only tuning RF model, stop process here


  } else {
    # use already tuned models.
    # should have the same varaibles in the follwing process
    # 'rf_model' should be the list of 2 rf results (Landsat and Sentinel)
    # in case only one model is available, use NULL for missing model

    cat(catTime(), "NOT run RF model building. Use the tuned RF models...")
    res_rf_l8 <- list(rf_model[[1]], NA, NA)
    res_rf_s1 <- list(rf_model[[2]], NA, NA)

  }


  # predict and map the result ---
  ## determin split n of blocks ---
  if(!is.null(ls_l8)){
    r_i <- raster(ls_l8[[1]])
  }else{
    r_i <- raster(ls_s1[[1]])
  }

  n_xys <- ncell(r_i)
  if(n_xys < 3e5){
    nx <- ny <- 1
  }else{
    nx <- (dim(r_i)[2] - 1) %/% 300 + 1
    ny <- (dim(r_i)[1] - 1) %/% 300 + 1
  }
  xys <- gridExtent(r_i, nx, ny)


  ## split into tiles ---
  cat(catTime(), "Predict and map disturbance by ", nrow(xys), " blocks...")
  for(j in 1:dim(xys)[1]){

    e1 <- extent(as.numeric(xys[j,]))

    # L8 ---
    if(!is.null(ls_l8)){
      cat(catTime(), " extract values from L8 data...")
      ## extract values from each image --
      col_names_l8 <- paste0("B", 2:7)
      dt_l8 <- extractParallelCrop(ls_l8, l8_doys, e1, col_names = col_names_l8, max_cores)

      # implement CCDC for each location (ID) ---
      id_uniq <- sort(unique(dt_l8$ID))

      ## extract values from each image ---
      n_cl <- min(max_cores, detectCores())
      cluster <- makeCluster(n_cl)
      registerDoParallel(cluster)
      cat(catTime(), " fit harmonic models to L8 data...")
      ## implement at each location ---
      dt_l8_c <-
        foreach(i = 1:length(id_uniq), .packages = c("raster", "dplyr"), .combine = rbind,
                .export = c("ccdcTimeSeries", "init_rirls", "rirls", "ccdc", "init_ccdc", "rmse")) %dopar% {
          ### extract data for ID ---
          dt_l8_i <- dt_l8[dt_l8$ID == id_uniq[i],]
          if(dim(dt_l8_i)[1] < 12) return(NULL)
          ### CCDC ---
          dt_l8_ccdc_i <- ccdcTimeSeries(dt_l8_i, VI = VI, startDOY = startDOY, fillNA = T, init_rirls,
                                         rirls, ccdc, init_ccdc, rmse)
          return(dt_l8_ccdc_i)
        }
      stopCluster(cluster)


      # assign elevation and slope ---
      if(!is.null(ls_dem)){
        dt_topo <- extractTopoCrop(ls_dem, e1, col_names = c("elevation", "slope"))
        dt_l8_c$elevation <- dt_topo$elevation[match(dt_l8_c$ID, dt_topo$ID)]
        dt_l8_c$slope     <- dt_topo$slope[match(dt_l8_c$ID, dt_topo$ID)]
      }


      # limit the period & select cols ---
      which_extract_l8 <- which((dt_l8_c$DOY >= startDOY)&(dt_l8_c$DOY <= endDOY))

      dt_l8_c <- dt_l8_c[which_extract_l8, ]
      dt_l8_c <- dt_l8_c[complete.cases(dt_l8_c),]

      # select cols used for RF model ---
      index <- sub("[[:digit:]]+", "", VI)
      colUse_l8 <- c(paste0(rep(index, each = 4), paste0("_coef", 1:4)),
                     paste0(index, "_RMSE"),
                     index)
      if(!is.null(ls_dem)) colUse_l8 <- c(colUse_l8, "slope", "elevation")


      # predict based on RF result ---
      cat(catTime(), " predict disturbance probability of L8 data...")
      dt_l8_c$predict <- predict(res_rf_l8[[1]], dt_l8_c[,colUse_l8], type = "prob")[,"1"]
    }


    # S1 ---
    if(!is.null(ls_s1)){
      cat(catTime(), " extract values from S1 data...")
      ## extract values from each image --
      col_names_s1 <- c("VV", "VH")
      if(raster::nlayers(stack(ls_s1[1])) == 1) col_names_s1 <- c("VV") # temporally. use VV if only 1 band available
      dt_s1 <- extractParallelCrop(ls_s1, s1_doys, e1, col_names = col_names_s1, max_cores)
      # implement CCDC for each location (ID) ---
      id_uniq <- sort(unique(dt_s1$ID))

      ## extract values from each image ---
      n_cl <- min(max_cores, detectCores())
      cluster <- makeCluster(n_cl)
      registerDoParallel(cluster)
      cat(catTime(), " fit harmonic models to S1 data...")
      ## implement at each location ---
      dt_s1_c <-
        foreach(i = 1:length(id_uniq), .packages = c("raster", "dplyr"), .combine = rbind,
                .export = c("ccdcTimeSeriesSAR", "init_rirls", "rirls", "ccdc", "init_ccdc", "rmse")) %dopar% {
          ### extract data for ID ---
          dt_s1_i <- dt_s1[dt_s1$ID == id_uniq[i],]
          if(dim(dt_s1_i)[1] < 12) return(NULL)
          ### CCDC ---
          dt_s1_ccdc_i <- ccdcTimeSeriesSAR(dt_s1_i, VI = col_names_s1, startDOY = startDOY, init_rirls,
                                            rirls, ccdc, init_ccdc, rmse)
          return(dt_s1_ccdc_i)
        }
      stopCluster(cluster)


      # assign elevation and slope ---
      if(!is.null(ls_dem)){
        dt_topo <- extractTopoCrop(ls_dem, e1, col_names = c("elevation", "slope"))
        dt_s1_c$elevation <- dt_topo$elevation[match(dt_s1_c$ID, dt_topo$ID)]
        dt_s1_c$slope     <- dt_topo$slope[match(dt_s1_c$ID, dt_topo$ID)]
      }

      # limit the period & select cols ---
      which_extract_s1 <- which((dt_s1_c$DOY >= startDOY)&(dt_s1_c$DOY <= endDOY))

      dt_s1_c <- dt_s1_c[which_extract_s1, ]
      dt_s1_c <- dt_s1_c[complete.cases(dt_s1_c),]

      # select cols used for RF model ---
      index <- col_names_s1
      colUse_s1 <- c(paste0(rep(index, each = 4), paste0("_coef", 1:4)),
                     paste0(index, "_RMSE"),
                     index)
      if(!is.null(ls_dem)) colUse_s1 <- c(colUse_s1, "slope", "elevation")

      # predict based on RF result ---
      cat(catTime(), " predict disturbance probability of S1 data...")
      dt_s1_c$predict <- predict(res_rf_s1[[1]], dt_s1_c[,colUse_s1], type = "prob")[,"1"]
    }


    ## detect disturbance ---
    cat(catTime(), " detect disturbance at the each pixel location...\n")
    if(!exists("dt_l8_c")) {
      dt_l8_c <- NULL
      ids <- sort(unique(dt_s1_c$ID))
    }
    if(!exists("dt_s1_c")) {
      dt_s1_c <- NULL
      ids <- sort(unique(dt_l8_c$ID))
    }
    if(exists("dt_l8_c") & exists("dt_s1_c")) ids <- sort(unique(c(dt_l8_c$ID, dt_s1_c$ID)))

    res_detect <- judgeProbSeriesDetect(dt_l8_c, dt_s1_c, ids = ids,
                                       threshold = threshold, seq_l = c(T,T,T), lag_m = 10, savename = NULL)


    # set as tile image ---
    if(!is.null(ls_l8) & !is.null(ls_s1)) val_mix <- res_detect$estDateMix
    if(!is.null(ls_l8) &  is.null(ls_s1)) val_mix <- res_detect$estDateL8
    if( is.null(ls_l8) & !is.null(ls_s1)) val_mix <- res_detect$estDateS1

    val_mix[val_mix == 0] <- NA
    r_crop <- crop(r_i, e1)
    r_crop <- setValues(r_crop, val_mix)


    # save as temp tile image ---
    cat(catTime(), " save temp", j, "of", nrow(xys), "blocks...")
    dir_temp <- file.path(dir_save, "temp_tiles")
    if(!dir.exists(dir_temp)) dir.create(dir_temp)
    name_temp <- paste0(dir_temp, "/temp_tile_", j, ".tif")
    writeRaster(r_crop, name_temp, options = c("COMPRESS=DEFLATE"), overwrite = T)
  }

  # combine all temp tile files ---
  cat(catTime(), "Combine splited blocks...\n")
  ls_blocks  <- list.files(dir_temp, "temp_tile_.+tif$", full.names = T)
  name_merge <- paste0(dir_save, "/disturbance_detection_L8S1.tif")
  r_merge    <- mosaic_rasters(ls_blocks, name_merge, output_Raster = T, co = c("COMPRESS=DEFLATE"))

  # apply mmu if needed
  if(!is.null(mmu)){
    cat(catTime(), "Apply MMU of", mmu, "pixels...\n")
    r_merge_mmu <- mapPatchMMU(name_merge, mmu = mmu, band = 1)
    name_mmu <- gsub("\\.tif", "_mmu.tif", name_merge)
    writeRaster(r_merge_mmu, name_mmu, options = c("COMPRESS=DEFLATE"), overwrite = T)
  }

  return(name_merge)
}


