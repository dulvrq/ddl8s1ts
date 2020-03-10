
#' Extract values from raster
#'
#' Extract values from L8, S1, or DEM at locations x and y.
#'
#' @param ls_l8 a file list of Landsat 8 or Sentinel-1. All files should have the same extent.
#' @param l8_doys a list of julian day of Landsat 8 or Sentinel-1 (numeric).
#' This should be the same order and length as \code{ls_l8}.
#' @param dt_ref a dataframe of reference data for RF. This should contain column 'x', 'y', and 'date'.
#' x and y should indicate locations in the same crs as Landsat or Sentinel. 'date' should indicate the timing of disturbance.
#' If there is no disturbance at that location, use NA.
#' @param col_names names of bands. If L8, this should be B2, B3...B7. If S1, this should be VV and VH.
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom foreach foreach
#' @importFrom raster stack extract
#' @importFrom data.table setnames
#' @importFrom dplyr mutate bind_cols %>% filter if_else
#'
#' @return
#' a dataframe of data values for each location.


extractParallel <- function(ls_l8, l8_doys, dt_ref, col_names = paste0("B", 2:7)){
  ## parallel ---
  n_cl <- detectCores()
  cluster <- makeCluster(n_cl)
  registerDoParallel(cluster)
  ## extract values from each image ---
  dt_l8 <-
    foreach(i = 1:length(ls_l8), .packages = c("raster", "dplyr"), .combine = rbind) %dopar% {
      ### stack ---
      stack_i <- stack(ls_l8[i])
      val_i <- raster::extract(stack_i, dt_ref[,c("x", "y")])
      ### tidy data ---
      dt_use <- dt_ref %>%
        mutate(DOY = l8_doys[i]) %>%
        bind_cols(data.frame(val_i) %>% data.table::setnames(col_names)) %>%
        filter(complete.cases(data.frame(val_i))) %>%
        mutate(manual = if_else(DOY_dis > DOY, 0, 1, missing = 0))

      return(dt_use)
    }
  stopCluster(cluster)

  ## reorder by date and ID ---
  dt_l8 <- dt_l8[order(dt_l8$DOY),]
  dt_l8 <- dt_l8[order(dt_l8$ID),]

  return(dt_l8)
}


#' Extract values from raster (DEM)
#'
#' @param ls_dem a file name of DEM with 2 bands (elevation and slope in this order).
#' @param dt_ref a dataframe of reference data for RF. This should contain column 'x', 'y', and 'date'.
#' x and y should indicate locations in the same crs as Landsat or Sentinel. 'date' should indicate the timing of disturbance.
#' If there is no disturbance at that location, use NA.
#' @param col_names names of bands. Should be elevation and slope.
#'
#' @importFrom raster stack extract
#' @importFrom dplyr bind_cols %>% filter
#'
#' @return
#' a dataframe of data values for each location.

extractTopo <- function(ls_dem, dt_ref, col_names = c("elevation", "slope")){
  ## extract values from each image ---
  ### stack ---
  r_i <- stack(ls_dem)
  val_i <- raster::extract(r_i, dt_ref[,c("x", "y")])
  ### tidy data ---
  dt_use <- dt_ref %>%
    bind_cols(data.frame(val_i) %>% data.table::setnames(col_names)) %>%
    filter(complete.cases(data.frame(val_i)))

  return(dt_use)
}
