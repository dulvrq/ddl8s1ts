
#' Extract values from raster and crop by extent
#'
#' @param ls_l8 a file list of Landsat 8 or Sentinel-1. All files should have the same extent.
#' @param l8_doys a list of julian day of Landsat 8 or Sentinel-1 (numeric).
#' This should be the same order and length as \code{ls_l8}.
#' @param e1 an extent of cropping.
#' @param col_names names of bands. If L8, this should be B2, B3...B7. If S1, this should be VV and VH.
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom foreach foreach
#' @importFrom raster stack extract crop getValues
#' @importFrom data.table setnames
#' @importFrom dplyr mutate bind_cols %>% filter
#'
#' @return
#' a dataframe of data values for each location.

extractParallelCrop <- function(ls_l8, l8_doys, e1, col_names = paste0("B", 2:7)){
  ## parallel ---
  n_cl <- detectCores()
  cluster <- makeCluster(n_cl)
  registerDoParallel(cluster)
  ## extract values from each image ---
  dt_l8 <-
    foreach(i = 1:length(ls_l8), .packages = c("raster", "dplyr"), .combine = rbind) %dopar% {
      ### stack ---
      stack_i <- stack(ls_l8[i])
      stack_crop <- crop(stack_i, e1)
      val_i <- getValues(stack_crop)
      ### tidy data ---
      dt_use <- data.frame(ID = 1:dim(val_i)[1]) %>%
        mutate(DOY = l8_doys[i]) %>%
        bind_cols(data.frame(val_i) %>% data.table::setnames(col_names)) %>%
        filter(complete.cases(data.frame(val_i)))

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
#' @param e1 an extent of cropping.
#' @param col_names names of bands. Should be elevation and slope.
#'
#' @importFrom raster stack crop getValues
#' @importFrom dplyr %>% bind_cols filter
#'
#' @return
#' a dataframe of data values for each location.

extractTopoCrop <- function(ls_dem, e1, col_names = c("elevation", "slope")){
  ## extract values from each image ---
  ### stack ---
  r_i <- stack(ls_dem)
  stack_crop <- crop(r_i, e1)
  val_i <- getValues(stack_crop)
  ### tidy data ---
  dt_use <-  data.frame(ID = 1:dim(val_i)[1]) %>%
    bind_cols(data.frame(val_i) %>% data.table::setnames(col_names)) %>%
    filter(complete.cases(data.frame(val_i)))

  return(dt_use)
}
