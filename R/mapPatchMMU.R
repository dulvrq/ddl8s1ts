
#' apply MMU for disturbance map
#'
#' Remove a small patches using a minimum mapping unit.
#'
#' @param dir_map a filename of disturbance map.
#' @param mmu a minimum mapping unit in pixels.
#' @param band a band for applying mmu.
#'
#' @importFrom raster raster getValues setValues clump
#'
#' @return a disturbance map with mmu. A new disturbance map is also generated.
#'
mapPatchMMU <- function(dir_map, mmu = 3, band = 1){
  map <- raster(dir_map, band = band)
  val_map <- val_replace <- getValues(map)
  map <- setValues(map, val_map)
  val_replace[which(val_map > 0)] <- 1
  map_rep <- setValues(map, val_replace)
  # make patch
  patch <- clump(map_rep, gaps = F)

  if(mmu > 1){
    # summary stats of patch
    tap_n <- tapply(patch@data@values, patch@data@values, length)
    id_remove <- as.numeric(names(which(tap_n < mmu)))
    # remove patches with less pixels than MMU
    which_remove <- which(patch@data@values %in% id_remove) # can be done with 'match'
    map@data@values[which_remove] <- NA
  }

  return(map)
}

