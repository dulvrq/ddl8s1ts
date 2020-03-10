
#' Generate grid tiles
#'
#' @param r_use raster object to extract extent.
#' @param nx a number of split in x.
#' @param ny a number of split in y.
#' @param buffer_xy a buffer for overlap between tiles.
#'
#' @importFrom raster extent
#'
#' @return a data frame (xmin, xmax, ymin, ymax) of grids split into tiles.
#'

gridExtent <- function(r_use, nx = 10, ny = 10, buffer_xy = 0){
  ext_r <- extent(r_use)
  x_along <- seq(ext_r[1], ext_r[2], (ext_r[2] - ext_r[1]) / nx)
  y_along <- seq(ext_r[3], ext_r[4], (ext_r[4] - ext_r[3]) / ny)

  xy_along <- expand.grid(head(x_along, nx), head(y_along, ny))
  xy_grid  <- data.frame(xmin = xy_along[, 1] - buffer_xy,
                         xmax = xy_along[, 1] + (ext_r[2] - ext_r[1]) / nx + buffer_xy,
                         ymin = xy_along[, 2] - buffer_xy,
                         ymax = xy_along[, 2] + (ext_r[4] - ext_r[3]) / ny + buffer_xy)
  return(xy_grid)
}
