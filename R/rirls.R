

#' Calculate RMSE
#'
rmse <- function(x, y) sqrt(mean((x-y)^2, na.rm = T))


#' Robust Iteratively Reweighted Least Squares (RIRLS) to remove outliers.
#'
#' Implement RIRLS to remove outliers from time series.
#' This is optional.
#'
#' @param doy a list of julain date
#' @param val values corresponding to \code{doy}.
#' @param fit NULL or julian dates to be fit. This change return of the function.
#' @param period A maximum julian date to limit fitting (optional). Use NULL to ignore.
#' @return If fit = NULL, coefficients of RIRLS implementation. Otherwise fitted values using coefficients.
#'
#'
rirls <- function(doy, val, fit = NULL, period = NULL){
  if(!is.null(period)){
    which_use <- which(doy <= period)
    doy <- doy[which_use]
    val <- val[which_use]
  }

  Nyr <- max(doy) - min(doy)
  res <- glm(val ~ cos(2*pi*doy) + sin(2*pi*doy) + cos(2*pi*doy/Nyr) + sin(2*pi*doy/Nyr))

  coef <- as.numeric(res$coefficients)

  if(is.null(fit)){
    return(coef)
  }else{
    value <- coef[1] + coef[2]*cos(2*pi*fit) + coef[3]*sin(2*pi*fit) +
      coef[4]*cos(2*pi*fit/Nyr) + coef[5]*sin(2*pi*fit/Nyr)
    return(value)
  }
}


#' Initial implementation of RIRLS
#'
#' Initial implementation of RIRLS before fitting CCDC (for Landsat only)
#'
#' Band values (surface reflectance) must be multiplied by 10000 (range from 0 to 10000).
#'
#' @param doys a list of julain date
#' @param val values corresponding to \code{doy}.
#' @param wch_strt numeric. Designate which is the start date of fitting in \code{doys}.
#' @param wch_end numeric. Designate which is the end date of fitting in \code{doys}.
#' @param index cols for used.
#' @param rirls (function)
#' @return a dataframe masked with NA
#'

init_rirls <- function(doys, val, wch_strt, wch_end, index, rirls){
  ## values are from original article by Zhu and Woodcock (2014) ---
  rirls_b3 <- rirls(doys, val$B3, doys[wch_strt:wch_end], doys[wch_end]) # band 2 of ETM+ in original
  rirls_b6 <- rirls(doys, val$B6, doys[wch_strt:wch_end], doys[wch_end]) # band 5 of ETM+ in original

  obs_omit <- (val$B3[wch_strt:wch_end] - rirls_b3 > 400)|(val$B6[wch_strt:wch_end] - rirls_b6 < -400)

  if(sum(!obs_omit) > 3){ # prevent omission of all values
    val$rirls_flag[wch_strt:wch_end][obs_omit] <- NA
    val[is.na(val$rirls_flag), index] <- NA
  }

  val$rirls_flag <- 0 # initialize for reuse

  return(val)
}

