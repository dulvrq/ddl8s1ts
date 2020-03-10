

#' initial implementation of CCDC
#'
#' @param doys a list of julain date
#' @param vals values corresponding to \code{doys}.
#' @param index cols for using ccdc.
#' @param col_adj adjustment for cols to insert results into dataframe.
#' @param wch_start numeric. Designate which is the start date of fitting in \code{doys}.
#' @param wch_end numeric. Designate which is the end date of fitting in \code{doys}.
#' @param ccdc (function)
#' @param rmse (function)
#'
#' @return a data frame with the result of initial CCDC implementation.
#'

init_ccdc <- function(doys, vals, index, col_adj, wch_start, wch_end, ccdc, rmse){
  wch <- wch_start:wch_end

  ccdc_k <- ccdc(doys[wch], vals[wch, index], NULL, NULL)
  vals[wch, col_adj] <- rep(ccdc_k, each = length(wch))

  pred_val <- ccdc(doys[wch], vals[wch, index], doys[wch], NULL)
  for(l in 1:length(wch)){
    vals[wch[l], paste0(index, "_RMSE")] <- rmse(vals[wch_start:wch[l], index], pred_val[1:l])
  }

  return(vals)
}

#' Implementation of simple CCDC
#'
#' Implementation of a simple CCDC harmonic regression. See original article by Zhu et al. (2014) for details.
#'
#' @param doy a list of julain date
#' @param val values corresponding to \code{doy}.
#' @param fit NULL or julian dates to be fit. This change return of the function.
#' @param period A maximum julian date to limit fitting (optional). Use NULL to ignore.
#'
#' @return If fit = NULL, coefficients of CCDC implementation. Otherwise fitted values using coefficients.
#'

## implement CCDC model fitting ---
ccdc <- function(doy, val, fit = NULL, period = NULL){
  ### adjust period to fit harmonic model ---
  if(!is.null(period)){
    which_use <- which(doy <= period)
    doy <- doy[which_use]
    val <- val[which_use]
  }
  ### fit model ---
  res <- glm(val ~ cos(2*pi*doy) + sin(2*pi*doy) + doy)
  coef <- as.numeric(res$coefficients)

  if(is.null(fit)){
    return(coef)
  }else{
    value <- coef[1] + coef[2]*cos(2*pi*fit) + coef[3]*sin(2*pi*fit) + coef[4]*fit
    return(value)
  }
}

