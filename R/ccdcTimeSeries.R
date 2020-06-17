

#' Implementation of CCDC for time series data
#'
#' Implementation of CCDC fitting and variable deriviation for time series Landsat data.
#'
#' Main function for CCDC implementation. For Sentinel-1 data, use \code{\link{ccdcTimeSeriesSAR}}.
#'
#' @param values a dataframe with time series Landsat data.
#' @param VI character. A spectral index to fit harmonic model.
#' @param startDOY numeric. A julian day of starting CCD change detection.
#' @param fillNA logical. If TRUE, fill NA values in time seires with original values.
#' @param init_rirls (function)
#' @param rirls (function)
#' @param ccdc (function)
#' @param init_ccdc (function)
#' @param rmse (function)
#'
#' @return a dataframe with CCDC results and potential disturbance flags.
#'
#'

ccdcTimeSeries <- function(values, VI = "NBR", startDOY = 2016, fillNA = F, init_rirls,
                           rirls, ccdc, init_ccdc, rmse){


  # coefficient ---
  coef_TCB5 <- c(0.2043, 0.4158, 0.5524, 0.5741, 0.3124, 0.2303)
  coef_TCG5 <- c(-0.1603, -0.2819, -0.4934, 0.7940, -0.0002, -0.1446)  # TC coef for landsat5
  coef_TCW5 <- c(0.0315, 0.2021, 0.3102, 0.1594, -0.6806, -0.6109)

  coef_TCB7 <- c(0.3561, 0.3972, 0.3904, 0.6966, 0.2286, 0.1596)
  coef_TCG7 <- c(-0.3344, -0.3544, -0.4556, 0.6966, -0.0242, -0.2630)  # TC coef for landsat7
  coef_TCW7 <- c(0.2626, 0.2141, 0.0926, 0.0656, -0.7629, -0.5388)

  coef_TCB8 <- c(0.3029, 0.2786, 0.4733, 0.5599, 0.5080, 0.1872)
  coef_TCG8 <- c(-0.2941, -0.2430, -0.5424, 0.7276, 0.0713, -0.1608)  # TC coef for landsat8
  coef_TCW8 <- c(0.1511, 0.1973, 0.3283, 0.3407, -0.7117, -0.4559)

  # add columns ---
  Values <- values # copy
  names.val <- c(sub("[[:digit:]]+", "", VI),
                 paste0(rep(sub("[[:digit:]]+", "", VI), each = 4), paste0("_coef", 1:4)),
                 paste0(sub("[[:digit:]]+", "", VI), "_RMSE"),
                 "rirls_flag",
                 paste0("disturbance_flag_", sub("[[:digit:]]+", "", VI))
  )

  Values[,names.val] <- 0
  doys <- Values[,"DOY"]

  # calculate VI values ---
  l8 <- Values[,c("B2", "B3", "B4", "B5", "B6", "B7")]
  for(k in 1:length(VI)){
    index <- VI[k]

    switch(index,
           "TCB5" = Res <- colSums(t(l8)*coef_TCB5) ,
           "TCG5" = Res <- colSums(t(l8)*coef_TCG5) ,
           "TCW5" = Res <- colSums(t(l8)*coef_TCW5) ,
           "TCA5" = Res <- atan(colSums(t(l8)*coef_TCG5) / colSums(t(l8)*coef_TCB5)) ,

           "TCB7" = Res <- colSums(t(l8)*coef_TCB7) ,
           "TCG7" = Res <- colSums(t(l8)*coef_TCG7) ,
           "TCW7" = Res <- colSums(t(l8)*coef_TCW7) ,
           "TCA7" = Res <- atan(colSums(t(l8)*coef_TCG7) / colSums(t(l8)*coef_TCB7)) ,

           "TCB8" = Res <- colSums(t(l8)*coef_TCB8) ,
           "TCG8" = Res <- colSums(t(l8)*coef_TCG8) ,
           "TCW8" = Res <- colSums(t(l8)*coef_TCW8) ,
           "TCA8" = Res <- atan(colSums(t(l8)*coef_TCG8) / colSums(t(l8)*coef_TCB8)) ,

           "NBR"  = Res <- (l8[,4] - l8[,6]) / (l8[,4] + l8[,6]),
           "NDVI" = Res <- (l8[,4] - l8[,3]) / (l8[,4] + l8[,3]),
           "NDMI" = Res <- (l8[,4] - l8[,5]) / (l8[,4] + l8[,5]),
           stop("Set Valid index")
    )
    Values[,sub("[[:digit:]]+", "", index)] <- Res
  }


  # RIRLS implementation ---
  # at least 12 observation is needed to start detect change
  # here observation before 2016/01/01 was used for model initialization.
  strt <- min(which(doys >= startDOY))
  Values <- init_rirls(doys, Values, 1, strt-3, sub("[[:digit:]]+", "", VI), rirls)

  # estimate coef for initial sequence ---
  for(k in 1:length(VI)){
    index <- sub("[[:digit:]]+", "", VI[k])
    Values <- init_ccdc(doys, Values, index,
                        col_adj = dim(values)[2]+length(VI)+seq(4*(k-1)+1,4*(k-1)+4,1),
                        wch_start = 1,
                        wch_end = strt-3, ccdc, rmse)
  }


  # start CCDC implementation ---
  for(i in strt:length(doys)){

    # _ccdc model ---
    for(k in 1:length(VI)){
      index <- sub("[[:digit:]]+", "", VI[k])
      col_coef <- dim(values)[2]+length(VI)+seq(4*(k-1)+1,4*(k-1)+4,1)

      ### _check the disturbance flag ---
      # flag=0: implement CCDC
      # flag=1: just after time series change
      # flag=2: less than 12 observation
      # flag=3: need model initializaion(RIRLS)
      if(Values[i, paste0("disturbance_flag_", index)] == 1) next
      if(Values[i, paste0("disturbance_flag_", index)] == 2) next

      ### _check flag for DOY of model initialization---
      disFlag_k <- which( Values[,paste0("disturbance_flag_",index)] == 1)
      if(length(disFlag_k) > 0){
        ini_k <- max(disFlag_k)
      }else{
        ini_k <- 1
      }

      ### _model initialization if needed ---
      if(Values[i, paste0("disturbance_flag_", index)] == 3){
        # _RIRLS
        if(i != length(doys))
          Values <- init_rirls(doys, Values, ini_k, i-3, index, rirls)
        # _model initilaization
        Values <- init_ccdc(doys, Values, index,
                            col_adj = dim(values)[2]+length(VI)+seq(4*(k-1)+1,4*(k-1)+4,1),
                            wch_start = ini_k,
                            wch_end = i-3, ccdc, rmse)
      }

      ### _ccdc coefficients ---
      ccdc_k <- ccdc(doys[ini_k:(i-3)], Values[ini_k:(i-3), index], NULL, NULL)
      Values[i-3, col_coef] <- ccdc_k

      ### _ccdc RMSE ---
      pred_val <- ccdc(doys[ini_k:(i-3)], Values[ini_k:(i-3), index], doys[ini_k:i], NULL)
      pred_val_head <- pred_val[1:(length(pred_val)-3)]
      pred_val_tail <- pred_val[(length(pred_val)-2):length(pred_val)]

      rmse_k <- rmse(Values[ini_k:(i-3), index], pred_val_head) # calculate RMSE
      Values[i-3, paste0(index, "_RMSE")] <- rmse_k

      obs3rmse <- sum(abs(pred_val_tail - Values[(i-2):i, index]) > rmse_k * 3, na.rm = T)

      ### _model initialization flag (if 3 observation > 3*RMSE) ---
      if(obs3rmse == 3){
        init.s <- i-2
        init.e <- i-2+18-1  # original = 15
        if(init.e > dim(Values)[1]) init.e <- dim(Values)[1]

        Values[init.s,                paste0("disturbance_flag_", index)] <- 1
        Values[(init.s+1):(init.e-1), paste0("disturbance_flag_", index)] <- 2
        Values[init.e,                paste0("disturbance_flag_", index)] <- 3
      }
    }
  }

  if(fillNA){
    # fill the NA values in L8 data
    l8 <- Values[,c("B2", "B3", "B4", "B5", "B6", "B7")]
    for(k in 1:length(VI)){
      index <- VI[k]

      switch(index,
             "TCB5" = Res <- colSums(t(l8)*coef_TCB5) ,
             "TCG5" = Res <- colSums(t(l8)*coef_TCG5) ,
             "TCW5" = Res <- colSums(t(l8)*coef_TCW5) ,
             "TCA5" = Res <- atan(colSums(t(l8)*coef_TCG5) / colSums(t(l8)*coef_TCB5)) ,

             "TCB7" = Res <- colSums(t(l8)*coef_TCB7) ,
             "TCG7" = Res <- colSums(t(l8)*coef_TCG7) ,
             "TCW7" = Res <- colSums(t(l8)*coef_TCW7) ,
             "TCA7" = Res <- atan(colSums(t(l8)*coef_TCG7) / colSums(t(l8)*coef_TCB7)) ,

             "TCB8" = Res <- colSums(t(l8)*coef_TCB8) ,
             "TCG8" = Res <- colSums(t(l8)*coef_TCG8) ,
             "TCW8" = Res <- colSums(t(l8)*coef_TCW8) ,
             "TCA8" = Res <- atan(colSums(t(l8)*coef_TCG8) / colSums(t(l8)*coef_TCB8)) ,

             "NBR"  = Res <- (l8[,4] - l8[,6]) / (l8[,4] + l8[,6]),
             "NDVI" = Res <- (l8[,4] - l8[,3]) / (l8[,4] + l8[,3]),
             "NDMI" = Res <- (l8[,4] - l8[,5]) / (l8[,4] + l8[,5]),
             stop("Set Valid index")
      )
      Values[,sub("[[:digit:]]+", "", index)] <- Res
    }
  }

  return(Values)

}

#' Implementation of CCDC for time series data (Sentinel-1)
#'
#' Implementation of CCDC fitting and variable deriviation for time series Sentinel-1 data.
#'
#' Main function for CCDC implementation. For Landsat 8 data, use \code{\link{ccdcTimeSeries}}.
#'
#' @param values a dataframe with time series Landsat data.
#' @param VI character. A spectral index to fit harmonic model.
#' @param startDOY numeric. A julian day of starting CCD change detection.
#' @param init_rirls (function)
#' @param rirls (function)
#' @param ccdc (function)
#' @param init_ccdc (function)
#' @param rmse (function)
#'
#' @return a dataframe with CCDC results and potential disturbance flags.
#'
#' @include rirls.R
#' @include ccdc.R

ccdcTimeSeriesSAR <- function(values, VI = c("VV","VH"), startDOY = 2016, init_rirls,
                              rirls, ccdc, init_ccdc, rmse){

  # add columns ---
  Values <- values # copy
  names.val <- c(paste0(rep(VI, each = 4), paste0("_coef", 1:4)),
                 paste0(VI, "_RMSE"),
                 "rirls_flag",
                 paste0("disturbance_flag_", VI)
  )

  Values[,names.val] <- 0
  doys <- Values[,"DOY"]


  # RIRLS implementation ---
  # at least 12 observation is needed to start detect change
  # here observation before 2016/01/01 was used for model initialization.
  strt <- min(which(doys >= startDOY))
  #Values <- init_rirls(doys, Values, 1, strt-3, sub("[[:digit:]]+", "", VI))

  # estimate coef for initial sequence ---
  for(k in 1:length(VI)){
    index <- VI[k]
    Values <- init_ccdc(doys, Values, index,
                        col_adj = dim(values)[2]+seq(4*(k-1)+1,4*(k-1)+4,1),
                        wch_start = 1,
                        wch_end = strt-3, ccdc, rmse)
  }


  # start CCDC implementation ---
  for(i in strt:length(doys)){

    # _ccdc model ---
    for(k in 1:length(VI)){
      index <- VI[k]
      col_coef <- dim(values)[2]+seq(4*(k-1)+1,4*(k-1)+4,1)

      ### _check the disturbance flag ---
      # flag=0: implement CCDC
      # flag=1: just after time series change
      # flag=2: less than 12 observation
      # flag=3: need model initializaion(RIRLS)
      if(Values[i, paste0("disturbance_flag_", index)] == 1) next
      if(Values[i, paste0("disturbance_flag_", index)] == 2) next

      ### _check flag for DOY of model initialization---
      disFlag_k <- which( Values[,paste0("disturbance_flag_",index)] == 1)
      if(length(disFlag_k) > 0){
        ini_k <- max(disFlag_k)
      }else{
        ini_k <- 1
      }

      ### _model initialization if needed ---
      if(Values[i, paste0("disturbance_flag_", index)] == 3){
        # _RIRLS
        #Values <- init_rirls(doys, Values, ini_k, i-3, index)
        # _model initilaization
        Values <- init_ccdc(doys, Values, index,
                            col_adj = col_coef,
                            wch_start = ini_k,
                            wch_end = i-3, ccdc, rmse)
      }

      ### _ccdc coefficients ---
      ccdc_k <- ccdc(doys[ini_k:(i-3)], Values[ini_k:(i-3), index], NULL, NULL)
      Values[i-3, col_coef] <- ccdc_k

      ### _ccdc RMSE ---
      pred_val <- ccdc(doys[ini_k:(i-3)], Values[ini_k:(i-3), index], doys[ini_k:i], NULL)
      pred_val_head <- pred_val[1:(length(pred_val)-3)]
      pred_val_tail <- pred_val[(length(pred_val)-2):length(pred_val)]

      rmse_k <- rmse(Values[ini_k:(i-3), index], pred_val_head) # calculate RMSE
      Values[i-3, paste0(index, "_RMSE")] <- rmse_k

      obs3rmse <- sum(abs(pred_val_tail - Values[(i-2):i, index]) > rmse_k * 3, na.rm = T)

      ### _model initialization flag (if 3 observation > 3*RMSE) ---
      if(obs3rmse == 3){
        init.s <- i-2
        init.e <- i-2+15-1
        if(init.e > dim(Values)[1]) init.e <- dim(Values)[1]

        Values[init.s,                paste0("disturbance_flag_", index)] <- 1
        Values[(init.s+1):(init.e-1), paste0("disturbance_flag_", index)] <- 2
        Values[init.e,                paste0("disturbance_flag_", index)] <- 3
      }
    }
  }

  return(Values)

}
