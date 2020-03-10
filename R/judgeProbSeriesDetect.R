
#' detect consective anomalies---
#'
#' @param x a vector of logical.
#' @param pattern a pattern of matching vector. Should be a vector of ligical.

matchVector<-function(x, pattern){
  n_pattern <- length(pattern)

  for(i in 1:n_pattern){
    if(i==1){
      l <- which(x == pattern[i])
    }else{
      mtch_i <- x[i:length(x)]
      wch_match <- which(mtch_i == pattern[i])
      l <- intersect(l, wch_match)
    }
  }
  return(l)
}


#' Detect disturbance in time series
#'
#' Detect disturbance in time series of CCDC and RF results. If \code{dt_l8_use} = NULL, then only S1 data is used for
#' disturbance detection.  If \code{dt_s1_use} = NULL, then only L8 data is used for disturbance detection.
#'
#' @param dt_l8_use a dataframe form Landsat 8 (or NULL). Should have 'ID', 'DOY', 'predict' columns.
#' @param dt_s1_use a dataframe form Sentinel-1 (or NULL). Should have 'ID', 'DOY', 'predict' columns.
#' @param ids vector. a list of IDs to be mapped.
#' @param threshold a threshold to identify disturbance.
#' @param seq_l a consective vector to judsge disturbance.
#' @param lag_m (Deprecated). lag date for identifying disturbance.
#' @param savename filenames for saving plots. If NULL, no plots.
#'
#' @importFrom ggplot2 ggplot geom_line geom_point geom_vline scale_x_continuous scale_y_continuous ggtitle xlab theme
#'
#' @return a dataframe with disturbance detection.
#'

judgeProbSeriesDetect <- function(dt_l8_use, dt_s1_use, ids,
                                  threshold, seq_l, lag_m, savename = NULL){

  res_judge <- data.frame(ID = ids,
                          estDateL8 = 0,
                          estDateS1 = 0,
                          estDateMix = 0)

  if(!is.null(dt_l8_use)){
    dt_l8_ci <- split(dt_l8_use, dt_l8_use$ID)
  }
  if(!is.null(dt_s1_use)){
    dt_s1_ci <- split(dt_s1_use, dt_s1_use$ID)
  }

  pb <- txtProgressBar(1, length(ids), style = 3)
  for(i in 1:length(ids)){
    ### extract target ID trajectory ---
    id_i <- ids[i]

    ## L8 ---
    if(!is.null(dt_l8_use)){
      # each ID ---
      dt_temp_l8 <- dt_l8_ci[[i]]
      dt_temp_l8$Sen <- "L8"

      # judge prediction by points (L8) ---
      ProbSeries <- dt_temp_l8$predict
      DisSeries <- ProbSeries >= threshold
      DisSequence <- matchVector(DisSeries, seq_l)
      isDis <- (length(DisSequence) > 0)

      if(isDis){
        DOY_i <- dt_temp_l8$DOY[min(DisSequence)]
        res_judge$estDateL8[i] <- DOY_i
      }else{

      }
    }
    ## S1 ---
    if(!is.null(dt_s1_use)){
      # each ID ---
      dt_temp_s1 <- dt_s1_ci[[i]]
      dt_temp_s1$Sen <- "S1"

      # judge prediction by points (S1) ---
      ProbSeries <- dt_temp_s1$predict
      DisSeries <- ProbSeries >= threshold
      DisSequence <- matchVector(DisSeries, seq_l)
      isDis <- length(DisSequence) > 0

      if(isDis){
        DOY_i <- dt_temp_s1$DOY[min(DisSequence)]
        res_judge$estDateS1[i] <- DOY_i
      }else{

      }
    }

    ## L8 & S1 ---
    if(!is.null(dt_l8_use) & !is.null(dt_s1_use)){
      ### prepare mix L8 and S1 disturbance prob---
      dt_temp_l8s1 <- rbind(dt_temp_l8[,c("DOY", "predict", "Sen")],
                            dt_temp_s1[,c("DOY", "predict", "Sen")])
      dt_temp_l8s1 <- dt_temp_l8s1[order(dt_temp_l8s1$DOY),]
      dt_temp_l8s1$predict_mix <- dt_temp_l8s1$predict
      dt_temp_l8s1$predict_sen <- dt_temp_l8s1$Sen

      #  use either of predicted value that is close to 0 or 1
      for(j in 2:dim(dt_temp_l8s1)[1]){
        temp_sen <- dt_temp_l8s1$Sen[j]
        if(temp_sen == dt_temp_l8s1$predict_sen[j-1]) next

        temp_mix <- dt_temp_l8s1$predict_mix[j-1]
        temp_pre <- dt_temp_l8s1$predict[j]
        if(temp_mix < 0.5) temp_mix <- 1 - temp_mix
        if(temp_pre < 0.5) temp_pre <- 1 - temp_pre

        if(temp_mix > temp_pre){
          dt_temp_l8s1$predict_mix[j] <- dt_temp_l8s1$predict_mix[j-1]
          dt_temp_l8s1$predict_sen[j] <- dt_temp_l8s1$predict_sen[j-1]
        }else{
          dt_temp_l8s1$predict_mix[j] <- dt_temp_l8s1$predict[j] # this makes no change
        }
      }

      # judge prediction by points (mix) ---
      ProbSeries <- dt_temp_l8s1$predict_mix
      DisSeries <- ProbSeries >= threshold
      DisSequence <- matchVector(DisSeries, seq_l)
      isDis <- length(DisSequence) > 0

      if(isDis){
        DOY_i <- dt_temp_l8s1$DOY[min(DisSequence)]
        res_judge$estDateMix[i] <- DOY_i
      }else{

      }

    }


    ### plot disturbance prob for S1 & L8 (optional)---
    if(!is.null(savename)){
      if(i == 1) pdf(savename, width = 8, height = 6)
      ggplot + theme_set(theme_bw(base_size = 12))

      a <- ggplot(dt_temp_l8s1, aes(x = DOY)) +
        geom_line(aes(y = predict_mix), colour = "grey70", size = 2) +
        geom_line(aes(y = predict, group = Sen, colour = Sen), size = 0.2) +
        geom_point(aes(y = predict, colour = Sen), size = 2) +
        geom_vline(xintercept = DisDate, colour = "blue", linetype = "dashed") +
        scale_x_continuous(breaks = seq(2016, 2018, 0.5), limits = c(2016, 2018)) +
        scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
        ggtitle(paste0("ID: ", id_use[id_i])) +
        xlab("Year") + ylab("Disturbance probability") +
        theme(plot.title = element_text(hjust = 0.5))

      print(a)

      if(i == length(which_test)) dev.off()
    }

    setTxtProgressBar(pb, i)
    if(i %% 1000 == 0) invisible({gc();gc()})
  }

  return(res_judge)
}

