
#' Build random forest model
#'
#' Build a RF model using caret and save result in a subfolder.
#'
#' @param dt_train dataframe. Training data for RF. Reference col should be 'manual'.
#' @param dt_test dataframe. Test data for RF. This should have the same columns as \code{dt_train}.
#' @param dir_rf a directory to save the result of RF. If NULL, a subfolder will be automatically generated.
#' @param name_rdata filefame for saving the RF result. If NULL, name will be automatically generated.
#' @param do_prallel logical. If TRUE, use parallel process to tune RF.
#' @param max_cores numeric. Maximum numbers of cores used for parallel processing.
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom randomForest randomForest
#' @importFrom foreach foreach
#' @importFrom caret trainControl train confusionMatrix varImp
#'
#' @return a list of the results of RF modeling. Rdata file is saved in a subfolder.
#'
implementRF <- function(dt_train, dt_test, dir_rf = NULL, name_rdata = NULL, do_prallel = T, max_cores){

  # save dirs ---
  if(is.null(dir_rf)) dir_rf <- file.path(getwd(), "rf_result")
  if(!dir.exists(dir_rf)) dir.create(dir_rf)

  # prep manual --
  dt_train$manual <- as.factor(dt_train$manual)
  dt_test$manual  <- as.factor(dt_test$manual)
  dt_train <- na.omit(dt_train)
  dt_test  <- na.omit(dt_test)


  if(do_prallel){
    n_cl <- min(max_cores, detectCores())
    cl <- makeCluster(n_cl)
    registerDoParallel(cl)
  }

  # tuning ---
  fitControl <- trainControl(method = "cv", number = 10)   # 10-fold cross validation

  result_tune <- train(manual ~.,
                       data = dt_train,
                       method = "rf",
                       tuneLength = 10,
                       trControl = fitControl,
                       importance = T)
  if(do_prallel) stopCluster(cl)

  # accuracy ---
  mat_conf <- confusionMatrix(data = predict(result_tune, dt_test[,-1]),
                              reference = dt_test[,1])

  # importance and selected mtry ---
  result_imp <- varImp(result_tune, scale = F)
  tunes <- c(result_tune$bestTune[1,1])

  # save RF data ---
  if(!is.null(name_rdata)){
    save(list = c("dt_train", "dt_test", "result_tune", "mat_conf", "result_imp", "tunes"),
         file = file.path(dir_rf, name_rdata))
  }
  # sink RF result ---
  name_sink <- file.path(dir_rf, gsub("Rdata", "txt", name_rdata))
  sinkRFresult(result_tune, mat_conf, name_sink)

  # return
  list_res <- list(result_tune, mat_conf, tunes)
  names(list_res) <- c("result_tune", "mat_conf", "tunes")

  return(list_res)
}

#' sink the result of RF
#'
#' @param rf_tune tune result of RF.
#' @param mat_conf confusion matrix.
#' @param filename a filename to save.
#'
sinkRFresult <- function(rf_tune, mat_conf, filename){
  sink(file = filename)

  cat("Date: ", as.character(Sys.time()))
  cat("\nAccuracy: "); print(mat_conf)
  cat("\n\nModel:"); print(rf_tune)
  cat("\nImportance: "); print(varImp(rf_tune))

  sink()
}


