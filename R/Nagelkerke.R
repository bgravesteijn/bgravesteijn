
#'Calculate Nagerkere R2 for fitted models
#'
#' @param fit a (list of) fitted model(s), which need to be compared to the null model
#' @param nullmod a null model, to which the fitted model(s) need to be compared
#' @param data the dataset on which the model is fitted
#'
#' @return the Nagelkerke R2, as a numer between 0 and 1
Nagelkerke<-function(fit=NULL, nullmod=NULL, data=NULL){
  if (length(x)==1){
    CoxSnell <- 1 - exp(-(2 / nrow(data)) * (logLik(fit) - logLik(nullmod)))
    R2 <- CoxSnell / (1 - exp((2 * nrow(data) ^ (-1)) * logLik(nullmod)))
  }else{
    R2v<-rep(NA, length(fit))
    for(i in 1:length(fit)){
      CoxSnell <- 1 - exp(-(2 / nrow(data) * (logLik(fit[[i]]) - logLik(nullmod))))
      R2v[i] <- CoxSnell / (1 - exp((2 * nrow(data) ^ (-1)) * logLik(nullmod)))
    }
    R2<-mean(R2v)
  }
  return(R2)
}

