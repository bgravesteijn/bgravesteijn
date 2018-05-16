#' Make a ROC plot, with AUC
#'
#' @param x A fitted logistic regression model
#' @param obs.v A vector of the observed outcome variable: (0,1)
#' @return Returns a ROC plot containing AUC
roc.model<-function(x, obs.v){
  library(ROCR)
  obs.df<-data.frame(obs.v)
  obs.df$pred<-predict(x, type="response")
  pred<-prediction(obs.df$pred, obs.df$obs.v)
  perf<-performance(pred, measure = "tpr", "fpr")
  auc<-round(performance(pred, measure="auc")@y.values[[1]],2)
  plot(perf)
  text(x=0.75, y=0.25 ,paste("AUC =", auc))
}
