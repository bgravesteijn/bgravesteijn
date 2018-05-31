#' Calculate the median odds ratio and CI using bootstrap
#'
#' @param formula a character containing the formula for the glmer function
#' @param data a (list of) data.frame(s) or data.table(s) containing the data on which the model was fitted
#' @param grp.var.t the grouping variable, e.g.: "country"
MOR.ci.calc <- function (formula=NULL, data=NULL, grp.var.t=NULL){

   fit.mor<-function(x, indices){
    fml<-as.formula(formula)
    fit<-lme4::glmer(formula=fml, data=x[indices,], family="binomial")
    var<-data.table::data.table(data.frame(lme4::VarCorr(fit)))[grp==grp.var.t,]$vcov
    mor<-bgravesteijn::MORcalc(my.var = var)
    return(mor)
   }

  n.fit<-length(fit)

 if(n.fit==1){
  MORS<-boot::boot(data = data, statistic = fit.mor,R = 1000)
  MORS<-MORS$t
 }else{
   mors.list<-vector("list", n.fit)
  for (i in 1:n.fit){
    mors<-boot::boot(data=data[[i]], statistic=fit.mor, R=1000)
    mors.list[[i]]<-mors$t
  }
   MORS<-unlist(mors.list)
 }
 MOR.final<-paste("MOR =",quantile(MORS, probs=0.5),
                  "(95% CI:", quantile(MORS, probs=0.025), "-", quantile(MORS, probs=0.0975))
 return(MOR.final)
}
