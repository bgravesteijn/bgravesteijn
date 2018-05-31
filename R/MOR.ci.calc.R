#' Calculate the median odds ratio and CI using bootstrap
#'
#' @param fit a (list of) fitted model(s), fitted with glmer
#' @param data a (list of) data.frame(s) or data.table(s) containing the data on which the model was fitted
MOR.ci.calc <- function (fit=NULL, data=NULL){
 n.fit<-length(fit)
 if(n.fit==1){
 fit.mor<-function(x){
   fml<-fit$call[[2]]
   fit<-lme4::glmer(formula=fml, data=x, family="binomial")
   vc<-lmer::VarCorr(fit)
   var<-print(vc, comp=("Variance"))
   mor<-bgravesteijn::MORcalc(my.var = var)
   return(MOR)
   }
 MORS<-boot::boot(data, fit.mor,1000)$t
 }else{
  print("Work in progress...")
 }
 MOR.final<-paste("MOR =",quantile(MORS, probs=0.5),
                  "(95% CI:", quantile(MORS, probs=0.025), "-", quantile(MORS, probs=0.0975))
 return(MOR.final)
}
