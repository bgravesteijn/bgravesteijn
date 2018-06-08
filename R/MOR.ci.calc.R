#' Calculate the median odds ratio and CI using bootstrap
#'
#' @param formula a character containing the formula for the glmer function
#' @param data a (list of) data.frame(s) or data.table(s), or a MIDS object
#'  containing the data on which the model was fitted
#' @param grp.var.t the grouping variable of interest, e.g.: "country"
#' @param cluster.var cluster variable to use for clustered bootstrap sampling (vector)
#' @param set.subset a logical vector, indicating whether a subject is included or not in the analysis
#' @param rep number of bootstrap replicates should be used 
#' @description  Calculate the median odds ratio and CI using bootstrap (1000 replicate samples). Note: this takes a long time... So maybe plan a Titanic night with your colleagues while using this function.
MOR.ci.calc <- function (formula=NULL, data=NULL, grp.var.t=NULL, cluster.var=NULL, set.subset=NULL, reps=1000){
   fit.mor<-function(x){
    fml<-as.formula(formula)
    fit<-lme4::glmer(formula=fml, data=x, subset=set.subset, family="binomial")
    var<-data.table::data.table(data.frame(lme4::VarCorr(fit)))[grp==grp.var.t,]$vcov
    mor<-bgravesteijn::MORcalc(my.var = var)
    return(mor)
   }
   
   if(class(data)=="mids"){
     data.l<-vector("list", data$m)
     for(i in 1:data$m){
       data.l[[i]]<-mice::complete(data, i)
     }
     data<-data.l
   }

  n.fit<-ifelse(class(data)=="list",length(data),1)

 if(n.fit==1){
  MORS<-cluster.boot(data = data, statistic = fit.mor, cluster.v=cluster.var,R = reps)
 }else{
   mors.list<-vector("list", n.fit)
  for (i in 1:n.fit){
    mors.list[[i]]<-cluster.boot(data=data[[i]], statistic=fit.mor, cluster.v = cluster.var, R=reps)
  }
   MORS<-unlist(mors.list)
 }
 MOR.final<-paste("MOR =",quantile(MORS, probs=0.5),
                  "(95% CI:", quantile(MORS, probs=0.025), "-", quantile(MORS, probs=0.975), ")")
 return(MOR.final)
}
