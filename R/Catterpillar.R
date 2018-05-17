

#' Make a caterpillar plot for random effects
#'
#' @param x Alist of fitted models (for multiple imputed dataset fitted models)
#' @param fitter A character indicating random effects fitting formula
#' @param grp.var.t = character indicating type of grouping variable, e.g.: "country"
#' @param MOR = logical indicating whether the median odds ratio should be printed, default is TRUE
#' @return Returns a plot containing mean and 95% CI, with median odds ratio
catterpillar<-function(x = NULL, fitter= NULL, grp.var.t = NULL, MOR=TRUE){
  #x = list of fitted models (for multiple imputed dataset fitted models)
  #fitter = character indicating random effects fitting formula
  #grp.var.t = character indicating type of grouping variable, e.g.: "country"
  #MOR = logical indicating whether the median odds ratio should be printed, default is TRUE

  library(lme4)
  library(ordinal)
  library(ggplot2)

  MORcalc <- function(my.var, digits = 2)
  { # MOR arguments: my.var = variance associated with level 2 clustering variable
    # digits = number of decimal places to which MOR value will be rounded.

    Median.OR <- round(exp(sqrt(2*my.var)*qnorm(.75)), digits)
    return(Median.OR)}

  length_list<-length(x)

  if(length_list==1){
    length_grp<-length(data.frame(lme4::ranef(x, which=grp.var.t))$grp)
  }else{
    length_grp<-length(data.frame(lme4::ranef(x[[i]], which=grp.var.t))$grp)
  }


  #grp = vector of the grouping variable, e.g. countries/schools...
  #re = vector of the random effects for each value of the grp variable
  #condvar = vector of the conditional variances of the random effects of the grouping variable
  #xlab = character indicating type of grouping variable, e.g.: "country"
  MOR<-rep(NA, length(x))

  if (fitter=="glmer"){
    if(length_list==1){
      plot.df<-data.frame(lme4::ranef(x, which=grp.var.t, condVar=TRUE))
      plot.df$lo<-plot.df$condval-1.96*plot.df$condsd
      plot.df$hi<-plot.df$condval+1.96*plot.df$condsd
      MOR<-MORcalc(data.frame(lme4::VarCorr(x))$sdcor^2)
    }else{
      condval<-matrix(data = rep(NA, length_grp*length_list), ncol = length_list)
      sd<-matrix(data = rep(NA, length_grp*length_list), ncol = length_list)
      for (i in 1:length_list){
        condval[,i]<-data.frame(lme4::ranef(x[[i]], which=grp.var.t, condVar=TRUE))$condval
        sd[,i]<-data.frame(lme4::ranef(x[[i]], which=grp.var.t, condVar=TRUE))$condsd
        MOR[i]<-MORcalc(data.frame(lme4::VarCorr(x[[i]]))$sdcor^2)
      }
      plot.df<-data.frame(lme4::ranef(x[[1]], which=grp.var.t, condVar=TRUE))
      plot.df$condval<-apply(condval, 1, mean)
      plot.df$condsd<-apply(sd, 1, mean)
      plot.df$lo<-plot.df$condval-1.96*plot.df$condsd
      plot.df$hi<-plot.df$condval+1.96*plot.df$condsd
      MOR<-mean(MOR)
    }
  }else{
    if(fitter%in%c("clmm2", "clmm")){
      if(length_list==1){
        plot.df<-data.frame(ordinal::ranef(x, which=grp.var.t, condVar=TRUE))
        plot.df$lo<-plot.df$condval-1.96*plot.df$condsd
        plot.df$hi<-plot.df$condval+1.96*plot.df$condsd
        MOR<-MORcalc(data.frame(ordinal::VarCorr(x))$sdcor^2)
      }else{
        condval<-matrix(data = rep(NA, length_grp*length_list), ncol = length_list)
        sd<-matrix(data = rep(NA, length_grp*length_list), ncol = length_list)
        for (i in 1:ncol(condval)){
          condval[,i]<-data.frame(ordinal::ranef(x[[i]], which=grp.var.t, condVar=TRUE))$condval
          sd[,i]<-data.frame(ordinal::ranef(x[[i]], which=grp.var.t, condVar=TRUE))$condsd
          MOR[i]<-MORcalc(data.frame(ordinal::VarCorr(x[[i]]))$sdcor^2)
        }
        plot.df<-ordinal::ranef(x[[1]], which=grp.var.t, condVar=TRUE)
        plot.df$condval<-apply(condval, 1, mean)
        plot.df$condsd<-apply(sd, 1, mean)
        plot.df$lo<-plot.df$condval-1.96*plot.df$condsd
        plot.df$hi<-plot.df$condval+1.96*plot.df$condsd
        MOR<-mean(MOR)
      }
    }else{
      print("This fitter is not yet available in this package")
    }
  }

  plot.df$grp <- factor(plot.df$grp, levels = plot.df$grp[order(plot.df$condval)])
  #catterpillar plot
  if(MOR==TRUE){
    ggplot(data = plot.df, aes(x=grp, y=condval, ymin=lo, ymax=hi))+geom_pointrange()+
      geom_hline(yintercept=0, linetype=2)+coord_flip()+xlab(label = grp.var.t)+
      scale_y_continuous(name = "Log odds")+
      annotate("text", x=3, y=2, label=paste("MOR =", MOR))
  }else{
    ggplot(data = plot.df, aes(x=grp, y=condval, ymin=lo, ymax=hi))+geom_pointrange()+
      geom_hline(yintercept=0, linetype=2)+coord_flip()+xlab(label = grp.var.t)+
      scale_y_continuous(name = "Log odds")
  }

}
