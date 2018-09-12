#' Instrumental variable analysis for multiple imputed datasets.
#'
#' @param instrument A character value of the name of the instrumental variable (any type)
#' @param therapy A character value of the name of the binary therapy
#' @param outcome A character value of the name of the outcome variable, both binary as ordinal.
#' @param mids A mids object, of the multiple imputed dataset which we will use for the IV analysis.
#' @return Returns a vector with the effect estimate with 95% CI.
#' @description Using an instrumental variable approach, calculate the effect estimate for a therapy. As instrument, any type of instrument may be used (continuous/ordinal/binary/categorical). The therapy must be a binary variable, and the outcome may be both binary (logistic regression) as ordinal (proportional odds regression). Make sure the assumptions are correct, before using this funciton!
IV.calc <- function(instrument=NULL, therapy=NULL, outcome=NULL, mids=NULL){

  calc.IV.lrm <- function(data, indices){
    d <- data[indices,]
    subset.cer <- d$center%in%centersubset & !is.na(d$fa.derv.GOSE)

    #stage1
    fml.st1 <- as.formula(paste(therapy, "~", instrument))
    st1 <- glm(fml.st1, data = d ,subset = subset.cer, family="binomial")

    #stage 2
    st2 <- lrm(d[subset.cer, outcome] ~ plogis(st1$linear.predictors))
    beta <- coef(st2)["st1"]
    return(beta)
  }

  results.list <- vector("list", mids$m)
  for(i in 1:dti.mice$m){
    bootres           <- boot::boot(data = complete(mids,i), statistic = calc.IV.lrm, R = 1000)
    results.list[[i]] <- bootres$t
  }

  ivres        <- quantile(unlist(results.list), probs = c(0.5,0.025,0.975))
  names(ivres) <- c("beta", "lo", "hi")
  return(ivres)
}
