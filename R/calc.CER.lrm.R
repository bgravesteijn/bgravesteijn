#' Calculate effect sizes for CER analyses
#'
#' @param mids.data A MIDS: a MICE imputed dataset
#' @param orig.data The original data.frame/data.table which was imputed
#' @param instrument A character of the name of the instrumental variable
#' @param therapy A character of the name of the therapy
#' @param outcome A character of the outcome parameter
#' @param subset.cer A logical vector (with length equal to the number of cases), which indicates included and excluded cases
#' @return Returns a vector of the estimate with 95% CI of the treatment effect, including the R2 of the stage 1 regression.
#' @description This function uses two-stage least squares regression to estimate a causal treatment effect with an instrumental variable. For info or explanation, please contact: b.gravesteijn@erasmusmc.nl
calc.CER.lrm <- function(mids.data = NULL, orig.data=NULL,
                         instrument = NULL, therapy= NULL, outcome=NULL,
                         subset.cer=NULL){
##stage 1##
fml.st1 <- as.formula(paste(therapy, "~", instrument))
st1 <- fit.mult.impute(fml.st1, xtrans=mids.data, data = dt.ih,
                       fitter = lrm, subset = subset.cer, pr = FALSE)

print(paste("R2 stage 1 regression =", round(st1$stats["R2"],3)))

##stage 2##
out.2   <- vector(mode = "list", length = 5)

for(i in 1:5){
  fit  <- lrm(complete(mids.data,i)[subset.cer,][,outcome]~st1$linear.predictors)
  beta <- fit$coefficients["st1"]
  se   <- sqrt(diag(fit$var))["st1"]
  out.2[[i]] <- c(beta, se)
}

#extract estimates
betas <- c(out.2[[1]][1],out.2[[2]][1],out.2[[3]][1],out.2[[4]][1],out.2[[5]][1])
ses   <- c(out.2[[1]][2],out.2[[2]][2],out.2[[3]][2],out.2[[4]][2],out.2[[5]][2])

#pool estimate
est <- mean(betas)

#pool se using rubin's rules
m <- 5
varwithin <- (1/m)*sum(ses^2)
varbetwen <- (1/(m-1))*sum((betas-est)^2)
setot <- sqrt(varwithin+varbetwen*(1+(1/m)))

#print result out
print(paste(outcome))
print(paste("OR =", round(exp(est),3)))
print(paste("95% CI =", round(exp(est-1.96*setot),3), "-", round(exp(est+1.96*setot),3)))

#make vector of effect size
res        <- c(exp(est), exp(est-1.96*setot), exp(est+1.96*setot))
names(res) <- c("OR", "2.5%", "97.5%")
return(res)
}
