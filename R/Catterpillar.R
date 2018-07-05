

#' Make a caterpillar plot for random effects
#'
#' @param x Alist of fitted models (for multiple imputed dataset fitted models)
#' @param fitter A character indicating random effects fitting formula
#' @param grp.var.t character indicating type of grouping variable, e.g.: "country"
#' @param printMOR logical indicating whether the median odds ratio should be printed, default is TRUE
#' @param xMOR x coordinate of the text of MOR (default is 3)
#' @param yMOR y coordinate of the text of MOR (default is 2)
#' @param col.group vector of binary (0 or 1) grouping variable
#' @return Returns a plot containing mean and 95% CI, with median odds ratio, using Rubin's rules to pool the variances
catterpillar<-function (x = NULL, fitter = NULL, grp.var.t = NULL, plotMOR = TRUE,
                        plotlabels = TRUE, xMOR = 3, yMOR = 2, set.MOR = NULL, col.group=NULL)
{
  library(ggplot2)
  library(mitools)
  MORcalc <- function(my.var, digits = 2) {
    Median.OR <- round(exp(sqrt(2 * my.var) * qnorm(0.75)),
                       digits)
    return(Median.OR)
  }
  length_list <- length(x)
  MOR <- rep(NA, length(x))
  if (fitter == "glmer") {
    if (length_list == 1) {
      length_grp <- length(data.frame(lme4::ranef(x, which = grp.var.t))$grp)
    }
    else {
      length_grp <- length(data.frame(lme4::ranef(x[[1]],
                                                  which = grp.var.t))$grp)
    }
    if (length_list == 1) {
      re.no <- which(data.frame(lme4::VarCorr(x, which = grp.var.t))$grp ==
                       grp.var.t)
      plot.df <- data.frame(lme4::ranef(x, which = grp.var.t,
                                        condVar = TRUE))
      plot.df$lo <- plot.df$condval - 1.96 * plot.df$condsd
      plot.df$hi <- plot.df$condval + 1.96 * plot.df$condsd
      MOR <- round(MORcalc(data.frame(lme4::VarCorr(x))[re.no,
                                                        ]$sdcor^2), 2)
    }
    else {
      re.no <- which(data.frame(lme4::VarCorr(x[[1]], which = grp.var.t))$grp ==
                       grp.var.t)
      condval <- matrix(data = rep(NA, length_grp * length_list),
                        ncol = length_list)
      sd <- matrix(data = rep(NA, length_grp * length_list),
                   ncol = length_list)
      for (i in 1:length_list) {
        condval[, i] <- data.frame(lme4::ranef(x[[i]],
                                               which = grp.var.t, condVar = TRUE))$condval
        sd[, i] <- data.frame(lme4::ranef(x[[i]], which = grp.var.t,
                                          condVar = TRUE))$condsd
        MOR[i] <- MORcalc(data.frame(lme4::VarCorr(x[[i]]))[re.no,
                                                            ]$sdcor^2)
      }
      plot.df <- data.frame(lme4::ranef(x[[1]], which = grp.var.t,
                                        condVar = TRUE))
      plot.df$condval <- apply(condval, 1, mean)
      extra.var <- ((length_list + 1)/(length_list * (length_list -
                                                        1))) * (apply(apply(condval, 2, function(x) {
                                                          (x - plot.df$condval)^2
                                                        }), 1, sum))
      plot.df$condsd <- apply(sd, 1, mean) + extra.var
      plot.df$lo <- plot.df$condval - 1.96 * plot.df$condsd
      plot.df$hi <- plot.df$condval + 1.96 * plot.df$condsd
      MOR <- round(mean(MOR), 2)
    }
  }
  else {
    if (fitter == "clmm") {
      if (length_list == 1) {
        length_grp <- nrow(data.frame(ordinal::ranef(x)))
      }
      else {
        length_grp <- nrow(data.frame(ordinal::ranef(x[[1]])))
      }
      if (length_list == 1) {
        plot.df <- data.frame(condval = data.frame(ordinal::ranef(x))$X.Intercept.)
        plot.df$condsd <- sqrt(data.frame(ordinal::condVar(x)))$X.Intercept.
        plot.df$lo <- plot.df$condval - 1.96 * plot.df$condsd
        plot.df$hi <- plot.df$condval + 1.96 * plot.df$condsd
        MOR <- round(MORcalc(data.frame(ordinal::VarCorr(x))$X.Intercept.),
                     2)
      }
      else {
        condval <- matrix(data = rep(NA, length_grp *
                                       length_list), ncol = length_list)
        sd <- matrix(data = rep(NA, length_grp * length_list),
                     ncol = length_list)
        for (i in 1:ncol(condval)) {
          condval[, i] <- data.frame(ordinal::ranef(x[[i]]))$X.Intercept.
          sd[, i] <- sqrt(data.frame(ordinal::condVar(x[[i]])))$X.Intercept.
          MOR[i] <- MORcalc(data.frame(ordinal::VarCorr(x[[i]]))$X.Intercept.)
        }
        plot.df <- data.frame(condval = data.frame(ordinal::ranef(x[[1]],
                                                                  condVar = TRUE))$X.Intercept.)
        plot.df$condval <- apply(condval, 1, mean)
        extra.var <- ((length_list + 1)/(length_list *
                                           (length_list - 1))) * (apply(apply(condval,
                                                                              2, function(x) {
                                                                                (x - plot.df$condval)^2
                                                                              }), 1, sum))
        plot.df$condsd <- apply(sd, 1, mean) + extra.var
        plot.df$lo <- plot.df$condval - 1.96 * plot.df$condsd
        plot.df$hi <- plot.df$condval + 1.96 * plot.df$condsd
        plot.df$grp <- rownames(plot.df)
        MOR <- round(mean(MOR), 2)
      }
    }
    else {
      if (fitter == "clmm2") {
        if (length_list == 1) {
          length_grp <- length(x$ranef)
        }
        else {
          length_grp <- length(x[[1]]$ranef)
        }
        if (length_list == 1) {
          plot.df <- data.frame(condval = x$ranef, condsd = sqrt(x$condVar))
          plot.df$lo <- plot.df$condval - 1.96 * plot.df$condsd
          plot.df$hi <- plot.df$condval + 1.96 * plot.df$condsd
          MOR <- round(MORcalc(x$stDev^2), 2)
          plot.df$grp <- rownames(plot.df)
        }
        else {
          condval <- matrix(data = rep(NA, length_grp *
                                         length_list), ncol = length_list)
          sd <- condval
          for (i in 1:ncol(sd)) {
            condval[, i] <- x[[i]]$ranef
            sd[, i] <- sqrt(x[[i]]$condVar)
            MOR[i] <- MORcalc(x[[i]]$stDev^2)
          }
          plot.df <- data.frame(condval = x[[1]]$ranef,
                                condsd = x[[1]]$condVar)
          plot.df$condval <- apply(condval, 1, mean)
          extra.var <- ((length_list + 1)/(length_list *
                                             (length_list - 1))) * (apply(apply(condval,
                                                                                2, function(x) {
                                                                                  (x - plot.df$condval)^2
                                                                                }), 1, sum))
          plot.df$condsd <- apply(sd, 1, mean) + extra.var
          plot.df$lo <- plot.df$condval - 1.96 * plot.df$condsd
          plot.df$hi <- plot.df$condval + 1.96 * plot.df$condsd
          MOR <- round(mean(MOR), 2)
          plot.df$grp <- rownames(plot.df)
        }
      }
      else {
        print("This fitter is not yet available in this package")
      }
    }
  }
  plot.df$grp <- factor(plot.df$grp, levels = plot.df$grp[order(plot.df$condval)])
  if (!is.null(set.MOR)) {
    MOR <- set.MOR
  }
  plot.df$col.group<-col.group
  if(col.group=NULL){
    plot <- ggplot(data = plot.df, aes(x = grp, y = condval, ymin = lo, ymax = hi))
    + geom_pointrange() + geom_hline(yintercept = 0, scale_y_continuous(name = "Log odds"))
  }else{
    plot <- ggplot(data = plot.df, aes(x = grp, y = condval, ymin = lo, ymax = hi, col=col.group))
    + geom_pointrange() + geom_hline(yintercept = 0, scale_y_continuous(name = "Log odds"))
  }

    if (plotMOR) {
    if (plotlabels) {
      plot + annotate("text", x = yMOR, y = xMOR, label = paste("MOR =",
                                                                MOR))
    }
    else {
      plot + annotate("text", x = yMOR, y = xMOR, label = paste("MOR =",
                                                                MOR)) + theme(axis.text.y = element_blank())
    }
  }
  else {
    if (plotlabels) {
      plot
    }
    else {
      plot + theme(axis.text.y = element_blank())
    }
  }
}

