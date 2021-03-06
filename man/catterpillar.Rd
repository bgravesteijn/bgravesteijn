% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Catterpillar.R
\name{catterpillar}
\alias{catterpillar}
\title{Make a caterpillar plot for random effects}
\usage{
catterpillar(x = NULL, fitter = NULL, grp.var.t = NULL,
  plotMOR = TRUE, plotlabels = TRUE, xMOR = 3, yMOR = 2,
  set.MOR = NULL, col.group = NULL)
}
\arguments{
\item{x}{Alist of fitted models (for multiple imputed dataset fitted models)}

\item{fitter}{A character indicating random effects fitting formula}

\item{grp.var.t}{character indicating type of grouping variable, e.g.: "country"}

\item{xMOR}{x coordinate of the text of MOR (default is 3)}

\item{yMOR}{y coordinate of the text of MOR (default is 2)}

\item{col.group}{vector of binary (0 or 1) grouping variable}

\item{printMOR}{logical indicating whether the median odds ratio should be printed, default is TRUE}
}
\value{
Returns a plot containing mean and 95% CI, with median odds ratio, using Rubin's rules to pool the variances
}
\description{
Make a caterpillar plot for random effects
}
