#' Calculate the median odds ratio
#'
#' @param my.var = variance associated with level 2 clustering variable
#' @param digits = number of decimal places to which MOR value will be rounded.
MORcalc <- function(my.var, digits = 2)
{ # MOR arguments: my.var = variance associated with level 2 clustering variable
  # digits = number of decimal places to which MOR value will be rounded.

  Median.OR <- round(exp(sqrt(2*my.var)*qnorm(.75)), digits)
  return(Median.OR)}
