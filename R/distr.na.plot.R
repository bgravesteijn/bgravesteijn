
#' Distribution and NA plot
#'
#' @param x An object with data
#' @return Returns a plot showing distribution and missingness of variables
#' @description This plot is copyright of Nicole Erler. Follow her course at NIHES: missing values.
distr.na<-function(x=NULL){
  nc <- max(5, ceiling(sqrt(ncol(x))))
  nr <- ceiling(ncol(x) / nc)
  par(mfrow = c(nr, nc), mgp = c(2, 0.6, 0), mar = c(2, 3, 3, 0.5))
  for (i in 1:ncol(x)) {
    x<-data.frame(x)
    if (is.numeric(x[, i])) {
      hist(x[, i], nclass = 50, xlab = "",
           main = paste0(names(x[i]), " (",
                         round(mean(is.na(x[, i])) * 100, 2), "% NA)")
      )
    } else {
      barplot(table(x[, i]), ylab = "Frequency",
              main = paste0(names(x[i]), " (",
                            round(mean(is.na(x[, i])) * 100, 2), "% NA)"))
    }
  }

}


