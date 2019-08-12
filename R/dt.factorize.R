dt.factorize <- function(data=NULL, facvars=NULL){
  data <- data.frame(data)
  for(i in facvars){
    df[,i] <- factor(df[,i])
  }
  df <- data.table(df)
  return(df)
}
