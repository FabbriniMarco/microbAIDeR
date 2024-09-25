se.median <- function(x){
  x = x[!is.na(x)]
  se.med <- 1.2533 * IQR(x) / sqrt(length(x))
  return(se.med)
}