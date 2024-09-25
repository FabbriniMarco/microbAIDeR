sem <- function(x) {
  x = x[!is.na(x)]
  n <- length(x)
  se <- sd(x) / sqrt(n)
  return(se)
}
