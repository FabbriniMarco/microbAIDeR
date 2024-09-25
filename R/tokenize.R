tokenize <- function(x) {
  if(!is.factor(x)){
    stop("Input object must be a factor")
  } else {
    tok <- c()
    for (len in 1:length(x)) {
      tok = c(tok, which(levels(x) == x[len]) )
    }
    return(tok)
  }
}