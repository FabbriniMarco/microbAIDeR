one_hot_encode <- function(x){
  if(!is.factor(x)){
    stop("Input object must be a factor")
  } else {
    encoded = data.frame()
    for (len in 1:length(x)){
      encoded = rbind(encoded, rep(0, nlevels(x)))
      encoded[ nrow(encoded) , which(levels(x) == x[len]) ] = 1
    }
    colnames(encoded) = levels(x)
    return(encoded)
  }
}