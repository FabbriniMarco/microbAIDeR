scaling.manual <- function(x, range.min, range.max) {
  return(range.min + (x-min(x))*(range.max-range.min)/(max(x)-min(x)))
}