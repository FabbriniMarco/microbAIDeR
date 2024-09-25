inc <- function(x, increment = 1) {
  eval.parent(substitute(x <- x + increment))
}