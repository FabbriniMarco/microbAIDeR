sigFunction <- function(x, trends=TRUE) {
  if (!is.na(x)) {
    if (x <= 0.0001) {
      "****"
    } else if (x <= 0.001) {
      "***"
    } else if (x <= 0.01) {
      "**"
    } else if (x <= 0.05) {
      "*"
    } else if (trends & x <= 0.1) {
      "°"
    } else {
      NA
    }
  } else {
    NA
  }
}