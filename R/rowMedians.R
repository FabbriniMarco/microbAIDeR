rowMedians <- function(df) {
  apply(df, 1, median, na.rm = TRUE)
}