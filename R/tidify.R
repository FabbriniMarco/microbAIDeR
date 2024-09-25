tidify <- function(x, as.numeric = TRUE, colNames = c("rowname", "colname", "value")){
  check_and_load_package("tidyr")
  check_and_load_package("tibble")
  if (length(colNames) != 3) {
    stop("colNames must be a vector of length 3, e.g., c('rowname', 'colname', 'value').")
  }
  y <- x %>%
    tibble::rownames_to_column() %>%
    tidyr::gather(colname, value, -rowname)
  if(as.numeric){
    y$value = as.numeric(y$value)
    if (any(is.na(y$value))) {
      message("Some values could not be converted to numeric. They have been replaced with NA.")
    }
  }
  colnames(y) = colNames
  return(y)
}