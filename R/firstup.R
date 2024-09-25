firstup <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  return(string)
}
