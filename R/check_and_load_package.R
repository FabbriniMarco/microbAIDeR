check_and_load_package <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    stop(paste("The required package", package_name, "is not installed. Please install it using install.packages('", package_name, "') and try again.", sep = ""), "\nHave you run the microbAIDeR::microbAIDeR_install_dependancies() function?\n")
  } else {
    library(package_name, character.only = TRUE)
  }
}