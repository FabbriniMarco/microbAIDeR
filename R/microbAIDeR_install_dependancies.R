microbAIDeR_install_dependancies <- function(lib = .libPaths()[1]){
  cran_packages <- c('devtools', 'parallel', 'doParallel', 'openxlsx', 
                     'tidyr', 'tibble', 'dplyr', 'ggplot2', 'ggsignif', 
                     'gridExtra', 'tidyverse', 'stringr', 'vegan', 
                     'cowplot', 'gtable', 'grid', 'ggrepel', 'scales', 'svglite')
  install_dependencies <- function(packages, lib) {
    for (package in packages) {
      if (!requireNamespace(package, quietly = TRUE)) {
        install.packages(package, lib = lib)
      } else {
        message(paste("Package", package, "is already installed."))
      }
    }
  }
  
  install_dependencies(cran_packages, lib)
  if (!requireNamespace("pairwiseAdonis", quietly = TRUE)) {
    tryCatch({
      devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis", lib = lib, upgrade = "never")
    }, error = function(e) {
      message("Error in installing pairwiseAdonis: ", e)
    })
  } else {
    message("pairwiseAdonis is already installed.")
  }
}