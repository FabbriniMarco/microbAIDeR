amp_up_models <- function(core_fraction = 0.8){
  check_and_load_package("parallel")
  check_and_load_package("doParallel")
  no_cores <- parallel::detectCores(logical=TRUE) * core_fraction
  cluster <- parallel::makePSOCKcluster(no_cores)
  doParallel::registerDoParallel(cluster)
  cat("Model amped and ready to go with:", no_cores, "cores. \n")
}