distr.stats <- function(x) {
  x <- x[!is.na(x) & !is.nan(x)]
  n <- length(x)
  avg <- mean(x)
  sd_avg <- sd(x)
  se_avg <- sd_avg / sqrt(n)
  
  med <- median(x)
  iqr <- IQR(x)
  se_med <- 1.2533 * iqr / sqrt(n)
  
  list(mean = avg, SD = sd_avg, SEM = se_avg, median = med, iqr = iqr, SEMed = se_med)
}