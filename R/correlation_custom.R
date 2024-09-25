correlation_custom <- function(feature_data, parameter_data, is.binary.parameter = FALSE, p.adjust.method = "fdr",
                               uncorr.p.threshold = 0.05, corr.p.threshold = 0.05, min.corr.threshold = 0.3,
                               force.numeric = TRUE, corr.method = "spearman"){
  check_and_load_package("dplyr")
  if ( !all(rownames(feature_data) == rownames(parameter_data)) )
  {
    stop("Rownames in both datasets must match (sample names)")
  }
  
  force_numeric <- function(data) {
    is_num <- sapply(data, class)
    if (!all(is_num %in% c("numeric", "integer"))) {
      if (force.numeric) {
        data <- as.data.frame(lapply(data, as.numeric))
        is_num <- sapply(data, class)
        if (!all(is_num %in% c("numeric", "integer"))) {
          stop("Non-numeric data detected in your dataset. Please check your data.")
        }
        message("Forced non-numeric data to numeric as force.numeric was TRUE. Check your data.")
      } else {
        stop("Non-numeric data detected in your dataset. Check your data.")
      }
    }
    return(data)
  }
  feature_data <- force_numeric(feature_data)
  parameter_data <- force_numeric(parameter_data)
  
  if (is.binary.parameter) {
    corr.method <- "pearson"
    message("Using Point-Biserial correlations (Pearson's special case) for binary data.")
  } else {
    message(paste("Using ", corr.method, " correlations."))
  }
  
  pb <- txtProgressBar(min = 0, max = ncol(feature_data), style = 3, width = 50, char = "#")
  cortab <- matrix(NA, nrow = ncol(feature_data), ncol = ncol(parameter_data))
  ptab <- matrix(NA, nrow = ncol(feature_data), ncol = ncol(parameter_data))
  
  for ( taxn in 1:ncol(feature_data))
  {
    setTxtProgressBar(pb, taxn)
    for (paramn in 1:ncol(parameter_data))
    {
      distrib = feature_data[,taxn]
      cats = parameter_data[,paramn]
      if (sum(!is.na(distrib)) > 3 & sum(!is.na(cats)) > 3 & length(unique(cats)) > 1) {
        cortemp <- suppressWarnings(cor.test(distrib, cats, method = corr.method))
        cortab[taxn, paramn] <- cortemp$estimate
        ptab[taxn, paramn] <- cortemp$p.value
      }
    }
  }
  
  rownames(cortab) <- colnames(feature_data)
  colnames(cortab) <- colnames(parameter_data)
  rownames(ptab) <- colnames(feature_data)
  colnames(ptab) <- colnames(parameter_data)
  # Remove columns/rows with all NAs
  cortab <- cortab[, colSums(is.na(ptab)) < nrow(ptab)]
  ptab <- ptab[, colSums(is.na(ptab)) < nrow(ptab)]
  cortab <- cortab[rowSums(is.na(ptab)) < ncol(ptab), ]
  ptab <- ptab[rowSums(is.na(ptab)) < ncol(ptab), ]
  # Substitute NA and NaNs
  ptab[is.na(ptab)] <- 1
  cortab[is.na(cortab)] <- 0
  ptab[is.nan.data.frame(ptab)]  <-  1
  cortab[is.nan.data.frame(cortab)]  <-  0
  if ( all(dim(cortab) < c(2,2)) )
  {
    return( list(uncorr_correlation=as.data.frame(cortab), uncorr_pvalues = as.data.frame(ptab)) )
    stop("No valid correlations. All correlations returned NAs. Check your data.")
  }
  uncorr_cortab = cortab
  uncorr_ptab = ptab
  # Filter on uncorrected p-values
  nosig_tax <- rowSums(ptab < uncorr.p.threshold) == 0
  ptab <- ptab[!nosig_tax, , drop = FALSE]
  cortab <- cortab[!nosig_tax, , drop = FALSE]
  nosig_param <- colSums(ptab < uncorr.p.threshold) == 0
  ptab <- ptab[, !nosig_param, drop = FALSE]
  cortab <- cortab[, !nosig_param, drop = FALSE]
  if ( all(dim(cortab) < c(2,2)) )
  {
    return( list(uncorr_correlation=cortab, uncorr_pvalues = ptab) )
    stop("No significant correlations even without correcting pvals.")
  }
  # Filter on corrected p-values
  p2 = ptab
  c2 = cortab
  for (i in 1:ncol(p2)){ p2[,i] = p.adjust(p2[,i], method = p.adjust.method, n = nrow(ptab) )}
  nosig_tax <- rowSums(p2 < corr.p.threshold) == 0
  if (sum(!nosig_tax) == 0) {
    return( list(uncorr_correlation=as.data.frame(cortab), uncorr_pvalues = as.data.frame(ptab)) )
    stop("No significant correlation correcting pvalues")
  }
  p2 <- p2[!nosig_tax, , drop = FALSE]
  c2 <- c2[!nosig_tax, , drop = FALSE]
  
  nosig_param <- colSums(p2 < corr.p.threshold) == 0
  if (sum(!nosig_param) == 0) {
    return( list(uncorr_correlation=as.data.frame(cortab), uncorr_pvalues = as.data.frame(ptab)) )
    stop("No significant correlation correcting pvalues")
  }
  p2 <- p2[, !nosig_param, drop = FALSE]
  c2 <- c2[, !nosig_param, drop = FALSE]
  
  if (ncol(c2) < 1 || nrow(c2) < 1) {
    return(list(
      uncorr_correlation = cortab,
      uncorr_pvalues = ptab,
      corrected_correlation = c2,
      corrected_pvalues = p2
    ))
    message( sum(ptab < 0.05) , " correlations found with uncorrected p-values")
    message( sum(p2 < 0.05) , " correlations found with ", p.adjust.method, " adjusted p-values")
    stop()
  }
  
  # Keep only correlation stronger than abs(0.3)
  ctemp <- c2
  ctemp[p2 >= corr.p.threshold] <- NA
  nostrong_tax <- rowSums(abs(ctemp) >= min.corr.threshold, na.rm = TRUE) == 0
  if (sum(!nostrong_tax) == 0) {
    return(list(
      uncorr_correlation = cortab,
      uncorr_pvalues = ptab,
      corrected_correlation = c2,
      corrected_pvalues = p2
    ))
    stop("No correlation stronger than ", min.corr.threshold, " detected")
  }
  p3 <- p2[!nostrong_tax, , drop = FALSE]
  c3 <- c2[!nostrong_tax, , drop = FALSE]
  ctemp <- ctemp[!nostrong_tax, , drop = FALSE]
  
  nostrong_param <- colSums(abs(ctemp) >= min.corr.threshold, na.rm = TRUE) == 0
  if (sum(!nostrong_param) == 0) {
    return(list(
      uncorr_correlation = cortab,
      uncorr_pvalues = ptab,
      corrected_correlation = c2,
      corrected_pvalues = p2
    ))
    stop("No correlation stronger than ", min.corr.threshold, " detected")
  }
  p3 <- p3[, !nostrong_param, drop = FALSE]
  c3 <- c3[, !nostrong_param, drop = FALSE]
  
  cor_out = list(uncorr_correlation = as.data.frame(uncorr_cortab), uncorr_pvalues = as.data.frame(uncorr_ptab),
                 corrected_correlation = as.data.frame(c2), corrected_pvalues = as.data.frame(p2),
                 corr_strong_correlation = as.data.frame(c3), corr_strong_pvalues = as.data.frame(p3),
                 call = as.data.frame(rbind(
                   if(is.binary.parameter){paste("Point Biserial correlation")}else{paste(corr.method , "correlation")},
                   paste("Uncorrected p-values threshold", uncorr.p.threshold),
                   paste("Multiple comparison p.adjust method", p.adjust.method),
                   paste("Corrected p-values threshold", corr.p.threshold)
                 ))
  )
  
  message( "\n", sum(ptab < 0.05) , " correlations found with uncorrected p-values")
  message( sum(p2 < 0.05) , " correlations found with ", p.adjust.method, " adjusted p-values")
  message( sum(p3 < 0.05) , " correlations found with ", p.adjust.method, " adjusted p-values and at least ", min.corr.threshold, " cor" )
  message("Computing completed.")
  return(cor_out)
}