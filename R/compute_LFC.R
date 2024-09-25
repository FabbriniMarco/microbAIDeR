compute_LFC <- function(data, group, w1=0.5, w2=0.5){
  if ( !is.data.frame(data) ){stop("Object provided to data is not class data.frame")}
  if(class(group) != "factor" ){ stop("The supplied group is not a factor.")}
  if ( any(is.na(group))){
    data = data[,!is.na(group)]
    group = group[!is.na(group)]
  }
  data = data[!rowSums(data) %in% 0,]
  if( ncol(data) == 0 ){ stop("No samples remaining after filtering NAs from the grouping factor. Check your grouping factor") }
  if( nrow(data) == 0 ){ stop("No features remaining after filtering features with 0 abundance across all samples. Check your data") }
  if( nlevels(group) != 2 ){ stop("The grouping factor must have exactly two groups.")}
  cat("\nFormula for Log-Fold Change (LFC):\n\n")
  cat(" LFC = log ( ( w1 * median(group1) + (1 - w1) * mean(group1) )\n")
  cat("               --------------------------------------------    )\n")
  cat("             ( w2 * median(group2) + (1 - w2) * mean(group2) )\n\n")
  group1 <- levels(group)[1]
  group2 <- levels(group)[2]
  LFC_result <- data.frame(taxon = rownames(data), LFC = rep(NA, nrow(data)))
  for (i in 1:nrow(data)) {
    values_group1 <- as.numeric(data[i, group == group1])
    M1 <- median(values_group1, na.rm = TRUE)
    m1 <- mean(values_group1, na.rm = TRUE)
    values_group2 <- as.numeric(data[i, group == group2])
    M2 <- median(values_group2, na.rm = TRUE)
    m2 <- mean(values_group2, na.rm = TRUE)
    group1_stat <- w1 * M1 + (1 - w1) * m1
    group2_stat <- w2 * M2 + (1 - w2) * m2
    # Ensure no division by zero
    if (group1_stat == 0 | group2_stat == 0) {
      LFC_result$LFC[i] <- NA
    } else {
      LFC_result$LFC[i] <- log(group1_stat / group2_stat)
    }
  }
  
  return(LFC_result)
}