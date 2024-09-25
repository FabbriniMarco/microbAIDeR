process_pairwise_wilcoxon <- function(testing, onlysig = TRUE) {
  significant_comparisons <- data.frame(Comparison = character(), Wilcox.pvalue = numeric(), stringsAsFactors = FALSE)
  process_comparison <- function(line, column) {
    pval <- testing$p.value[line, column]
    if (!is.na(pval) && (!onlysig || pval <= 0.1)) {
      comparison <- paste(rownames(testing$p.value)[line], "vs", colnames(testing$p.value)[column])
      return(c(comparison, pval))
    }
    return(NULL)
  }
  for (line in seq_len(nrow(testing$p.value))) {
    for (column in seq_len(ncol(testing$p.value))) {
      result <- process_comparison(line, column)
      if (!is.null(result)) {
        significant_comparisons <- rbind(significant_comparisons, setNames(as.list(result), colnames(significant_comparisons)))
      }
    }
  }
  if (nrow(significant_comparisons) == 0) {
    return(NA)
  }
  significant_comparisons$Wilcox.pvalue = round(as.numeric(significant_comparisons$Wilcox.pvalue), 5)
  return(significant_comparisons)
}