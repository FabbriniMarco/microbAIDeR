compute_wilcoxon_and_plot <- function(data, group, taxlevel, save.path = getwd(), color.grouping, comparison.list = NULL, p.adjust.method = "fdr", trends = TRUE, plot.not.sig = FALSE,
                                      paired = FALSE, nrow.graph = 2, ncol.graph = 2, width.graph = 4.5, height.graph = 3.5, horiz = FALSE, ggplot.margins = c(.18, .18, .18, .6),
                                      box.lwd = 0.4, jitter.pch = 21, jitter.stroke = 0.15, jitter.size = 0.7, jitter.color = "grey22",
                                      signif.step.increase = 0.12, signif.text.size = 3, signif.line.size = 0.4, contrast.color = "ivory1",
                                      text.x.size = 6, text.y.size = 6, text.y.title.size = 8, smoothing = FALSE, smoothing.lwd = 1, smoothing.color = "darkred", smoothing.se = FALSE, 
                                      smoothing.method="loess", additional.params = NULL, align.legend = FALSE, plot.order = "kruskal", pattern.fill = FALSE, pattern = "stripe",
                                      pattern.angle = 45, pattern.alpha = 0.4, pattern.density = 0.1, pattern.spacing = 0.05) {
  check_and_load_package("dplyr")
  check_and_load_package("ggplot2")
  check_and_load_package("ggsignif")
  check_and_load_package("gridExtra")
  check_and_load_package("tidyverse")
  check_and_load_package("tidyr")
  check_and_load_package("openxlsx")
  if (pattern.fill) check_and_load_package("ggpattern")
  
  microbAIDeR::mkdir(save.path)
  if (any(is.na(group))) {
    data <- data[, !is.na(group)]
    group <- group[!is.na(group)]
  }
  data <- data[rowSums(data) != 0, ]
  
  if (ncol(data) == 0) stop("No samples remaining after filtering NAs from the grouping factor. Check your grouping factor")
  if (nrow(data) == 0) stop("No features remaining after filtering features with 0 abundance across all samples. Check your data")
  if (!p.adjust.method %in% p.adjust.methods) stop("Supplied an invalid value to p.adjust.method parameter.")
  if (!plot.order %in% c("kruskal", "rownames")) stop("Invalid value supplied to plot.order parameter. Can either be 'kruskal' or 'rownames'")
  
  call.print = as.data.frame(rbind(p.adjust.method, taxlevel, save.path, color.grouping = paste(color.grouping, collapse = ", "), smoothing.color, contrast.color,
                                   comparison.list = paste(sapply(comparison.list, function(x) paste(x, collapse = " vs ")), collapse=","),
                                   nrow.graph , ncol.graph , width.graph , height.graph , ggplot2.margins = paste(ggplot.margins, collapse = ", ") ,
                                   box.lwd , jitter.pch , jitter.stroke , jitter.size , jitter.color ,
                                   signif.step.increase,  signif.text.size ,
                                   text.x.size , text.y.size ,  text.y.title.size, additional.params=paste(additional.params, collapse=", ")
  ))
  output_excel_path = paste(save.path, "/", taxlevel, "_wilcoxon_pairwise.xlsx", sep="")
  write.xlsx( call.print , output_excel_path, sheetName = "call", colNames = FALSE, rowNames = TRUE)
  
  if (plot.order == "kruskal") {
    kw <- apply(data, 1, function(x) kruskal.test(x, group)$p.value)
    data <- data[order(ifelse(is.na(kw), 1, kw)), ]
  }
  gvec <- list()
  gvec_corr <- list()
  p_threshold <- if (trends) 0.1 else 0.05
  if (is.null(comparison.list)) comparison.list <- combn(levels(group), 2, simplify = FALSE)
  ptab <- data.frame(matrix(nrow = nrow(data), ncol = length(comparison.list) + 4 * nlevels(group)))
  rownames(ptab) <- rownames(data)
  colnames(ptab) <- c(
    sapply(comparison.list, function(x) paste(x, collapse = " vs ")),
    unlist(apply(expand.grid(c("Mean", "SEM", "Median", "SEMedian"), levels(group)), 1, function(x) paste(x[1], x[2], sep = " ")))
  )
  
  
  generate_plot <- function(tempdf, taxa, onlysig) {
    plot <- ggplot(tempdf, aes(x = grouping_column, y = abundances)) +
      { if (pattern.fill) geom_boxplot_pattern(aes(pattern_fill = grouping_column), pattern = pattern, outlier.shape = NA, lwd = box.lwd, pattern_angle = pattern.angle, pattern_alpha = pattern.alpha, pattern_density = pattern.density, pattern_spacing = pattern.spacing, pattern_colour = color.grouping) else geom_boxplot(aes(fill = grouping_column), outlier.shape = NA, lwd = box.lwd) } +
      geom_jitter(aes(color = grouping_column), width = 0.1, height = 0, pch = jitter.pch, color = jitter.color, stroke = jitter.stroke, size = jitter.size) +
      { if (!is.null(onlysig)) geom_signif(comparisons = strsplit(onlysig$V1, " vs "), annotations = sapply(onlysig$V2, sigFunction, trends = trends), size = signif.line.size, color = "grey22", textsize = signif.text.size, step_increase = signif.step.increase, tip_length = 0.02, margin_top = 0.05, vjust = 0.5) } +
      { if (smoothing) geom_smooth(aes(x = microbAIDeR::tokenize(grouping_column)), method = smoothing.method, lwd = smoothing.lwd, color = contrast.color, se = smoothing.se) } +
      scale_fill_manual(values = color.grouping) +
      labs(x = "", y = taxa) +
      { if (horiz) coord_flip() } +
      theme(legend.position = 'none', plot.margin = unit(ggplot.margins, "cm"), axis.text.x = element_text(size = text.x.size), axis.text.y = element_text(size = text.y.size), axis.title.y = element_text(size = text.y.title.size)) +
      { if (!is.null(additional.params)) additional.params }
    
    if (align.legend) plot <- align_legend(plot)
    return(plot)
  }
  
  
  for (taxa in rownames(data)) {
    tempdf <- data.frame(abundances = as.numeric(data[taxa, ]), grouping_column = group)
    wtests <- sapply(comparison.list, function(confr) suppressWarnings(wilcox.test(tempdf$abundances[tempdf$grouping_column == confr[1]], tempdf$abundances[tempdf$grouping_column == confr[2]], paired = paired)$p.value))
    wtests[is.na(wtests) | is.infinite(wtests)] <- 1
    
    if (any(wtests <= p_threshold)){
      onlysig = data.frame(V1 = unlist(lapply(comparison.list[wtests <= p_threshold], paste, collapse = " vs ")), V2 = round(wtests[wtests <= p_threshold], 5))
    } else { onlysig = NULL }
    
    ptab[taxa, 1:length(comparison.list)] <- signif(wtests, 3)
    for (lvl in levels(group)) {
      ptab[taxa, paste("Mean", lvl)] <- signif(mean(tempdf$abundances[tempdf$grouping_column == lvl], na.rm = TRUE), 3)
      ptab[taxa, paste("SEM", lvl)] <- signif(microbAIDeR::sem(tempdf$abundances[tempdf$grouping_column == lvl]), 3)
      ptab[taxa, paste("Median", lvl)] <- signif(median(tempdf$abundances[tempdf$grouping_column == lvl], na.rm = TRUE), 3)
      ptab[taxa, paste("SEMedian", lvl)] <- signif(microbAIDeR::se.median(tempdf$abundances[tempdf$grouping_column == lvl]), 3)
    }
    
    if( is.null(onlysig) )
    { if ( plot.not.sig ) { gvec <- c(gvec, list(generate_plot(tempdf, taxa, onlysig))) } else { next }
    } else { gvec <- c(gvec, list(generate_plot(tempdf, taxa, onlysig)))  }
  }
  file_uncorrected = paste0(save.path, "/", taxlevel, "_wilcox_uncorrected.pdf")
  message("Saving uncorrected boxplots to ", file_uncorrected)
  suppressWarnings(ggsave(
    filename = file_uncorrected,
    plot = marrangeGrob(gvec, nrow = nrow.graph, ncol = ncol.graph, top = NULL,
                        layout_matrix = matrix(seq_len(nrow.graph * ncol.graph), nrow = nrow.graph, ncol = ncol.graph, byrow = TRUE)),
    width = width.graph, height = height.graph, dpi = 330
  ))
  message("Adding uncorrected pvalues to ", output_excel_path)
  addSheet(path = output_excel_path, sheet.name = "uncorrected",addition = ptab,col.save = TRUE,row.save = TRUE)
  
  
  clean_ptab <- function(x, comparison.list, p_threshold, digits = 3) {
    for (col in seq_along(comparison.list)) {
      x[x[, col] > p_threshold, col] <- ""
      x[, col] <- ifelse(x[, col] != "", 
                         signif(as.numeric(x[, col]), digits = digits), 
                         "")
    }
    return(x)
  }
  
  
  cleaned_ptab = clean_ptab(ptab, comparison.list = comparison.list, p_threshold = p_threshold)
  addSheet( path= output_excel_path, sheet.name = "uncorrected_clean", addition = cleaned_ptab, col.save=TRUE, row.save=TRUE)
  
  # If pvalue correction is requested
  if (p.adjust.method != "none") {
    corrected_ptab <- ptab
    for (col in 1:length(comparison.list)) {
      corrected_ptab[, col] <- p.adjust(corrected_ptab[, col], method = p.adjust.method)
    }
    
    for (taxa in rownames(data)) {
      tempdf <- data.frame(abundances = as.numeric(data[taxa, ]), grouping_column = group)
      wtests_corrected <- as.numeric(corrected_ptab[taxa,1:length(comparison.list)])
      
      if ( any(wtests_corrected <= p_threshold) ){
        onlysig_corrected = data.frame(V1 = unlist(lapply(comparison.list[wtests_corrected <= p_threshold], paste, collapse = " vs ")), V2 = round(wtests_corrected[wtests_corrected <= p_threshold], 5))
      } else { onlysig_corrected = NULL }
      
      if( is.null(onlysig_corrected) )
      { if ( plot.not.sig ) { gvec_corr <- c(gvec_corr, list(generate_plot(tempdf, taxa, onlysig_corrected))) } else { next }
      } else { gvec_corr <- c(gvec_corr, list(generate_plot(tempdf, taxa, onlysig_corrected)))  }
      
    }
    file_corrected = paste0(save.path, "/", taxlevel, "_wilcox_", p.adjust.method,".pdf")
    message("Saving corrected boxplots to ", file_corrected)
    suppressWarnings(ggsave(
      filename = file_corrected,
      plot = marrangeGrob(gvec_corr, nrow = nrow.graph, ncol = ncol.graph, top = NULL,
                          layout_matrix = matrix(seq_len(nrow.graph * ncol.graph), nrow = nrow.graph, ncol = ncol.graph, byrow = TRUE)),
      width = width.graph, height = height.graph, dpi = 330
    ))
    message("Adding corrected pvalues to ", output_excel_path)
    addSheet(path = output_excel_path, sheet.name = p.adjust.method ,addition = corrected_ptab ,col.save = TRUE,row.save = TRUE)
    cleaned_corrected_ptab = clean_ptab(corrected_ptab, comparison.list = comparison.list, p_threshold = p_threshold)
    addSheet( path= output_excel_path, sheet.name = paste(p.adjust.method, "clean", sep="_"), addition = cleaned_corrected_ptab, col.save=TRUE, row.save=TRUE)
  }
  message("Computing completed.")
}
