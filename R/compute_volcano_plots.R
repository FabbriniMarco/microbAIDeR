compute_volcano_plots <- function(data, group, taxlevel, save.path = getwd(), p.adjust.method = "fdr", w1 = 0.5, w2 = 0.5, trends=TRUE, sigcolors=c("ivory3", "orange3", "darkred"),
                                  nrow.graph = 1, ncol.graph = 2, width.graph = 20, height.graph = 8, ggplot.margins = c(1,1,1,.5),
                                  label.text.size = 2.8, label.padding=4,
                                  arrows=TRUE, arrow.length=1.5,arrow.lwd=1.8,arrow.origin.offset=0.5,
                                  arrow.text.cex=.7,arrow.text.color="grey22",
                                  text.x.size = 8, text.y.size = 7,  text.y.title.size = 9,
                                  additional.params = NULL, align.legend = FALSE){
  check_and_load_package("dplyr")
  check_and_load_package("ggplot2")
  check_and_load_package("gridExtra")
  check_and_load_package("tidyverse")
  check_and_load_package("tidyr")
  check_and_load_package("openxlsx")
  check_and_load_package("grid")
  check_and_load_package("ggrepel")
  check_and_load_package("cowplot")
  
  output_xlsx = paste(save.path, "/", taxlevel, "_volcano.xlsx", sep="")
  output_pdf = paste(save.path, "/", taxlevel, "_volcano.pdf", sep="")
  
  microbAIDeR::mkdir(save.path)
  
  if ( any(is.na(group))){
    data = data[,!is.na(group)]
    group = group[!is.na(group)]
  }
  
  data = data[!rowSums(data) %in% 0,]
  
  if( ncol(data) == 0 ){ stop("No samples remaining after filtering NAs from the grouping factor. Check your grouping factor") }
  if( nrow(data) == 0 ){ stop("No features remaining after filtering features with 0 abundance across all samples. Check your data") }
  if(!p.adjust.method %in% p.adjust.methods){ stop("Supplied an invalid value to p.adjust.method parameter.")}
  if(class(group) != "factor" ){ stop("The supplied group is not a factor.")}
  if( nlevels(group) != 2 ){ stop("The grouping factor must have exactly two groups.")}
  
  wt <- rep(NA, nrow(data))
  for( i in 1:nrow(data))
  {
    wt[i] = suppressWarnings(wilcox.test( as.numeric(data[i,]) ~ group )$p.value)
  }
  wt[is.nan(wt)] = 1
  wt[is.na(wt)] = 1
  data = data[order(wt),]
  
  if( isTRUE(trends) )
  {
    p_threshold <- 0.1
  } else if ( isFALSE(trends) ) {
    p_threshold <- 0.05
  } else {
    stop("Wrong value supplied to the 'trends' argument")
  }
  
  gvec <- vector("list",length=0 )
  iter1 = 1
  
  if(length(sigcolors) != 3)
  {
    message("WARNING: You have not supplied 3 colors to sigcolors. Reverting to default values")
    sigcolors=c("ivory3", "orange3", "darkred")
  }
  
  call.print = as.data.frame(rbind(p.adjust.method, taxlevel, save.path,
                                   sigcolors = paste(sigcolors, collapse = ", "),
                                   nrow.graph , ncol.graph , width.graph , height.graph , ggplot2.margins = paste(ggplot.margins, collapse = ", ") ,
                                   label.text.size, label.padding,arrow.length,arrow.lwd,arrow.origin.offset,
                                   arrow.text.cex,arrow.text.color,
                                   text.x.size , text.y.size ,  text.y.title.size, additional.params=paste(additional.params, collapse=", ")
  ))
  openxlsx::write.xlsx( call.print , output_xlsx, sheetName = "call", colNames = FALSE, rowNames = TRUE)
  
  ptab = as.data.frame(
    matrix(
      nrow = nrow(data),
      ncol = 9
    )
  )
  rownames(ptab) = rownames(data)
  colnames(ptab) = c(
    paste(levels(group)[1], levels(group)[2], sep=" vs "),
    unlist(apply(expand.grid(c("Mean", "SEM", "Median", "SEMedian"), levels(group)), 1, function(x) paste(x[1], x[2], sep = " ")))
  )
  # Get pvalues into ptab and save
  message("Saving stats into ", output_xlsx)
  ptab[,1] = wt
  ptab[,2] = signif(as.numeric(rowMeans(data[,group == levels(group)[1]])), 2)
  ptab[,3] = signif(as.numeric(sapply( as.data.frame(t(data[,group == levels(group)[1]])) , sem)), 2)
  ptab[,4] = signif(as.numeric(rowMedians(data[,group == levels(group)[1]])), 2)
  ptab[,5] = signif(as.numeric(sapply( as.data.frame(t(data[,group == levels(group)[1]])) , se.median)), 2)
  ptab[,6] = signif(as.numeric(rowMeans(data[,group == levels(group)[2]])), 2)
  ptab[,7] = signif(as.numeric(sapply( as.data.frame(t(data[,group == levels(group)[2]])) , sem)), 2)
  ptab[,8] = signif(as.numeric(rowMedians(data[,group == levels(group)[2]])), 2)
  ptab[,9] = signif(as.numeric(sapply( as.data.frame(t(data[,group == levels(group)[2]])) , se.median)), 2)
  
  addSheet(path = output_xlsx, sheet.name = "uncorrected",addition = ptab, col.save = TRUE, row.save = TRUE)
  if( p.adjust.method != "none")
  {
    ptab_corr <- ptab
    ptab_corr[,1] = p.adjust(ptab_corr[,1], method = p.adjust.method)
    addSheet(path =output_xlsx, sheet.name = p.adjust.method,addition = ptab_corr, col.save = TRUE, row.save = TRUE)
  }
  
  # Plot volcanos for uncorrected P-values
  message("Terraforming volcano plots..")
  tidata = as.data.frame(
    matrix(
      nrow=nrow(data),
      ncol = 7,
      dimnames = list(rownames(data), c("pval", "corr", "LFC", "Median_Relabb", "Mean_Relabb", "Color_uncorrected", "Color_corr"))))
  for (i in 1:nrow(tidata)){tidata[i,1] <- wt[i]}
  tidata[,2] <- p.adjust(tidata[,1], method = p.adjust.method)
  
  LFC_tab <- compute_LFC(data = data, group = group, w1 = w1, w2 = w2)
  LFC_tab = LFC_tab[!is.na(LFC_tab$LFC),]
  tidata = tidata[LFC_tab$taxon,]
  tidata$LFC = LFC_tab$LFC
  for(i in 1:nrow(tidata))
  {
    if ( tidata$LFC[i] < 0 )
    {
      tidata$Median_Relabb[i] = median( as.numeric(data[rownames(tidata)[i], group == levels(group)[1]]), na.rm = TRUE )
      tidata$Mean_Relabb[i] = mean( as.numeric(data[rownames(tidata)[i], group == levels(group)[1]]), na.rm = TRUE )
    } else {
      tidata$Median_Relabb[i] = median( as.numeric(data[rownames(tidata)[i], group == levels(group)[2]]), na.rm = TRUE )
      tidata$Mean_Relabb[i] = mean( as.numeric(data[rownames(tidata)[i], group == levels(group)[1]]), na.rm = TRUE )
    }
    if( tidata$pval[i] <= 0.05 ){
      tidata$Color_uncorrected[i] = sigcolors[3]
    } else if ( isTRUE(trends) & tidata$pval[i] <= 0.1){
      tidata$Color_uncorrected[i] = sigcolors[2]
    } else if ( isFALSE(trends) & tidata$pval[i] <= 0.1){
      tidata$Color_uncorrected[i] = sigcolors[1]
    } else { tidata$Color_uncorrected[i] = sigcolors[1] }
    # Corrected ones
    if( tidata$corr[i] <= 0.05 ){
      tidata$Color_corr[i] = sigcolors[3]
    } else if ( isTRUE(trends) & tidata$corr[i] <= 0.1){
      tidata$Color_corr[i] = sigcolors[2]
    } else if ( isFALSE(trends) & tidata$corr[i] <= 0.1){
      tidata$Color_corr[i] = sigcolors[1]
    } else { tidata$Color_corr[i] = sigcolors[1] }
  }
  tidata$Plotlabel_uncorrected = rownames(tidata)
  tidata$Plotlabel_uncorrected[tidata$pval > p_threshold]=NA
  tidata$Plotlabel_corr = rownames(tidata)
  tidata$Plotlabel_corr[tidata$corr > p_threshold]=NA
  
  
  create_volcano <- function(tidata, x_var, y_var, color_var, label_var, title_y, size_var, additional.params = NULL, arrows = FALSE) {
    if( size_var == "Median_Relabb" ) text_label = "Median Rel. Ab." else text_label = "Avg. Rel. Ab."
    
    graphy <- ggplot(data = tidata, aes_string(x = x_var, y = -log2(tidata[[y_var]]), label = label_var, size = size_var, color = color_var)) +
      geom_point(alpha = 0.8) +
      geom_vline(xintercept = 0, col = "grey80") +
      geom_hline(yintercept = -log2(0.05), col = "grey20", linetype = "dashed") +
      scale_color_identity() +
      labs(x = "LFC", size = text_label, y = title_y) +
      scale_y_continuous(limits = c(0, max(-log2(tidata[[y_var]])) + 0.5)) +
      scale_x_continuous(limits = c(min(tidata$LFC) - 0.5, max(tidata$LFC) + 0.5)) +
      geom_text_repel(show.legend = FALSE, size = label.text.size, point.padding = label.padding, aes(fontface = "italic")) +
      theme(legend.key = element_rect(fill = NA), plot.margin = unit(ggplot.margins, "cm"),
            axis.text.x = element_text(size = text.x.size),
            axis.text.y = element_text(size = text.y.size), axis.title.y = element_text(size = text.y.title.size))
    
    if (arrows) {
      graphy <- graphy +
        annotation_custom(grob = linesGrob(arrow = arrow(type = "open", ends = "first", length = unit(arrow.length, "mm")),
                                           gp = gpar(col = "grey22", lwd = arrow.lwd)), xmax = -arrow.origin.offset, xmin = min(tidata$LFC) / 1.5, ymin = 0, ymax = 0) +
        annotation_custom(grob = linesGrob(arrow = arrow(type = "open", ends = "last", length = unit(arrow.length, "mm")),
                                           gp = gpar(col = "grey22", lwd = arrow.lwd)), xmax = arrow.origin.offset, xmin = max(tidata$LFC) / 1.5, ymin = 0, ymax = 0) +
        annotation_custom(grob = grid::textGrob(label = levels(group)[2], hjust = 1.05, gp = gpar(col = arrow.text.color, cex = arrow.text.cex)),
                          xmin = min(tidata$LFC) / 1.5, xmax = min(tidata$LFC) / 1.5, ymin = 0, ymax = 0) +
        annotation_custom(grob = grid::textGrob(label = levels(group)[1], hjust = -0.05, gp = gpar(col = arrow.text.color, cex = arrow.text.cex)),
                          xmin = max(tidata$LFC) / 1.5, xmax = max(tidata$LFC) / 1.5, ymin = 0, ymax = 0)
    }
    
    if (!is.null(additional.params)) {
      graphy <- graphy + additional.params
    }
    if ( align.legend ){ graphy <- cowplot::ggdraw(align_legend(graphy)) }
    graphy
  }
  
  
  add_plot_to_gvec <- function(graphy, gvec) {
    gvec = c(gvec, list(graphy))
    return(gvec)
  }
  
  
  # Main plotting
  graphy <- create_volcano(tidata, "LFC", "pval", "Color_uncorrected", "Plotlabel_uncorrected", "-log2(p)", "Median_Relabb", additional.params, arrows)
  gvec <- add_plot_to_gvec(graphy, gvec)
  if (p.adjust.method != "none") {
    graphy <- create_volcano(tidata, "LFC", "corr", "Color_corr", "Plotlabel_corr", paste("-log2(", p.adjust.method, ")", sep = ""), "Median_Relabb", additional.params, arrows)
    gvec <- add_plot_to_gvec(graphy, gvec)
  }
  
  graphy <- create_volcano(tidata, "LFC", "pval", "Color_uncorrected", "Plotlabel_uncorrected", "-log2(p)", "Mean_Relabb", additional.params, arrows)
  gvec <- add_plot_to_gvec(graphy, gvec)
  if (p.adjust.method != "none") {
    graphy <- create_volcano(tidata, "LFC", "corr", "Color_corr", "Plotlabel_corr", paste("-log2(", p.adjust.method, ")", sep = ""), "Mean_Relabb", additional.params, arrows)
    gvec <- add_plot_to_gvec(graphy, gvec)
  }
  
  # Save plots and excels
  message("Saving plots into ", output_pdf)
  suppressWarnings(ggsave(
    filename = output_pdf,
    plot = marrangeGrob(gvec, nrow=nrow.graph, ncol=ncol.graph, top=NULL ,
                        layout_matrix = matrix(seq_len(nrow.graph*ncol.graph),
                                               nrow = nrow.graph,
                                               ncol = ncol.graph , byrow = TRUE) ),
    width = width.graph, height = height.graph, dpi = 330))
}