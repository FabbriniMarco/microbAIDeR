compute_beta_diversity <- function(beta_metrics = c("braycurtis", "jaccard", "unweighted_unifrac", "weighted_unifrac"),
                                   save.path = getwd(), adonis_n_perm = 9999, beta.folder.path = paste(getwd(), "RESULTS/beta_diversity", sep="/"), manual.beta.path = NULL, sample.list, mds = c(1,2), group, color.grouping ,
                                   alpha.points = 1, spiders = FALSE, spiders.lwd = 1, ellipses = FALSE, ellipse.focus = FALSE,ellipse.fill = FALSE, ellipse.conf = 0.95, ellipse.alpha = 0.45 , ellipse.lwd = 1, manual.bordercol = NULL,
                                   svg.width = 5, svg.height = 3,
                                   cex.points = 1, adonis.p.adjust = "fdr", wilcoxon.p.adjust = "fdr",
                                   nrow.graph = 2, ncol.graph = 2, width.graph = 6, height.graph = 4, ggplot.margins = c(0.15, 0.15, 0.15, 0.6),
                                   text.x.size = 8, text.y.size = 7,  text.title.size = 9,
                                   additional.params=NULL, use.ggplot = TRUE, additional.params.beta=NULL){
  check_and_load_package("openxlsx")
  check_and_load_package("dplyr")
  check_and_load_package("tidyr")
  check_and_load_package("ggplot2")
  check_and_load_package("vegan")
  check_and_load_package("gridExtra")
  check_and_load_package("pairwiseAdonis")
  check_and_load_package("stringr")
  check_and_load_package("svglite")
  
  if ( any(is.na(group))){ message("WARNING: You have NAs in your group, removing such samples from the beta distance matrices. BE SURE not to have NAs in the sample.list") }
  if ( length(group) != length(sample.list) ){ stop("ERROR: The grouping factor provided has a different length compared to the sample.list") }
  if ( !is.factor(group) ){ stop("Input group provided is not a factor") }
  
  microbAIDeR::mkdir(save.path)
  call.print = as.data.frame(rbind( paste(beta_metrics, collapse=", "), color.grouping = paste(color.grouping, collapse = ", "),
                                    save.path, adonis_n_perm , beta.folder.path ,  mds = paste(mds, collapse=", ") ,
                                    alpha.points, spiders, ellipses , ellipse.focus , ellipse.conf , ellipse.alpha , 
                                    svg.width , svg.height ,
                                    cex.points , adonis.p.adjust, wilcoxon.p.adjust, additional.params=paste(additional.params, collapse=", ")
  ))
  adonis_excel = paste(save.path, "/beta_adonis.xlsx", sep="")
  variance_excel = paste(save.path, "/beta_variance.xlsx", sep="")
  openxlsx::write.xlsx( call.print , adonis_excel, sheetName = "call", colNames = FALSE, rowNames = TRUE)
  openxlsx::write.xlsx( call.print , variance_excel, sheetName = "call", colNames = FALSE, rowNames = TRUE)
  
  iter1 = 1
  gvec <- vector("list",length=0 )
  for ( metrics in beta_metrics)
  {
    message("Running ", toupper(metrics), " beta diversity analyses")
    if (is.null(manual.beta.path)) {
      beta = read.delim( paste(beta.folder.path , metrics, "distance-matrix.tsv", sep="/"), header=T, row.names=1, sep="\t")
    } else if( !is.null(manual.beta.path)) {
      beta = read.delim( manual.beta.path, header=T, row.names=1, sep="\t")
    }
    if (!all(sample.list %in% rownames(beta)) | !all(sample.list %in% colnames(beta))){ stop("The provided sample list and beta diversity matrices sample do not correspond") }
    beta = beta[ sample.list, sample.list ]
    if ( any(is.na(group)) ){
      beta = beta[!is.na(group), !is.na(group)]
      group_beta = group[!is.na(group)]
    } else { group_beta = group }
    metric_beta <- capscale(as.dist(beta)~1)
    coord_beta = as.data.frame( scores(metric_beta, display="sites", choices=mds ) )
    eighenvalues = metric_beta$CA$eig
    sig1 = round(as.numeric( eighenvalues[mds[1]]*100/sum(eighenvalues) ), 2)
    sig2 = round(as.numeric( eighenvalues[mds[2]]*100/sum(eighenvalues) ), 2)
    color_beta = rep(NA, length(group_beta))
    for (x in 1:length(group_beta)){color_beta[x] = color.grouping[which(levels(group_beta) == as.character(group_beta)[x] )]}
    if(any(is.na(color_beta))){
      message("Errors in producing the color vector. Check that your group is a factor and nlevels(group) == length(color.grouping)")
      message("WARNING:Procededing anyway without colors")
      color_beta[is.na(color_beta)] = "black"
    }
    colnames(coord_beta) = c("AX1", "AX2")
    if( !is.numeric(alpha.points) ){alpha.points = 1}
    if( alpha.points > 1 ){ alpha.points = 1}
    if( alpha.points < 1 ){ color_beta = scales::alpha(color_beta, alpha.points) }
    
    if( !use.ggplot ){
      plot_beta <- function(save.path, metrics, plot_type, metric_beta, coord_beta, color_beta, group_beta, mds, sig1, sig2,
                            color.grouping, cex.points, spiders, ellipses, ellipse.focus, 
                            ellipse.fill, svg.width, svg.height, spiders.lwd, ellipse.conf, 
                            ellipse.lwd, ellipse.alpha, manual.bordercol) {
        file_name <- switch(plot_type,
                            "base" = paste0(save.path, "/", metrics, ".svg"),
                            "spiders" = paste0(save.path, "/", metrics, "_spiders.svg"),
                            "ellipse" = paste0(save.path, "/", metrics, "_ellipse.svg"),
                            "ellipse_focus" = paste0(save.path, "/", metrics, "_ellipse_focus.svg"),
                            "ellipse_fill" = paste0(save.path, "/", metrics, "_ellipse_fill.svg"))
        svglite(file_name, width = svg.width, height = svg.height)
        plot(metric_beta, type = "n", main = microbAIDeR::firstup(metrics), choices = mds,
             xlab = paste("MDS", mds[1], " - [", sig1, "%]", sep = ""),
             ylab = paste("MDS", mds[2], " - [", sig2, "%]", sep = ""))
        with(coord_beta, points(AX1, AX2, pch = 21, bg = color_beta, cex = cex.points, col = NULL))
        if (plot_type == "spiders") {
          for (groupiter in levels(group_beta)) {
            ordispider(metric_beta, group_beta, display = "sites", spiders = "centroid", 
                       show.groups = groupiter, 
                       col = if (is.null(manual.bordercol)) color.grouping[which(levels(group_beta) == groupiter)] else manual.bordercol, 
                       lwd = spiders.lwd)
          }
        } else if (plot_type == "ellipse" || plot_type == "ellipse_focus" || plot_type == "ellipse_fill") {
          for (groupiter in levels(group_beta)) {
            ordiellipse(metric_beta, group_beta, display = "sites", kind = "se", conf = ellipse.conf, 
                        draw = if (plot_type == "ellipse") "lines" else "polygon", 
                        show.groups = groupiter, 
                        col = color.grouping[which(levels(group_beta) == groupiter)], 
                        border = if (is.null(manual.bordercol)) color.grouping[which(levels(group_beta) == groupiter)] else manual.bordercol, 
                        lwd = ellipse.lwd, alpha = if (plot_type == "ellipse_fill" || plot_type == "ellipse_focus") ellipse.alpha else 0.50)
          }
        }
        legend("topright", legend = levels(group_beta), pt.cex = 1, fill = color.grouping, horiz = FALSE,
               x.intersp = 0.4, y.intersp = 0.7, bty = "n", cex = 0.8)
        graphics.off()
      }
      # Base plot
      plot_beta(save.path, metrics, "base", metric_beta, coord_beta, color_beta, group_beta, mds, sig1, sig2, color.grouping, cex.points, )
      # Spider plot
      if (spiders) {
        plot_beta(save.path, metrics, "spiders", metric_beta, coord_beta, color_beta, group_beta, mds, sig1, sig2, color.grouping, cex.points, spiders.lwd = spiders.lwd, spiders = TRUE)
      }
      # Ellipse plot
      if (ellipses) {
        plot_beta(save.path, metrics, "ellipse", metric_beta, coord_beta, color_beta, group_beta, mds, sig1, sig2, color.grouping, cex.points, ellipse.conf = ellipse.conf, ellipse.lwd = ellipse.lwd)
      }
      # Ellipse focus plot
      if (ellipse.focus) {
        plot_beta(save.path, metrics, "ellipse_focus", metric_beta, coord_beta, color_beta, group_beta, mds, sig1, sig2, color.grouping, cex.points/4, ellipse.conf = ellipse.conf, ellipse.lwd = ellipse.lwd, ellipse.alpha = ellipse.alpha)
      }
      # Ellipse fill plot
      if (ellipse.fill) {
        plot_beta(save.path, metrics, "ellipse_fill", metric_beta, coord_beta, color_beta, group_beta, mds, sig1, sig2, color.grouping, cex.points, ellipse.conf = ellipse.conf, ellipse.lwd = ellipse.lwd, ellipse.alpha = ellipse.alpha)
      }
    } else {
      # Plot the PCoA using ggplot with ellipses (if requested)
      beta.pcoa <- cmdscale(as.dist(beta), eig=TRUE, k=max(mds))
      eigenvalues <- beta.pcoa$eig
      explained_variance <- eigenvalues / sum(eigenvalues)
      explained_variance_PC1 <- explained_variance[1]
      explained_variance_PC2 <- explained_variance[2]
      pcoa_df <- data.frame(PC1 = beta.pcoa$points[,mds[1]], PC2 = beta.pcoa$points[,mds[2]], Group = group)
      ellipse_data <- pcoa_df %>% group_by(Group) %>% do(calc_ellipse(., level=ellipse.conf))
      if(ellipse.lwd > 2){message("WARNING: Setting an ellipse.lwd value over 2 with ggplot results in very thick ellipse borders")}
      
      graphy <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
        geom_point(pch = 21, aes(fill = Group), color = "grey10", stroke = 0.25, size = cex.points) +
        labs(x = paste("MDS1 - ", round(explained_variance_PC1*100, 2), "%", sep=""),
             y = paste("MDS2 - ", round(explained_variance_PC2*100, 2), "%", sep=""),
             fill = "", color = "", title = firstup(metrics)) +
        scale_color_manual(values = color.grouping) +
        scale_fill_manual(values = color.grouping) +
        coord_fixed() +
        theme_minimal() +
        {if (!is.null(additional.params.beta)) additional.params.beta}
      svglite(paste0(save.path, "/", metrics, "_ggplot.svg") , width = svg.width, height = svg.height)
      print(graphy)
      graphics.off()
      
      if (ellipses) {
        graphy <- graphy +
          geom_polygon(data = ellipse_data, aes(x = PC1, y = PC2, group = Group), fill = NA, color = "grey10", linewidth = ellipse.lwd) +
          geom_polygon(data = ellipse_data, aes(x = PC1, y = PC2, fill = Group, color = Group), alpha = ellipse.alpha)
        svglite(paste0(save.path, "/", metrics, "_ggplot_ellipses.svg") , width = svg.width, height = svg.height)
        print(graphy)
        graphics.off()
      }
      
      if (spiders) {
        centroids <- pcoa_df %>%
          group_by(Group) %>%
          summarize(PC1_centroid = mean(PC1), PC2_centroid = mean(PC2))
        pcoa_df <- pcoa_df %>%
          left_join(centroids, by = "Group")
        
        graphy <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
          geom_segment(aes(xend = PC1_centroid, yend = PC2_centroid, color = Group), size = 0.5) +
          geom_point(pch = 21, aes(fill = Group), color = "grey10", stroke = 0.25, size = cex.points) +
          labs(x = paste("MDS1 - ", round(explained_variance_PC1*100, 2), "%", sep=""),
               y = paste("MDS2 - ", round(explained_variance_PC2*100, 2), "%", sep=""),
               fill = "", color = "", title = firstup(metrics)) +
          scale_color_manual(values = color.grouping) +
          scale_fill_manual(values = color.grouping) +
          coord_fixed() +
          theme_minimal() +
          {if (!is.null(additional.params.beta)) additional.params.beta}
        
        svglite(paste0(save.path, "/", metrics, "_ggplot_spiders.svg") , width = svg.width, height = svg.height)
        print(graphy)
        graphics.off()
      }
      
    }
    # Adonis
    message("Computing pairwise Adonis for: ", metrics)
    adone = pairwise.adonis(beta, group_beta, perm = adonis_n_perm, p.adjust.m = adonis.p.adjust)
    trends = TRUE #Display trends in pairwise adonis p-values table
    adone$uncorrected_pvalues <- sapply(adone$p.value, sigFunction, trends = trends)
    adone[,adonis.p.adjust] <- sapply(adone$p.adjusted, sigFunction, trends = trends)
    addSheet( path=adonis_excel, sheet.name = metrics, addition = adone, col.save=TRUE, row.save=FALSE)
    # Calculate and plot the intra-group variances
    compute_intra_group_variance <- function(beta, group_beta, color_grouping, save.path, metrics, text.x.size, text.y.size, text.title.size, ggplot.margins, additional.params = NULL, wilcoxon.p.adjust = "fdr") {
      message("Computing intra-group variance")
      beta_half <- beta
      beta_half[lower.tri(beta_half, diag = TRUE)] <- NA
      beta_half <- as.data.frame(beta_half)
      values_beta_var <- c()
      group_beta_var <- c()
      for (z in levels(group_beta))
      {
        if ( sum(group_beta==z) > 1 ){
          tempvalues = as.numeric(unlist(beta_half[ group_beta==z, group_beta==z]))
          tempvalues = tempvalues[!is.na(tempvalues)]
          tempvalues = tempvalues[!is.nan(tempvalues)]
        } else { tempvalues = 0 }
        values_beta_var = c(values_beta_var, tempvalues )
        group_beta_var = c(group_beta_var, rep(z, length(tempvalues)) )
      }
      plot_beta_variance <- data.frame(
        values_beta_var = as.numeric(values_beta_var),
        group_beta_var = factor(group_beta_var, levels = levels(group_beta))
      )
      # Group variance calculation
      var_df <- plot_beta_variance %>%
        group_by(group_beta_var) %>%
        summarise(Variance = var(values_beta_var))
      theme_params = theme(
        legend.position = "none",
        plot.margin = unit(ggplot.margins, "cm"),
        axis.text.x = element_text(size = text.x.size),
        axis.text.y = element_text(size = text.y.size),
        axis.title.y = element_text(size = text.title.size)
      )
      # Plot variance bar chart
      gvec <- list()
      graphy <- ggplot(var_df, aes(x = group_beta_var, y = Variance, fill = group_beta_var)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = color.grouping) +
        labs(x = "", y = "Intra-Group Variance", title = firstup(metrics)) +
        theme_params +
        {if (!is.null(additional.params)) additional.params}
      gvec[[1]] <- graphy
      # Compute pairwise distances, ranges, and IQR
      pairwise_distances <- lapply(levels(group_beta), function(gg) {
        group_data <- plot_beta_variance[plot_beta_variance$group_beta_var == gg, "values_beta_var"]
        as.vector(dist(group_data))
      })
      ranges <- sapply(pairwise_distances, function(distances) diff(range(distances)))
      iqr <- sapply(pairwise_distances, IQR)
      # Combine range and IQR results
      comparison_result <- data.frame(Range = ranges, IQR = iqr, row.names = levels(group_beta))
      comparison_result <- tidify(comparison_result, colNames=c("Group", "Measure", "Value"))
      comparison_result$Group = factor(comparison_result$Group, levels = levels(group))
      # Plot IQR
      graphy_iqr <- ggplot(comparison_result[comparison_result$Measure == "IQR", ], aes(x = Group, y = Value, fill = Group)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = color.grouping) +
        labs(title = firstup(metrics), x = "", y = "IQR") +
        theme_params +
        {if (!is.null(additional.params)) additional.params}
      gvec[[2]] <- graphy_iqr
      # Plot Range
      graphy_range <- ggplot(comparison_result[comparison_result$Measure == "Range", ], aes(x = Group, y = Value, fill = Group)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = color.grouping) +
        labs(title = firstup(metrics), x = "", y = "Range") +
        theme_params +
        {if (!is.null(additional.params)) additional.params}
      gvec[[3]] <- graphy_range
      # Kruskal test and pairwise Wilcoxon tests
      ktest <- kruskal.test(plot_beta_variance$values_beta_var, plot_beta_variance$group_beta_var)$p.value
      wtest_raw <- process_pairwise_wilcoxon(
        suppressWarnings(pairwise.wilcox.test(plot_beta_variance$values_beta_var, plot_beta_variance$group_beta_var, p.adjust.method = "none")),
        onlysig = FALSE
      )
      wtest_corrected <- process_pairwise_wilcoxon(
        suppressWarnings(pairwise.wilcox.test(plot_beta_variance$values_beta_var, plot_beta_variance$group_beta_var, p.adjust.method = wilcoxon.p.adjust)),
        onlysig = FALSE
      )
      comparisondf <- data.frame(
        Comparison = wtest_raw$Comparison,
        uncorrected_p = signif(wtest_raw$Wilcox.pvalue, 3),
        wilcox_p_adjust = signif(wtest_corrected$Wilcox.pvalue, 3),
        Kruskal_test_p = rep("", nrow(wtest_raw))
      )
      comparisondf[1, "Kruskal_test_p"] <- signif(ktest, 3)
      var_df = var_df %>%
        dplyr::rename(Group = group_beta_var) %>%
        tidyr::pivot_longer(cols = Variance, names_to = "Measure", values_to = "Value")
      var_out = bind_rows(var_df, comparison_result)
      var_out$Value = signif(var_out$Value, 3)
      return(list(gvec = gvec, comparisondf = comparisondf, var_df = var_out ) )
    }
    result <- compute_intra_group_variance(beta, group_beta, color_grouping, save.path, metrics, text.x.size, text.y.size, text.title.size, ggplot.margins, additional.params, wilcoxon.p.adjust)
    # Export the data in the Excel files
    addSheet( path=variance_excel, sheet.name = paste( gsub("weighted", "w", metrics), "Var", sep="-" ), addition = result$var_df, col.save=TRUE, row.save=FALSE)
    addSheet( path=variance_excel, sheet.name = paste( gsub("weighted", "w", metrics), "p", sep="-" ), addition = result$comparisondf, col.save=TRUE, row.save=FALSE)
    
    message("Saving plots to ", save.path, "/", metrics, "_intragroup_variance.pdf", sep="")
    suppressWarnings(ggsave(
      filename = paste(save.path, "/", metrics, "_intragroup_variance.pdf", sep="") ,
      plot = marrangeGrob(result$gvec, nrow=nrow.graph, ncol=ncol.graph, top=NULL ,
                          layout_matrix = matrix(seq_len(nrow.graph*ncol.graph), nrow = nrow.graph, ncol = ncol.graph , byrow = TRUE) ),
      width = width.graph, height = height.graph, dpi = 330
    ))
  }
  message("Computing completed.")
}