microbAIDeR_install_dependancies <- function(){
   if(!require(devtools)){
      install.packages('devtools')
   }
   if(!require(pairwiseAdonis)){
      devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
   }
   if(!require(parallel)){
      install.packages("parallel")
   }
   if(!require(doParallel)){
      install.packages("doParallel")
   }
   if(!require(openxlsx)){
      install.packages("openxlsx")
   }
   if(!require(tidyr)){
      install.packages("tidyr")
   }
   if(!require(tibble)){
      install.packages("tibble")
   }
   if(!require(dplyr)){
      install.packages("dplyr")
   }
   if(!require(ggplot2)){
      install.packages("ggplot2")
   }
   if(!require(ggsignif)){
      install.packages("ggsignif")
   }
   if(!require(gridExtra)){
      install.packages("gridExtra")
   }
   if(!require(tidyverse)){
      install.packages("tidyverse")
   }
   if(!require(stringr)){
      install.packages("stringr")
   }
   if(!require(vegan)){
      install.packages("vegan")
   }
   if(!require(cowplot)){
      install.packages("cowplot")
   }
   if(!require(gtable)){
      install.packages("gtable")
   }
   if(!require(grid)){
      install.packages("grid")
   }
   if(!require(ggrepel)){
      install.packages("ggrepel")
   }
}


check_and_load_package <- function(package_name) {
   if (!requireNamespace(package_name, quietly = TRUE)) {
      stop(paste("The required package", package_name, "is not installed. Please install it using install.packages('", package_name, "') and try again.", sep = ""), "\nHave you run the microbAIDeR::microbAIDeR_install_dependancies() function?\n")
   } else {
      library(package_name, character.only = TRUE)
   }
}


tokenize <- function(x) {
   if(!is.factor(x)){
      print(paste("Input object must be a factor"))
   } else {
      tok <- c()
      for (len in 1:length(x)) {
         tok = c(tok, which(levels(x) == x[len]) )
      }
      return(tok)
   }
}


one_hot_encode <- function(x){
   if(!is.factor(x)){
      print(paste("Input object must be a factor"))
   } else {
         encoded = data.frame()
      for (len in 1:length(x)){
         encoded = rbind(encoded, rep(0, nlevels(x)))
         encoded[ nrow(encoded) , which(levels(x) == x[len]) ] = 1
      }
      colnames(encoded) = levels(x)
      return(encoded)
   }
}


amp_up_models <- function(core_fraction = 0.8){
   check_and_load_package("parallel")
   check_and_load_package("doParallel")
   no_cores <- parallel::detectCores() * core_fraction
   cluster <- parallel::makePSOCKcluster(no_cores)
   doParallel::registerDoParallel(cluster)
   cat("Model amped and ready to go with:", no_cores, "cores. \n")
}


tidify <- function(x, as.numeric = TRUE){
   check_and_load_package("tidyr")
   check_and_load_package("tibble")
   y <- x %>%
      tibble::rownames_to_column() %>%
      tidyr::gather(colname, value, -rowname)
   if(as.numeric){y$value = as.numeric(y$value)}
   return(y)
}


sigFunction <- function(x, trends) {
   if (!is.na(x)) {
      if (x <= 0.0001) {
         "****"
      } else if (x <= 0.001) {
         "***"
      } else if (x <= 0.01) {
         "**"
      } else if (x <= 0.05) {
         "*"
      } else if (trends & x <= 0.1) {
         "°"
      } else {
         NA
      }
   } else {
      NA
   }
}


addSheet <- function(path,sheet.name, addition, col.save=TRUE, row.save=TRUE, overwrite = TRUE){
   check_and_load_package("openxlsx")
   wb <- openxlsx::loadWorkbook(path)
   if ( sheet.name %in% openxlsx::getSheetNames(path) & overwrite ){ openxlsx::removeWorksheet(wb, sheet.name) }
   openxlsx::addWorksheet(wb,sheet.name)
   openxlsx::writeData(wb,sheet.name, as.data.frame(addition), colNames = col.save, rowNames = row.save )
   openxlsx::saveWorkbook(wb,path,overwrite = TRUE)
}


process_pairwise_wilcoxon <- function(testing, onlysig = TRUE){
   if (onlysig)
   {
      significant_comparisons = as.data.frame(matrix(ncol=2))
      colnames(significant_comparisons) = c("Comparison", "Wilcox.pvalue")
      insert = 1
      if ( any(testing$p.value <= 0.1, na.rm=TRUE) )
      {
         for ( line in 1:nrow(testing$p.value) )
         {
            if ( any(testing$p.value[line,] <= 0.1, na.rm=TRUE) )
            {
               for ( column in 1:length(testing$p.value[line,]) )
               {
                  if ( !is.na(testing$p.value[line, column]) & testing$p.value[line, column] <= 0.1  )
                  {
                     significant_comparisons[insert, 1] = paste( rownames(testing$p.value)[line], "vs", colnames(testing$p.value)[column] )
                     significant_comparisons[insert, 2] = as.numeric( testing$p.value[line,column] )
                     microbAIDeR::inc(insert, increment = 1)
                  }
               }
            }
         }
      }
   } else if (!onlysig) {
      significant_comparisons = as.data.frame(matrix(ncol=2))
      colnames(significant_comparisons) = c("Comparison", "Wilcox.pvalue")
      insert = 1
      for ( line in 1:nrow(testing$p.value) )
      {
         for ( column in 1:length(testing$p.value[line,]) )
         {
            if ( !is.na(testing$p.value[line, column]) )
            {
               significant_comparisons[insert, 1] = paste( rownames(testing$p.value)[line], "vs", colnames(testing$p.value)[column] )
               significant_comparisons[insert, 2] = as.numeric( testing$p.value[line,column] )
               microbAIDeR::inc(insert, increment = 1)
            }
         }
      }
   } else { significant_comparisons = NA }

   return(significant_comparisons)
}



is.nan.data.frame <- function(x)
{do.call(cbind, lapply(x, is.nan))}


distr.stats <- function(x) {
   x = x[!is.na(x)]
   x = x[!is.nan(x)]

   avg <- mean(x)
   sd.avg <- sd(x)
   se.avg <- sd.avg / sqrt(length(x))

   med <- median(x)
   iqr <- IQR(x)
   se.med <- 1.2533 * iqr / sqrt(length(x))

   output = list(avg, sd.avg, se.avg, med, iqr, se.med)
   names(output) = c("mean", "SD", "SEM", "median", "iqr", "SEMed")

   return(output)

}


sem <- function(x) {
   x = x[!is.na(x)]
   n <- length(x)
   se <- sd(x) / sqrt(n)
   return(se)
}


se.median <- function(x){
   x = x[!is.na(x)]
   se.med <- 1.2533 * iqr / sqrt(length(x))
   return(se.med)
}



firstup <- function(string) {
   substr(string, 1, 1) <- toupper(substr(string, 1, 1))
   return(string)
}


normalize <- function(x) {
   return ((x - min(x)) / (max(x) - min(x))) }


inc <- function(x, increment = 1) {
   eval.parent(substitute(x <- x + increment))
}


align_legend <- function(p, hjust = 0.5)
{
   check_and_load_package("cowplot")
   check_and_load_package("gtable")
   # extract legend
   g <- cowplot::plot_to_gtable(p)
   grobs <- g$grobs
   if ( any(sapply(grobs, function(x) "guide-box" %in% x$name)) )
   {
      legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
      legend <- grobs[[legend_index]]

      # extract guides table
      guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")

      # there can be multiple guides within one legend box
      for (gi in guides_index) {
         guides <- legend$grobs[[gi]]

         # add extra column for spacing
         # guides$width[5] is the extra spacing from the end of the legend text
         # to the end of the legend title. If we instead distribute it by `hjust:(1-hjust)` on
         # both sides, we get an aligned legend
         spacing <- guides$width[5]
         guides <- gtable::gtable_add_cols(guides, hjust*spacing, 1)
         guides$widths[6] <- (1-hjust)*spacing
         title_index <- guides$layout$name == "title"
         guides$layout$l[title_index] <- 2

         # reconstruct guides and write back
         legend$grobs[[gi]] <- guides
      }

      # reconstruct legend and write back
      g$grobs[[legend_index]] <- legend
      g
   }

}




compute_wilcoxon_and_plot <- function(data, group, taxlevel, save.path = getwd(), color.grouping, comparison.list = NULL, p.adjust.method = "fdr", trends = TRUE, plot.not.sig = FALSE,
                                      paired = FALSE, nrow.graph = 2, ncol.graph = 2, width.graph = 6, height.graph = 4, horiz = FALSE, ggplot.margins = c(.18, .18, .18, .6),
                                      box.lwd = 0.4, jitter.pch = 21, jitter.stroke = 0.15, jitter.size = 0.7, jitter.color = "grey22",
                                      signif.step.increase = 0.12, signif.text.size = 3.5, signif.line.size = 0.4, contrast.color = "ivory1",
                                      text.x.size = 8, text.y.size = 7,  text.y.title.size = 9,
                                      smoothing = FALSE, smoothing.lwd = 1, smoothing.color = "darkred", smoothing.se = FALSE, smoothing.method="loess",
                                      additional.params = NULL, align.legend = FALSE){
   check_and_load_package("dplyr")
   check_and_load_package("ggplot2")
   check_and_load_package("ggsignif")
   check_and_load_package("gridExtra")
   check_and_load_package("tidyverse")
   check_and_load_package("tidyr")
   check_and_load_package("openxlsx")


   if(!dir.exists(save.path)){dir.create(save.path)}

   if ( any(is.na(group))){
      data = data[,!is.na(group)]
      group = group[!is.na(group)]
   }

   data = data[!rowSums(data) %in% 0,]

   if( ncol(data) == 0 ){ stop("No samples remaining after filtering NAs from the grouping factor. Check your grouping factor") }
   if( nrow(data) == 0 ){ stop("No features remaining after filtering features with 0 abundance across all samples. Check your data") }
   if(!p.adjust.method %in% p.adjust.methods){ stop("Supplied an invalid value to p.adjust.method parameter.")}

   kw <- rep(NA, nrow(data))
   for( i in 1:nrow(data))
   {
      kw[i] = kruskal.test( as.numeric(data[i,]) , group )$p.value
   }
   kw[is.nan(kw)] = 1
   kw[is.na(kw)] = 1
   data = data[order(kw),]

   gvec <- vector("list",length=0 )
   gvec_corr <- vector("list",length=0 )
   iter1 = 1
   iter2 = 1

   if( isTRUE(trends) )
   {
      p_threshold <- 0.1
   } else if ( isFALSE(trends) ) {
      p_threshold <- 0.05
   } else {
      stop("Wrong value supplied to the 'trends' argument")
   }

   if ( class(comparison.list) != "list" & !is.null(comparison.list) )
   {
      cat("Wrong data structure passed to the item 'comparison.list'. Please ensure to supply a correctly formatted list or set it to NULL")
   }

   call.print = as.data.frame(rbind(p.adjust.method, taxlevel, save.path, color.grouping = paste(color.grouping, collapse = ", "), smoothing.color, contrast.color,
                                    comparison.list = paste(sapply(comparison.list, function(x) paste(x, collapse = " vs ")), collapse=","),
                                    nrow.graph , ncol.graph , width.graph , height.graph , ggplot2.margins = paste(ggplot.margins, collapse = ", ") ,
                                    box.lwd , jitter.pch , jitter.stroke , jitter.size , jitter.color ,
                                    signif.step.increase,  signif.text.size ,
                                    text.x.size , text.y.size ,  text.y.title.size, additional.params=paste(additional.params, collapse=", ")
   ))

   write.xlsx( call.print , paste(save.path, "/", taxlevel, "_wilcoxon_pairwise.xlsx", sep=""), sheetName = "call",
               colNames = FALSE, rowNames = TRUE)

   if( is.null(comparison.list) )
   {
      comparison.list <- combn(levels(group), 2, simplify = FALSE)
   }

   ptab = as.data.frame(
      matrix(
         nrow = nrow(data),
         ncol = length(comparison.list) + 2*nlevels(group)
      )
   )
   rownames(ptab) = rownames(data)

   collapse_elements <- function(x, collapse = " vs ") {
      paste(x, collapse = collapse)
   }

   colnames(ptab) = c(
      sapply(comparison.list, collapse_elements),
      unlist(apply(expand.grid(c("Mean", "SEM"), levels(group)), 1, function(x) paste(x[1], x[2], sep = " ")))
   )

   # Compute all (pairwise)-Wilcoxon tests and plot uncorrected boxplots
   for (taxa in rownames(data) )
   {
      tempdf = as.data.frame(t(data)) %>% select(abundances = all_of(taxa)) %>% mutate(grouping_column = group )
      wtests = as.data.frame(matrix(ncol = 2))
      for ( iterconfr in c(1:length(comparison.list)))
      {
         confr = comparison.list[[iterconfr]]
         wtests[iterconfr , 1] = paste(confr, collapse = " vs ")
         wtests[iterconfr , 2] = suppressWarnings(wilcox.test( tempdf$abundances[tempdf$grouping_column == confr[1]] ,tempdf$abundances[tempdf$grouping_column == confr[2]], paired = paired )$p.value)
      }
      wtests$V2[is.nan(wtests$V2)] = 1
      wtests$V2[is.na(wtests$V2)] = 1
      wtests$V2[wtests$V2 %in% c(Inf, -Inf)] = 1
      if (any(wtests$V2 <= p_threshold)) {onlysig <- wtests[wtests$V2 <= p_threshold, ]} else {onlysig <- NULL}
      ############ Uncorrected boxplots
      if ( !is.null(onlysig) | plot.not.sig ) {
         suppressWarnings(graphy <- ggplot( data=tempdf, aes(x=grouping_column, y=as.numeric(abundances)) )+
                             geom_boxplot( aes(fill=grouping_column),outlier.shape = NA, lwd = box.lwd) +
                             geom_jitter( aes(color=grouping_column), width = 0.1, height=0, pch = jitter.pch, color=jitter.color, stroke = jitter.stroke, size = jitter.size)+ {
                                if (!is.null(onlysig))
                                   geom_signif(comparisons = strsplit(onlysig$V1, split = " vs "),
                                               annotations = sapply(onlysig$V2, sigFunction, trends = trends),
                                               size = signif.line.size,
                                               tip_length = 0.02, color = "grey22", textsize = signif.text.size,
                                               margin_top = 0.05, vjust = 0.5, step_increase = signif.step.increase)
                             } +
                             {if (smoothing)
                                geom_smooth(aes(x=microbAIDeR::tokenize(grouping_column)), method=smoothing.method, lwd=smoothing.lwd, color=contrast.color, se=smoothing.se, formula = 'y ~ x' )
                             }+
                             {if (smoothing)
                                geom_smooth(aes(x=microbAIDeR::tokenize(grouping_column)), method=smoothing.method, lwd=(0.5*smoothing.lwd), color=smoothing.color, se=smoothing.se, formula = 'y ~ x' )
                             }+
                             scale_fill_manual(values=color.grouping)+
                             labs(x="", y=taxa ) +
                             {if (horiz)
                                coord_flip()
                             }+
                             theme(legend.position = 'none', plot.margin = unit(ggplot.margins, "cm"),
                                   axis.text.x = element_text(size = text.x.size),
                                   axis.text.y = element_text(size = text.y.size), axis.title.y = element_text(size = text.y.title.size)
                             )+
                             { if( !is.null(additional.params)) additional.params}

         )
         if ( align.legend ){ graphy <- align_legend(graphy) }
         gvec[[iter1]] = graphy
         microbAIDeR::inc(iter1)
      }
      # Save the pvalues inside ptab
      ptab[taxa,wtests$V1] = wtests$V2
      # Add mean relabb and SEM for each taxa in each group
      for( iterlevel in levels(group))
      {
         ptab[taxa,paste("Mean", iterlevel)] = signif(mean(tempdf$abundances[tempdf$grouping_column == iterlevel]), 3)
         ptab[taxa,paste("SEM", iterlevel)] = signif(microbAIDeR::sem(tempdf$abundances[tempdf$grouping_column == iterlevel]),3)
      }
   }
   # Correct the pvalues if requested
   if (p.adjust.method != "none" & p.adjust.method %in% p.adjust.methods) {
      print(paste("Correcting pvalues with", p.adjust.method))
      corrected_ptab = ptab
      for ( itercol in 1:length(comparison.list) )
      {
         corrected_ptab[,itercol] = p.adjust(corrected_ptab[,itercol], method = p.adjust.method)
      }
      for (taxa in rownames(data) )
      {
         tempdf = as.data.frame(t(data)) %>% select(abundances = all_of(taxa)) %>% mutate(grouping_column = group )
         if ( length(comparison.list) == 1 )
         {
            wtests_corrected = as.data.frame(cbind(colnames(corrected_ptab)[1], corrected_ptab[taxa,1]))
            wtests_corrected$V2 = as.numeric(wtests_corrected$V2)
         } else {
            wtests_corrected = tidify(corrected_ptab[taxa,1:length(comparison.list)])[,c("colname", "value")]
            colnames(wtests_corrected) = colnames(wtests)
         }
         wtests_corrected$V2[is.nan(wtests_corrected$V2)] = 1
         wtests_corrected$V2[is.na(wtests_corrected$V2)] = 1
         wtests_corrected$V2[wtests_corrected$V2 %in% c(Inf, -Inf)] = 1
         if (any(wtests_corrected$V2 <= p_threshold)) {onlysig_corrected <- wtests_corrected[wtests_corrected$V2 <= p_threshold, ]} else {onlysig_corrected <- NULL}
         ############ Uncorrected boxplots
         if ( !is.null(onlysig_corrected) | plot.not.sig ) {
            suppressWarnings(graphy <- ggplot( data=tempdf, aes(x=grouping_column, y=as.numeric(abundances)) )+
                                geom_boxplot( aes(fill=grouping_column),outlier.shape = NA, lwd = box.lwd) +
                                geom_jitter( aes(color=grouping_column), width = 0.1, height=0, pch = jitter.pch, color=jitter.color, stroke = jitter.stroke, size = jitter.size)+ {
                                   if (!is.null(onlysig_corrected))
                                      geom_signif(comparisons = strsplit(onlysig_corrected$V1, split = " vs "),
                                                  annotations = sapply(onlysig_corrected$V2, sigFunction, trends = trends),
                                                  size = signif.line.size,
                                                  tip_length = 0.02, color = "grey22", textsize = signif.text.size,
                                                  margin_top = 0.05, vjust = 0.5, step_increase = signif.step.increase)
                                } +
                                {if (smoothing)
                                   geom_smooth(aes(x=microbAIDeR::tokenize(grouping_column)), method=smoothing.method, lwd=smoothing.lwd, color=contrast.color, se=smoothing.se, formula = 'y ~ x' )
                                }+
                                {if (smoothing)
                                   geom_smooth(aes(x=microbAIDeR::tokenize(grouping_column)), method=smoothing.method, lwd=(0.5*smoothing.lwd), color=smoothing.color, se=smoothing.se, formula = 'y ~ x' )
                                }+
                                scale_fill_manual(values=color.grouping)+
                                labs(x="", y=taxa ) +
                                {if (horiz)
                                   coord_flip()
                                }+
                                theme(legend.position = 'none', plot.margin = unit(ggplot.margins, "cm"),
                                      axis.text.x = element_text(size = text.x.size),
                                      axis.text.y = element_text(size = text.y.size), axis.title.y = element_text(size = text.y.title.size)
                                )+
                                { if( !is.null(additional.params)) additional.params}

            )
            if ( align.legend ){ graphy <- align_legend(graphy) }
            gvec_corr[[iter2]] = graphy
            microbAIDeR::inc(iter2)
         }
      }
   }
   ### Save pvalues
   print(paste("Saving P-values to ", save.path, "/", taxlevel, "_wilcoxon_pairwise.xlsx", sep=""))
   addSheet( path=paste(save.path, "/", taxlevel, "_wilcoxon_pairwise.xlsx", sep=""), sheet.name = "uncorrected",
             addition = ptab, col.save=TRUE, row.save=TRUE)
   p_value_columns <- which(colnames(ptab) %in% sapply(comparison.list, collapse_elements) )
   for (col in p_value_columns) {
      ptab[ ptab[,col] > p_threshold, col] <- ""
   }
   for (col in p_value_columns) {
      ptab[,col] <- ifelse(ptab[,col] != "", signif( as.numeric(ptab[,col]), digits = 3), "")
   }
   addSheet( path=paste(save.path, "/", taxlevel, "_wilcoxon_pairwise.xlsx", sep=""), sheet.name = "uncorrected_clean",
             addition = ptab, col.save=TRUE, row.save=TRUE)
   if (p.adjust.method != "none") {
      addSheet(path = paste(save.path, "/", taxlevel, "_wilcoxon_pairwise.xlsx", sep=""), sheet.name = p.adjust.method, addition = corrected_ptab, col.save = T, row.save = T)
      p_value_columns <- which(colnames(corrected_ptab) %in% sapply(comparison.list, collapse_elements) )
      for (col in p_value_columns) {
         corrected_ptab[ corrected_ptab[,col] > p_threshold, col] <- ""
      }
      for (col in p_value_columns) {
         corrected_ptab[,col] <- ifelse(corrected_ptab[,col] != "", signif( as.numeric(corrected_ptab[,col]), digits = 3), "")
      }
      addSheet(path = paste(save.path, "/", taxlevel, "_wilcoxon_pairwise.xlsx", sep=""), sheet.name = paste(p.adjust.method, "clean", sep="_"), addition = corrected_ptab, col.save = T, row.save = T)
   }
   ### Save uncorrected plots
   if (smoothing){
      print(paste("Saving plots to ", save.path, "/", taxlevel, "_boxplot_wilcox_uncorrected_smooth.pdf", sep=""))
      suppressWarnings(ggsave(
         filename = paste(save.path, "/", taxlevel, "_boxplot_wilcox_uncorrected_smooth.pdf", sep="") ,
         plot = marrangeGrob(gvec, nrow=nrow.graph, ncol=ncol.graph, top=NULL,
                             layout_matrix = matrix(seq_len(nrow.graph*ncol.graph), nrow = nrow.graph, ncol = ncol.graph , byrow = TRUE) ),
         width = width.graph, height = height.graph, dpi = 330
      ))
   } else {
      print(paste("Saving plots to ", save.path, "/", taxlevel, "_boxplot_wilcox_uncorrected.pdf", sep=""))
      suppressWarnings(ggsave(
         filename = paste(save.path, "/", taxlevel, "_boxplot_wilcox_uncorrected.pdf", sep="") ,
         plot = marrangeGrob(gvec, nrow=nrow.graph, ncol=ncol.graph, top=NULL ,
                             layout_matrix = matrix(seq_len(nrow.graph*ncol.graph), nrow = nrow.graph, ncol = ncol.graph , byrow = TRUE) ),
         width = width.graph, height = height.graph, dpi = 330
      ))
   }
   ### Save corrected plots
   if ( p.adjust.method != 'none' )
   {
      if (smoothing){
         print(paste("Saving plots to ", save.path, "/", taxlevel, "_boxplot_wilcox_", p.adjust.method, "_smooth.pdf", sep=""))
         suppressWarnings(ggsave(
            filename = paste(save.path, "/", taxlevel, "_boxplot_wilcox_", p.adjust.method, "_smooth.pdf", sep=""),
            plot = marrangeGrob(gvec_corr, nrow=nrow.graph, ncol=ncol.graph, top=NULL ,
                                layout_matrix = matrix(seq_len(nrow.graph*ncol.graph), nrow = nrow.graph, ncol = ncol.graph , byrow = TRUE) ),
            width = width.graph, height = height.graph, dpi = 330
         ))
      } else {
         print(paste("Saving plots to ", save.path, "/", taxlevel, "_boxplot_wilcox_", p.adjust.method, ".pdf", sep=""))
         suppressWarnings(ggsave(
            filename = paste(save.path, "/", taxlevel, "_boxplot_wilcox_", p.adjust.method, ".pdf", sep=""),
            plot = marrangeGrob(gvec_corr, nrow=nrow.graph, ncol=ncol.graph, top=NULL ,
                                layout_matrix = matrix(seq_len(nrow.graph*ncol.graph), nrow = nrow.graph, ncol = ncol.graph , byrow = TRUE) ),
            width = width.graph, height = height.graph, dpi = 330
         ))
      }
   }

}


compute_beta_diversity <- function(beta_metrics = c("braycurtis", "jaccard", "unweighted_unifrac", "weighted_unifrac"),
                                   save.path = getwd(), adonis_n_perm = 9999, beta.folder.path = paste(getwd(), "RESULTS/beta_diversity", sep="/"), manual.beta.path = NULL, sample.list, mds = c(1,2), group, color.grouping ,
                                   spiders = FALSE, spiders.lwd = 1.5, ellipses = FALSE, ellipse.focus = FALSE,ellipse.fill = FALSE, ellipse.conf = 0.95, ellipse.alpha = 0.75 , ellipse.lwd = 2.5, manual.bordercol = NULL,
                                   svg.width = 7, svg.height = 5,
                                   cex.points = 1.2, adonis.p.adjust = "fdr", wilcoxon.p.adjust = "fdr",
                                   nrow.graph = 2, ncol.graph = 2, width.graph = 6, height.graph = 4, ggplot.margins = c(0.15, 0.15, 0.15, 0.6),
                                   text.x.size = 8, text.y.size = 7,  text.title.size = 9,
                                   additional.params=NULL){
   check_and_load_package("openxlsx")
   check_and_load_package("dplyr")
   check_and_load_package("tidyr")
   check_and_load_package("ggplot2")
   check_and_load_package("vegan")
   check_and_load_package("gridExtra")
   check_and_load_package("pairwiseAdonis")
   check_and_load_package("stringr")

   if ( any(is.na(group))){print(paste("WARNING: You have NAs in your group, removing such samples from the beta distance matrices. BE SURE not to have NAs in the sample.list"))}
   if ( length(group) != length(sample.list)){
      print(paste("WARNING: The grouping factor provided has a different length compared to the sample.list"))
      stop()
   }

   if(!dir.exists(save.path)){dir.create(save.path)}

   call.print = as.data.frame(rbind( paste(beta_metrics, collapse=", "), color.grouping = paste(color.grouping, collapse = ", "),
                                     save.path, adonis_n_perm , beta.folder.path ,  mds = paste(mds, collapse=", ") ,
                                     spiders, ellipses , ellipse.focus , ellipse.conf , ellipse.alpha , 
                                     svg.width , svg.height ,
                                     cex.points , adonis.p.adjust, wilcoxon.p.adjust, additional.params=paste(additional.params, collapse=", ")
   ))

   openxlsx::write.xlsx( call.print , paste(save.path, "/beta_adonis.xlsx", sep=""), sheetName = "call",
               colNames = FALSE, rowNames = TRUE)
   openxlsx::write.xlsx( call.print , paste(save.path, "/beta_variance.xlsx", sep=""), sheetName = "call",
               colNames = FALSE, rowNames = TRUE)

   iter1 = 1
   gvec <- vector("list",length=0 )

   for ( metrics in beta_metrics)
   {
      print(paste("Running", toupper(metrics), "beta diversity analyses"))
      if (is.null(manual.beta.path)) {
         beta = read.delim( paste(beta.folder.path , metrics, "distance-matrix.tsv", sep="/"), header=T, row.names=1, sep="\t")
      } else if( !is.null(manual.beta.path)) {
         beta = read.delim( manual.beta.path, header=T, row.names=1, sep="\t")
      }
      beta = beta[ sample.list, sample.list ]
      if ( any(is.na(group))){
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
         print(paste("Errors in producing the color vector. Check that your group is a factor and nlevels(group) == length(color.grouping)"))
         print(paste("Procededing anyway without colors"))
         color_beta[is.na(color_beta)] = "black"
      }
      colnames(coord_beta) = c("AX1", "AX2")
      # Base plot
      svg( paste(save.path,"/", metrics, ".svg", sep=""), width=svg.width, height = svg.height)
      plot(metric_beta, type="n", main=metrics, choices=mds, xlab=paste("MDS", mds[1]," - [" ,sig1, "%]", sep="") , ylab=paste("MDS", mds[2]," - [" ,sig2, "%]", sep="") )
      with(coord_beta, points(AX1, AX2, pch=21, bg=color_beta, cex=cex.points, col=NULL))
      legend("topright", legend=levels(group_beta), pt.cex=1, fill=color.grouping, horiz=F, x.intersp = 0.4, y.intersp = 0.7, bty="n", cex=0.8)
      graphics.off()
      #Spiders
      if(spiders){
         svg( paste(save.path, "/", metrics, "_spiders.svg", sep=""), width=svg.width, height = svg.height)
         plot(metric_beta, type="n", main=metrics, choices=mds, xlab=paste("MDS", mds[1]," - [" ,sig1, "%]", sep="") , ylab=paste("MDS", mds[2]," - [" ,sig2, "%]", sep="") )
         with(coord_beta, points(AX1, AX2, pch=21, bg=color_beta, cex=(cex.points/1.5), col=NULL))
         for ( groupiter in levels(group_beta))
         {
            ordispider(metric_beta, group_beta, display="sites", spiders="centroid", show.groups=groupiter, col= if ( is.null(manual.bordercol) ){ color.grouping[which(levels(group_beta) == groupiter)] } else { manual.bordercol } , lwd=spiders.lwd, lty=1)
         }
         legend("topright", legend=levels(group_beta), pt.cex=1, fill=color.grouping, horiz=F, x.intersp = 0.4, y.intersp = 0.7, bty="n", cex=0.8)
         graphics.off()
      }
      # Ellipse plot
      if(ellipses){
         svg( paste(save.path, "/", metrics, "_ellipse.svg", sep=""), width=svg.width, height = svg.height)
         plot(metric_beta, type="n", main=metrics, choices=mds, xlab=paste("MDS", mds[1]," - [" ,sig1, "%]", sep="") , ylab=paste("MDS", mds[2]," - [" ,sig2, "%]", sep="") )
         with(coord_beta, points(AX1, AX2, pch=21, bg=color_beta, cex=(cex.points/2), col=NULL))
         for ( groupiter in levels(group_beta))
         {
            ordiellipse(metric_beta, group_beta, display="sites", kind=c("se"), conf=ellipse.conf, draw=c("lines"), show.groups=groupiter, col= if ( is.null(manual.bordercol) ){ color.grouping[which(levels(group_beta) == groupiter)] } else { manual.bordercol } , lwd=ellipse.lwd, lty=1, alpha=50)
         }
         legend("topright", legend=levels(group_beta), pt.cex=1, fill=color.grouping, horiz=F, x.intersp = 0.4, y.intersp = 0.7, bty="n", cex=0.8)
         graphics.off()
      }
      # Ellipse focus
      if (ellipse.focus){
         svg( paste(save.path, "/", metrics, "_ellipse_focus.svg", sep=""), width=svg.width, height = svg.height)
         plot(metric_beta, type="n", main=metrics, choices=mds, xlab=paste("MDS", mds[1]," - [" ,sig1, "%]", sep="") , ylab=paste("MDS", mds[2]," - [" ,sig2, "%]", sep="") )
         with(coord_beta, points(AX1, AX2, pch=21, bg=color_beta, cex=(cex.points/4), col=NULL))
         for ( groupiter in levels(group_beta))
         {
            ordiellipse(metric_beta, group_beta, display="sites", kind=c("se"), conf=ellipse.conf, draw=c("polygon"), show.groups=groupiter, col= color.grouping[which(levels(group_beta) == groupiter)], border = if ( is.null(manual.bordercol) ){ color.grouping[which(levels(group_beta) == groupiter)] } else { manual.bordercol }, lwd=ellipse.lwd, lty=1, alpha=ellipse.alpha)
         }
         legend("topright", legend=levels(group_beta), pt.cex=1, fill=color.grouping, horiz=F, x.intersp = 0.4, y.intersp = 0.7, bty="n", cex=0.8)
         graphics.off()
      }
      # Ellipse fill
      if (ellipse.fill){
         svg( paste(save.path, "/", metrics, "_ellipse_fill.svg", sep=""), width=svg.width, height = svg.height)
         plot(metric_beta, type="n", main=metrics, choices=mds, xlab=paste("MDS", mds[1]," - [" ,sig1, "%]", sep="") , ylab=paste("MDS", mds[2]," - [" ,sig2, "%]", sep="") )
         with(coord_beta, points(AX1, AX2, pch=21, bg=color_beta, cex=cex.points, col=NULL))
         for ( groupiter in levels(group_beta))
         {
            ordiellipse(metric_beta, group_beta, display="sites", kind=c("se"), conf=ellipse.conf, draw=c("polygon"), show.groups=groupiter, col= color.grouping[which(levels(group_beta) == groupiter)], border = if ( is.null(manual.bordercol) ){ color.grouping[which(levels(group_beta) == groupiter)] } else { manual.bordercol }, lwd=1, lty=1, alpha=ellipse.alpha)
         }
         legend("topright", legend=levels(group_beta), pt.cex=1, fill=color.grouping, horiz=F, x.intersp = 0.4, y.intersp = 0.7, bty="n", cex=0.8)
         graphics.off()
      }
      # ADONIS
      print(paste("Computing pairwise Adonis for", metrics))
      adone = pairwise.adonis(beta, group_beta, perm = adonis_n_perm, p.adjust.m = adonis.p.adjust)
      trends = TRUE #Display trends in pairwise adonis p-values table
      adone$uncorrected_pvalues <- sapply(adone$p.value, sigFunction, trends = trends)
      adone[,adonis.p.adjust] <- sapply(adone$p.adjusted, sigFunction, trends = trends)

      addSheet( path=paste(save.path, "/beta_adonis.xlsx", sep=""), sheet.name = metrics,
                addition = adone[,c(1,6,9, 7,10)], col.save=TRUE, row.save=FALSE)

      # Calculate and plot the intra-group variances
      print(paste("Computing intra-group variance"))
      beta_half = beta
      beta_half[lower.tri(beta_half, diag=T)] <- NA
      beta_half = as.data.frame(beta_half)
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
      plot_beta_variance = as.data.frame(cbind(values_beta_var, group_beta_var))
      plot_beta_variance$values_beta_var = as.numeric( plot_beta_variance$values_beta_var )
      plot_beta_variance$group_beta_var = factor(plot_beta_variance$group_beta_var, levels=levels(group_beta))
      var_df <- plot_beta_variance %>%
         group_by(group_beta_var) %>%
         summarise(GroupVariance = var(values_beta_var))
      colnames(var_df) = c("Group", "Variance")
      ## Plot variance
      graphy <- ggplot(var_df, aes(x = Group, y = Variance, fill=Group)) +
         geom_bar(stat = "identity", color="black", lwd = 0.4) +
         scale_fill_manual(values=color.grouping)+
         labs(x = "", y = "Intra-Group Variance", title = metrics)+
         theme(legend.position = "none", plot.margin = unit(ggplot.margins, "cm"),
               axis.text.x = element_text(size = text.x.size),
               axis.text.y = element_text(size = text.y.size), axis.title.y = element_text(size = text.title.size)
         )+
         { if( !is.null(additional.params)) additional.params}

      gvec[[iter1]] = graphy
      inc(iter1, increment = 1)
      # IQR and Kruskal tests
      pairwise_distances <- list()
      for (gg in levels(group_beta)) {
         group_data <- plot_beta_variance[group_beta_var == gg, "values_beta_var"]
         distances <- dist(group_data)
         pairwise_distances[[as.character(gg)]] <- as.vector(distances)
      }
      ranges <- sapply(pairwise_distances, function(distances) diff(range(distances)))
      iqr <- sapply(pairwise_distances, IQR)
      comparison_result <- data.frame(
         Range = ranges,
         IQR = iqr
      )
      comparison_result <- tidify(comparison_result)
      #Plot IQR
      graphy <- ggplot(comparison_result[comparison_result$colname=="IQR",], aes(x = rowname, y = value, fill=rowname)) +
         geom_bar(stat = "identity") +
         scale_fill_manual(values=color.grouping)+
         labs(title = metrics,
              x = "Group", y = "IQR") +
         theme(legend.position = "none", plot.margin = unit(ggplot.margins, "cm"),
               axis.text.x = element_text(size = text.x.size),
               axis.text.y = element_text(size = text.y.size), axis.title.y = element_text(size = text.title.size)
         )+
         { if( !is.null(additional.params)) additional.params}
      gvec[[iter1]] = graphy
      inc(iter1, increment = 1)
      #Plot Ranges
      graphy <- ggplot(comparison_result[comparison_result$colname=="Range",], aes(x = rowname, y = value, fill=rowname)) +
         geom_bar(stat = "identity") +
         scale_fill_manual(values=color.grouping)+
         labs(title = metrics,
              x = "Group", y = "Range") +
         theme(legend.position = "none", plot.margin = unit(ggplot.margins, "cm"),
               axis.text.x = element_text(size = text.x.size),
               axis.text.y = element_text(size = text.y.size), axis.title.y = element_text(size = text.title.size)
         )+
         { if( !is.null(additional.params)) additional.params}
      gvec[[iter1]] = graphy
      inc(iter1, increment = 1)
      ktest = kruskal.test(plot_beta_variance$values_beta_var, plot_beta_variance$group_beta_var)$p.value
      wtest_raw = process_pairwise_wilcoxon(
         suppressWarnings(
            pairwise.wilcox.test(plot_beta_variance$values_beta_var, plot_beta_variance$group_beta_var, p.adjust.method = "none") ),
         onlysig = FALSE)
      wtest_corrected = process_pairwise_wilcoxon(
         suppressWarnings(
            pairwise.wilcox.test(plot_beta_variance$values_beta_var, plot_beta_variance$group_beta_var, p.adjust.method = wilcoxon.p.adjust)),
         onlysig = FALSE)
      comparisondf = data.frame(
         Comparison = wtest_raw$Comparison,
         uncorrected_p = signif(wtest_raw$Wilcox.pvalue,3),
         wilcox.p.adjust = signif(wtest_corrected$Wilcox.pvalue,3),
         Kruskal_test_p = ""
      )
      comparisondf[1, "Kruskal_test_p"] <- signif(ktest,3)
      # Export the data in the Excel files
      colnames(var_df) = c("rowname", "value")
      var_df$colname = "Variance"
      addSheet( path=paste(save.path, "/beta_variance.xlsx", sep=""), sheet.name = paste( gsub("weighted", "w", metrics), "Var", sep="-" ),
                addition = as.data.frame(rbind(var_df, comparison_result)), col.save=TRUE, row.save=FALSE)
      addSheet( path=paste(save.path, "/beta_variance.xlsx", sep=""), sheet.name = paste( gsub("weighted", "w", metrics), "p", sep="-" ),
                addition = comparisondf, col.save=TRUE, row.save=FALSE)
   }
   print(paste("Saving plots to ", save.path, "/intragroup_variance.pdf", sep=""))
   suppressWarnings(ggsave(
      filename = paste(save.path, "/intragroup_variance.pdf", sep="") ,
      plot = marrangeGrob(gvec, nrow=nrow.graph, ncol=ncol.graph, top=NULL ,
                          layout_matrix = matrix(seq_len(nrow.graph*ncol.graph), nrow = nrow.graph, ncol = ncol.graph , byrow = TRUE) ),
      width = width.graph, height = height.graph, dpi = 330
   ))
}



flattenCorrMatrix <- function(cormat, pmat) {
   ut <- upper.tri(cormat)
   df = data.frame(
         row = rownames(cormat)[row(cormat)[ut]],
         column = rownames(cormat)[col(cormat)[ut]],
         cor  =(cormat)[ut],
         p = pmat[ut]
      )
   return(df)
}


correlation_custom <- function(feature_data, parameter_data, is.binary.parameter = FALSE, p.adjust.method = "fdr",
                               uncorr.p.threshold = 0.05, corr.p.threshold = 0.05, min.corr.threshold = 0.3,
                               force.numeric = TRUE, corr.method = "spearman"){
   check_and_load_package("dplyr")
   # Check if data sets are matched
   if ( !all(rownames(feature_data) == rownames(parameter_data)) )
   {
      stop("Rownames in both datasets must match (sample names)")
   }

   # Initializes the progress bar
   pb <- txtProgressBar(min = 1,      # Minimum value of the progress bar
                        max = ncol(feature_data), # Maximum value of the progress bar
                        style = 3,    # Progress bar style (also available style = 1 and style = 2)
                        width = 50,   # Progress bar width. Defaults to getOption("width")
                        char = "#")   # Character used to create the bar

   # Check if parameter data is numeric
   is_num <- c()
   for (i in 1:ncol(parameter_data)){is_num <- c(is_num, class(parameter_data[,i]))}

   if( !all(is_num %in% c("numeric", "integer")) )
   {
      if ( force.numeric )
      {
         for (i in 1:ncol(parameter_data)){ parameter_data[,i] = as.numeric(parameter_data[,i])}
         is_num <- c()
         for (i in 1:ncol(parameter_data)){is_num <- c(is_num, class(parameter_data[,i]))}
         if( !all(is_num %in% c("numeric", "integer")) )
         {
            print(paste("Check columns", which(! is_num %in% c("numeric", "integer") ), "of your parameter data"))
            stop("Detected non-numeric data in your parameter data")
         } else if ( all(is_num %in% c("numeric", "integer")) )
         {
            print(paste("Found some non-numeric values in your parameter data, coerced as.numeric because force.numeric was TRUE"))
            print(paste("I strongly suggest to have a look at your parameter data"))
         }
      } else if ( !force.numeric )
      {
         print(paste("Check columns", which(! is_num %in% c("numeric", "integer") ), "of your parameter data"))
         stop("Detected non-numeric data in your parameter data")
      }
   }

   # Check if feature data is numeric
   is_num <- c()
   for (i in 1:ncol(feature_data)){is_num <- c(is_num, class(feature_data[,i]))}

   if( !all(is_num %in% c("numeric", "integer")) )
   {
      if ( force.numeric )
      {
         for (i in 1:ncol(feature_data)){ feature_data[,i] = as.numeric(feature_data[,i])}
         is_num <- c()
         for (i in 1:ncol(feature_data)){is_num <- c(is_num, class(feature_data[,i]))}
         if( !all(is_num %in% c("numeric", "integer")) )
         {
            print(paste("Check columns", which(! is_num %in% c("numeric", "integer") ), "of your parameter data"))
            stop("Detected non-numeric data in your parameter data")
         } else if ( all(is_num %in% c("numeric", "integer")) )
         {
            print(paste("Found some non-numeric values in your parameter data, coerced as.numeric because force.numeric was TRUE"))
            print(paste("I strongly suggest to have a look at your parameter data"))
         }
      } else if ( !force.numeric )
      {
         print(paste("Check columns", which(! is_num %in% c("numeric", "integer") ), "of your parameter data"))
         stop("Detected non-numeric data in your parameter data")
      }
   }

   if( is.binary.parameter )
   {
      corr.method = "pearson"
      print(paste("Proceeding with Point-Biserial correlations (Pearson's special case)"))
   } else { print(paste("Proceeding with", corr.method, "correlations"))}


   # Compute correlations
   cortab <- as.data.frame(matrix(nrow = ncol(feature_data), ncol = ncol(parameter_data)))
   ptab <- as.data.frame(matrix(nrow = ncol(feature_data), ncol = ncol(parameter_data)))

   for ( taxn in 1:ncol(feature_data))
   {
      setTxtProgressBar(pb, taxn)
      for (paramn in 1:ncol(parameter_data))
      {
         distrib = feature_data[,taxn]
         cats = parameter_data[,paramn]
         if ( sum(!is.na(distrib)) > 3 & sum(!is.na(cats)) > 3 )
         {
            if ( length(unique(cats[!is.na(cats)])) != 1 )
            {
               cortemp = suppressWarnings( cor.test(distrib, cats, method = corr.method) )
               ptab[taxn,paramn] = cortemp$p.value
               cortab[taxn,paramn] = cortemp$estimate
            }
         }
      }
   }

   rownames(cortab) = colnames(feature_data)
   colnames(cortab) = colnames(parameter_data)
   rownames(ptab) = colnames(feature_data)
   colnames(ptab) = colnames(parameter_data)

   # Filter empty columns and rows
   torem <- c()
   for( i in 1:ncol(ptab)){ if(all(is.na(ptab[,i]))){ torem <- c(torem, i)} }
   if( !is.null(torem))
   {
      cortab = cortab[,-torem]
      ptab = ptab[,-torem]
   }
   torem <- c()
   for( i in 1:nrow(ptab)){ if(all(is.na(ptab[i,]))){ torem <- c(torem, i)} }
   if( !is.null(torem))
   {
      cortab = cortab[-torem,]
      ptab = ptab[-torem,]
   }
   ptab[is.na(ptab)] = 1
   ptab[is.nan.data.frame(ptab)] = 1
   cortab[is.na(cortab)] = 0
   cortab[is.nan.data.frame(cortab)] = 0
   if ( all(dim(cortab) < c(2,2)) )
   {
      stop("No correlations values were computed correcting. All correlations returned NAs. Check your data")
      return( list(uncorr_correlation=cortab, uncorr_pvalues = ptab) )
   }
   # Keep only uncorrected significant comparisons
   nosig_tax = rowSums(ptab < uncorr.p.threshold) == 0
   ptab = ptab %>% dplyr::filter( !nosig_tax )
   cortab = cortab %>% dplyr::filter( !nosig_tax )
   nosig_param = colSums(ptab < uncorr.p.threshold) == 0
   ptab = ptab %>% dplyr::select( names(nosig_param)[!nosig_param]  )
   cortab = cortab %>% dplyr::select( names(nosig_param)[!nosig_param] )
   if ( all(dim(cortab) < c(2,2)) )
   {
      print("Done!")
      stop("No significant correlations even without correcting pvals.")
      return( list(uncorr_correlation=cortab, uncorr_pvalues = ptab) )
   }
   # Correct pvalues and keep only significant
   p2 = ptab
   c2 = cortab
   for (i in 1:ncol(p2)){ p2[,i] = p.adjust(p2[,i], method = p.adjust.method, n = nrow(ptab) )}
   nosig_tax = rowSums(p2 < corr.p.threshold) == 0
   if( sum(!nosig_tax) == 0 ){print(paste("No significant correlation correcting pvalues")) ; return( list(uncorr_correlation=cortab, uncorr_pvalues = ptab) ) ; stop()}
   p2 = p2 %>% dplyr::filter( !nosig_tax )
   c2 = c2 %>% dplyr::filter( !nosig_tax )
   nosig_param = colSums(p2 < corr.p.threshold) == 0
   if( sum(!nosig_param) == 0 ){print(paste("No significant correlation correcting pvalues")) ; return( list(uncorr_correlation=cortab, uncorr_pvalues = ptab) ) ; stop()}
   p2 = p2 %>% dplyr::select( names(nosig_param)[!nosig_param] )
   c2 = c2 %>% dplyr::select( names(nosig_param)[!nosig_param] )
   if ( any(dim(c2) == 1) )
   {
      print(paste("Done!"))
      print(paste( sum(ptab < 0.05) , "correlations found with uncorrected p-values"))
      print(paste( sum(p2 < 0.05) , "correlations found with", p.adjust.method, "adjusted p-values"))
      return( list(uncorr_correlation=cortab, uncorr_pvalues = ptab,
                   corrected_correlation = c2, corrected_pvalues = p2) )
      stop()
   }
   # Keep only correlation stronger than abs(0.3)
   p3 = p2
   c3 = c2
   ctemp = c3
   ctemp[p3 >= corr.p.threshold] = NA
   nostrong_tax = rowSums( abs(ctemp) >= min.corr.threshold, na.rm = TRUE ) == 0
   p3 = p3 %>% dplyr::filter( !nostrong_tax )
   c3 = c3 %>% dplyr::filter( !nostrong_tax )
   ctemp = ctemp %>% dplyr::filter( !nostrong_tax )
   nostrong_param = colSums( abs(ctemp) >= min.corr.threshold, na.rm = TRUE ) == 0
   p3 = p3  %>% dplyr::select(names(nostrong_param)[!nostrong_param] )
   c3 = c3 %>% dplyr::select( names(nostrong_param)[!nostrong_param] )


   cor_out = list(uncorr_correlation = cortab, uncorr_pvalues = ptab,
                  corrected_correlation = c2, corrected_pvalues = p2,
                  corr_strong_correlation = c3, corr_strong_pvalues = p3,
                  call = as.data.frame(rbind(
                     if(is.binary.parameter){paste("Point Biserial correlation")}else{paste(corr.method , "correlation")},
                     paste("Uncorrected p-values threshold", uncorr.p.threshold),
                     paste("Multiple comparison p.adjust method", p.adjust.method),
                     paste("Corrected p-values threshold", corr.p.threshold)
                  ))
   )

   print(paste("Done!"))
   print(paste( sum(ptab < 0.05) , "correlations found with uncorrected p-values"))
   print(paste( sum(p2 < 0.05) , "correlations found with", p.adjust.method, "adjusted p-values"))
   print(paste( sum(p3 < 0.05) , "correlations found with", p.adjust.method, "adjusted p-values and at least", min.corr.threshold, "cor" ))

   return(cor_out)


}


scaling.manual <- function(x, range.min, range.max) {
   return(range.min + (x-min(x))*(range.max-range.min)/(max(x)-min(x)))
}


compute_volcano_plots <- function(data, group, taxlevel, save.path = getwd(), p.adjust.method = "fdr", trends=TRUE, sigcolors=c("ivory3", "orange3", "darkred"),
                                  paired = FALSE, nrow.graph = 1, ncol.graph = 2, width.graph = 20, height.graph = 8, ggplot.margins = c(1,1,1,.5),
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

   if(!dir.exists(save.path)){dir.create(save.path)}

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
      wt[i] = suppressWarnings(wilcox.test( as.numeric(data[i,]) ~ group, paired = paired )$p.value)
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
      print(paste("You have not supplied 3 colors to sigcolors. Reverting to default values"))
      sigcolors=c("ivory3", "orange3", "darkred")
   }

   call.print = as.data.frame(rbind(p.adjust.method, taxlevel, save.path,
                                    sigcolors = paste(sigcolors, collapse = ", "),
                                    nrow.graph , ncol.graph , width.graph , height.graph , ggplot2.margins = paste(ggplot.margins, collapse = ", ") ,
                                    paired, label.text.size, label.padding,arrow.length,arrow.lwd,arrow.origin.offset,
                                    arrow.text.cex,arrow.text.color,
                                    text.x.size , text.y.size ,  text.y.title.size, additional.params=paste(additional.params, collapse=", ")
   ))
   write.xlsx( call.print , paste(save.path, "/", taxlevel, "_volcano.xlsx", sep=""), sheetName = "call",
               colNames = FALSE, rowNames = TRUE)

   ptab = as.data.frame(
      matrix(
         nrow = nrow(data),
         ncol = 5
      )
   )
   rownames(ptab) = rownames(data)
   colnames(ptab) = c(
      paste(levels(group)[1], levels(group)[2], sep=" vs "),
      unlist(apply(expand.grid(c("Mean", "SEM"), levels(group)), 1, function(x) paste(x[1], x[2], sep = " ")))
   )
   if (p.adjust.method != "none") {
      ptab_corr <- ptab
   }
   # Plot volcanos for uncorrected P-values
   print(paste("Terraforming volcano plots.."))
   tidata = as.data.frame(
      matrix(
         nrow=nrow(data),
         ncol = 4,
         dimnames = list(rownames(data), c("pval", "corr", "Group1", "Group2"))))
   for (i in 1:nrow(tidata)){tidata[i,1] <- wt[i]}
   tidata[,2] <- p.adjust(tidata[,1], method = p.adjust.method)
   tidata$Group1 = rowMeans(data[group==levels(group)[1]])
   tidata$Group2 = rowMeans(data[group==levels(group)[2]])
   tidata$LFC <- log((tidata$Group1+0.001)/(tidata$Group2+0.001))
   tidata$Color_uncorrected = NA
   tidata$Color_corr = NA
   for (i in 1:nrow(tidata))
   {
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
   tidata$Relabb <- as.numeric(rowMeans(data))

   graphy <- ggplot(data = tidata, aes(x = LFC, y = -log2(pval), label = Plotlabel_uncorrected, size = Relabb, color = Color_uncorrected)) +
      geom_point(alpha = 0.8) +
      geom_vline(xintercept = 0, col = "grey80") +
      geom_hline(yintercept = -log2(0.05), col = "grey20", linetype = "dashed") +
      scale_color_identity() +
      labs(x="Log-Fold change", size="Rel. Ab.", y="-log2(p)")+
      scale_y_continuous(limits = c(0, max(-log2(tidata$pval))+.5 ) )+
      scale_x_continuous(limits = c( min(tidata$LFC)-.5, max(tidata$LFC)+.5 ))+
      geom_text_repel(show.legend=F, size = label.text.size, point.padding = label.padding, aes(fontface = "italic"))+
      theme( legend.key = element_rect(fill = NA),
             plot.margin = unit(ggplot.margins, "cm"),
             axis.text.x = element_text(size = text.x.size),
             axis.text.y = element_text(size = text.y.size), axis.title.y = element_text(size = text.y.title.size)
      )

   if(arrows)
   {
      graphy <- graphy +
         #coord_cartesian(clip="off") +
         annotation_custom(grob = linesGrob(arrow=arrow(type="open", ends="first", length=unit(arrow.length,"mm")),
                                            gp=gpar(col="grey22", lwd=arrow.lwd)), xmax = -arrow.origin.offset, xmin = min(tidata$LFC)/1.5,
                           ymin = 0, ymax = 0 )+
         annotation_custom(grob = linesGrob(arrow=arrow(type="open", ends="last", length=unit(arrow.length,"mm")),
                                            gp=gpar(col="grey22", lwd=arrow.lwd)), xmax = arrow.origin.offset, xmin = max(tidata$LFC)/1.5,
                           ymin = 0, ymax = 0 ) +
         annotation_custom(grob = grid::textGrob(label = levels(group)[1], hjust=1.05, gp=gpar(col=arrow.text.color, cex=arrow.text.cex)),
                           xmin = min(tidata$LFC)/1.5, xmax = min(tidata$LFC)/1.5, ymin = 0, ymax = 0) +
         annotation_custom(grob = grid::textGrob(label = levels(group)[2], hjust=-0.05, gp=gpar(col=arrow.text.color, cex=arrow.text.cex)),
                           xmin = max(tidata$LFC)/1.5, xmax = max(tidata$LFC)/1.5, ymin = 0, ymax = 0)
   }
   if ( align.legend ){ graphy <- align_legend(graphy) }
   gvec[[iter1]] <- graphy
   inc(iter1, increment = 1)
   # Corrected pvalues
   if(p.adjust.method != "none" )
   {
      print(paste("Generating volcano plot for corrected P-values as well..."))
      graphy <- ggplot(data = tidata, aes(x = LFC, y = -log2(corr), label = Plotlabel_corr, size = Relabb, color = Color_corr)) +
      geom_point(alpha = 0.8) +
      geom_vline(xintercept = 0, col = "grey80") +
      geom_hline(yintercept = -log2(0.05), col = "grey20", linetype = "dashed") +
      scale_color_identity() +
      labs(x="Log-Fold change", size="Rel. Ab.", y=paste("-log2(", p.adjust.method, ")", sep=""))+
      scale_y_continuous(limits = c(0, max(-log2(tidata$corr))+.5 ) )+
      scale_x_continuous(limits = c( min(tidata$LFC)-.5, max(tidata$LFC)+.5 ))+
      geom_text_repel(show_guide=F, size = label.text.size, point.padding = label.padding, aes(fontface = "italic"))+
      theme( legend.key = element_rect(fill = NA), plot.margin = unit(ggplot.margins, "cm"),
             axis.text.x = element_text(size = text.x.size),
             axis.text.y = element_text(size = text.y.size), axis.title.y = element_text(size = text.y.title.size)
             )

      if(arrows)
      {
         graphy <- graphy +
            #coord_cartesian(clip="off") +
            annotation_custom(grob = linesGrob(arrow=arrow(type="open", ends="first", length=unit(arrow.length,"mm")),
                                               gp=gpar(col="grey22", lwd=arrow.lwd)), xmax = -arrow.origin.offset, xmin = min(tidata$LFC)/1.5,
                              ymin = 0, ymax = 0 )+
            annotation_custom(grob = linesGrob(arrow=arrow(type="open", ends="last", length=unit(arrow.length,"mm")),
                                               gp=gpar(col="grey22", lwd=arrow.lwd)), xmax = arrow.origin.offset, xmin = max(tidata$LFC)/1.5,
                              ymin = 0, ymax = 0 ) +
            annotation_custom(grob = grid::textGrob(label = levels(group)[1], hjust=1.05, gp=gpar(col=arrow.text.color, cex=arrow.text.cex)),
                              xmin = min(tidata$LFC)/1.5, xmax = min(tidata$LFC)/1.5, ymin = 0, ymax = 0) +
            annotation_custom(grob = grid::textGrob(label = levels(group)[2], hjust=-0.05, gp=gpar(col=arrow.text.color, cex=arrow.text.cex)),
                              xmin = max(tidata$LFC)/1.5, xmax = max(tidata$LFC)/1.5, ymin = 0, ymax = 0)
      }
      if ( align.legend ){ graphy <- align_legend(graphy) }
      gvec[[iter1]] <- graphy
   }
   # Get pvalues into ptab and save
   print(paste("Saving stats into ", save.path, "/", taxlevel, "_volcano.xlsx", sep=""))
   ptab[,1] = tidata$pval
   ptab[,2] = tidata$Group1
   ptab[,3] = sapply( as.data.frame(t(data[group==levels(group)[1]])) , sem)
   ptab[,4] = tidata$Group2
   ptab[,5] = sapply( as.data.frame(t(data[group==levels(group)[2]])) , sem)
   addSheet(path = paste(save.path, "/", taxlevel, "_volcano.xlsx", sep=""), sheet.name = "uncorrected",addition = ptab,
            col.save = TRUE, row.save = TRUE)
   if( p.adjust.method != "none")
   {
      ptab_corr[,1] = tidata$corr
      ptab_corr[,2] = tidata$Group1
      ptab_corr[,3] = sapply( as.data.frame(t(data[group==levels(group)[1]])) , sem)
      ptab_corr[,4] = tidata$Group2
      ptab_corr[,5] = sapply( as.data.frame(t(data[group==levels(group)[2]])) , sem)
      addSheet(path = paste(save.path, "/", taxlevel, "_volcano.xlsx", sep=""), sheet.name = p.adjust.method,addition = ptab_corr,
               col.save = TRUE, row.save = TRUE)
   }
   # Save plots and excels
   print(paste("Saving plots into ", save.path, "/", taxlevel, "_volcano.xlsx", sep=""))
   suppressWarnings(ggsave(
      filename = paste(save.path, "/", taxlevel, "_volcano.pdf", sep=""),
      plot = marrangeGrob(gvec, nrow=nrow.graph, ncol=ncol.graph, top=NULL ,
                          layout_matrix = matrix(seq_len(nrow.graph*ncol.graph),
                                                 nrow = nrow.graph,
                                                 ncol = ncol.graph , byrow = TRUE) ),
      width = width.graph, height = height.graph, dpi = 330))
}




