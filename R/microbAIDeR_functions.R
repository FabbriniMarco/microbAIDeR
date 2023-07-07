microbAIDeR_install_dependancies <- function(){
   install.packages('devtools')
   devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
   install.packages("parallel")
   install.packages("doParallel")
   install.packages("openxlsx")
   install.packages("tidyr")
   install.packages("tibble")
   install.packages("dplyr")
   install.packages("ggplot2")
   install.packages("ggsignif")
   install.packages("gridExtra")
   install.packages("tidyverse")
   install.packages("stringr")
   install.packages("vegan")
   install.packages("cowplot")
   install.packages("gtable")
}


tokenize <- function(x) {
   if(!is.factor(x)){
      print(paste("Input object must be a factor"))
   }
   tok <- c()
   for (len in 1:length(x)) {
      tok = c(tok, which(levels(x) == x[len]) )
   }
   return(tok)
}


one_hot_encode <- function(x){
   if(!is.factor(x)){
      print(paste("Input object must be a factor"))
   }
   encoded = data.frame()
   for (len in 1:length(x)){
      encoded = rbind(encoded, rep(0, nlevels(x)))
      encoded[ nrow(encoded) , which(levels(x) == x[len]) ] = 1
   }
   colnames(encoded) = levels(x)
   return(encoded)
}


amp_up_models <- function(core_fraction = 0.8){
   require(parallel)
   require(doParallel)
   no_cores <- parallel::detectCores() * core_fraction
   #Leave one core available for Operating system
   cluster <- makePSOCKcluster(no_cores)
   registerDoParallel(cluster)
   cat("Model amped and ready to go with:", no_cores, "cores. \n")
}


tidify <- function(x, as.numeric = TRUE){
   require(tidyr)
   require(tibble)
   y <- x %>%
      tibble::rownames_to_column() %>%
      tidyr::gather(colname, value, -rowname)
   if(as.numeric){y$value = as.numeric(y$value)}
   return(y)
}


sigFunc = function(x){
   if ( !is.na(x)){
      if(x <= 0.0001){"****"}
      else if(x <= 0.001){"***"}
      else if(x <= 0.01){"**"}
      else if(x <= 0.05){"*"}
      else if(x <= 0.1){"°"}
      else{NA}
   } else {NA}
}


addSheet <- function(path,sheet.name, addition, col.save=TRUE, row.save=TRUE, overwrite = TRUE){
   require(openxlsx)
   wb <- loadWorkbook(path)
   if ( sheet.name %in% getSheetNames(path) & overwrite ){ removeWorksheet(wb, sheet.name) }
   addWorksheet(wb,sheet.name)
   writeData(wb,sheet.name, as.data.frame(addition), colNames = col.save, rowNames = row.save )
   saveWorkbook(wb,path,overwrite = TRUE)
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
                     insert = insert +1
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
               insert = insert +1
            }
         }
      }
   } else { significant_comparisons = NA }

   return(significant_comparisons)
}


compute_wilcoxon_and_plot <- function(data, p.adjust.method = "fdr", group, taxlevel, save.path = getwd(), color.grouping , plot.not.sig = FALSE,
                                      nrow.graph = 2, ncol.graph = 2, width.graph = 6, height.graph = 4, horiz = FALSE, ggplot.margins = c(.18, .18, .18, .6),
                                      box.lwd = 0.4, jitter.pch = 21, jitter.stroke = 0.15, jitter.size = 0.7, jitter.color = "grey22",
                                      signif.step.increase = 0.12, signif.text.size = 3.5, signif.line.size = 0.4, contrast.color = "ivory1",
                                      text.x.size = 8, text.y.size = 7,  text.y.title.size = 9,
                                      smoothing = FALSE, smoothing.lwd = 1, smoothing.color = "darkred", smoothing.se = FALSE, smoothing.method="loess",
                                      additional.params = NULL, align.legend = FALSE){
   require(dplyr)
   require(ggplot2)
   require(ggsignif)
   require(gridExtra)
   require(tidyverse)
   require(tidyr)
   require(openxlsx)
   require(pairwiseAdonis)
   require(stringr)

   if(!dir.exists(save.path)){dir.create(save.path)}

   if ( any(is.na(group))){
      data = data[,!is.na(group)]
      group = group[!is.na(group)]
   }

   data = data[!rowSums(data) %in% 0,]

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

   call.print = as.data.frame(rbind(p.adjust.method, taxlevel, save.path, color.grouping = paste(color.grouping, collapse = ", "), smoothing.color, contrast.color,
                                    nrow.graph , ncol.graph , width.graph , height.graph , ggplot2.margins = paste(ggplot.margins, collapse = ", ") ,
                                    box.lwd , jitter.pch , jitter.stroke , jitter.size , jitter.color ,
                                    signif.step.increase,  signif.text.size ,
                                    text.x.size , text.y.size ,  text.y.title.size, additional.params=paste(additional.params, collapse=", ")
   ))

   write.xlsx( call.print , paste(save.path, "/", taxlevel, "_wilcoxon_pairwise.xlsx", sep=""), sheetName = "call",
               colNames = FALSE, rowNames = TRUE)

   wilcox_uncorrected_list = vector("list", 0)

   # Compute all (pairwise)-Wilcoxon tests and plot uncorrected boxplots
   for (taxa in rownames(data) )
   {
      tempdf = as.data.frame(t(data)) %>% select(abundances = all_of(taxa)) %>% mutate(grouping_column = group )
      suppressWarnings(testing <-  pairwise.wilcox.test(tempdf$abundances, tempdf$grouping_column, p.adjust.method = "none"))
      wilcox_uncorrected_list[[taxa]] = testing
      onlysig = process_pairwise_wilcoxon(testing, onlysig = TRUE)
      wilcox_list = as.list(as.data.frame(t(str_split_fixed(onlysig$Comparison, " vs ",2))))
      ############ Uncorrected
      if ( !all(is.na(onlysig)) | plot.not.sig ) {
         suppressWarnings(graphy <- ggplot( data=tempdf, aes(x=grouping_column, y=as.numeric(abundances)) )+
                             geom_boxplot( aes(fill=grouping_column),outlier.shape = NA, lwd = 0.4) +
                             geom_jitter( aes(color=grouping_column), width = 0.1, height=0, pch = jitter.pch, color=jitter.color, stroke = jitter.stroke, size = jitter.size)+
                             geom_signif( comparisons = wilcox_list  ,
                                          map_signif_level=sigFunc,
                                          test = "wilcox.test", size = signif.line.size,
                                          tip_length = 0.02, color = "grey22", textsize = signif.text.size, margin_top = 0.05, vjust = .5, step_increase = signif.step.increase)+
                             {if (smoothing)
                                geom_smooth(aes(x=tokenize(grouping_column)), method=smoothing.method, lwd=smoothing.lwd, color=contrast.color, se=smoothing.se, formula = 'y ~ x' )
                             }+
                             {if (smoothing)
                                geom_smooth(aes(x=tokenize(grouping_column)), method=smoothing.method, lwd=(0.5*smoothing.lwd), color=smoothing.color, se=smoothing.se, formula = 'y ~ x' )
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
         iter1 = iter1 + 1
      }
   }

   wilcox_processed_uncorrected = lapply(wilcox_uncorrected_list, process_pairwise_wilcoxon, onlysig = FALSE)
   extracted_pvalues <- data.frame(
      name = rep( names(wilcox_processed_uncorrected), as.numeric(lapply(wilcox_processed_uncorrected, nrow)) ),
      Comparison = unlist(lapply(wilcox_processed_uncorrected, "[[", "Comparison")) ,
      Wilcox.pvalue = unlist(lapply(wilcox_processed_uncorrected, "[[", "Wilcox.pvalue")),
      row.names=NULL
   )
   uncorrected_p_list = extracted_pvalues

   # Gather and correct all pvalues, if requested
   if ( p.adjust.method != 'none' )
   {
      extracted_pvalues$Wilcox.pvalue = p.adjust(extracted_pvalues$Wilcox.pvalue, method = p.adjust.method)
      wilcox_corrected_list = split(extracted_pvalues[,c(2,3)], f=extracted_pvalues$name)

      for (taxa in rownames(data))
      {
         tempdf = as.data.frame(t(data)) %>% select(abundances = all_of(taxa)) %>% mutate(grouping_column = group )
         testing_corrected = wilcox_corrected_list[taxa][[1]]
         if ( any(testing_corrected$Wilcox.pvalue<0.05, na.rm=TRUE) ){
            keep = testing_corrected$Wilcox.pvalue<0.05
            keep[is.na(keep)]=FALSE
            onlysig_corr = testing_corrected[keep,]
            wilcox_list = as.list(as.data.frame(t(str_split_fixed(onlysig_corr$Comparison, " vs ",2))))
         } else {
            onlysig_corr = data.frame(Comparison=NA, Wilcox.pvalue=NA)
            wilcox_list = NA}

         #
         if ( !all(is.na(onlysig_corr)) | plot.not.sig ) {
            suppressWarnings(graphy <- ggplot( data=tempdf, aes(x=grouping_column, y=as.numeric(abundances)) )+
                                geom_boxplot( aes(fill=grouping_column),outlier.shape = NA, lwd = 0.4) +
                                geom_jitter( aes(color=grouping_column), width = 0.1, height=0, pch = jitter.pch, color=jitter.color, stroke = jitter.stroke, size = jitter.size)+
                                geom_signif( comparisons = wilcox_list  ,
                                             annotations = sapply(onlysig_corr$Wilcox.pvalue, sigFunc), size = signif.line.size,
                                             tip_length = 0.02, color = "grey22", textsize = signif.text.size, margin_top = 0.05, vjust = .5, step_increase = signif.step.increase) +
                                {if (smoothing)
                                   geom_smooth(aes(x=tokenize(grouping_column)), method=smoothing.method, lwd=smoothing.lwd, color=contrast.color, se=smoothing.se, formula = 'y ~ x' )
                                }+
                                {if (smoothing)
                                   geom_smooth(aes(x=tokenize(grouping_column)), method=smoothing.method, lwd=(0.5*smoothing.lwd), color=smoothing.color, se=smoothing.se, formula = 'y ~ x' )
                                }+
                                scale_fill_manual(values=color.grouping)+
                                labs(x="", y=taxa ) +
                                {if (horiz)
                                   coord_flip()
                                }+
                                theme(legend.position = 'none', plot.margin = unit(ggplot.margins, "cm"),
                                      axis.text.x = element_text(size = text.x.size),
                                      axis.text.y = element_text(size = text.y.size), axis.title.y = element_text(size = text.y.title.size)
                                ) +
                                { if( !is.null(additional.params)) additional.params}
            )
            if ( align.legend ){ graphy <- align_legend(graphy) }
            gvec_corr[[iter2]] = graphy
            iter2 = iter2 + 1
         }
      }
   }


   print(paste("Saving P-values to ", save.path, "/", taxlevel, "_wilcoxon_pairwise.xlsx", sep=""))
   addSheet( path=paste(save.path, "/", taxlevel, "_wilcoxon_pairwise.xlsx", sep=""), sheet.name = paste(taxlevel),
             addition = uncorrected_p_list, col.save=T, row.save=FALSE)
   if ( p.adjust.method != 'none' )
   {addSheet( path=paste(save.path, "/", taxlevel, "_wilcoxon_pairwise.xlsx", sep=""), sheet.name = paste(taxlevel, p.adjust.method, sep="_"),
              addition = extracted_pvalues, col.save=T, row.save=FALSE) }

   ### Saving uncorrected
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

   ### Saving corrected
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
                                   cex.plot = 1.3, adonis.p.adjust = "fdr", wilcoxon.p.adjust = "fdr",

                                   nrow.graph = 2, ncol.graph = 2, width.graph = 6, height.graph = 4, ggplot.margins = c(0.15, 0.15, 0.15, 0.6),
                                   box.lwd = 0.4, jitter.pch = 21, jitter.stroke = 0.15, jitter.size = 0.7, jitter.color = "grey22",
                                   signif.step.increase = 0.12, signif.text.size = 3.5, signif.line.size = 0.4, horiz = FALSE,
                                   text.x.size = 8, text.y.size = 7,  text.title.size = 9, plot.not.sig = FALSE,
                                   smoothing = FALSE, smoothing.lwd = 1, smoothing.color = "darkred", contrast.color="ivory1", smoothing.se = FALSE, smoothing.method = "loess",
                                   additional.params=NULL, align.legend = FALSE){

   require(vegan)
   require(gridExtra)
   require(pairwiseAdonis)
   require(stringr)

   if ( any(is.na(group))){print(paste("WARNING: You have NAs in your group, removing such samples from the beta distance matrices. BE SURE not to have NAs in the sample.list"))}
   if ( length(group) != length(sample.list)){
      print(paste("WARNING: The grouping factor provided has a different length compared to the sample.list"))
      stop()
   }

   if(!dir.exists(save.path)){dir.create(save.path)}

   call.print = as.data.frame(rbind( paste(beta_metrics, collapse=", "), color.grouping = paste(color.grouping, collapse = ", "),
                                     save.path, adonis_n_perm , beta.folder.path ,  mds = paste(mds, collapse=", ") ,
                                     spiders, ellipses , ellipse.focus , ellipse.conf , ellipse.alpha , smoothing.color, contrast.color,
                                     svg.width , svg.height ,
                                     cex.plot , adonis.p.adjust, wilcoxon.p.adjust, additional.params=paste(additional.params, collapse=", ")
   ))

   write.xlsx( call.print , paste(save.path, "/beta_adonis.xlsx", sep=""), sheetName = "call",
               colNames = FALSE, rowNames = TRUE)
   write.xlsx( call.print , paste(save.path, "/beta_wilcoxon_variance.xlsx", sep=""), sheetName = "call",
               colNames = FALSE, rowNames = TRUE)

   iter1 = 1
   iter2 = 1

   gvec <- vector("list",length=0 )
   gvec_corr <- vector("list",length=0 )

   for ( metrics in beta_metrics)
   {
      print(paste("Running", toupper(metrics), "beta diversity analyses"))
      if (is.null(manual.beta.path))
      {
         beta = read.delim( paste(beta.folder.path , metrics, "distance-matrix.tsv", sep="/"), header=T, row.names=1, sep="\t")
      } else if( !is.null(manual.beta.path))
      {
         beta = read.delim( manual.beta.path, header=T, row.names=1, sep="\t")
      }
      beta = beta[ sample.list, sample.list ]
      if ( any(is.na(group)))
      {
         beta = beta[!is.na(group), !is.na(group)]
         group_beta = group[!is.na(group)]
      } else { group_beta = group }
      metric_beta <- capscale(as.dist(beta)~1)
      coord_beta = as.data.frame( scores(metric_beta, display="sites", choices=mds ) )
      eighenvalues = metric_beta$CA$eig
      sig1 = round(as.numeric( eighenvalues[mds[1]]*100/sum(eighenvalues) ), 2)
      sig2 = round(as.numeric( eighenvalues[mds[2]]*100/sum(eighenvalues) ), 2)
      color_beta = rep(NA, length(group_beta))
      for (x in 1:length(group_beta)){color_beta[x] = color.grouping[which(levels(group_beta) == group_beta[x] )]}
      if(any(is.na(color_beta))){
         print(paste("Errors in producing the color vector. Check that your group is a factor and nlevels(group) == length(color.grouping)"))
         color_beta[is.na(color_beta)] = "black"
      }
      colnames(coord_beta) = c("AX1", "AX2")
      # Base plot
      svg( paste(save.path,"/beta_", metrics, ".svg", sep=""), width=svg.width, height = svg.height)
      plot(metric_beta, type="n", main=metrics, choices=mds, xlab=paste("MDS", mds[1]," - [" ,sig1, "%]", sep="") , ylab=paste("MDS", mds[2]," - [" ,sig2, "%]", sep="") )
      with(coord_beta, points(AX1, AX2, pch=21, bg=color_beta, cex=cex.plot, col=NULL))
      legend("topright", legend=levels(group_beta), pt.cex=1, fill=color.grouping, horiz=F, x.intersp = 0.4, y.intersp = 0.7, bty="n", cex=0.8)
      graphics.off()
      #Spiders
      if(spiders){
         svg( paste(save.path, "/beta_", metrics, "_spiders.svg", sep=""), width=svg.width, height = svg.height)
         plot(metric_beta, type="n", main=metrics, choices=mds, xlab=paste("MDS", mds[1]," - [" ,sig1, "%]", sep="") , ylab=paste("MDS", mds[2]," - [" ,sig2, "%]", sep="") )
         with(coord_beta, points(AX1, AX2, pch=21, bg=color_beta, cex=cex.plot, col=NULL))
         for ( groupiter in levels(group_beta))
         {
            ordispider(metric_beta, group_beta, display="sites", spiders="centroid", show.groups=groupiter, col= if ( is.null(manual.bordercol) ){ color.grouping[which(levels(group_beta) == groupiter)] } else { manual.bordercol } , lwd=spiders.lwd, lty=1)
         }
         legend("topright", legend=levels(group_beta), pt.cex=1, fill=color.grouping, horiz=F, x.intersp = 0.4, y.intersp = 0.7, bty="n", cex=0.8)
         graphics.off()
      }
      # Ellipse plot
      if(ellipses){
         svg( paste(save.path, "/beta_", metrics, "_ellipse.svg", sep=""), width=svg.width, height = svg.height)
         plot(metric_beta, type="n", main=metrics, choices=mds, xlab=paste("MDS", mds[1]," - [" ,sig1, "%]", sep="") , ylab=paste("MDS", mds[2]," - [" ,sig2, "%]", sep="") )
         with(coord_beta, points(AX1, AX2, pch=21, bg=color_beta, cex=cex.plot, col=NULL))
         for ( groupiter in levels(group_beta))
         {
            ordiellipse(metric_beta, group_beta, display="sites", kind=c("se"), conf=ellipse.conf, draw=c("lines"), show.groups=groupiter, col= if ( is.null(manual.bordercol) ){ color.grouping[which(levels(group_beta) == groupiter)] } else { manual.bordercol } , lwd=ellipse.lwd, lty=1, alpha=50)
         }
         legend("topright", legend=levels(group_beta), pt.cex=1, fill=color.grouping, horiz=F, x.intersp = 0.4, y.intersp = 0.7, bty="n", cex=0.8)
         graphics.off()
      }
      # Ellipse focus
      if (ellipse.focus){
         svg( paste(save.path, "/beta_", metrics, "_ellipse_focus.svg", sep=""), width=svg.width, height = svg.height)
         plot(metric_beta, type="n", main=metrics, choices=mds, xlab=paste("MDS", mds[1]," - [" ,sig1, "%]", sep="") , ylab=paste("MDS", mds[2]," - [" ,sig2, "%]", sep="") )
         with(coord_beta, points(AX1, AX2, pch=21, bg=color_beta, cex=(cex.plot/4), col=NULL))
         for ( groupiter in levels(group_beta))
         {
            ordiellipse(metric_beta, group_beta, display="sites", kind=c("se"), conf=ellipse.conf, draw=c("polygon"), show.groups=groupiter, col= color.grouping[which(levels(group_beta) == groupiter)], border = if ( is.null(manual.bordercol) ){ color.grouping[which(levels(group_beta) == groupiter)] } else { manual.bordercol }, lwd=ellipse.lwd, lty=1, alpha=ellipse.alpha)
         }
         legend("topright", legend=levels(group_beta), pt.cex=1, fill=color.grouping, horiz=F, x.intersp = 0.4, y.intersp = 0.7, bty="n", cex=0.8)
         graphics.off()
      }
      # Ellipse fill
      if (ellipse.fill){
         svg( paste(save.path, "/beta_", metrics, "_ellipse_fill.svg", sep=""), width=svg.width, height = svg.height)
         plot(metric_beta, type="n", main=metrics, choices=mds, xlab=paste("MDS", mds[1]," - [" ,sig1, "%]", sep="") , ylab=paste("MDS", mds[2]," - [" ,sig2, "%]", sep="") )
         with(coord_beta, points(AX1, AX2, pch=21, bg=color_beta, cex=cex.plot, col=NULL))
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
      adone$uncorrected_pvalues <- sapply(adone$p.value, sigFunc)
      adone[,adonis.p.adjust] <- sapply(adone$p.adjusted, sigFunc)

      addSheet( path=paste(save.path, "/beta_adonis.xlsx", sep=""), sheet.name = metrics,
                addition = adone[,c(1,6,9, 7,10)], col.save=TRUE, row.save=FALSE)

      # Plot the intra-group inter-individual variance in terms of beta with boxplots
      print(paste("Computing intra-individual variance and Wilcoxon tests"))
      beta_half = beta
      beta_half[lower.tri(beta_half, diag=T)] <- NA
      beta_half = as.data.frame(beta_half)
      values_beta_var <- c()
      group_beta_var <- c()
      for (z in levels(group_beta))
      {
         if ( sum(group_beta==z) > 1 )
         {  tempvalues = as.numeric(unlist(beta_half[ group_beta==z, group_beta==z]))  }
         else { tempvalues = 0 }
         if( any(is.nan(tempvalues)) ) { tempvalues = tempvalues[ - is.nan(tempvalues) ] }
         values_beta_var = c(values_beta_var, tempvalues )
         group_beta_var = c(group_beta_var, rep(z, length(tempvalues)) )
      }
      plot_beta_variance = as.data.frame(cbind(values_beta_var, group_beta_var))
      plot_beta_variance$values_beta_var = as.numeric( plot_beta_variance$values_beta_var )
      plot_beta_variance$group_beta_var = factor(plot_beta_variance$group_beta_var, levels=levels(group_beta))
      ## Plot variance
      ############ Uncorrected
      suppressWarnings(testing <-  pairwise.wilcox.test(plot_beta_variance$values_beta_var, plot_beta_variance$group_beta_var, p.adjust.method = "none"))
      onlysig = process_pairwise_wilcoxon(testing, onlysig = TRUE)
      wilcox_list = as.list(as.data.frame(t(str_split_fixed(onlysig$Comparison, " vs ",2))))
      if ( !all(is.na(onlysig)) | plot.not.sig ) {
         suppressWarnings(graphy <- ggplot( data=plot_beta_variance, aes(x=group_beta_var, y=as.numeric(values_beta_var )) )+
                             geom_boxplot( aes(fill=group_beta_var),outlier.shape = NA, lwd = 0.4) +
                             geom_jitter( aes(color=group_beta_var), width = 0.1, height=0, pch = jitter.pch, color=jitter.color, stroke = jitter.stroke, size = jitter.size)+
                             geom_signif( comparisons = wilcox_list  ,
                                          map_signif_level=sigFunc,
                                          test = "wilcox.test",size = signif.line.size,
                                          tip_length = 0.02, color = "grey22", textsize = signif.text.size, margin_top = 0.05, vjust = .5, step_increase = signif.step.increase)+
                             {if (smoothing)
                                geom_smooth(aes(x=tokenize(group_beta_var)), method=smoothing.method, lwd=smoothing.lwd, color=contrast.color, se=smoothing.se, formula = 'y ~ x' )
                             }+
                             {if (smoothing)
                                geom_smooth(aes(x=tokenize(group_beta_var)), method=smoothing.method, lwd=(0.5*smoothing.lwd), color=smoothing.color, se=smoothing.se, formula = 'y ~ x' )
                             }+
                             scale_fill_manual(values=color.grouping)+
                             labs(x="", y=metrics ) +
                             {if (horiz)
                                coord_flip()
                             }+
                             theme(legend.position = 'none', plot.margin = unit(ggplot.margins, "cm"),
                                   axis.text.x = element_text(size = text.x.size),
                                   axis.text.y = element_text(size = text.y.size), axis.title.y = element_text(size = text.title.size)
                             ) +
                             { if( !is.null(additional.params)) additional.params}
         )
         if ( align.legend ){ graphy <- align_legend(graphy) }
         gvec[[iter1]] = graphy
         iter1 = iter1 + 1
         onlysig$Wilcox.pvalue = signif(onlysig$Wilcox.pvalue, 5)
         addSheet( path=paste(save.path, "/beta_wilcoxon_variance.xlsx", sep=""), sheet.name = paste( gsub("weighted", "w", metrics), "unc"),
                   addition = onlysig, col.save=TRUE, row.save=FALSE)
         addSheet( path=paste(save.path, "/beta_wilcoxon_variance.xlsx", sep=""), sheet.name = paste( gsub("weighted", "w", metrics), "unc_all", sep="_"),
                   addition = process_pairwise_wilcoxon(testing, onlysig = FALSE), col.save=TRUE, row.save=FALSE)
      }
      ############# Corrected pvals
      suppressWarnings(testing <-  pairwise.wilcox.test(plot_beta_variance$values_beta_var, plot_beta_variance$group_beta_var, p.adjust.method = wilcoxon.p.adjust))
      onlysig_corr = process_pairwise_wilcoxon(testing)
      wilcox_list = as.list(as.data.frame(t(str_split_fixed(onlysig_corr$Comparison, " vs ",2))))
      if ( !all(is.na(onlysig_corr)) | plot.not.sig ) {
         suppressWarnings(graphy <- ggplot( data=plot_beta_variance, aes(x=group_beta_var, y=as.numeric(values_beta_var )) )+
                             geom_boxplot( aes(fill=group_beta_var),outlier.shape = NA, lwd = 0.4) +
                             geom_jitter( aes(color=group_beta_var), width = 0.1, height=0, pch = jitter.pch, color=jitter.color, stroke = jitter.stroke, size = jitter.size)+
                             geom_signif( comparisons = wilcox_list  ,
                                          map_signif_level=sigFunc,
                                          test = "wilcox.test",size = signif.line.size,
                                          tip_length = 0.02, color = "grey22", textsize = signif.text.size, margin_top = 0.05, vjust = .5, step_increase = signif.step.increase)+
                             {if (smoothing)
                                geom_smooth(aes(x=tokenize(group_beta_var)), method=smoothing.method, lwd=smoothing.lwd, color=contrast.color, se=smoothing.se, formula = 'y ~ x' )
                             }+
                             {if (smoothing)
                                geom_smooth(aes(x=tokenize(group_beta_var)), method=smoothing.method, lwd=(0.5*smoothing.lwd), color=smoothing.color, se=smoothing.se, formula = 'y ~ x' )
                             }+
                             scale_fill_manual(values=color.grouping)+
                             labs(x="", y=metrics ) +
                             {if (horiz)
                                coord_flip()
                             }+
                             theme(legend.position = 'none', plot.margin = unit(ggplot.margins, "cm"),
                                   axis.text.x = element_text(size = text.x.size),
                                   axis.text.y = element_text(size = text.y.size), axis.title.y = element_text(size = text.title.size)
                             ) +
                             { if( !is.null(additional.params)) additional.params}
         )
         if ( align.legend ){ graphy <- align_legend(graphy) }
         gvec_corr[[iter2]] = graphy
         iter2 = iter2 + 1
         onlysig_corr$Wilcox.pvalue = signif(onlysig_corr$Wilcox.pvalue, 5)
         addSheet( path=paste(save.path, "/beta_wilcoxon_variance.xlsx", sep=""), sheet.name = paste( gsub("weighted", "w", metrics), gsub("bonferroni", "bonf", wilcoxon.p.adjust), sep="_"),
                   addition = onlysig_corr, col.save=TRUE, row.save=FALSE)
         addSheet( path=paste(save.path, "/beta_wilcoxon_variance.xlsx", sep=""), sheet.name = paste( gsub("weighted", "w", metrics) , gsub("bonferroni", "bonf", wilcoxon.p.adjust), "all", sep="_"),
                   addition = process_pairwise_wilcoxon(testing, onlysig = FALSE), col.save=TRUE, row.save=FALSE)
      }
   }

   ### Save uncorrected
   if (smoothing){
      print(paste("Saving plots to ", save.path, "/beta_intragroup_variance_smooth.pdf", sep=""))
      suppressWarnings(ggsave(
         filename = paste(save.path, "/beta_intragroup_variance_smooth.pdf", sep="") ,
         plot = marrangeGrob(gvec, nrow=nrow.graph, ncol=ncol.graph, top=NULL ,
                             layout_matrix = matrix(seq_len(nrow.graph*ncol.graph), nrow = nrow.graph, ncol = ncol.graph , byrow = TRUE) ),
         width = width.graph, height = height.graph, dpi = 330
      ))
   } else {
      print(paste("Saving plots to ", save.path, "/beta_intragroup_variance.pdf", sep=""))
      suppressWarnings(ggsave(
         filename = paste(save.path, "/beta_intragroup_variance.pdf", sep="") ,
         plot = marrangeGrob(gvec, nrow=nrow.graph, ncol=ncol.graph, top=NULL ,
                             layout_matrix = matrix(seq_len(nrow.graph*ncol.graph), nrow = nrow.graph, ncol = ncol.graph , byrow = TRUE) ),
         width = width.graph, height = height.graph, dpi = 330
      ))
   }
   ### Save corrected
   if (smoothing){
      print(paste("Saving plots to ", save.path, "/beta_", wilcoxon.p.adjust, "_intragroup_variance_smooth.pdf", sep=""))
      suppressWarnings(ggsave(
         filename = paste(save.path, "/beta_", wilcoxon.p.adjust, "_intragroup_variance_smooth.pdf", sep=""),
         plot = marrangeGrob(gvec_corr, nrow=nrow.graph, ncol=ncol.graph, top=NULL ,
                             layout_matrix = matrix(seq_len(nrow.graph*ncol.graph), nrow = nrow.graph, ncol = ncol.graph , byrow = TRUE) ),
         width = width.graph, height = height.graph, dpi = 330
      ))
   } else {
      print(paste("Saving plots to ", save.path, "/beta_", wilcoxon.p.adjust, "_intragroup_variance.pdf", sep=""))
      suppressWarnings(ggsave(
         filename = paste(save.path, "/beta_", wilcoxon.p.adjust, "_intragroup_variance.pdf", sep=""),
         plot = marrangeGrob(gvec_corr, nrow=nrow.graph, ncol=ncol.graph, top=NULL ,
                             layout_matrix = matrix(seq_len(nrow.graph*ncol.graph), nrow = nrow.graph, ncol = ncol.graph , byrow = TRUE) ),
         width = width.graph, height = height.graph, dpi = 330
      ))
   }

}


firstup <- function(string) {
   substr(string, 1, 1) <- toupper(substr(string, 1, 1))
   return(string)
}


normalize <- function(x) {
   return ((x - min(x)) / (max(x) - min(x))) }


inc <- function(x, increment = 1) {
   eval.parent(substitute(x <- x + increment)) }


align_legend <- function(p, hjust = 0.5)
{
   require(cowplot)
   require(gtable)
   # extract legend
   g <- cowplot::plot_to_gtable(p)
   grobs <- g$grobs
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


se.avg <- function(x){
   x = x[!is.na(x)]
   sd(x) / sqrt(length(x))
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
   require(dplyr)
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

