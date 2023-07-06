# microbAIDeR
An ensemble of function for easier anc quicker preliminary microbiome analyses

microbAIDeR is a R-package made for easy preliminar microbiome analyses starting from QIIME2 outputs.

# Installation

```R
install.packages('devtools')
devtools::install_github("https://github.com/FabbriniMarco/microbAIDeR.git")
library(microbAIDeR)
microbAIDeR_install_dependancies() #To install all the packages that are necessary
```

Depends: devtools, pairwiseAdonis, parallel, doParallel, openxlsx, tidyr, tibble, dplyr, ggplot2, ggsignif, gridExtra, tidyverse, stringr, vegan, cowplot, gtable


# EXAMPLES


# Alpha diversity boxplots

```R
alpha <- read.delim("alpha.tsv", header = TRUE, row.names = 1)
all( rownames(alpha) == rownames(map_file) )
grouping = factor( map_file$grouping_factor, levels = c("Control", "T0", "T1") )
color_grouping = c("violet", "red", "green") #One color for each group

if(!dir.exists("Alphadiv")){dir.create("Alphadiv")}

compute_wilcoxon_and_plot( data = as.data.frame(t(alpha)), 
                           plot.not.sig = TRUE, 
                           p.adjust.method = "fdr", 
                           group = grouping, 
                           color.grouping = color_grouping,
                           taxlevel = "Alpha", 
                           save.path = "Alphadiv", 
                           signif.step.increase = 0.15, 
                           signif.text.size = 3 )
```

# Composition boxplots
```R
L2 = read.delim("L2_cleaned_table.tsv", header=T, row.names=1, sep='\t')
# L2 is the otu table, filtered and with cleared rownams. Taxa on the rows, samples on the columns 
all( colnames(L2) %in% rownames(map_file) )

if(!dir.exists("Composition")){dir.create("Composition")}

compute_wilcoxon_and_plot( data = L2, 
                           p.adjust.method = "fdr", 
                           group = grouping, 
                           plot.not.sig = TRUE,
                           taxlevel = "L2", 
                           save.path = "Composition", 
                           color.grouping = color_grouping )
# It works at L2 as at other levels, just input an otu table with clear taxa names in the rownames, a grouping factor matching the samples in the colnames and a color vector for each group
```


# Beta diversity PcoA 
```R
if(!dir.exists("Betadiv")){dir.create("Betadiv")}
compute_beta_diversity( beta.folder.path="beta_diversity", 
                        save.path="Betadiv", 
                        beta_metrics = c("unweighted_unifrac", "weighted_unifrac", "braycurtis", "jaccard"),
                        ellipses = TRUE, 
                        ellipse.focus = TRUE, 
                        ellipse.fill = TRUE, 
                        ellipse.alpha = 40,
                        ellipse.conf = .95, 
                        color.grouping = color_grouping, 
                        group = grouping, 
                        adonis_n_perm = 9999,
                        sample.list = rownames(map_file) )
```
Feel free to explore the other functions, some of them are very useful!
