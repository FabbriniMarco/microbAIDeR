\name{compute_wilcoxon_and_plot}
\alias{compute_wilcoxon_and_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quickly generate boxplots with tested pairwise Wilcoxon test significance between groups starting from a data.frame object
}
\description{
This function handles non-tidy data.frames (e.g. OTU tables or alpha diversity tables) and produces one boxplot per each parameters found in the rownames, testing the pairwise wilcoxon differences between the observation found in the colnames of the input data.frame, grouped by a grouping factor
}
\usage{
compute_wilcoxon_and_plot(data, group, taxlevel,save.path = getwd(), color.grouping, comparison.list=NULL,
    p.adjust.method = "fdr", trends = TRUE, plot.not.sig = FALSE, paired=FALSE, 
    nrow.graph = 2, ncol.graph = 2, width.graph = 4.5, height.graph = 3.5, 
    horiz = FALSE, ggplot.margins = c(0.18, 0.18, 0.18, 0.6), box.lwd = 0.4, 
    jitter.pch = 21, jitter.stroke = 0.15, jitter.size = 0.7, jitter.color = "grey22", 
    signif.step.increase = 0.12, signif.text.size = 3, signif.line.size = 0.4, 
    text.x.size = 6, text.y.size = 6, text.y.title.size = 8, 
    smoothing = FALSE, smoothing.lwd = 1, smoothing.color = "darkred", 
    smoothing.se = FALSE, smoothing.method = "loess", contrast.color = "ivory1",
    additional.params = NULL, align.legend = FALSE, 
    plot.order = "kruskal", pattern.fill = FALSE, pattern = "stripe",
    pattern.angle = 45, pattern.alpha = 0.4, pattern.density = 0.1, pattern.spacing = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A \bold{data.frame} feature table (OTU table or transposed alpha diversity table) with parameters in the rows (e.g., taxa or alpha diversity measures) and observations (i.e., samples) in the columns. For each parameter in the rows, a boxplot will be plotted, and pairwise Wilcoxon test differences assessed between grouped observations.
}
  \item{group}{
A \bold{factor} of \emph{length(group) == ncol(data)}. The level order would be the one appearing in the X-axis of the plotted boxplots.
}
  \item{taxlevel}{
A prefix to append to file names, to pinpoint for exampl at which level the analysis was computed (e.g., "L2", "L5", "L6") or a name of choice. \cr
Output file names start with this text. Using short names is recommended to avoid potential issues with file name or sheet name length.
}
  \item{save.path}{
A full \bold{PATH}, or the expansion of the current working directory into the folder where all the produced files will be saved (Default = current working directory).
}
  \item{color.grouping}{
A \bold{vector} of \emph{length(color.grouping) == nlevels(group)} specifying for each group a color for the boxplot fill (either color names or RGB HEX codes).
}
  \item{comparison.list}{
A \bold{list} where each element consists of two characters corresponding to two levels of the grouping factor. This list specifies the custom pairwise comparisons that the function will execute instead of evaluating all possible pairwise combinations of the grouping factor's levels. For instance, if the grouping factor has levels A, B, C, and D, a comparison list \code{list(c("A", "B"), c("C", "D"))} would instruct the function to only compare A to B and C to D, ignoring other combinations. If not specified, the default value is NULL, and the function will test all possible pairwise combinations of the levels of the factor supplied to the \emph{group} argument.
}
  \item{p.adjust.method}{
Method to adjust pairwise Wilcoxon test P-values (can be "none", "fdr", "bonferroni", etc... check \emph{p.adjust.methods} to see all the possible methods to be used)\cr
P-values are first generated for feature iterating through pairwise comparison with the wilcox.test(). Once \bold{all} the features have been tested, the P-values are corrected for each comparison using the desired method. \cr
}
  \item{trends}{
A \bold{logical} (Default = TRUE), specifying whether to represent in significance bars and output tables also trends (P-value <= 0.1) with the "°" symbol. If set to FALSE, the P-value threshold is set to 0.05 for significance (P-value <= 0.05) and trends are not taken in consideration.
}
  \item{plot.not.sig}{
A \bold{logical} (Default = FALSE), specifying whether to plot also boxplots with no significant pairwise Wilcox test differences.
}
  \item{paired}{
A \bold{logical} (Default = FALSE), specifying whether to perform paired Wilcoxon tests (signed-rank tests) in case of longitudinal studies with paired observations (e.g. matched pairs or before/after samples, or repeated measures).
}
  \item{nrow.graph}{
Number of rows of the output .pdf file (Default = 2). Used to group multiple boxplots in the same sheet.
}
  \item{ncol.graph}{
Number of columns of the output .pdf file (Default = 2). Used to group multiple boxplots in the same sheet.
}
  \item{width.graph}{
A \bold{numerical} value specifying the width of each boxplot (Default = 4.5).
}
  \item{height.graph}{
A \bold{numerical} value specifying the height of each boxplot (Default = 3.5).
}
  \item{horiz}{
A \bold{logical} (Default = FALSE), specifying whether to plot horizontal boxplots.
}
  \item{ggplot.margins}{
Concatenation of four \bold{numeric} values passed to \emph{ggplot} for handling the margins of each boxplot (Default = c(0.18, 0.18, 0.18, 0.6)).
}
  \item{box.lwd}{
A \bold{numeric} value to adjust line width of boxplots (Default = 0.4).
}
  \item{jitter.pch}{
A \bold{numeric} value setting the shape of the jittered points over the boxplot (Default = 21).
}
  \item{jitter.stroke}{
A \bold{numeric} value setting the width of each jitter point border (Default = 0.15). Set both this and \emph{jitter.size} to 0 in order to avoid plotting jitter over boxplots.
}
  \item{jitter.size}{
A \bold{numeric} value setting the diameter of each jitter point (Default = 0.7). Set both this and \emph{jitter.stroke} to 0 in order to avoid plotting jitter over boxplots.
}
  \item{jitter.color}{
A \bold{character} specifying the color of each jitter border (Default = "grey22")
}
  \item{signif.step.increase}{
A \bold{numeric} value handling the distance between each sig bar plotted by ggsignif (Default = 0.12).
}
  \item{signif.text.size}{
A \bold{numeric} value setting the font size of statistical significance symbols over sig bars (Default = 3).
}
  \item{signif.line.size}{
A \bold{numeric} value setting the width of the lines of the significance bracket (Default = 0.4).
}
  \item{text.x.size}{
A \bold{numeric} value setting X-axis units font size for boxplots (Default = 6).
}
  \item{text.y.size}{
A \bold{numeric} value setting Y-axis units font size for boxplots (Default = 6).
}
  \item{text.y.title.size}{
A \bold{numeric} value setting title font size of Y axis title for boxplots (Default = 8).
}
  \item{smoothing}{
A \bold{logical} describing whether to plot the smoothed conditional means with \emph{geom_smooth()} over the boxplots (Default = FALSE).
}
  \item{smoothing.lwd}{
If \emph{smoothing = TRUE}, this set the line width of the smoothing line (Default = 1).
}
  \item{smoothing.color}{
If \emph{smoothing = TRUE},this set the line color of the smoothing line (Default = "darkred").
}
  \item{smoothing.se}{
If \emph{smoothing = TRUE}, this \bold{logical} specify whether to plot the confidence interval shadow (Default = FALSE).
}
  \item{smoothing.method}{
If \emph{smoothing = TRUE}, this parameter handles the method used to compute the smoothing (Default = "loess").\cr Accepts a character vector, e.g. "lm", "glm", "gam", "loess" or a function, e.g. MASS::rlm or mgcv::gam, stats::lm, or stats::loess
}
  \item{contrast.color}{
If \emph{smoothing = TRUE}, this \bold{character} specify the color used to contrast the borders of \bold{smoothing.color} (Default = "ivory1").
}
  \item{additional.params}{
An object of class \bold{gg} (a ggplot object, function, or parameter) to be passed at the end of each ggplot code in the script. \cr
This is useful to override parameters in the plot as well as bypass further parameters as panel background, themes presets, etc.. (Default = NULL). \cr
Example #1: additional.params = theme(panel.background = element_rect(fill="ivory1"))\cr
Example #2: additional.params = theme_bw()+theme(legend.position="none")
  }
  \item{align.legend}{
A \bold{logical} stating whether to center align the legend title to the legend elements in all boxplots. Only applicable if \code{legend.position}!="none" is passed via \bold{additional.params}, otherwise legend is not plotted as default (Default = FALSE).
  }
  \item{pattern.fill}{
A \bold{logical} (Default = FALSE), specifying whether to fill the boxplot with patterns instead of continuous colors. If set to TRUE, the boxplot will be filled using the ggpattern package, which allows for a variety of pattern styles. This is useful for visual differentiation, especially in black-and-white prints or when color is not sufficient.
}
  \item{pattern}{
A \bold{character} specifying the type of pattern to be used in the ggpattern filling of the boxplots (Default = "stripe"). Some available patterns include "stripe", "crosshatch", "circle", and "none". This parameter determines the visual style of the patterns inside the boxplots when \code{pattern.fill = TRUE}.
}
  \item{pattern.angle}{
A \bold{numeric} value specifying the angle of the pattern lines (Default = 45). This is particularly relevant for line-based patterns like "stripe" or "crosshatch" and can help differentiate boxplots when multiple patterns are used.
}
  \item{pattern.density}{
A \bold{numeric} value (Default = 0.1) that controls the density of the pattern, determining how closely the lines or shapes are spaced in the boxplot pattern fill. A higher value results in more closely spaced patterns.
}
  \item{pattern.spacing}{
A \bold{numeric} value (Default = 0.05) that controls the spacing between pattern elements, affecting how much space is between the individual lines or shapes in the pattern. This is particularly useful for fine-tuning the appearance of patterns to achieve visual clarity in the plot.
}
  \item{plot.order}{
A \bold{character} (Default = "kruskal"), specifying the order in which boxplots are arranged in the output PDF. This parameter can be set to "kruskal" to arrange boxplots based on the significance of the Kruskal-Wallis test, where the most significant differences between groups appear first, or "rownames" to maintain the order of the rows in the input \bold{data} object.
}

}
\details{
This function is designed for visualizing group differences using boxplots, with statistical significance indicated by significance bars derived from pairwise Wilcoxon tests. The function generates both visual and tabular outputs that are useful for preliminary exploratory analysis of group differences in microbiome or other biological datasets.\cr

This function takes as input only \bold{data.frame} objects, so consider testing the \code{class()} of your data before running the script. \cr
Mandatory items are: \bold{data}, \bold{group}, \bold{taxlevel}, \bold{color.grouping}\cr
The resulting output in the \bold{save.path} folder consist of two \emph{.pdf} files and a \emph{.xlsx} file.\cr
Note that this function performs all the possible pairwise Wilcox test comparisons between groups. If you have a lot of timepoints and feel like some comparisons are not noteworthy or biologically senseless, consider this analysis only as preliminar and evaluate your group comparisons manually and carefully.
}
\value{
\item{'taxlevel'_boxplot_wilcox_uncorrected.pdf }{
A 330dpi PDF file with boxplots arranged in a \bold{nrow.graph} row and \bold{ncol.graph} column matrix \cr
A single boxplot is created for each row of the starting \bold{data}, and samples are grouped according to \bold{group} levels.
If \bold{plot.not.sig = FALSE} only boxplots with \emph{at least} one significant pairwise Wilcoxon comparison considering uncorrected P-values are plotted.
}
\item{'taxlevel'_boxplot_wilcox_'p.adjust.method'.pdf }{
A 330dpi PDF file with boxplots arranged in a \bold{nrow.graph} row and \bold{ncol.graph} column matrix \cr
A single boxplot is created for each row of the starting \bold{data}, and samples are grouped according to \bold{group} levels.
If \bold{plot.not.sig = FALSE} only boxplots with \emph{at least} one significant pairwise Wilcoxon comparison according to corrected P-values are plotted.
}
\item{'taxlevel'_wilcoxon_pairwise.xlsx }{
An Excel file with three sheets:\cr
1. The \bold{call} sheet contains the parameters passed to the function. It's useful for future edits on the code to remember which colours have been used in the plot, as well as the other graphical parameters and correction method.\cr
2. The second sheet takes the name from the \bold{taxlevel} parameter and contains all the \bold{uncorrected} pairwise Wilcoxon test P-values, as well as the average and SEM abundances of each taxa in each group.
3. Same as the second sheed but with blank cells for the non-significant P-values according to \emph{trends} .
4. Same as the second, but with \bold{p.adjust.method} corrected P-values. The correction is made for each pairwise comparison individually accounting for \emph{n} multiple comparisons, where \emph{n} equalt the number of taxa supplied. Consider therefore that submitting an input data with very few rows will provide quite light P-values adjustments.
5. Same as the fourth, but again with blank cells where the P-values are not significant.
}
}
\author{
Marco Fabbrini - m.fabbrini@unibo.it
}
\examples{
env <- read.delim("env.tsv", header=T, row.names = 1, sep='\t')
env$Group = factor(env$Group, levels = c("Group 1", "Group 2"))
genera <- read.delim("genera_example.tsv", header=T, row.names = 1, sep='\t')
all( colnames(genera) == rownames(env) )

compute_wilcoxon_and_plot( data = genera, 
                           group = env$Group, 
                           p.adjust.method = "fdr", 
                           taxlevel = "L6", save.path = "Results/Compositional",
                           color.grouping = c("navajowhite", "salmon"), plot.not.sig = TRUE )
# Setting plot.not.sig = TRUE to include also non-significant boxplots in the output
                           
# If you want to experiment with more fancy plots

compute_wilcoxon_and_plot( data = genera, 
                           group = env$Group, 
                           p.adjust.method = "fdr", 
                           taxlevel = "L6_custom", save.path = "Results/Compositional",
                           color.grouping = c("navajowhite", "salmon"), plot.not.sig = FALSE,
                           nrow.graph = 2, ncol.graph = 3, jitter.pch = 24,
                           additional.params = theme(panel.background = element_blank(), 
                                                     axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                                                     axis.title.y = element_text(size = 5, face = "italic"),
                                                     panel.grid.major = element_line(colour = "grey90", linewidth = 0.38), 
                                                     legend.position = 'right', legend.key.size = unit(0.2, "cm"), 
                                                     legend.text = element_text(size = 5), legend.title = element_blank(), 
                                                     legend.background = element_blank(), legend.box.background = element_blank(),
                                                     panel.grid.minor = element_line(colour = "grey95", linewidth = 0.08) )
)
                           
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
