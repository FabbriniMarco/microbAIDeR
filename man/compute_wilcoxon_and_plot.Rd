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
compute_wilcoxon_and_plot(data, group, taxlevel,save.path = getwd(), color.grouping, \cr
    p.adjust.method = "fdr", plot.not.sig = FALSE, \cr
    nrow.graph = 2, ncol.graph = 2, width.graph = 6, height.graph = 4, \cr
    horiz = FALSE, ggplot.margins = c(0.18, 0.18, 0.18, 0.6), box.lwd = 0.4, \cr
    jitter.pch = 21, jitter.stroke = 0.15, jitter.size = 0.7, jitter.color = "grey22", \cr
    signif.step.increase = 0.12, signif.text.size = 3.5, signif.line.size = 0.4, \cr
    text.x.size = 8, text.y.size = 7, text.y.title.size = 9, \cr
    smoothing = FALSE, smoothing.lwd = 1, smoothing.color = "darkred", \cr
    smoothing.se = FALSE, smoothing.method = "loess", contrast.color = "ivory1",
    additional.params = NULL, align.legend = FALSE)
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
Output file names start with this text.
}
  \item{save.path}{
A full \bold{PATH}, or the expansion of the current working directory into the folder where all the produced files will be saved (Default = current working directory).
}
  \item{color.grouping}{
A \bold{vector} of \emph{length(color.grouping) == nlevels(group)} specifying for each group a color for the boxplot fill (either color names or RGB HEX codes).
}
  \item{p.adjust.method}{
Method to adjust pairwise Wilcoxon test P-values (can be "none", "fdr", "bonferroni", etc... check p.adjust.methods to see all the possible methods to be used)\cr
P-values are first generated for each feature with the pairwise.wilcox.test() function without multiple comparison correction. Once \bold{all} the features have been tested, the P-values are corrected using the desired method considering \bold{all} the comparisons computed. \cr

Example: imagine feeding an OTU table with 30 taxa. For each taxa a pairwise.wilcox.test() is computed between all the groups specified with the group parameter (let's say 3 groups). Once all the possible pairwise comparisons are computed for all the 30 taxa (90 comparisons in total), the resulting P-values list is corrected for 90 items.
}
  \item{plot.not.sig}{
A \bold{logical} (Default = FALSE), specifying whether to plot also boxplots with no significant pairwise Wilcox test differences.
}
  \item{nrow.graph}{
Number of rows of the output .pdf file (Default = 2). Used to group multiple boxplots in the same sheet.
}
  \item{ncol.graph}{
Number of columns of the output .pdf file (Default = 2). Used to group multiple boxplots in the same sheet.
}
  \item{width.graph}{
A \bold{numerical} value specifying the width of each boxplot (Default = 6).
}
  \item{height.graph}{
A \bold{numerical} value specifying the height of each boxplot (Default = 4).
}
  \item{horiz}{
A \bold{logical} (Default = FALSE), specifying whether to plot horizontal boxplots.
}
  \item{ggplot.margins}{
Concatenation of four \bold{numeric} values passed to \emph{ggplot} for handling the margins of each boxplot (Default = c(0.15, 0.15, 0.15, 0.6)).
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
A \bold{numeric} value setting the font size of statistical significance symbols over sig bars (Default = 3.5).
}
  \item{signif.line.size}{
A \bold{numeric} value setting the width of the lines of the significance bracket (Default = 0.4).
}
  \item{text.x.size}{
A \bold{numeric} value setting X-axis units font size for boxplots (Default = 8).
}
  \item{text.y.size}{
A \bold{numeric} value setting Y-axis units font size for boxplots (Default = 7).
}
  \item{text.y.title.size}{
A \bold{numeric} value setting title font size of Y axis title for boxplots (Default = 9).
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
This is useful to override parameters in the plot as well as bypass further parameters as panel background, themes presets, etc.. (Defauult = NULL). \cr
Example #1: additional.params = theme(panel.background = element_rect(fill="ivory1"))\cr
Example #2: additional.params = theme_bw()+theme(legend.position="none")
  }
  \item{align.legend}{
A \bold{logical} stating whether to center align the legend title to the legend elements in all boxplots. Only applicable if \emph{legend.position}!="none" is passed via \bold{additional.params}, otherwise legend is not plotted as default (Default = FALSE).
  }
}
\details{
This function takes as input only \bold{data.frame} objects, so consider testing the \emph{class()} of your data before running the script. \cr
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
2. The second sheet takes the name from the \bold{taxlevel} parameter and contains all the \bold{uncorrecred} pairwise Wilcoxon test P-values.
3. The third sheet takes the name from the \bold{taxlevel} parameter and reports the \bold{p.adjust.method} as well. This sheet contains all the \bold{correcred} pairwise Wilcoxon test P-values (corrected for all the \emph{n} comparisons made within all the observations).
}
}
\author{
Marco Fabbrini - m.fabbrini@unibo.it
}
\examples{

# Make sure to check in which directory you are working
getwd()
[1] "C:/Users/.../my_working_directory/"

# If you want to store your output in a subfolder of your currently working directory, first create it
# For example:
if( !dir.exists("Plots")){ dir.create("Plots") }


########## Example 1 : Processing OTU tables

# Phylum2 is the OTU table at L2 (clean phyla names on the rows, samples on the columns)
# Group is a factor for grouping the samples so that:  length(group) == ncol(phylum2) & is.factor(group)
# Make sure to have the directories created for the save path before running the script
# the color_grouping vector must contain a matched color for each level of group:  length(color_grouping) == nlevels(group)
# MAKE SURE TO MATCH your otu_table_L?.tsv table's sample names with the sample names of yout env or map file used to generate the group factor variable!!!!

> compute_wilcoxon_and_plot( data = phylum2, p.adjust.method = "bonferroni", group = group,
                           taxlevel = "L2", save.path = "Plots", color.grouping = color_grouping )

[1] "Saving P-values to C:/Users/.../my_working_directory/Plots/L2_wilcoxon_pairwise.xlsx"
[1] "Saving plots to C:/Users/.../my_working_directory/Plots/L2_boxplot_wilcox_uncorrected.pdf"
[1] "Saving plots to C:/Users/.../my_working_directory/Plots/L2_boxplot_wilcox_fdr.pdf"


########## Example 2 : Processing alpha diversity tables

# alpha is the alpha diversity table, with samples on the rows and alpha metrics on the columns
# Each sample is assigned to a group, and the group variable is a factor:  nrow(alpha) == length(group) & is.factor(group)
# Each level of the group factor is assigned to a color through the color_grouping variable:  length(color_grouping) == nlevels(group)
# MAKE SURE TO MATCH your alpha.tsv table's sample names with the sample names of yout env or map file used to generate the group factor variable!!!!

compute_wilcoxon_and_plot( data = as.data.frame(t(alpha)), p.adjust.method = "bonferroni", group = group, color.grouping = color_grouping,
                           taxlevel = "alpha", save.path = "Plots", signif.step.increase = 0.15, signif.text.size = 3 )

[1] "Saving P-values to C:/Users/.../my_working_directory/Plots/alpha_wilcoxon_pairwise.xlsx"
[1] "Saving plots to C:/Users/.../my_working_directory/Plots/alpha_boxplot_wilcox_uncorrected.pdf"
[1] "Saving plots to C:/Users/.../my_working_directory/Plots/alpha_boxplot_wilcox_fdr.pdf"

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
