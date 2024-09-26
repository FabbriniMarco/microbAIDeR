\name{compute_volcano_plots}
\alias{compute_volcano_plots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Volcano Plots with Wilcoxon Test Significance Between Two Groups
}
\description{
This function handles non-tidy data.frames (e.g. OTU tables) and produces a scatter volcano plot, marking differences between the observation found in the colnames of the input data.frame, grouped by a two-level grouping factor. \cr
The tendency of each taxa is rendered in Log-Fold Change, considering both average and median values of the distribution, with customizable control by the user via the \code{w1} and \code{w2}  parameter according to the LFC formula:\cr\cr
\code{
\verb{            }(w1 * median(group1) + (1 - w1) * mean(group1))\cr
LFC = log ( ----------------------------------------------    )\cr
\verb{            }(w2 * median(group2) + (1 - \code{w2} ) * mean(group2))
}
}
\usage{
compute_volcano_plots(data, group, taxlevel, save.path = getwd(), p.adjust.method = "fdr", 
    w1 = 0.5, w2 = 0.5, trends=TRUE, sigcolors=c("ivory3", "orange3", "darkred"),
    nrow.graph = 1, ncol.graph = 2, width.graph = 20, height.graph = 8, 
    ggplot.margins = c(1, 1, 1, .5), label.text.size = 2.8, label.padding=4,
    arrows=TRUE, arrow.length=1.5,arrow.lwd=1.8,arrow.origin.offset=0.5,
    arrow.text.cex=.7,arrow.text.color="grey22",
    text.x.size = 8, text.y.size = 7,  text.y.title.size = 9,
    additional.params = NULL, align.legend = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A \bold{data.frame} feature table (OTU table) with parameters in the rows (e.g., taxa) and observations (i.e., samples) in the columns. For each parameter in the rows, a volcano plot will be plotted, and pairwise Wilcoxon test differences assessed between grouped observations.
}
  \item{group}{
A \bold{factor} of \emph{length(group) == ncol(data)}. The level order would be the one appearing in the X-axis of the plotted volcano plots.
}
  \item{taxlevel}{
A prefix to append to file names, to pinpoint for exampl at which level the analysis was computed (e.g., "L2", "L5", "L6") or a name of choice. \cr
Output file names start with this text. Use shorts names to avoid problems with the length of file and excel sheet names.
}
  \item{save.path}{
A full \bold{PATH}, or the expansion of the current working directory into the folder where all the produced files will be saved (Default = current working directory).
}
  \item{p.adjust.method}{
Method to adjust pairwise Wilcoxon test P-values (can be "none", "fdr", "bonferroni", etc... check p.adjust.methods to see all the possible methods to be used)\cr
P-values are computed for each feature and the correction is applied to all P-values all together.
}
  \item{w1}{
A \bold{numeric} value setting the \code{w1} coefficient for the LFC calculation.
}
  \item{w2}{
A \bold{numeric} value setting the \code{w2}  coefficient for the LFC calculation.
}
  \item{trends}{
A \bold{logical} (Default = TRUE), specifying whether to consider also trends (P-value <= 0.1). In this way trends will be colored differently (Default = "orange3") from non significant points. If set to FALSE, the P-value threshold is set to 0.05 for significance (P-value <= 0.05) and trends are plotted as non significant colored according to \bold{sigcolors}.
}
  \item{sigcolors}{
A character vector with three colors used to represent non-significant points, trends, and significant points, respectively (Default = c("ivory3", "orange3", "darkred")). If \code{trends = FALSE}, only the first and third colors are used.
}
  \item{nrow.graph}{
Number of rows of the output .pdf file (Default = 1). Used to group the  volcano plots in the same sheet.
}
  \item{ncol.graph}{
Number of columns of the output .pdf file (Default = 2). Used to group the volcano plots in the same sheet.
}
  \item{width.graph}{
A \bold{numerical} value specifying the width of each volcano plot (Default = 20).
}
  \item{height.graph}{
A \bold{numerical} value specifying the height of each volcano plot (Default = 8).
}
  \item{ggplot.margins}{
Concatenation of four \bold{numeric} values passed to \emph{ggplot} for handling the margins of each volcano plot (Default = c(1, 1, 1, .5)).
}
  \item{label.text.size}{
A \bold{numerical} value setting the size of the labels text size for the significant (depending on \emph{trends}) data points (Default = 2.8).
}
  \item{label.padding}{
A \bold{numerical} value setting the point.padding of the geom_text_repel() function (Default = 4). (See \code{?geom_text_repel} for further details)
}
  \item{arrows}{
A \bold{logical} setting whether to plot or not the two lateral arrows pointing the direction of fold change towards the two groups (Default = TRUE).
}
  \item{arrow.length}{
If \code{arrows == TRUE} the length (size) of the arrowhead (Default = 1.5).
}
  \item{arrow.lwd}{
If \code{arrows == TRUE} the line width of the arrow and arrowhead (Default = 1.8).
}
  \item{arrow.origin.offset}{
If \code{arrows == TRUE} the horizontal offset of the arrows in respect to the \code{x=0} position (Default = 0.5).
}
  \item{arrow.text.cex}{
If \code{arrows == TRUE} the size of the text at the extremity of the arrows (Default = 0.7).
}
  \item{arrow.text.color}{
If \code{arrows == TRUE} the color of the text at the extremity of the arrows (Default = "grey22").
}
  \item{text.x.size}{
A \bold{numeric} value setting X-axis units font size for volcano plots (Default = 8).
}
  \item{text.y.size}{
A \bold{numeric} value setting Y-axis units font size for volcano plots (Default = 7).
}
  \item{text.y.title.size}{
A \bold{numeric} value setting title font size of Y axis title for volcano plots (Default = 9).
}
  \item{additional.params}{
An object of class \bold{gg} (a ggplot object, function, or parameter) to be passed at the end of each ggplot code in the script. \cr
This is useful to override parameters in the plot as well as bypass further parameters as panel background, themes presets, etc.. (Default = NULL). \cr
Example #1: additional.params = theme(panel.background = element_rect(fill="ivory1"))\cr
Example #2: additional.params = theme_bw()+theme(legend.position="none")
  }
  \item{align.legend}{
A \bold{logical} stating whether to center align the legend title to the legend elements in all volcano plots. Only applicable if \code{legend.position}!="none" is passed via \bold{additional.params}, otherwise legend is not plotted as default (Default = FALSE).
  }
}
\details{
This function takes as input only \bold{data.frame} objects, so consider testing the \code{class()} of your data before running the script. \cr
Mandatory items are: \bold{data}, \bold{group}, \bold{taxlevel}, \cr
The resulting output in the \bold{save.path} folder consist of one \emph{.pdf} file and a \emph{.xlsx} file.\cr
\bold{w1 and \code{w2}  coefficients}\cr
Adjusting the values of the \code{w1} and \code{w2} coefficients allows control over the influence of the median and mean in the LFC calculation.\cr
\itemize{
  \item \code{w1 = 1} or \code{w2 = 1}: gives full weight to the median, exclusively using the median of the group’s values in the calculation.
  \item \code{w1 = 0} or \code{w2 = 0}: gives full weight to the mean, exclusively using the mean.
  \item Intermediate values (e.g., \code{w1 = 0.5} or \code{w2 = 0.5}): balance the influence of both the median and the mean, assigning equal importance to these two measures of central tendency.
}
In datasets with a large number of zero values (as is common in microbiome data), the median may provide a more robust measure, as it is less influenced by outliers and zeros.\cr
The mean tends to smooth the data but may be biased in zero-inflated datasets where certain taxa have high relative abundances, skewing the average.\cr
For skewed data or when outliers are present, setting \code{w1} and \code{w2} closer to 1 (favoring the median) may lead to more reliable results. For normally distributed data, favoring the mean by setting \code{w1} and \code{w2} closer to 0 may be appropriate.\cr

}
\value{
\item{`'taxlevel'_volcano.pdf`}{
A PDF file at 330dpi containing volcano plots arranged in a matrix of \code{nrow.graph} rows and \code{ncol.graph} columns.
}
\item{`'taxlevel'_volcano.xlsx`}{
An Excel file with three sheets:\cr
1. \bold{call} sheet: contains the parameters passed to the function for reproducibility.\cr
2. \bold{uncorrected}: Contains uncorrected pairwise Wilcoxon test p-values along with average and SEM abundances of each taxa.\cr
3. \bold{corrected}: Same as the second sheet but with \code{p.adjust.method} corrected p-values.
}
}
\author{
Marco Fabbrini, Gabriele Conti
}
\examples{
env <- read.delim("example data/env.tsv", header=T, row.names = 1, sep='\t')
env$Group = factor(env$Group, levels = c("Group 1", "Group 2"))

genera <- read.delim("example data/genera_example.tsv", header=T, row.names = 1, sep='\t')
all( colnames(genera) == rownames(env) )

compute_volcano_plots(data = genera,
                      group = env$Group,
                      taxlevel = "L6", 
                      save.path = "Results/Compositional")
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
