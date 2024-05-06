\name{compute_volcano_plots}
\alias{compute_volcano_plots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quickly generate volcano plots with tested pairwise Wilcoxon test significance between two groups (healthy/disease or treatment/control) starting from a data.frame object
}
\description{
This function handles non-tidy data.frames (e.g. OTU tables) and produces a scatter volcano plot, marking differences between the observation found in the colnames of the input data.frame, grouped by a two-level grouping factor
}
\usage{
compute_volcano_plots(data, group, taxlevel, save.path = getwd(), p.adjust.method = "fdr", \cr
    trends=TRUE, sigcolors=c("ivory3", "orange3", "darkred"),\cr
    paired = FALSE, nrow.graph = 1, ncol.graph = 2, width.graph = 20, height.graph = 8, \cr
    ggplot.margins = c(1, 1, 1, .5), label.text.size = 2.8, label.padding=4,\cr
    arrows=TRUE, arrow.length=1.5,arrow.lwd=1.8,arrow.origin.offset=0.5,\cr
    arrow.text.cex=.7,arrow.text.color="grey22",\cr
    text.x.size = 8, text.y.size = 7,  text.y.title.size = 9,\cr
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
  \item{trends}{
A \bold{logical} (Default = TRUE), specifying whether to consider also trends (P-value <= 0.1). In this way trends will be colored differently (Default = "orange3") from non significant points. If set to FALSE, the P-value threshold is set to 0.05 for significance (P-value <= 0.05) and trends are plotted as non significant colored according to \bold{sigcolors}.
}
  \item{sigcolors}{
A three element \emph{character} vector containing the colors for the non significant, trends and significant points (Default = c("ivory3", "orange3", "darkred")). If trends is set to FALSE, only the first and third colors will be used, however the vector must always contain three colors.
}
  \item{paired}{
A \bold{logical} (Default = FALSE), specifying whether to perform paired Wilcoxon tests (signed-rank tests) in case of paired observations (e.g. matched pairs or before/after samples, or repeated measures).
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
If \code{arrows == TRUE} the size of the text at the extremity of the arrows (Default = 7).
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
}
\value{
\item{'taxlevel'_volcano.pdf }{
A 330dpi PDF file with volcano plots arranged in a \bold{nrow.graph} row and \bold{ncol.graph} column matrix \cr
A single volcano plot is always created with uncorrected Wilcoxon P-values. Dots are colored according to P-value significance and the diameter is representative of the average relative abundance of each feature in the entire data.
If \code{p.adjust.method != "none"} an additional volcano plot is created with P-values corrected with the specified method.
}
\item{'taxlevel'_volcano.xlsx }{
An Excel file with three sheets:\cr
1. The \bold{call} sheet contains the parameters passed to the function. It's useful for future edits on the code to remember which colours have been used in the plot, as well as the other graphical parameters and correction method.\cr
2. The second sheet contains all the \bold{uncorrecred} pairwise Wilcoxon test P-values, as well as the average and SEM abundances of each taxa in each group.
4. Same as the second, but with \bold{p.adjust.method} corrected P-values.
}
}
\author{
Gabriele Conti - gabriele.conti12@unibo.it ; Marco Fabbrini - m.fabbrini@unibo.it
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

> compute_volcano_plots( data = phylum2, p.adjust.method = "bonferroni", group = group,
                           taxlevel = "L2", save.path = "Plots" )

[1] "Terraforming volcano plots.."
[1] "Generating volcano plot for corrected P-values as well..."
[1] "Saving stats into C:/Users/.../my_working_directory/Plots/L2_volcano.xlsx"
[1] "Saving plots into C:/Users/.../my_working_directory/Plots/L2_volcano.pdf"



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
