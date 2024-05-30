\name{compute_beta_diversity}
\alias{compute_beta_diversity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quickly compute MDS analyses starting from beta diversity matrices output of the QIIME2 pipeline (ideally qiime2_multipipeline)
}
\description{
This function automatically computes beta diversity analyses performind MDS, plotting PCoA, testing inter-group differences with pairwise-PERMANOVA and intra-group variance.
}
\usage{
compute_beta_diversity(beta.folder.path = paste(getwd(), "RESULTS/beta_diversity", sep="/"), \cr
        manual.beta.path = NULL, save.path = getwd(),  \cr
        group, color.grouping, sample.list, \cr
        beta_metrics = c("braycurtis", "jaccard", "unweighted_unifrac", "weighted_unifrac"), \cr
        mds = c(1, 2), adonis_n_perm = 9999, \cr
        adonis.p.adjust = "fdr", wilcoxon.p.adjust = "fdr" \cr
        plot.not.sig = FALSE, alpha.points=1 \cr
        spiders = FALSE, spiders.lwd = 1.5, ellipses = FALSE, ellipse.focus = FALSE, ellipse.fill = FALSE, \cr
        ellipse.lwd = 2.5, ellipse.conf = 0.95, ellipse.alpha = 0.75, manual.bordercol = NULL,\cr
        svg.width = 7, svg.height = 5, \cr
        cex.points = 1.2, width.graph = 6, height.graph = 4, \cr
        ggplot.margins = c(0.1, 0.1, 0.1, 0.5), \cr
        nrow.graph = 2, ncol.graph = 2, \cr
        text.x.size = 8, text.y.size = 7, text.title.size = 9, \cr
        additional.params = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta.folder.path}{
The directory where the beta diversity matrices are stores. \cr

Please note that the script requires inside the \bold{beta.folder.path} a folder for each metric used, named according to the \bold{beta_metrics} parameter and inside of each folder a \emph{distance-matrix.tsv} file for each metric. \cr

Please also note that this hierarchy is the default configuration of the \emph{qiime2_multipipeline} script outputs. By default, the \bold{beta.folder.path} is set to  ".RESULTS/beta_analyses" in the current working directory. If you downloaded the RESULTS folder from the \emph{qiime2_multipipeline} outputs in your current working directory, you can avoid specifying this parameter.
}
  \item{manual.beta.path}{
A \bold{character} string pointing to the PATH to a TSV distance matrix to use \bold{instead} of the \bold{beta.folder.path}. The matrix must be delimited by tabs and have samples in the first row and first column. The function is going to read this table importing header as \emph{colnames} and the first column as \emph{rownames}. Remind to \bold{adjust} the beta_metrics accordingly, specifying a single metrics.
}
  \item{save.path}{
A \bold{character} string pointing to the output folder where all the graphs and output files will be stored (Default = current working directory).
}
  \item{group}{
A \bold{factor} of \emph{length(group) == length(sample.list)}. The level order would be the one appearing in the X-axis of the plotted barplots for intra-group distances.
}
  \item{color.grouping}{
A \bold{vector} of \emph{length(color.grouping) == nlevels(group)} specifying for each group a color for the PCoA points fill and intra-group beta barplot fill (either color names or RGB HEX codes).
}
  \item{sample.list}{
A \bold{character} vector containing the sample names that we want to analyze. \cr

The order of the samples must match the meaning of the group factor, as it will be used to group samples in categories.\cr

Note that providing a smaller \bold{sample.list} compared to the samples used for the generation beta diversity matrix would simply result in the analysis of the subset of the given samples. This can be useful to perform sub-group beta diversities from the same beta matrix.
}
  \item{beta_metrics}{
A \bold{character} vector of beta diversity metrics that have been computed. These names \bold{must} match exactly the folder names of each metric, contained in the \bold{beta.folder.path} folder.This is the default configuration of the beta diversity output of the qiime2_multipipeline. If you used the qiime2_multipipeline command, you can leave unaltered this parameter (Default = c("braycurtis", "jaccard", "unweighted_unifrac", "weighted_unifrac"))
}
  \item{mds}{
A length == 2 \bold{numeric} vector containing the number of the principal component axes we want to use for the analyses (Default = c(1,2), meaning that "MDS1" and "MDS2" axes will be considered)
}
  \item{adonis_n_perm}{
A \bold{numeric} value of permutations for each of the pairwise Adonis tests that will be performed (Default = 9999).
}
  \item{adonis.p.adjust}{
A \bold{character} value pointing to the P-value correction method for the pairwise Adonis testing (Default = "fdr"). Can be any value from p.adjust.methods.\cr

Please note that in case of a two-group beta diversity analysis, the reported corrected P-values would be the same of the uncorrected ones.
}
  \item{wilcoxon.p.adjust}{
A \bold{character} value pointing to the P-value correction method for the pairwise Wilcoxon testing for the intra-group variance (Default = "fdr"). Can be any value from p.adjust.methods.\cr
}
  \item{alpha.points}{
A \bold{numeric} value handling points fill transparency in plots (Default = 1 - fully opaque) \cr
}
  \item{spiders}{
A \bold{logical} (Default = FALSE), specifying whether to produce an additional PCoA plot with spiders.
}
  \item{spiders.lwd}{
A \bold{numeric} (Default = 1.5), specifying the line width for PCoA spiders.
}
  \item{ellipses}{
A \bold{logical} (Default = FALSE), specifying whether to produce an additional PCoA plot with ellipses (just the border, no fill). Note that ellipses are computed with the "se" method at \bold{ellipse.conf} confidence level.
}
  \item{ellipse.focus}{
A \bold{logical} (Default = FALSE), specifying whether to produce an additional PCoA plot with filled ellipses and much smaller points (ideal for huge datasets (> 800 points). Note that ellipses are computed with the "se" method at \bold{ellipse.conf} confidence level.
}
  \item{ellipse.fill}{
A \bold{logical} (Default = FALSE), specifying whether to produce an additional PCoA plot with filled ellipses and normal-size data points. Note that ellipses are computed with the "se" method at \bold{ellipse.conf} confidence level.
}
  \item{ellipse.lwd}{
A \bold{numeric} (Default = 1.5), specifying the line width for PCoA ellipses and ellipses.focus (not affecting ellipse.fill).
}
  \item{ellipse.conf}{
A \bold{numeric} value passing ellipse confidence level for the eventual ellipses plot (Default = 0.95).
}
  \item{ellipse.alpha}{
A \bold{numeric} value handling ellipse fill transparency in plots (Default = 0.75).
}
  \item{manual.bordercol}{
If NULL, ellipse border colors are assigned from the \emph{color.grouping} variable. If a single \bold{character} value is specified, all ellipses border are colored accordingly (Default = NULL).
}
  \item{svg.width}{
\bold{numeric} value of the width for output PCoA plots in scalable vector graphics format (Default = 7).
}
  \item{svg.height}{
\bold{numeric} value of the height for output PCoA plots in scalable vector graphics format (Default = 5).
}
  \item{cex.points}{
A \bold{numeric} value passing the size of PCoA data points.
}
  \item{width.graph}{
A \bold{numerical} value specifying the width of intra-group beta barplots (Default = 6).
}
  \item{height.graph}{
A \bold{numerical} value specifying the height of intra-group beta barplots (Default = 4).
}
  \item{ggplot.margins}{
Concatenation of four \bold{numeric} values passed to \emph{ggplot} for handling the margins of each intra-group beta barplot (Default = c(0.1, 0.1, 0.1, 0.5)).
}
  \item{nrow.graph}{
Number of rows of the output .pdf file (Default = 2). Used to group multiple barplots in the same sheet.
}
  \item{ncol.graph}{
Number of columns of the output .pdf file (Default = 2). Used to group multiple barplots in the same sheet.
}
  \item{text.x.size}{
A \bold{numeric} value setting X-axis font size for barplots (Default = 8).
}
  \item{text.y.size}{
A \bold{numeric} value setting Y-axis font size for barplots (Default = 7).
}
  \item{text.title.size}{
A \bold{numeric} value setting title font size for barplots (Default = 9).
}
  \item{additional.params}{
An object of class \bold{gg} (a ggplot object, function, or parameter) to be passed at the end of each ggplot code in the script. \cr
This is useful to override parameters in the plot as well as bypass further parameters as panel background, themes presets, etc.. (Default = NULL). \cr
Example: additional.params = theme(panel.background = element_rect(fill="ivory1"))
  }
}
\details{
If additional distance metrics were computed by the user, the script can perform the analyses with no problems. Simply place the beta diversity squared matrix named \emph{distance-matrix.tsv} inside a folder named according to the string passed to \bold{beta_metrics}.

}
\value{
The script produces up to 4 PCoA SVG plots for each metrics, as well as one PDF files containing intra-group beta barplots, an Excel file with all pairwise Adonis P-values and an excel file with intra-group variance P-values and IQR values:

\item{'beta_metrics'*.svg }{
For each metric found in the \bold{beta.folder.path} and passed to the \bold{beta_metrics} parameter, at least a PCoA plot is produced in SVG format. For each TRUE logical passed to one of the \bold{ellipses}, \bold{ellipses.focus} or \bold{ellipses.fill} parameters, an additional SVG plot is produced and named accordingly.
}
\item{beta_intragroup_variance.pdf }{
A 330dpi PDF file with intra-group beta distances barplots arranged in a \bold{nrow.graph} row and \bold{ncol.graph} column matrix \cr
Three barplots are computed for each metrics (variance, IQR, Range).
}
\item{beta_adonis.xlsx }{
An Excel file with several sheets:\cr
1. The \bold{call} sheet contains the parameters passed to the function. It's useful for future edits on the code to remember which colours have been used in the plot, as well as the other graphical parameters and correction method.\cr
\emph{n}. Each of the other sheets takes the name of a metric and reports the Adonis P-values for each pairwise comparison. For multi-group analyses, P-values are adjusted according to the number of pairwise comparisons computed for each metric individually.
}
\item{beta_variance.xlsx }{
An Excel file with several sheets:\cr
1. The \bold{call} sheet contains the parameters passed to the function. It's useful for future edits on the code to remember which colours have been used in the plot, as well as the other graphical parameters and correction method.\cr
2. \emph{metric}-Var. Contains the Variance, IQR and Range values plotted in the barplots
3. \emph{metric}-p. Contains the Wilcoxon and Kruskal P-values of the comparison between the group's distances distributions.
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



# In the beta_diversity folder there is a folder for each metric
## named "braycurtis", "jaccard", "unweighted_unifrac", "weighted_unifrac"
# Inside of each folder there must be a distance-matrix.tsv file
# This is the default configuration for the qiime2_multipipeline script output, so just point to the RESULT/beta_diversity folder for the beta.folder.path

# The sample.list variable needs to contain the sample names, matched with the group variable!!!

compute_beta_diversity( beta.folder.path="C:/Users/.../my_working_directory/RESULTS/beta_diversity", save.path="Plots",
                        ellipses = TRUE,
                        ellipse.focus = TRUE,
                        ellipse.conf = .99,
                        color.grouping = color_grouping,
                        group = group,
                        sample.list = mapfile$Samplecodes)


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
