\name{correlation_custom}
\alias{correlation_custom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make correlation easier for continuous-vs-continuous distribution and continuous-vs-binary ones.
}
\description{
Starting from distribution of values easily compute correlation and filter retaining only significant comparisons. Enables p-value filtering, multiple comparison correction and minimum correlation coefficient filtering. The output includes comparisons and raw p-values, comparisons with adjusted p-values and comparisons with correlation coefficient over a certain threshold and the matching adjusted p-values.
}
\usage{
correlation_custom(feature_data, parameter_data, is.binary.parameter = FALSE, \cr
         p.adjust.method = "fdr", uncorr.p.threshold = 0.05, corr.p.threshold = 0.05, \cr
         min.corr.threshold = 0.3, force.numeric = TRUE, corr.method = "spearman")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{feature_data}{
A \bold{data.frame} object with the feature data, with samples on the rows and features on the columns. For example a L6 otu table.
}
  \item{parameter_data}{
A \bold{data.frame} object with the parameters data, with samples on the rows and parameters on the columns. For example a table of clinical blood parameters.
}
  \item{is.binary.parameter}{
A \bold{logical} determining whether the \emph{parameter_data} is binary. In case of continuous-vs-binary correlation a Point-Biserial correlation is implemented.
}
  \item{p.adjust.method}{
A \bold{character} specifying the correction method for p-values. Regardless of this parameter, you will always have the raw uncorrected p-values. Check \emph{p.adjust.methods} for the available methods that can be used. Default = "fdr".
}
  \item{uncorr.p.threshold}{
A \bold{numeric} value that represent the threshold for filtering p-values. Only features or parameters with at least one raw uncorrected p-values strictly below this threshold will be kept in the output. Default = 0.05.
}
  \item{corr.p.threshold}{
A \bold{numeric} value that represent the threshold for filtering corrected p-values. Only features or parameters with at least one multiple-comparison adjusted p-values strictly below this threshold will be kept in the corresponding output. Default = 0.05.
}
  \item{min.corr.threshold}{
A \bold{numeric} value that represent the threshold for filtering the correlation coefficient. Only features or parameters with at least one comparison with the absolute value of the correlation coefficient greater than this threshold will be kept in the corresponding output. This correlation coefficient filtering is applied only to the correlation table after adjusting p-values. Default = 0.3.
}
  \item{force.numeric}{
A \bold{logical} stating whether to force the input data to be numeric. Default = TRUE.
}
  \item{corr.method}{
A \bold{character} pointing the correlation method to be used. In case of binary data, the \bold{is.binary.parameter = TRUE} parameter will override any value supplied to this option and a Point-Biserial correlation will be used. If \bold{is.binary.parameter = FALSE} "pearson", "kendall", or "spearman" can be used.
}

}
\details{
This function first computes the correlations bewteen each feature and all the parameters. Then the paired resulting tables are filtered to remove those features and those parameters that did not show any significant (p < \emph{uncorr.p.threshold}) correlation. Then, the p-value table is adjusted for multiple comparison using \emph{p.adjust.method}. Both the correlation table and adjusted p-value table are then filtered again, removing features and parameters that did not show any significant (p < \emph{corr.p.threshold}) correlation. The last filter, on the correlation coefficient is applied on the resulting tables, discarding features and parameters that did not show any strong ( abs(cor) > \emph{min.corr.threshold}) correlation.
The output is returned as a list of 7 elements:
- A table of correlation coefficients between features and parameters that showed at least one significant (p < \emph{uncorr.p.threshold}) correlation;
- The corresponding uncorrected p-value table;
- A table of correlation coefficients between features and parameters that showed at least one significant (p < \emph{corr.p.threshold}) correlation;
- The corresponding adjusted p-value table;
- A table of correlation coefficients between features and parameters that showed at least one strong ( abs(cor) > \emph{min.corr.threshold}) correlation;
- The corresponding adjusted p-value table;
- A call argument reporting the parameters supplied to the function.
}
\value{
A \emph{list} object with 7 elements. In case the are no remaining features or parameters with passing-filter values, the elements will be less and a warning message is displayed.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Marco Fabbrini - m.fabbrini@unibo.it
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

rownames(microbiota_data) == rownames(parameter_data)

out_corr = correlation_custom( feature_data = microbiota_data, parameter_data = parameter_data ,
                     corr.method = "spearman", p.adjust.method = "fdr")

names(out_corr)
[1] "uncorr_correlation"      "uncorr_pvalues"          "corrected_correlation"   "corrected_pvalues"       "corr_strong_correlation" "corr_strong_pvalues"
[7] "call"

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
