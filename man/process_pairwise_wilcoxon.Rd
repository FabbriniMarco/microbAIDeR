\name{process_pairwise_wilcoxon}
\alias{process_pairwise_wilcoxon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract pairwise comparisons and P-values from the output of pairwise.wilcox.test() function
}
\description{
This function process a \bold{pairwise.htest} object resulting from the pairwise.wilcox.test() function into a 2-column \bold{data.frame} object reporting a single pairwise comparison for each line and the corresponding P-value.
}
\usage{
process_pairwise_wilcoxon(testing, onlysig = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{testing}{
A \bold{pairwise.htest} object
}
  \item{onlysig}{
A \bold{logical} (Default = TRUE) stating whether to keep only significant pairwise comparisons and trends (P-values <= 0.1) in the resulting data.frame. \cr
Setting this option to FALSE would result in the final data.frame having all the unique pairwise comparisons contained in the pairwise.htest object
}
}

\author{
Marco Fabbrini - m.fabbrini@unibo.it
}

\examples{
Consider having the variable a storing the result of the pairwise.wilcox.test() function

> a = pairwise.wilcox.test( as.numeric(table[1,]), group, method="bonferroni")

> class(a)
[1] "pairwise.htest"

> print(a)

  Pairwise comparisons using Wilcoxon rank sum test with continuity correction

data:  as.numeric(table[1, ]) and group

         Uninfected Mild    Moderate Severe
Mild     0.6230     -       -        -
Moderate 1.1e-08    < 2e-16 -        -
Severe   1.0e-08    < 2e-16 0.6230   -
Dead     0.0032     5.6e-07 0.2111   0.1266

P value adjustment method: holm

# Use the function only to retrieve significant comparisons.
> process_pairwise_wilcoxon(a, onlysig = TRUE)
              Comparison Wilcox.pvalue
1 Moderate vs Uninfected  1.121082e-08
2       Moderate vs Mild  4.277167e-22
3   Severe vs Uninfected  1.030190e-08
4         Severe vs Mild  1.105551e-21
5     Dead vs Uninfected  3.168674e-03
6           Dead vs Mild  5.629542e-07


# Keep also non-significant comparisons
> process_pairwise_wilcoxon(a, onlysig = FALSE)
               Comparison Wilcox.pvalue
1      Mild vs Uninfected  6.229920e-01
2  Moderate vs Uninfected  1.121082e-08
3        Moderate vs Mild  4.277167e-22
4    Severe vs Uninfected  1.030190e-08
5          Severe vs Mild  1.105551e-21
6      Severe vs Moderate  6.229920e-01
7      Dead vs Uninfected  3.168674e-03
8            Dead vs Mild  5.629542e-07
9        Dead vs Moderate  2.110720e-01
10         Dead vs Severe  1.266130e-01

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
