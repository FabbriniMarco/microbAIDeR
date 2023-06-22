\name{flattenCorrMatrix}
\alias{flattenCorrMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transforms a squared correlation matrix with a squared p-values matrix in a single \emph{data.frame} object
}
\description{
Transforms a squared correlation matrix with a squared p-values matrix in a single \emph{data.frame} object.
}
\usage{
flattenCorrMatrix(corrmat, pmat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{corrmat}{
A matrix output of a correlation for example made with \emph{rcorr} function. The matrix contains features on rows and columns and correlation coefficients numeric values.
}
  \item{pmat}{
A matrix output of a correlation for example made with \emph{rcorr} function. The matrix contains features on rows and columns and correlation p-values numeric values. The output is a tidy data.frame object. Particularly useful when dealing with correlation networks.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
library(Hmisc)
data(mtcars)

corr_output = rcorr( as.matrix(mtcars) )

head(corr_output$r)

            mpg        cyl       disp         hp       drat         wt        qsec         vs         am       gear       carb
mpg   1.0000000 -0.8521620 -0.8475514 -0.7761684  0.6811719 -0.8676594  0.41868403  0.6640389  0.5998324  0.4802848 -0.5509251
cyl  -0.8521620  1.0000000  0.9020329  0.8324475 -0.6999381  0.7824958 -0.59124207 -0.8108118 -0.5226070 -0.4926866  0.5269883
disp -0.8475514  0.9020329  1.0000000  0.7909486 -0.7102139  0.8879799 -0.43369788 -0.7104159 -0.5912270 -0.5555692  0.3949769
hp   -0.7761684  0.8324475  0.7909486  1.0000000 -0.4487591  0.6587479 -0.70822339 -0.7230967 -0.2432043 -0.1257043  0.7498125
drat  0.6811719 -0.6999381 -0.7102139 -0.4487591  1.0000000 -0.7124406  0.09120476  0.4402785  0.7127111  0.6996101 -0.0907898
wt   -0.8676594  0.7824958  0.8879799  0.6587479 -0.7124406  1.0000000 -0.17471588 -0.5549157 -0.6924953 -0.5832870  0.4276059

flattened = flattenCorrMatrix( corr_output$r, corr_output$P )

head(flattened)

   row column        cor            p
1  mpg    cyl -0.8521620 6.112688e-10
2  mpg   disp -0.8475514 9.380328e-10
3  cyl   disp  0.9020329 1.803002e-12
4  mpg     hp -0.7761684 1.787835e-07
5  cyl     hp  0.8324475 3.477861e-09
6 disp     hp  0.7909486 7.142679e-08

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
