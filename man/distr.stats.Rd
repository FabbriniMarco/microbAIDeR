\name{distr.stats}
\alias{distr.stats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quickly have a look at the distribution of your data
}
\description{
This function takes as input a \bold{numeric} vector and computes the mean, standard deviation, standard error of the mean, median, interquantile range (IQR) and the standard deviation of the median of the distribution.
}
\usage{
distr.stats(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \bold{numeric} vector. Can contain NAs.
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
Marco Fabbrini
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

distrib = runif(300, min = 0, max = 1)

distr.stats(distrib)

#$mean
[1] 0.5071967

$SD
[1] 0.2792346

$SEM
[1] 0.01612162

$median
[1] 0.5045009

$iqr
[1] 0.4405633

$SEMed
[1] 0.03187885

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
