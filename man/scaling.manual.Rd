\name{scaling.manual}
\alias{scaling.manual}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scale a vector of data to fit an interval of values
}
\description{
This function takes a numeric vector as input and scales the values to fit inside a provided range. \cr
This function is particularly useful to restrain data size for plotting purposes.
}
\usage{
scaling.manual(x, range.min, range.max)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \bold{numeric} vector of input values to be scaled
}
  \item{range.min}{
A \bold{numeric} value to set the minimum value of the scaling range.
}
  \item{range.max}{
A \bold{numeric} value to set the maximum value of the scaling range.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A \bold{numeric} vector is returned as output, with values ranging from \emph{range.min} to \emph{range.max}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Gabriele Conti
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

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
