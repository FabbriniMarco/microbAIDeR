\name{se.avg}
\alias{se.avg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the standard error of the mean
}
\description{
Computes the standard error of the mean
}
\usage{
se.avg(x, na.rm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \bold{numeric} vactor containing the values of the distribution we want to compute the standard error of the mean of. Can contain NAs, which will be removed.
}
  \item{na.rm}{
%%     ~~Describe \code{na.rm} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, na.rm = TRUE)
{
    if (na.rm) {
        x = x[!is.na(x)]
    }
    sd(x)/sqrt(length(x))
  }
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
