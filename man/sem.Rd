\name{sem}
\alias{sem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute standard error of the mean
}
\description{
Compute standard error of the mean (SEM)
}
\usage{
sem(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \bold{numeric} vector containing the distribution values. Can contain NAs, which will be removed.
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
distrib = runif(200, min = 3, max = 6)

sem(distrib)
[1] 0.06085993

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
