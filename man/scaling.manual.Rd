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
distrib = runif(20, min = 3, max = 6)
distrib

# [1] 5.846163 5.350342 4.624601 4.048783 3.281904 4.335436 3.680077 4.025879 3.828067 5.104560 3.540106 4.475808 3.441325 4.030834 5.735098
# [16] 4.463539 3.435960 4.212136 5.746419 4.392540

scaling.manual(distrib, range.min = 2, range.max = 33)
# [1] 33.000000 27.005882 18.232217 11.271001  2.000000 14.736424  6.813624 10.994116  8.602706 24.034570  5.121469 16.433417  3.927282
# [14] 11.054007 31.657307 16.285090  3.862423 13.245823 31.794160 15.426768

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
