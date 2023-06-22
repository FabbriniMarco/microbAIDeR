\name{normalize}
\alias{normalize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normalize a vector of value setting the highest to 1 and the other proportionally
}
\description{
This function takes as input a numeric vector and returns a numeric vector of the same length, scaling the data from 0 to 1 according to the highest value detected in the input vector.
}
\usage{
normalize(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \bold{numeric} vector
}
}

\author{
Marco Fabbrini - m.fabbrini@unibo.it
}

\examples{

> a <- runif(10, min = 0.2, max = 123)

> a
 [1] 60.333382 57.381081 38.889767  3.020177  3.543527 17.908551 48.452833 76.978163 84.090733 31.044141


> normalize(a)
 [1] 0.706954631 0.670538195 0.442449048 0.000000000 0.006455485 0.183647123 0.560408835 0.912266908 1.000000000 0.345673756

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
