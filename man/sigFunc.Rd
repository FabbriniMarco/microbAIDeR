\name{sigFunc}
\alias{sigFunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert a P-value into symbolic pinpoint
}
\description{
This function converts a single \bold{numeric} P-value into a character of encoded statistical significance.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sigFunc(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  A \bold{numeric} object
}
}
\details{
0.05 < values <= 0.1 : the function returns a ° to indicate a trend; \cr

0.01 < values <= 0.05 : the function returns a * ; \cr

0.001 < values <= 0.01 : the function returns  ** ; \cr

0.0001 < values <= 0.001 : the function returns  *** ; \cr

Values <= 0.0001 : the function returns  **** ; \cr

Otherwise the function returns NA. \cr

If a NA is passed to the function, a NA is given. \cr
}
\author{
Marco Fabbrini - m.fabbrini@unibo.it
}
\examples{
> sigFunc(0.3)
[1] NA

> sigFunc(0.1)
[1] "°"

> sigFunc(0.06)
[1] "°"

> sigFunc(0.05)
[1] "*"

> sigFunc(0.03)
[1] "*"

> sigFunc(0.01)
[1] "**"

> sigFunc(0.007)
[1] "**"

> sigFunc(0.001)
[1] "***"

> sigFunc(0.0005)
[1] "***"

> sigFunc(0.0001)
[1] "****"

> sigFunc(2.3e-14)
[1] "****"


Apply the conversion to a numeric vector of P-values
> p_list = c(0.3, 0.1, 0.06, 0.05, 0.03, 0.01, 0.005, 0.001, 0.0001, 1.3e-10)

> sapply(p_list, sigFunc)
[1] NA     "°"     "°"    "*"    "*"    "**"    "**"   "***"   "****"  "****"

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
