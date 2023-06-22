\name{tokenize}
\alias{tokenize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Tokenize: turn a factor into a indexed positional vector
%%  ~~function to do ... ~~
}
\description{
tokenize() is a function used to turn a \bold{factor} into a \bold{vector} of integers based on the index position of the corresponding level for each element in the factor.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tokenize(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  A \bold{factor} with defined levels. Each element in the factor would be converted in the corresponding index of the matching level
%%     ~~Describe \code{x} here~~
}
}

\author{
Marco Fabbrini - m.fabbrini@unibo.it
%%  ~~who you are~~
}

\examples{
> print(my_factor)
[1] Three Two   Two   Three Two   One   One   Three Two   Three
[11] One   One   Three One   Two
Levels: One Two Three

> tokenize(my_factor)
[1] 3 2 2 3 2 1 1 3 2 3 1 1 3 1 2


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
