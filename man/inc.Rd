\name{inc}
\alias{inc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Increment a variable by a chosen amount
}
\description{
This function takes a variable as input and re-assigns the value of such variable adding an increment.
Forget about reassigning a variable again in a for-loop!
}
\usage{
inc(x, increment = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \bold{numeric} or \bold{integer} variable already present in the environment.
}
  \item{increment}{
A \bold{numeric} value corresponding to the increase that the input variable \emph{x} will receive.
}
}

\author{
Marco Fabbrini - m.fabbrini@unibo.it
}
\examples{

> a = 1

> for ( cycles in 1:10){ inc(a, increment = 1) }
# This is equivalent to: for ( cycles in 1:10){ a = a + 1 }
> print(a)
[1] 11

> inc(a)
> print(a)
[1] 12

> inc(a, increment = 0.3)
> a
[1] 12.3


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
