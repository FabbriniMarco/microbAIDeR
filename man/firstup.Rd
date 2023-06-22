\name{firstup}
\alias{firstup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set the first letter of a string to MAIUSC
}
\description{
Simply take an input string and return with the first capital letter
}
\usage{
firstup(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \bold{character} element (not a vector, for vectors iterate using \emph{sapply}, or \emph{lapply} )
}
}

\author{
%%  ~~who you are~~
}
\note{
Marco Fabbrini - m.fabbrini@unibo.it
}

\examples{

> lowercase_text = c("my text")

> firstup("my text")
[1] "My text"

> some_more_text = c("this", "is", "some", "more", "text")

> sapply( some_more_text, firstup )
  this     is   some   more   text
"This"   "Is" "Some" "More" "Text"

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
