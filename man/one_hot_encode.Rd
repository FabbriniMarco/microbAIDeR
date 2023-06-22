\name{one_hot_encode}
\alias{one_hot_encode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{One-hot encode a factor returning a 0/1 dataframe
%%  ~~function to do ... ~~
}
\description{
The one_hot_encode() function takes a factor as the input argument and yields a \bold{data.frame} object with in each line a one-hot encoding presence of the corresponding factor element. Columns are ordered according to factor levels
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
one_hot_encode(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
A \bold{factor} with matching levels.
}
}

\author{
Marco Fabbrini - m.fabbrini@unibo.it
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
> print(my_factor)
[1] Three Two   Two   Three Two   One   One   Three Two   Three
[11] One   One   Three One   Two
Levels: One Two Three

> one_hot_encode(b)
  One Two Three
#1    0   0     1
#2    0   1     0
#3    0   1     0
#4    0   0     1
#5    0   1     0
#6    1   0     0
#7    1   0     0
#8    0   0     1
#9    0   1     0
#10   0   0     1
#11   1   0     0
#12   1   0     0
#13   0   0     1
#14   1   0     0
#15   0   1     0
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
