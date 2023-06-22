\name{is.nan.data.frame}
\alias{is.nan.data.frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extend the is.nan() function to data.frame objects
}
\description{
As \bold{is.na()} works on \bold{data.frame} objects, this function allows to extend \bold{is.nan()} function to \bold{data.frame}s.
}
\usage{
is.nan.data.frame(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \bold{data.frame} object
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function returns an object of class \bold{matrix} \bold{array} with logicals pointing to NaNs.
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
> a = c(2, 6, 4, NaN, 2, 1, 0)
> b = c(1, 7, 3, 56, NaN , 1, 2)
> mydf = data.frame(a, b)

> mydf
    a   b
1   2   1
2   6   7
3   4   3
4 NaN  56
5   2 NaN
6   1   1
7   0   2

> is.nan.data.frame(mydf)
         a     b
[1,] FALSE FALSE
[2,] FALSE FALSE
[3,] FALSE FALSE
[4,]  TRUE FALSE
[5,] FALSE  TRUE
[6,] FALSE FALSE
[7,] FALSE FALSE

> mydf[is.nan.data.frame(mydf)] = "HERE_IT_IS"
> mydf
           a          b
1          2          1
2          6          7
3          4          3
4 HERE_IT_IS         56
5          2 HERE_IT_IS
6          1          1
7          0          2

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
