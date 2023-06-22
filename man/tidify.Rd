\name{tidify}
\alias{tidify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quickly convert a tabular \bold{data.frame} object in tidy \bold{data.frame}
%%  ~~function to do ... ~~
}
\description{
In a \emph{tidy} \bold{data.frame} every variable is stored in its own column. Every observation is stored in its own row — that is, every row corresponds to a single case. Each value of a variable is stored in a cell of the table. This format is what you need to work with \emph{ggplot} and to keep your data more.. tidy!
Working with non-tidy data might be easier sometimes, and converting one into another might be messy. This function turns a canonical observation/samples data.frame into a tidy data.frame
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tidify(x, as.numeric = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \bold{data.frame} object with defined colnames and rownames
}
  \item{as.numeric}{
A \bold{logical} (Default = TRUE) specifying whether to convert to numeric the data values present in the data.frame. The resulting \emph{value} column of the tidy data frame would therefore be of class \bold{numeric}.
}
}

\author{
Marco Fabbrini - m.fabbrini@unibo.it
%%  ~~who you are~~
}

\examples{
> a = data.frame( A=c(1,3,4,5), B=c(2,3,4,6), C=c(2,3,5,7))
> rownames(a) = c("Line1", "Line2", "Line3", "Line4")

> print(a)
      A B C
Line1 1 2 2
Line2 3 3 3
Line3 4 4 5
Line4 5 6 7

> class(a)
[1] "data.frame"

> tidify(a, as.numeric = TRUE)
   rowname colname value
1    Line1       A     1
2    Line2       A     3
3    Line3       A     4
4    Line4       A     5
5    Line1       B     2
6    Line2       B     3
7    Line3       B     4
8    Line4       B     6
9    Line1       C     2
10   Line2       C     3
11   Line3       C     5
12   Line4       C     7


> class(tidify(a))
[1] "data.frame"
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
