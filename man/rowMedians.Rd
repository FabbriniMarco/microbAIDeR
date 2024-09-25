\name{rowMedians}
\alias{rowMedians}
\title{
Compute Row-Wise Medians for a Data Frame
}
\description{
This function calculates the median value for each row of a given \code{data.frame}, ignoring missing values (\code{NA}).
}
\usage{
rowMedians(df)
}
\arguments{
  \item{df}{
  A \bold{data.frame} or \bold{matrix} where the medians are computed for each row. The input must contain numeric values.
  }
}
\details{
The function uses the \code{apply()} function to compute the median for each row of the input \code{df}. The calculation is performed by applying the \code{median()} function row-wise, with the argument \code{na.rm = TRUE} to handle missing values (i.e., \code{NA}). This ensures that rows with missing values still return a median, ignoring those \code{NA}s.
}
\value{
A numeric vector containing the row-wise medians of the input \code{df}.
}
\examples{
# Create a data frame with numeric values
df <- data.frame(
  A = c(1, 2, 3, NA),
  B = c(4, 5, NA, 8),
  C = c(7, 8, 9, 10)
)

# Compute row-wise medians
row_medians <- rowMedians(df)
print(row_medians)

# Example output:
# [1] 4.0 5.0 6.0 9.0
}
\seealso{
\code{\link{apply}}, \code{\link{median}}
}
\author{
Marco Fabbrini
}
