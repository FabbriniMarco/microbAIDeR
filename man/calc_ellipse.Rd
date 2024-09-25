\name{calc_ellipse}
\alias{calc_ellipse}
\title{
Calculate Ellipse for Beta Diversity Coordinates
}
\description{
This function calculates the ellipse around points in a PCoA or MDS space based on a specified confidence level. It computes the ellipse coordinates for the principal coordinates (PC1 and PC2) from a beta diversity matrix.
}
\usage{
calc_ellipse(coord_beta, level = 0.99)
}
\arguments{
  \item{coord_beta}{
  A \bold{data.frame} containing at least three columns: PC1, PC2, and Group. The PC1 and PC2 columns represent the coordinates of the samples, while the Group column identifies the grouping of samples (e.g., by condition or treatment).
  }
  
  \item{level}{
  A \bold{numeric} value representing the confidence level for the ellipse calculation. Default is 0.99, corresponding to a 99\% confidence ellipse.
  }
}
\details{
This function calculates the parameters for an ellipse that encloses a specified proportion of the points (by default 99\%). The ellipse is centered at the mean of the provided coordinates and is scaled by the covariance matrix of the coordinates. The function returns a data frame with the coordinates for plotting the ellipse.
}
\value{
A \bold{data.frame} containing the coordinates of the ellipse, with columns:
  \item{PC1}{
  The x-axis coordinates of the ellipse.
  }
  \item{PC2}{
  The y-axis coordinates of the ellipse.
  }
  \item{Group}{
  The group identifier, which matches the unique value of the Group column from \emph{coord_beta}.
  }
}
\author{
Marco Fabbrini
}
\examples{
# Example usage with sample data

# Create example coordinate data
coord_beta <- data.frame(
  PC1 = rnorm(20, mean = 0, sd = 1),
  PC2 = rnorm(20, mean = 0, sd = 1),
  Group = rep("Treatment", 20)
)

# Calculate the ellipse
ellipse_data <- calc_ellipse(coord_beta, level = 0.95)

# Visualize the result with ggplot2
if(require(ggplot2)) {
  ggplot(coord_beta, aes(x = PC1, y = PC2)) +
    geom_point() +
    geom_path(data = ellipse_data, aes(x = PC1, y = PC2), color = "blue") +
    ggtitle("95\% Confidence Ellipse for Treatment Group")
}
}
\seealso{
\code{\link{ggplot2}}, \code{\link{stat_ellipse}}
}
\references{
%% Provide relevant references or further reading if needed
}
\keyword{ellipse}
