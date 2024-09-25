\name{compute_LFC}
\alias{compute_LFC}
\title{
Calculate Log-Fold Change (LFC) with Custom Weighting for Median and Mean
}
\description{
This function computes the Log-Fold Change (LFC) for each feature (taxon) in a \code{data.frame}, comparing two groups. The LFC formula allows users to customize the contribution of median and mean of the distributions in the two groups via the \code{w1} and \code{w2} parameters.
}
\usage{
compute_LFC(data, group, w1 = 0.5, w2 = 0.5)
}
\arguments{
  \item{data}{
  A \bold{data.frame} where each row represents a feature (e.g., a taxon or gene), and each column represents a sample. The data contains numeric values (e.g., relative abundances).
  }
  \item{group}{
  A \bold{factor} indicating the group assignment for each sample in \code{data}. This factor must contain exactly two levels, as the function computes LFC between two groups.
  }
  \item{w1}{
  A \bold{numeric} value between 0 and 1, representing the weight applied to the median of \code{group1} in the LFC calculation. The weight applied to the mean is \code{1 - w1}. Default is \code{w1 = 0.5}.
  }
  \item{w2}{
  A \bold{numeric} value between 0 and 1, representing the weight applied to the median of \code{group2} in the LFC calculation. The weight applied to the mean is \code{1 - w2}. Default is \code{w2 = 0.5}.
  }
}
\details{
This function computes Log-Fold Change (LFC) for each feature (row) in the provided data. LFC compares the values of the feature between two groups, with a customizable formula that uses both the median and mean of the feature values in each group. Users can control the relative contribution of the median and mean through the \code{w1} and \code{w2} parameters.

\bold{LFC formula:}

\code{
\verb{            }(w1 * median(group1) + (1 - w1) * mean(group1))\cr
LFC = log ( ----------------------------------------------    )\cr
\verb{            }(w2 * median(group2) + (1 - \code{w2} ) * mean(group2))
}

The function ensures no division by zero by returning \code{NA} for taxa where either group has zero values after weighting.
}
\value{
A \bold{data.frame} with two columns:
  \item{taxon}{
  The rownames of the input \code{data.frame}, representing the feature or taxon names.
  }
  \item{LFC}{
  The calculated Log-Fold Change (LFC) between the two groups. If division by zero occurs, the LFC is set to \code{NA}.
  }
}
\examples{
# Create example data with two groups
data <- data.frame(
  Taxon1 = c(10, 20, 30, 40),
  Taxon2 = c(5, 15, 25, 35),
  Taxon3 = c(2, 3, 4, 5)
)
group <- factor(c("Group1", "Group1", "Group2", "Group2"))

# Compute LFC with default weights (w1 = 0.5, w2 = 0.5)
LFC_result <- compute_LFC(data, group)

# Adjust weighting to favor the median for group 1 and the mean for group 2
LFC_result <- compute_LFC(data, group, w1 = 0.8, w2 = 0.2)

print(LFC_result)
}
\seealso{
\code{\link{log}}, \code{\link{mean}}, \code{\link{median}}
}
\author{
Marco Fabbrini
}
