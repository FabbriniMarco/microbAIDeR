\name{generate_random_colors}
\alias{generate_random_colors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Distinct Random Colors Using RColorBrewer
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function generates a vector of distinct colors by randomly sampling from the qualitative palettes available in the \code{RColorBrewer} package. The first color in the vector is specified by the \code{start.col} parameter, and the remaining colors are randomly selected. This is useful for creating color palettes for plots where distinct colors are required. \cr
If you call the function with \code{show.colors = TRUE} without assigning the output to a variable, a barplot of the generated colors will be automatically displayed.
}
\usage{
generate_random_colors(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An \bold{integer} specifying the total number of colors to generate (including the starting color).
  }
  \item{start.col}{
    A \bold{character} string specifying the color code for the first color in the palette. Default is \code{"grey90"}.
  }
  \item{seed}{
    An \bold{integer} used to set the random seed for reproducibility. Default is \code{1996}.
  }
  \item{show.colors}{
    A \bold{logical} indicating whether to display a barplot of the generated colors. Default is \code{FALSE}.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The function works by extracting all available colors from the qualitative palettes in \code{RColorBrewer}. It then randomly selects \code{x - 1} colors from this pool and combines them with the provided starting color. If the requested number of colors exceeds the available colors, a message is returned and \code{NA} is output.
}
\value{
%%  ~Describe the value returned
A \bold{character} vector containing \code{x} color codes.
}
\references{
%% ~put references to the literature/web site here ~
Fabbrini, M. & Conti, G. (2024). \emph{microbAIDeR - An ensemble of functions for easier and quicker preliminary microbiome analyses}. Available on Github: \url{https://github.com/FabbriniMarco/microbAIDeR/}. \cr
Brewer, C. A. (2002). \emph{ColorBrewer 2.0: Color Advice for Maps}. Retrieved from \url{https://colorbrewer2.org/}.
}
\author{
%%  ~~who you are~~
Gabriele Conti, Marco Fabbrini
}
\note{
This function depends on the \code{RColorBrewer} package. Make sure that the package is installed and loaded.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{barplot}}, \code{\link[RColorBrewer]{brewer.pal}}
}
\examples{
## Generate n_cols distinct colors and display them in a barplot:
n_cols <- 8
my_colors <- generate_random_colors(n_cols, start.col = "grey90", seed = 1996, show.colors = TRUE)
my_colors

# [1] "grey90"  "#B3DE69" "#377EB8" "#FC8D62" "#F781BF" "#FDB462" "#B15928" "#FFF2AE"
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
