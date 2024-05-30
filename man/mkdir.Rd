\name{mkdir}
\alias{mkdir}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a folder if it doesn't already exists
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This script takes a PATH as input and checks if the given folder already exists. If not, it creates it.
}
\usage{
mkdir(path)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
Full or local \bold{PATH} to the folder that we want to create
}
}
\details{

}

\author{
Marco Fabbrini - m.fabbrini@unibo.it
}

\examples{
mkdir("Plots")
mkdir("Plots/Subfolder of Plots")
mkdir("Plots/Subfolder of Plots/Additional subfolder")
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
