\name{mkdir}
\alias{mkdir}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a folder if it doesn't already exists
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This script takes a PATH as input and checks if the given folder (and eventual subfolders) already exist. If not, it creates it/them.
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
Marco Fabbrini 
}

\examples{
mkdir("Plots")
mkdir("Plots/Subfolder of Plots")
mkdir("Plots/Subfolder of Plots/Additional subfolder")

# Running directly the third command will automatically create the 'Plots', 'Subfolder of Plots' and 'Additional subfolder' folders (recursively)
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
