\name{amp_up_models}
\alias{amp_up_models}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Allow machine learning models to work in CPU multi-thread
%%  ~~function to do ... ~~
}
\description{
Machine learning models especially from the \emph{caret} package could benefit of muulti-threading performance boosting via the \bold{allowParallel = TRUE} parameter.
This function provides regulated access to computer threads fractioning all the available threads through the core_fraction parameter. This function requires the foreach and iterators packages to be installed.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
amp_up_models(core_fraction = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{core_fraction}{
          (Value from 0.1 to 1) Fraction of the available cores to be unlocked
}
}
\author{
Marco Fabbrini - m.fabbrini@unibo.it
%%  ~~who you are~~
}
\examples{
# Example on a 20 threads system to use 40\% of system horsepower
> amp_up_models(core_fraction = 0.4)
Loading the required package: foreach
Loading the required package: iterators
Model amped and ready to go with: 8 cores.
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
