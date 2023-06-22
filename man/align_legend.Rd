\name{align_legend}
\alias{align_legend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Align legend title to legend items in a ggplot2 plot
}
\description{
This function aligns the legend title to the legend keys of a ggplot2 object.
}
\usage{
ggdraw( align_legend(x) )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \bold{ggplot} with a plotted legend.
}
}

\author{
Marco Fabbrini - m.fabbrini@unibo.it
}

\examples{

data("mtcars")

a <- myscripts::tidify(mtcars[,1:2])


b <- ggplot(data=a, aes(x=rowname, y=colname, fill=value))+
    geom_tile()+
    scale_fill_distiller(palette = "OrRd", direction = 1)+
    labs(x="", y="", fill="This is really\na long title for \na legend!")

ggdraw( align_legend(b) )

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
