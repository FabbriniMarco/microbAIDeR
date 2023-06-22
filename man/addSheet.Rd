\name{addSheet}
\alias{addSheet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add an excel sheet to a previously created excel workbook - with openxlsx
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
addSheet(path, sheet.name, addition, col.save = TRUE, row.save = TRUE, overwrite = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
Full \bold{PATH} or short path to the excel workbook (*.xls or *.xlsx) to be edited
}
  \item{sheet.name}{
Name for the new sheet that is going to be created. Max 31 characters including space. Must NOT include special characters as /, \, ?, [ or ]
}
  \item{addition}{
A \bold{data.frame} object to be added to the sheeet
}
  \item{col.save}{
Logical specifying whether to save also data.frame \bold{colnames} in the excel sheet. Default = TRUE.
}
  \item{row.save}{
Logical specifying whether to save also data.frame \bold{rownames} in the excel sheet. Default = TRUE.
}
  \item{overwrite}{
Logical specifying whether to overwrite the sheet in the Excel file if there's already a sheet named as \bold{sheet.name}. Default = TRUE.
}
}
\details{
Requires the \emph{openxlsx} package to be installed.
}

\author{
Marco Fabbrini - m.fabbrini@unibo.it
}

\examples{
addSheet( "excel_file.xlsx", sheet.name = "My new sheet",
             addition = dataframe_to_add , col.save = TRUE, row.save = TRUE)
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
