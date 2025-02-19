<div align=center>
  <img width="196" alt="Circle logo" src="https://github.com/user-attachments/assets/9953b84a-12ac-4704-84c1-678bbc92eafd" height="186">
</div>


# microbAIDeR - An ensemble of functions for easier and quicker preliminary microbiome analyses.

microbAIDeR is a R-package made for easy preliminar microbiome analyses starting from typical 16S outputs, such as the one produced in the [16S-KrakenBracken-Pipeline](https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline). The package provides all-in-one functions for statistical analysis and graphs generation, as well as useful function for general data analyses in R.
This package requires RTools to be installed in order to compiled the required dependancies. Please follow the developer's instruction from the [Comprehensive R Archive Network](https://cran.r-project.org/) to install RTools on your platform.

<br>

# Installation of the package and the required dependancies

```R
install.packages('devtools')
devtools::install_github("https://github.com/FabbriniMarco/microbAIDeR.git")
library(microbAIDeR)
microbAIDeR_install_dependancies() #To install all the packages that are necessary
```

Depends: devtools, pairwiseAdonis, parallel, doParallel, openxlsx, tidyr, tibble, dplyr, ggplot2, ggsignif, gridExtra, tidyverse, stringr, vegan, cowplot, gtable, grid, ggpattern, svglite, RColorBrewer

<br>

# Examples

Usage example for this package are included in the [Wiki](https://github.com/FabbriniMarco/microbAIDeR/wiki) section.

<br>

# Citation

If you use this R package please cite this GitHub in your article:

```
@Manual{,
  title = {microbAIDeR: An ensemble of functions for easier and quicker preliminary microbiome analyses},
  author = {Marco Fabbrini, Gabriele Conti},
  year = {2025},
  note = {R package version 0.3.2},
  url = {https://github.com/FabbriniMarco/microbAIDeR.git},
}
```
Example: 
> Fabbrini, M. & Conti, G. (2025) microbAIDeR - An ensemble of functions for easier and quicker preliminary microbiome analyses (v0.3.2). Available on Github: https://github.com/FabbriniMarco/microbAIDeR/

