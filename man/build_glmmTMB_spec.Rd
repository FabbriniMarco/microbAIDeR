\name{build_glmmTMB_spec}
\alias{build_glmmTMB_spec}
\title{Assemble a \code{glmmTMB}-ready model specification from simple inputs}
\description{
Convenience constructor that builds a main model formula and optional
zero-inflation and dispersion formulas for \pkg{glmmTMB}, with common
random-effects patterns (subject intercepts/slopes, grouping intercepts,
and by-taxa intercepts for stacked features). Optionally validates that
referenced variables exist in a provided data frame.
}
\usage{
build_glmmTMB_spec(
  response,
  fixed_effects = NULL,
  random_effects = NULL,
  family = NULL,
  dispformula = TRUE,
  ziformula = NULL,
  var.timepoint.random.slope = NULL,
  var.subject.random.intercept = NULL,
  taxa_var,
  include_random_intercept_by_taxa = TRUE,
  data = NULL
)
}
\arguments{
  \item{response}{\strong{Character string}. Name of the response variable to
  appear on the left-hand side of the main formula (e.g., \code{"count"}).}

  \item{fixed_effects}{Character vector of fixed-effect terms for the right-hand
  side (RHS) of the main formula. Each element is pasted with \code{+}. Can
  contain interaction syntax, e.g., \code{c("PCR * timepoint", "age", "sex")}.
  If \code{NULL} or length 0, an intercept-only RHS (\code{"1"}) is used.}

  \item{random_effects}{Character vector of grouping variables for additional
  random intercepts, each translated to a term of the form \code{(1 | group)}.
  Default \code{NULL}.}

  \item{family}{A \pkg{glmmTMB}-compatible family object or identifier. This is
  returned unchanged for convenience but not used internally when building
  formulas.}

  \item{dispformula}{Controls the dispersion model formula. 
  \itemize{
    \item If \code{TRUE}, uses \code{~ 0 + taxa_var} when \code{taxa_var} is provided, otherwise \code{~ 1}.
    \item If \code{FALSE} or \code{NULL}, the dispersion model is omitted.
    \item If a \strong{character} or a \strong{formula}, it is converted to a formula
      (character strings may omit the leading \code{"~"}).
  }
  \strong{Microbiome context:}
  \code{dispformula = ~ 0 + Genus} fits a completely separate dispersion parameter
  for each genus. This is often desirable for compositional abundance data,
  because different taxa can have very different variance structures even after
  accounting for the mean.}

  \item{ziformula}{Controls the zero-inflation model formula.
  \itemize{
    \item If \code{TRUE}, uses \code{~ 0 + taxa_var} when \code{taxa_var} is provided, otherwise \code{~ 1}.
    \item If \code{FALSE} or \code{NULL}, the zero-inflation model is omitted.
    \item If a \strong{character} or a \strong{formula}, it is converted to a formula.
  }
  \strong{Microbiome context:}
  In stacked abundance tables (one row per sample–genus pair), it’s often important
  to allow each genus to have its own zero-inflation probability. Some genera are
  ubiquitous (rarely zero), while others are rare and frequently absent. There’s
  no single “baseline genus” to compare against.  
  Thus, \code{~ 0 + Genus} fits separate zero-inflation coefficients (one per genus),
  effectively giving each genus its own zero-inflation sub-model.}

  \item{var.timepoint.random.slope}{Optional character scalar naming a
  (typically temporal) covariate to include as a random slope within
  \code{var.subject.random.intercept}, producing a term
  \code{(1 + slope | subject)}. Default \code{NULL}.}

  \item{var.subject.random.intercept}{Optional character scalar naming the
  subject/cluster ID for a subject-level random intercept, producing a term
  \code{(1 | subject)} (or a correlated intercept–slope if
  \code{var.timepoint.random.slope} is provided). Default \code{NULL}.}

  \item{taxa_var}{\strong{Character string}. Name of the stacked feature / taxa
  column (e.g., \code{"Genus"} or \code{"Taxa"}). When provided, it enables
  by-taxa random intercepts (if requested) and affects the defaults for
  \code{ziformula} and \code{dispformula}.}

  \item{include_random_intercept_by_taxa}{Logical. If \code{TRUE} and
  \code{taxa_var} is provided, adds a by-taxa random intercept term
  \code{(1 | taxa_var)} to the main formula. Default \code{TRUE}.

  \strong{What it does:}
  Adds a random intercept per genus (e.g., \code{(1 | Genus)}), meaning each
  genus has its own baseline level of abundance, randomly drawn around a global
  mean.

  \strong{When to use:}
  \itemize{
    \item ✅ Use it when you have a stacked model across many genera (one row per sample–genus pair).
    \item Different genera have different typical abundances.
    \item You want to borrow strength across genera while allowing each genus its own baseline mean.
  }

  \strong{When not needed:}
  \itemize{
    \item If you are fitting one genus at a time (per-genus models).}
  }

  \item{data}{Optional \code{data.frame}. If supplied, the function will warn
  about variables referenced in the specification that are not found among
  \code{colnames(data)}.}
}
\details{
\subsection{Main formula construction.}
The fixed-effects RHS is created by pasting \code{fixed_effects} with \code{+}
(or \code{"1"} if empty). Random-effect terms are then appended:
\enumerate{
  \item Subject-level random intercept \code{(1 | subject)} if
  \code{var.subject.random.intercept} is given.
  \item Optional subject-level random slope \code{(1 + slope | subject)} if
  \code{var.timepoint.random.slope} is also given.
  \item Additional grouping intercepts \code{(1 | g)} for each element of
  \code{random_effects}.
  \item By-taxa intercept \code{(1 | taxa_var)} when \code{taxa_var} is
  provided and \code{include_random_intercept_by_taxa = TRUE}.
}
These pieces are combined into \code{response ~ <fixed and random terms>}.

\subsection{Zero-inflation and dispersion formulas.}
For \code{ziformula} and \code{dispformula}:
\itemize{
  \item \code{TRUE} → \code{~ 0 + taxa_var} if \code{taxa_var} is provided, otherwise \code{~ 1}.
  \item \code{FALSE} or \code{NULL} → omit the component.
  \item \code{character} or \code{formula} → converted to a formula.
}

\strong{Guidance for microbiome data:}
\itemize{
  \item \code{ziformula = ~ 0 + Genus} → Each genus gets its own zero-inflation probability.  
        This captures the idea that some genera are common and rarely zero, while others are
        sparse and often absent.
  \item \code{dispformula = ~ 0 + Genus} → Each genus gets its own dispersion parameter,
        which is often beneficial for compositional or overdispersed count data.
}

\subsection{Validation.}
If \code{data} is supplied, the function heuristically tokenizes the fixed RHS
and checks that referenced variable names (including random-effect groupings,
subject/slope vars, and \code{taxa_var}) appear as columns in \code{data}.
Missing variables trigger a warning.
}
\value{
A named \code{list} with components:
\itemize{
  \item \code{formula}: the main model \code{\link[stats]{formula}}.
  \item \code{ziformula}: a \code{formula} for zero inflation, or \code{NULL}.
  \item \code{dispformula}: a \code{formula} for dispersion, or \code{NULL}.
  \item \code{family}: the \code{family} argument passed in (returned unchanged).
}
This object can be supplied directly to \code{\link[glmmTMB]{glmmTMB}} via
\code{do.call} or by extracting components.
}
\seealso{
\code{\link[glmmTMB]{glmmTMB}}, \code{\link[stats]{as.formula}}
}
\examples{
## Not run: 
## Basic Poisson GLMM with fixed effects only
spec <- build_glmmTMB_spec(
  response   = "count",
  fixed_effects = c("PCR * timepoint", "age", "sex"),
  taxa_var   = "Genus",
  family     = poisson()
)

## Microbiome-style zero-inflation and dispersion per genus
spec2 <- build_glmmTMB_spec(
  response   = "count",
  fixed_effects = c("PCR", "timepoint", "age"),
  var.subject.random.intercept = "Subject_ID",
  var.timepoint.random.slope   = "timepoint",
  random_effects = c("Center", "Batch"),
  taxa_var   = "Genus",
  include_random_intercept_by_taxa = TRUE,
  ziformula  = TRUE,        # => ~ 0 + Genus
  dispformula = TRUE,       # => ~ 0 + Genus
  family     = nbinom2()
)

## If running per-genus models, you can disable the genus-level intercept
spec3 <- build_glmmTMB_spec(
  response   = "count",
  fixed_effects = "treatment + age",
  var.subject.random.intercept = "ID",
  taxa_var   = "Genus",
  include_random_intercept_by_taxa = FALSE,
  ziformula  = FALSE,
  dispformula = FALSE,
  family     = poisson()
)

## Other examples
spec4 <- build_glmmTMB_spec(
  response = "abundance",
  fixed_effects = c("group * Timepoint * taxa_name", "age", "gender"),
  var.subject.random.intercept = "Subject_ID",
  var.timepoint.random.slope = "Timepoint",
  taxa_var = "taxa_name",
  include_random_intercept_by_taxa = TRUE,
  dispformula = TRUE,  # -> ~ 0 + Genus
  ziformula = TRUE     # -> ~ 0 + Genus
)

## Fit with glmmTMB
## library(glmmTMB)
## fit <- glmmTMB(
##   formula     = spec2$formula,
##   ziformula   = spec2$ziformula,
##   dispformula = spec2$dispformula,
##   family      = spec2$family,
##   data        = your_data
## )
## End(Not run)
}


\author{
Marco Fabbrini - m.fabbrini@unibo.it
}
