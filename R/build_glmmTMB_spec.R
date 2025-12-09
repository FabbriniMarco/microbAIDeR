build_glmmTMB_spec <- function(
    response,
    fixed_effects = NULL,                  # character vector of fixed-effect terms (e.g., c("New_PCR * New_timepoint","age","gender"))
    random_effects = NULL,                 # character vector of grouping vars for random intercepts (e.g., c("Center","Batch"))
    family_model = NULL,                         # model family, returned for convenience
    dispformula = TRUE,                    # TRUE => ~ 0 + taxa_var (if provided) else ~ 1; FALSE/NULL => omit; character => custom
    ziformula = NULL,                      # TRUE => ~ 0 + taxa_var (if provided) else ~ 1; FALSE/NULL => omit; character => custom
    var.timepoint.random.slope = NULL,     # e.g., "New_timepoint"
    var.subject.random.intercept = NULL,   # e.g., "Subject_ID"
    taxa_var = NULL,                              # name of stacked feature column (e.g., "Genus" or "Taxa"); can be NULL if not stacked
    include_random_intercept_by_taxa = TRUE,
    data = NULL                            # optional: validate variable names presence
) {
  # --- helpers ---
  .to_formula <- function(x) {
    if (is.null(x) || isFALSE(x)) return(NULL)
    if (isTRUE(x)) stop("TRUE not allowed here; pass character or NULL.")
    if (inherits(x, "formula")) return(x)
    if (is.character(x) && length(x) == 1L) {
      if (grepl("^\\s*~", x)) stats::as.formula(x) else stats::as.formula(paste("~", x))
    } else stop("Unsupported type for formula input.")
  }
  
  # --- fixed effects ---
  fixed_rhs <- if (is.null(fixed_effects) || length(fixed_effects) == 0L) "1" else paste(fixed_effects, collapse = " + ")
  
  # --- random effects ---
  rand_terms <- character(0)
  
  # Subject random intercept (+ optional time slope)
  if (!is.null(var.subject.random.intercept)) {
    if (!is.null(var.timepoint.random.slope)) {
      rand_terms <- c(rand_terms, sprintf("(1 + %s | %s)", var.timepoint.random.slope, var.subject.random.intercept))
    } else {
      rand_terms <- c(rand_terms, sprintf("(1 | %s)", var.subject.random.intercept))
    }
  }
  
  # Additional random intercepts
  if (!is.null(random_effects) && length(random_effects)) {
    rand_terms <- c(rand_terms, sprintf("(1 | %s)", random_effects))
  }
  
  # By-taxa random intercept (only if stacked and requested)
  if (!is.null(taxa_var) && isTRUE(include_random_intercept_by_taxa)) {
    rand_terms <- c(rand_terms, sprintf("(1 | %s)", taxa_var))
  }
  
  # Combine RHS
  rhs_all <- fixed_rhs
  if (length(rand_terms)) rhs_all <- paste(rhs_all, paste(rand_terms, collapse = " + "), sep = " + ")
  
  # Full model formula
  main_formula <- stats::as.formula(paste(response, "~", rhs_all))
  
  # --- ziformula ---
  zi_formula <- NULL
  if (isTRUE(ziformula)) {
    zi_formula <- if (!is.null(taxa_var)) stats::as.formula(paste("~ 0 +", taxa_var)) else ~ 1
  } else if (!is.null(ziformula) && !isFALSE(ziformula)) {
    zi_formula <- .to_formula(ziformula)
  }
  
  # --- dispformula ---
  disp_formula <- NULL
  if (isTRUE(dispformula)) {
    disp_formula <- if (!is.null(taxa_var)) stats::as.formula(paste("~ 0 +", taxa_var)) else ~ 1
  } else if (!is.null(dispformula) && !isFALSE(dispformula)) {
    disp_formula <- .to_formula(dispformula)
  }
  
  # --- optional validation ---
  if (!is.null(data)) {
    # Rough parse of fixed terms into tokens to check presence
    fixed_tokens <- gsub("[*:/()|+~-]", " ", fixed_rhs)
    fixed_tokens <- unique(strsplit(fixed_tokens, "\\s+")[[1]])
    fixed_tokens <- fixed_tokens[fixed_tokens != "" & fixed_tokens != "1"]
    
    vars_to_check <- unique(c(
      fixed_tokens,
      random_effects,
      var.subject.random.intercept,
      var.timepoint.random.slope,
      taxa_var
    ))
    vars_to_check <- vars_to_check[!is.na(vars_to_check) & vars_to_check != ""]
    missing_vars <- setdiff(vars_to_check, colnames(data))
    if (length(missing_vars)) {
      warning("These variables were referenced but not found in 'data': ",
              paste(missing_vars, collapse = ", "))
    }
  }
  
  list(
    formula     = main_formula,
    ziformula   = zi_formula,
    dispformula = disp_formula,
    family_model      = family_model
  )
}
