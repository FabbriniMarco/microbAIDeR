compute_LMM_and_plot <- function(
    data,                         # features x samples (rows = features like taxa or alpha values, columns = samples)
    group,                        # name of all meta columns acting as grouping factors, length = ncol(data), must be class factor
    taxlevel,
    comparison.list = NULL,       # used to display only selective comparisons
    save.path = getwd(),
    color.grouping,
    p.adjust.method = "fdr",
    trends = TRUE,
    plot.not.sig = FALSE,
    paired = FALSE,
    mode = c("singlefeature","wholetable"),     # Fit one model per each feature (ideal for alpha diversity) or fit a single model with all the features (taxa study)
    data_type = "relabb",         # "relabb" or "counts"
    meta = NULL,                  # data.frame with covariates; rownames = sample ids (must match colnames(data))
    fixed_effects = NULL,         # vector with meta column names of fixed_effects
    random_effects = NULL,        # vector with meta column names of random_effects
    family_model_model = NULL ,               # model family_model used to fit the model in glmmTMB
    dispformula = NULL,
    ziformula = NULL,
    var.timepoint.random.slope = NULL,          # If specified, a random sloper for time is added to the model (1 + Time)
    var.subject.random.intercept = NULL,        # meta column names of subject ID
    include_random_intercept_by_taxa = TRUE,
    conflevel = 0.95,             # emmeans confidence level
    save_single_model_summaries = FALSE,
    save_singlemode_emmeans = FALSE,
    nrow.graph = 2, ncol.graph = 2, width.graph = 4.5, height.graph = 3.5, horiz = FALSE, ggplot.margins = c(.18, .18, .18, .6),
    box.lwd = 0.4, jitter.pch = 21, jitter.stroke = 0.15, jitter.size = 0.7, jitter.color = "grey22",
    signif.step.increase = 0.12, signif.text.size = 3, signif.line.size = 0.4, contrast.color = "ivory1",
    text.x.size = 6, text.y.size = 6, text.y.title.size = 8, smoothing = FALSE, smoothing.lwd = 1, smoothing.color = "darkred", smoothing.se = FALSE,
    smoothing.method="loess", additional.params = NULL, align.legend = FALSE, plot.order = "kruskal", pattern.fill = FALSE, pattern = "stripe",
    pattern.angle = 45, pattern.alpha = 0.4, pattern.density = 0.1, pattern.spacing = 0.05
) {
  check_and_load_package("dplyr")
  check_and_load_package("ggplot2")
  check_and_load_package("ggsignif")
  check_and_load_package("gridExtra")
  check_and_load_package("tidyverse")
  check_and_load_package("tidyr")
  check_and_load_package("tibble")
  check_and_load_package("openxlsx")
  check_and_load_package("glmmTMB")
  check_and_load_package("emmeans")

  if (pattern.fill) check_and_load_package("ggpattern")
  
  microbAIDeR::mkdir(save.path)
  if( is.null(group) ){stop("No group argument supplied")}
  if (is.null(meta)) stop("For LMM usage you must provide 'meta' with rows matching sample order.")
  for (iterc in length(group))
  {
    if(!is.factor(meta[,group[iterc]])){stop("The supplied group argument inside the meta table is not a factor.")}
  }
  
  data <- data[rowSums(data) != 0, , drop=FALSE]
  if (ncol(data) == 0) stop("No samples remaining after filtering NAs from the grouping factor. Check your grouping factor")
  if (nrow(data) == 0) stop("No features remaining after filtering features with 0 abundance across all samples. Check your data")
  if (!p.adjust.method %in% p.adjust.methods) stop("Invalid p.adjust.method.")
  if (!plot.order %in% c("kruskal", "rownames")) stop("Invalid plot.order.")
  if( is.null(family_model) & data_type == "relabb" ){family_model = glmmTMB::beta_family(link="logit")}
  if( is.null(family_model) & data_type == "counts" ){family_model = glmmTMB::betabinomial(link = "logit")}
  
  # Ensure meta alignment when used
  if (!all(colnames(data) %in% rownames(meta))) stop("All sample names (colnames of data) must be present as rownames in meta.")
  meta <- meta[colnames(data), , drop=FALSE]
  if( length(group) == 1 )
  {
    grouping_factor = meta[,group]
    contrast_df <- data.frame(
      contrast = paste(levels(meta[,group]), collapse = " - ")
    )
  } else {
    comb <- expand.grid(lapply(meta[group], levels))
    pairs <- t(combn(seq_len(nrow(comb)), 2))
    contrast_df <- data.frame(
      contrast = apply(pairs, 1, function(idx) {
        paste(apply(comb[idx, ], 1, paste, collapse = " "), collapse = " - ")
      })
    )
    grouping_factor = as.factor(apply(meta[, group, drop = FALSE], 1, paste, collapse = " "))
  }
  
  if (plot.order == "kruskal") {
    kw <- apply(data, 1, function(x) kruskal.test(x, grouping_factor)$p.value)
    data <- data[order(ifelse(is.na(kw), 1, kw)), , drop=FALSE]
  }
  
  if( is.null(comparison.list) ){
    contrast_of_interest = contrast_df$contrast
  } else {
    if( !all(unlist(comparison.list) %in% unlist(strsplit(contrast_df$contrast, split = " - "))) )
    {
      message("Values supplied to the comparison.list are not included in the grouping columns:")
      paste( unlist(comparison.list)[!unlist(comparison.list) %in% unlist(strsplit(contrast_df$contrast, split = " - "))], collapse = " // " )
      stop("Exiting")
    } else {
      contrast_of_interest = sapply(comparison.list, function(x) paste(x, collapse = " - "))
    }
  }
  
  ptab <- data.frame(matrix(nrow = nrow(data), ncol = length(contrast_of_interest) + 4 * nlevels(grouping_factor)))
  rownames(ptab) <- rownames(data)
  colnames(ptab) <- c(
    contrast_of_interest,
    unlist(apply(expand.grid(c("Mean", "SEM", "Median", "SEMedian"), levels(grouping_factor)), 1, function(x) paste(x[1], x[2], sep = " ")))
  )
  
  # call log
  call.print = as.data.frame(rbind(
    p.adjust.method, taxlevel, save.path,
    mode = mode, 
    fixed_effects = paste(fixed_effects, collapse = ', '),
    random_effects = paste(random_effects, collapse =', '),
    family_model = paste(family_model$family_model, 'family_model, link:', family_model$link), 
    dispformula = dispformula,
    ziformula = ziformula,
    var.timepoint.random.slope = var.timepoint.random.slope,
    var.subject.random.intercept = var.subject.random.intercept,
    include_random_intercept_by_taxa = include_random_intercept_by_taxa,
    conflevel = conflevel,
    color.grouping = paste(color.grouping, collapse = ", "),
    smoothing.color, contrast.color,
    comparison.list = paste(sapply(comparison.list, function(x) paste(x, collapse = " vs ")), collapse=","),
    nrow.graph , ncol.graph , width.graph ,
    ggplot2.margins = paste(ggplot.margins, collapse = ", "),
    box.lwd , jitter.pch , jitter.stroke , jitter.size , jitter.color ,
    signif.step.increase, signif.text.size ,
    text.x.size , text.y.size , text.y.title.size,
    additional.params=paste(additional.params, collapse=", ")
  ))
  output_excel_path = file.path(save.path, paste0(taxlevel, "_LMM.xlsx"))
  write.xlsx(call.print, output_excel_path, sheetName = "call", colNames = FALSE, rowNames = TRUE)
  addSheet(output_excel_path, sheet.name = "kruskal-test", addition = as.data.frame(sort(kw, decreasing = F)), col.save = F, row.save = T)
  
  # helpers
  p_threshold <- if (trends) 0.1 else 0.05
  to_prop <- function(x) x/100
  sv_squeeze <- function(y) {
    # Smithson & Verkuilen (2006) squeeze to (0,1)
    n <- length(y)
    (y*(n-1) + 0.5)/n
  }
  
  if ( max(data) > 1 & data_type == "relabb" )
  {
    data_prop <- apply(data, 2, to_prop)              # features x samples in [0,1]
  } else { data_prop = data }
  if ( family_model$family == "beta" )
  {
    data_model <- apply(data_prop, 1:2, sv_squeeze)    # remove exact 0/1
  } else { data_model = data_prop }
  
  generate_plot <- function(tempdf, taxa, onlysig) {
    p <- ggplot(tempdf, aes(x = grouping_column, y = abundances)) +
      { if (pattern.fill) ggpattern::geom_boxplot_pattern(aes(pattern_fill = grouping_column), pattern = pattern, outlier.shape = NA, lwd = box.lwd, pattern_angle = pattern.angle, pattern_alpha = pattern.alpha, pattern_density = pattern.density, pattern_spacing = pattern.spacing, pattern_colour = color.grouping) else geom_boxplot(aes(fill = grouping_column), outlier.shape = NA, lwd = box.lwd) } +
      geom_jitter(aes(color = grouping_column), width = 0.1, height = 0, pch = jitter.pch, color = jitter.color, stroke = jitter.stroke, size = jitter.size) +
      { if (!is.null(onlysig)) ggsignif::geom_signif(comparisons = strsplit(onlysig$V1, " - "), annotations = sapply(onlysig$V2, sigFunction, trends = trends), size = signif.line.size, color = "grey22", textsize = signif.text.size, step_increase = signif.step.increase, tip_length = 0.02, margin_top = 0.05, vjust = 0.5) } +
      { if (smoothing) geom_smooth(aes(x = microbAIDeR::tokenize(grouping_column)), method = smoothing.method, lwd = smoothing.lwd, color = contrast.color, se = smoothing.se) } +
      scale_fill_manual(values = color.grouping) +
      labs(x = "", y = taxa) +
      { if (horiz) coord_flip() } +
      theme(legend.position = 'none', plot.margin = unit(ggplot.margins, "cm"), axis.text.x = element_text(size = text.x.size), axis.text.y = element_text(size = text.y.size), axis.title.y = element_text(size = text.y.title.size)) +
      { if (!is.null(additional.params)) additional.params }
    if (align.legend) p <- align_legend(p)
    p
  }
  
  clean_fixed_effects_single <- function(fixed_effects, taxa_var = "taxa_name") {
    norm <- function(s) {
      s <- gsub("\\s+", " ", s)
      s <- gsub("\\s*([*:])\\s*", " \\1 ", s)  # normalize around * and :
      trimws(gsub("\\s+", " ", s))
    }
    strip_taxa <- function(term) {
      t0 <- norm(term)
      if (!grepl(paste0("\\b", taxa_var, "\\b"), t0)) return(t0)
      
      # If the whole term is just taxa_var (possibly with spaces), drop it
      if (grepl(paste0("^\\b", taxa_var, "\\b$"), t0)) return("")
      
      has_star  <- grepl("\\*", t0)
      has_colon <- grepl(":",  t0)
      
      # choose splitter by operator precedence in the term
      if (has_star) {
        parts <- strsplit(t0, "\\*", perl = TRUE)[[1]]
        parts <- trimws(parts)
        parts <- parts[parts != "" & parts != taxa_var]
        if (!length(parts)) return("")
        return(paste(parts, collapse = " * "))
      } else if (has_colon) {
        parts <- strsplit(t0, ":", fixed = TRUE)[[1]]
        parts <- trimws(parts)
        parts <- parts[parts != "" & parts != taxa_var]
        # keep interaction only if >=2 factors remain; else drop
        if (length(parts) >= 2) return(paste(parts, collapse = ":")) else return("")
      } else {
        # some other odd case that still includes taxa_var → drop
        return("")
      }
    }
    
    out <- vapply(fixed_effects, strip_taxa, character(1))
    out <- out[out != ""]
    unname(out)
  }

  
  # containers
  gvec <- list()
  gvec_corr <- list()

  
  # ======= PER-FEATURE MODE =======
  if ("singlefeature" %in% mode) {
    message("Computing per feature model")
    nfeat <- length(rownames(data))
    pb <- utils::txtProgressBar(min = 0, max = nfeat, style = 3)
    i <- 0
    
    if( any(grepl("taxa_name", fixed_effects)) )
    {
      fixed_effects_single = clean_fixed_effects_single(fixed_effects, taxa_var = "taxa_name")
    } else { fixed_effects_single = fixed_effects}
    
    mode1_formula <- build_glmmTMB_spec(
      response = "abundances",
      fixed_effects = fixed_effects_single,
      var.subject.random.intercept = var.subject.random.intercept,
      var.timepoint.random.slope = var.timepoint.random.slope,
      taxa_var = NULL,
      include_random_intercept_by_taxa = FALSE, # per-taxa model
      dispformula = dispformula,
      ziformula = ziformula,
      family_model = family_model
    )
    
    for (iter_feature in rownames(data)) {
      i <- i + 1
      # build per-feature data with covariates
      tempdf <- cbind.data.frame(
        abundances = as.numeric(data_model[iter_feature, ]),
        taxa_name = iter_feature,
        meta
      )
      # fit
      fit <- try(glmmTMB::glmmTMB(formula = mode1_formula$formula, 
                                  ziformula = mode1_formula$ziformula, 
                                  family = mode1_formula$family_model, 
                                  data = tempdf), silent = TRUE)
      if (inherits(fit, "try-error")) {
        # fallback: no ZI
        fit <- try(glmmTMB::glmmTMB(formula = mode1_formula$formula, 
                                    family = mode1_formula$family_model, 
                                    data = tempdf), silent = TRUE)
      }
      # pairwise p-values on group
      if (inherits(fit, "try-error")) {
        # if still failing, mark NA / skip plot
        pvals <- rep(1, length(contrast_of_interest))
      } else {
        emm <- try(emmeans::emmeans(fit, as.formula(paste("~", fixed_effects_single[1], collapse = " ")), level = conflevel), silent = TRUE)
        if (inherits(emm, "try-error")) {
          pvals <- rep(1, length(contrast_of_interest))
        } else {
          contr_table <- emmeans::contrast(emm, method = "pairwise")
          contr <- data.frame(
            row.names = summary(contr_table)$contrast,
            p = summary(contr_table)$p.value
          )
          pvals = contr[contrast_of_interest,]
          names(pvals) <- contrast_of_interest
        }
      }
      pvals[is.na(pvals) | is.infinite(pvals)] <- 1
      
      if ( isTRUE(save_single_model_summaries) )
      {
        microbAIDeR::mkdir(paste(save.path, "singlemode_summaries", sep = "/"))
        sink(paste(save.path, "/singlemode_summaries/", iter_feature, ".txt",  sep = ""))
        print(summary(fit))
        print(paste("\n\n##################################################################\n\n"))
        print(summary(emm))
        print(paste("\n\n##################################################################\n\n"))
        print(summary(contr_table))
        sink()
      }
      
      # fill table + summary stats
      ptab[iter_feature, 1:length(contrast_of_interest)] <- signif(pvals, 3)
      for (lvl in levels(grouping_factor)) {
        vv <- as.numeric(data[iter_feature, grouping_factor == lvl])
        ptab[iter_feature, paste("Mean", lvl)]     <- signif(mean(vv, na.rm=TRUE), 3)
        ptab[iter_feature, paste("SEM", lvl)]      <- signif(microbAIDeR::sem(vv), 3)
        ptab[iter_feature, paste("Median", lvl)]   <- signif(median(vv, na.rm=TRUE), 3)
        ptab[iter_feature, paste("SEMedian", lvl)] <- signif(microbAIDeR::se.median(vv), 3)
      }
      
      onlysig <- if (any(pvals <= p_threshold))
        data.frame(V1 = unlist(lapply(comparison.list[pvals <= p_threshold], paste, collapse = " - ")),
                   V2 = round(pvals[pvals <= p_threshold], 5)) else NULL
      
      # for plotting we need a simple df with abundances + grouping
      plotdf <- data.frame(abundances = as.numeric(data[iter_feature, ]), grouping_column = grouping_factor)
      if (is.null(onlysig)) { if (plot.not.sig) gvec <- c(gvec, list(generate_plot(plotdf, iter_feature, onlysig))) } else { gvec <- c(gvec, list(generate_plot(plotdf, iter_feature, onlysig))) }
      
      utils::setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # save uncorrected
    file_uncorrected = file.path(save.path, paste0(taxlevel, "_uncorrected.pdf"))
    message("Saving uncorrected boxplots to ", file_uncorrected)
    suppressWarnings(ggsave(
      filename = file_uncorrected,
      plot = marrangeGrob(gvec, nrow = nrow.graph, ncol = ncol.graph, top = NULL,
                          layout_matrix = matrix(seq_len(nrow.graph * ncol.graph), nrow = nrow.graph, ncol = ncol.graph, byrow = TRUE)),
      width = width.graph, height = height.graph, dpi = 330
    ))
    message("Adding uncorrected pvalues to ", output_excel_path)
    addSheet(path = output_excel_path, sheet.name = paste0("uncorrected"), addition = ptab, col.save = TRUE, row.save = TRUE)
    
    # cleaned
    clean_ptab <- function(x, contrast_of_interest, p_threshold, digits = 3) {
      for (col in seq_along(contrast_of_interest)) {
        x[x[, col] > p_threshold, col] <- ""
        x[, col] <- ifelse(x[, col] != "", signif(as.numeric(x[, col]), digits = digits), "")
      }
      return(x)
    }
    cleaned_ptab = clean_ptab(ptab, contrast_of_interest = contrast_of_interest, p_threshold = p_threshold)
    addSheet(path = output_excel_path, sheet.name = paste0("uncorrected_clean"), addition = cleaned_ptab, col.save=TRUE, row.save=TRUE)
    
    # correction
    if (p.adjust.method != "none") {
      corrected_ptab <- ptab
      for (col in 1:length(contrast_of_interest)) {
        corrected_ptab[, col] <- p.adjust(as.numeric(corrected_ptab[, col]), method = p.adjust.method)
      }
      addSheet(path = output_excel_path, sheet.name = paste0(p.adjust.method, ""), addition = corrected_ptab , col.save = TRUE, row.save = TRUE)
      cleaned_corrected_ptab = clean_ptab(corrected_ptab, contrast_of_interest = contrast_of_interest, p_threshold = p_threshold)
      addSheet(path= output_excel_path, sheet.name = paste(p.adjust.method, "clean", sep="_"), addition = cleaned_corrected_ptab, col.save=TRUE, row.save=TRUE)
      
      # plots with corrected p
      gvec_corr <- list()
      for (taxa in rownames(data)) {
        plotdf <- data.frame(abundances = as.numeric(data[taxa, ]), grouping_column = grouping_factor)
        wtests_corrected <- as.numeric(corrected_ptab[taxa, 1:nrow(contrast_df)])
        onlysig_corrected <- if ( any(wtests_corrected <= p_threshold) )
          data.frame(V1 = unlist(lapply(comparison.list[wtests_corrected <= p_threshold], paste, collapse = " - ")),
                     V2 = round(wtests_corrected[wtests_corrected <= p_threshold], 5)) else NULL
        if (is.null(onlysig_corrected)) { if (plot.not.sig) gvec_corr <- c(gvec_corr, list(generate_plot(plotdf, taxa, onlysig_corrected))) }
        else { gvec_corr <- c(gvec_corr, list(generate_plot(plotdf, taxa, onlysig_corrected))) }
      }
      file_corrected = file.path(save.path, paste0(taxlevel, "_", p.adjust.method,".pdf"))
      message("Saving corrected boxplots to ", file_corrected)
      suppressWarnings(ggsave(
        filename = file_corrected,
        plot = marrangeGrob(gvec_corr, nrow = nrow.graph, ncol = ncol.graph, top = NULL,
                            layout_matrix = matrix(seq_len(nrow.graph * ncol.graph), nrow = nrow.graph, ncol = ncol.graph, byrow = TRUE)),
        width = width.graph, height = height.graph, dpi = 330
      ))
    }
    message("Computing completed.")
    return(invisible(NULL))
  }
  
  # ======= TAXA MODE =======
  # Whole-table model once (all features), then feature selection, then refit
  if (mode == "wholetable") {
    message("Computing wholetable model")

    mode2_formula <- build_glmmTMB_spec(
      response = "abundances",
      fixed_effects = fixed_effects,
      var.subject.random.intercept = var.subject.random.intercept,
      var.timepoint.random.slope = var.timepoint.random.slope,
      taxa_var = "taxa_name",
      include_random_intercept_by_taxa = TRUE, # whole table model
      dispformula = dispformula,
      ziformula = ziformula,
      family_model = family_model
    )
    
    backup_mode2 = mode2_formula
    mode2_formula$formula = as.formula(
      "abundances ~ Group * Timepoint * taxa_name + age + gender + (1 + Timepoint | Subject_ID)"
    )
    mode2_formula$ziformula = as.formula("~ taxa_name + (1|Subject_ID)")
    mode2_formula$dispformula = as.formula("~ 0")
    
    wholedata <- microbAIDeR::tidify(as.data.frame(data_model), colNames = c("taxa_name", "placeholder_sampleid", "abundances"))
    meta_join <- meta %>% tibble::rownames_to_column(var = "placeholder_sampleid")
    wholedata <- wholedata %>%
      left_join(meta_join, by = "placeholder_sampleid")
    
    fit2 <- try(glmmTMB::glmmTMB(formula = mode2_formula$formula, 
                                ziformula = mode2_formula$ziformula, 
                                family = mode2_formula$family_model, 
                                data = wholedata), silent = TRUE)
    
    emm2 <- emmeans::emmeans(fit2, ~ taxa_name * (Group * Timepoint), level = conflevel)
    contr_table2 <- emmeans::contrast(emm2, method = "pairwise")
    
    contr2 <- contrast(emm2,
                      method = list(
                        "NEG T0 - NEG T1" = c(1, -1, 0, 0),
                        "NEG T0 - POS T0" = c(1, 0, -1, 0),
                        "POS T0 - POS T1" = c(0, 1, 0, -1),
                        "NEG T1 - POS T1" = c(0, 0, 1, -1)
                      ),
                      by = "taxa_name", adjust = "fdr"
    )
    as.data.frame(contr)
    
    
    
    
    
    
    if (inherits(fit, "try-error")) {
      # fallback: no ZI
      fit <- try(glmmTMB::glmmTMB(formula = mode1_formula$formula, 
                                  family = mode1_formula$family_model, 
                                  data = tempdf), silent = TRUE)
    }
    # pairwise p-values on group
    if (inherits(fit, "try-error")) {
      # if still failing, mark NA / skip plot
      pvals <- rep(1, length(contrast_of_interest))
    } else {
      emm <- try(emmeans::emmeans(fit, as.formula(paste("~", fixed_effects_single[1], collapse = " ")), level = conflevel), silent = TRUE)
      if (inherits(emm, "try-error")) {
        pvals <- rep(1, length(contrast_of_interest))
      } else {
        contr_table <- emmeans::contrast(emm, method = "pairwise")
        contr <- data.frame(
          row.names = summary(contr_table)$contrast,
          p = summary(contr_table)$p.value
        )
        pvals = contr[contrast_of_interest,]
        names(pvals) <- contrast_of_interest
      }
    }
    pvals[is.na(pvals) | is.infinite(pvals)] <- 1
    
    if ( isTRUE(save_single_model_summaries) )
    {
      microbAIDeR::mkdir(paste(save.path, "singlemode_summaries", sep = "/"))
      sink(paste(save.path, "/singlemode_summaries/", iter_feature, ".txt",  sep = ""))
      print(summary(fit))
      print(paste("\n\n##################################################################\n\n"))
      print(summary(emm))
      print(paste("\n\n##################################################################\n\n"))
      print(summary(contr_table))
      sink()
    }
    
    # fill table + summary stats
    ptab[iter_feature, 1:length(contrast_of_interest)] <- signif(pvals, 3)
    for (lvl in levels(grouping_factor)) {
      vv <- as.numeric(data[iter_feature, grouping_factor == lvl])
      ptab[iter_feature, paste("Mean", lvl)]     <- signif(mean(vv, na.rm=TRUE), 3)
      ptab[iter_feature, paste("SEM", lvl)]      <- signif(microbAIDeR::sem(vv), 3)
      ptab[iter_feature, paste("Median", lvl)]   <- signif(median(vv, na.rm=TRUE), 3)
      ptab[iter_feature, paste("SEMedian", lvl)] <- signif(microbAIDeR::se.median(vv), 3)
    }
    
    onlysig <- if (any(pvals <= p_threshold))
      data.frame(V1 = unlist(lapply(comparison.list[pvals <= p_threshold], paste, collapse = " - ")),
                 V2 = round(pvals[pvals <= p_threshold], 5)) else NULL
    
    # for plotting we need a simple df with abundances + grouping
    plotdf <- data.frame(abundances = as.numeric(data[iter_feature, ]), grouping_column = grouping_factor)
    if (is.null(onlysig)) { if (plot.not.sig) gvec <- c(gvec, list(generate_plot(plotdf, iter_feature, onlysig))) } else { gvec <- c(gvec, list(generate_plot(plotdf, iter_feature, onlysig))) }
    
    
    
    
    
    
    # 2) Produce boxplots using p-values from FIRST model — we’ll compute per-feature emmeans p for plotting labels
    # (ptab is filled with per-feature p-values vs pairwise group contrasts)
    for (taxa in rownames(data)) {
      plotdf <- data.frame(abundances = as.numeric(data[taxa, ]), grouping_column = group)
      # p via GLMM emmeans
      tmp <- cbind.data.frame(abundances = plotdf$abundances, meta)
      fml <- build_glmm_formula("abundances", fixed_effects, random_effects)
      fit <- try(glmmTMB::glmmTMB(formula = fml, ziformula = ziformula_glmm, family_model = family_model_glmm, data = tmp), silent=TRUE)
      pvals <- rep(1, nrow(contrast_df))
      if (!inherits(fit,"try-error")) {
        emm <- try(emmeans::emmeans(fit, ~ grouping_column), silent=TRUE)
        if (!inherits(emm,"try-error")) {
          contr <- emmeans::contrast(emm, method="pairwise")$contrasts
          name_map <- setNames(contr@grid$contrast, seq_len(nrow(contr)))
          p_series <- summary(contr)$p.value; names(p_series) <- name_map
          want_names <- sapply(comparison.list, function(v) paste(v, collapse="-"))
          pvals <- sapply(seq_along(comparison.list), function(i) {
            nm1 <- want_names[i]; nm2 <- paste(rev(strsplit(nm1, "-", fixed=TRUE)[[1]]), collapse="-")
            if (nm1 %in% names(p_series)) p_series[[nm1]] else if (nm2 %in% names(p_series)) p_series[[nm2]] else 1
          })
        }
      }
      pvals[is.na(pvals) | is.infinite(pvals)] <- 1
      
      ptab[taxa, 1:nrow(contrast_df)] <- signif(pvals, 3)
      for (lvl in levels(group)) {
        vv <- as.numeric(data[taxa, group == lvl])
        ptab[taxa, paste("Mean", lvl)]     <- signif(mean(vv, na.rm=TRUE), 3)
        ptab[taxa, paste("SEM", lvl)]      <- signif(microbAIDeR::sem(vv), 3)
        ptab[taxa, paste("Median", lvl)]   <- signif(median(vv, na.rm=TRUE), 3)
        ptab[taxa, paste("SEMedian", lvl)] <- signif(microbAIDeR::se.median(vv), 3)
      }
      onlysig <- if (any(pvals <= p_threshold))
        data.frame(V1 = unlist(lapply(comparison.list[pvals <= p_threshold], paste, collapse = " vs ")),
                   V2 = round(pvals[pvals <= p_threshold], 5)) else NULL
      gvec <- c(gvec, list(generate_plot(plotdf, taxa, onlysig)))
    }
    
    # Save first-pass plots and table
    file_uncorrected = file.path(save.path, paste0(taxlevel, "_taxa_first.pdf"))
    message("Saving first whole-table pass plots to ", file_uncorrected)
    suppressWarnings(ggsave(
      filename = file_uncorrected,
      plot = marrangeGrob(gvec, nrow = nrow.graph, ncol = ncol.graph, top = NULL,
                          layout_matrix = matrix(seq_len(nrow.graph * ncol.graph), nrow = nrow.graph, ncol = ncol.graph, byrow = TRUE)),
      width = width.graph, height = height.graph, dpi = 330
    ))
    addSheet(path = output_excel_path, sheet.name = "taxa_first_pass", addition = ptab, col.save = TRUE, row.save = TRUE)
    
    # 3) Refit second whole-table model on selected features (report selection list)
    addSheet(path = output_excel_path, sheet.name = "taxa_selected_features", addition = data.frame(feature=selected_features), col.save=TRUE, row.save=FALSE)
    
    # optional: a compact joint model summary on selected features (mvabund preferred)
    if (requireNamespace("mvabund", quietly=TRUE)) {
      Y2 <- t(as.matrix(data[selected_features, , drop=FALSE]))
      fx_terms <- c("grouping_column", fixed_effects)
      rhs <- paste(fx_terms[!is.na(fx_terms) & nzchar(fx_terms)], collapse = " + ")
      fml_mv2 <- as.formula(paste("Y2 ~", ifelse(nzchar(rhs), rhs, "1")))
      d <- data.frame(meta)
      fit_sel <- try(mvabund::manyglm(fml_mv2, data = d, family_model = "negative.binomial"), silent=TRUE)
      if (!inherits(fit_sel,"try-error")) {
        summ2 <- capture.output(suppressWarnings(mvabund::anova.manyglm(fit_sel, p.uni="adjusted")))
        writeLines(summ2, con = file.path(save.path, paste0(taxlevel,"_taxa_selected_manyglm.txt")))
      }
    }
    
    message("Computing completed.")
    return(invisible(NULL))
  }
}
