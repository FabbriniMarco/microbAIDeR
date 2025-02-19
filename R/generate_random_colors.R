generate_random_colors <- function(x, start.col = "grey90", seed = 1996, show.colors = FALSE) {
  check_and_load_package("RColorBrewer")
  set.seed(seed)
  n <- x-1
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  if (x > length(col_vector)) {
    message(paste0("Error: exceeding number of requested colors (max = ", length(col_vector), ")"))
    return(NA)
  }
  colori=c(start.col, sample(col_vector, n))
  if (show.colors == TRUE) {
    par(mgp = c(3, 0.2, 0))
    barplot(rep(1, x), col = colori, border = NA,
            main = paste(n_colors, "Distinct Random Colors"), 
            names.arg = 1:x, yaxt = "n", cex.names = 0.6)
  }
  return(colori)
} 