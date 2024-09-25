align_legend <- function(p, hjust = 0.5)
{
  check_and_load_package("cowplot")
  check_and_load_package("gtable")
  g <- cowplot::plot_to_gtable(p)
  grobs <- g$grobs
  if ( any(sapply(grobs, function(x) "guide-box" %in% x$name)) )
  {
    legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
    legend <- grobs[[legend_index]]
    
    guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")
    
    for (gi in guides_index) {
      guides <- legend$grobs[[gi]]
      # add extra column for spacing
      # guides$width[5] is the extra spacing from the end of the legend text
      # to the end of the legend title. If we instead distribute it by `hjust:(1-hjust)` on
      # both sides, we get an aligned legend
      spacing <- guides$width[5]
      guides <- gtable::gtable_add_cols(guides, hjust*spacing, 1)
      guides$widths[6] <- (1-hjust)*spacing
      title_index <- guides$layout$name == "title"
      guides$layout$l[title_index] <- 2
      
      legend$grobs[[gi]] <- guides
    }
    
    g$grobs[[legend_index]] <- legend
    g
  }
  
}