

theme_panel <- function(base_family = "Helvetica") {
  ggplot2::theme_classic(base_size = 7, base_family = base_family) %+replace%
    ggplot2::theme(
      # ---- Text ----
      text              = ggplot2::element_text(color = "black"),
      plot.title        = ggplot2::element_text(size = 7, face = "plain", hjust = 0.5, margin = ggplot2::margin(b = 1.5, unit = "mm")),
      plot.subtitle     = ggplot2::element_text(size = 6.5, hjust = 0.5, margin = ggplot2::margin(b = 1.5, unit = "mm")),
      plot.caption      = ggplot2::element_text(size = 6, hjust = 0.5, margin = ggplot2::margin(t = 1.5, unit = "mm")),
      
      axis.title        = ggplot2::element_text(size = 7),
      axis.title.x      = ggplot2::element_text(margin = ggplot2::margin(t = 1.2, unit = "mm")),
      axis.title.y      = ggplot2::element_text(margin = ggplot2::margin(r = 1.2, unit = "mm")),
      
      axis.text         = ggplot2::element_text(size = 6),
      axis.text.x       = ggplot2::element_text(margin = ggplot2::margin(t = 0, unit = "mm")),
      axis.text.y       = ggplot2::element_text(margin = ggplot2::margin(r = 0.6, unit = "mm")),
      
      # ---- Legend (compact) ----
      legend.title      = ggplot2::element_text(size = 7),
      legend.text       = ggplot2::element_text(size = 6),
      legend.key.size   = grid::unit(3, "mm"),
      legend.key.height = grid::unit(3, "mm"),
      legend.key.width  = grid::unit(3, "mm"),
      legend.spacing.x  = grid::unit(1, "mm"),
      legend.spacing.y  = grid::unit(1, "mm"),
      legend.box.spacing= grid::unit(1, "mm"),
      legend.margin     = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
      
      legend.position   = "right",
      legend.direction  = "vertical",
      legend.background = ggplot2::element_blank(),
      
      # ---- Panel sizing / margins ----
      plot.margin       = ggplot2::margin(2.0, 1.0, 2.0, 1.0, unit = "mm"),
      
      # ---- Strips (facet labels) ----
      strip.background  = ggplot2::element_blank(),
      strip.text        = ggplot2::element_text(size = 6.5, face = "plain", margin = ggplot2::margin(b = 1, t = 1, unit = "mm"))
    )
}
