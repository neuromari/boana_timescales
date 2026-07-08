################## B. pulchellus Timescales of Call Variability #################
#### author: M. Rodriguez-Santiago
#### last update: 6.17.2026
####
#### this script contains all functions necessary for analysis  
#### and color schemes for the boana timescales characterization manuscript
#### this scripts need only be called once and stored in environment
####
####
####
#
#### color palettes #####
#
## by timescale
#
timescales <- c("#99cc99","#cc99cc","#66cccc") #notes, calls, series
TimeSeries <- c("seriesdur" = "#2a5674",
                "seriesisi" = "#66cccc", 
                "seriesiso" = "#d1eeea")
TimeCalls <- c("callsdur" = "#9f839d", 
               "callsici" = "#cfa9d1", 
               "callsioi" = "#e4cee6")
TimeNotes <- c("notesdur_note1" = "#839c7e", "notesdur_note2" = "#839c7e",
               "notesini_note1" = "#99cc99", 
               "notesino_note1" = "#cde4cc")

CVColors <- c("#ccffff","#99cccc")  #CV within series
#
#
## by vocal feature 
#
SeriesTypeColors <- c("#d1eeea", "#85c4c9","#66cccc","#01665e") #solo two three four+
#
CallColors3 <- c("#778868","#de8a5a","black") #doublet squeak uc
CallColors2 <- c("#778868","#de8a5a") #doublet squeak
CallBinColors <- c("#e4cee6","#cfa9d1","#cc99cc","#9f839d") #three four five 
#
NoteColors <- c("#cde4cc","#839c7e","#de8a5a","black")
NoteBinColors <- c("#bcdab7","#67b461","#336633","#003300") #1 50 100 0 
#
#
#### functions ####
#
## sem
#
sem <- function(x) sd(x, na.rm=T)/sqrt(length(x))
#
## calculate mode
#
get_mode <- function(x) {
  if(length(x) == 0) return(NA)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
#
## function for patchwork plot
#
plot_seq <- function(data, show_y_axis = TRUE, bg_color = "white") {
  p <- ggplot(data,
              aes(x = grand_CVw, y = CV_ratio,
                  color = feature_label, shape = feature_label)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = bg_color, alpha = 0.3) +
    geom_hline(yintercept = 1, linetype = "dashed",
               color = "gray60", linewidth = 0.4) +
    geom_point(size = 8, alpha = 1) +
    scale_color_manual(values = label_colors_seq) +
    scale_shape_manual(values = label_shapes_seq) +
    xlim(0, 75) +
    ylim(0, 5) +
    xlab("CVw (%)") +
    coord_cartesian(clip = "off") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line        = element_line(colour = "black", size = 0.3),
      legend.position  = "none",
      axis.title.x     = element_text(size = 26),
      axis.text.x      = element_text(size = 20),
      axis.text.y      = element_text(size = 20),
      plot.margin      = margin(t = 10, r = 20, b = 10, l = 10)
    )
  
  if (show_y_axis) {
    p <- p + ylab("CVb / CVw") +
      theme(axis.title.y = element_text(size = 26))
  } else {
    p <- p + ylab(NULL) +
      theme(axis.title.y = element_blank(),
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y  = element_blank())
  }
  p
}
#
## function to find troughs
#
find_trough <- function(x, bw = "nrd0", n_grid = 2048) {
  x <- x[is.finite(x)]
  out <- list(trough = NA, mode1 = NA, mode2 = NA,
              n_modes = NA, bimodal = FALSE, sep_depth = NA)
  if (length(x) < 5 || diff(range(x)) == 0) return(out)
  
  d  <- density(x, bw = bw, n = n_grid)
  gx <- d$x; gy <- d$y
  
  # Local maxima (modes) and minima via slope sign changes on the density grid
  dy    <- diff(gy)
  ismax <- which(diff(sign(dy)) == -2) + 1
  out$n_modes <- length(ismax)
  
  if (length(ismax) < 2) {            # unimodal at this bandwidth -> no trough
    out$mode1 <- gx[which.max(gy)]
    return(out)
  }
  
  # Two tallest modes; trough = lowest-density point between them
  ord <- ismax[order(gy[ismax], decreasing = TRUE)]
  m   <- sort(ord[1:2])               # left-to-right
  seg <- m[1]:m[2]
  ti  <- seg[which.min(gy[seg])]
  
  out$trough    <- gx[ti]
  out$mode1     <- gx[m[1]]
  out$mode2     <- gx[m[2]]
  out$bimodal   <- TRUE
  out$sep_depth <- 1 - gy[ti] / min(gy[m[1]], gy[m[2]])  # 0 = shallow, ->1 = deep/clean
  out
}
#
## function to find trough boundaries
#
find_boundaries <- function(counts, breaks, trough) {
  n          <- length(counts)
  bin_center <- (breaks[1:n] + breaks[2:(n + 1)]) / 2
  
  if (is.na(trough)) {                 # unimodal: one distribution over all non-empty bins
    nz <- which(counts > 0)
    return(list(dist1_min = breaks[min(nz)], dist1_max = breaks[max(nz) + 1],
                dist2_min = NA, dist2_max = NA, trough = NA))
  }
  
  left  <- which(bin_center <  trough)
  right <- which(bin_center >= trough)
  ln <- left[counts[left]   > 0]
  rn <- right[counts[right] > 0]
  
  list(
    dist1_min = if (length(ln)) breaks[min(ln)]     else NA,
    dist1_max = if (length(ln)) breaks[max(ln) + 1] else NA,
    dist2_min = if (length(rn)) breaks[min(rn)]     else NA,
    dist2_max = if (length(rn)) breaks[max(rn) + 1] else NA,
    trough    = trough)
}
