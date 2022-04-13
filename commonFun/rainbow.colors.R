rainbow.colors <- function(n) {
  # Return a rainbox colors color-map (vector of colors)
  # Need to clip R 'rainbow' range at 11/15 to avoid red wraparound
  gamma=1
  black=0
  icolor=0
  colormap0 <- rev(rainbow(n, s=1, v=1, start=0, end=11/15, alpha=gamma))
  
  # black should be supported on other schemes too...
  # (it is not documented for this reason.)
  if ({black}) {
    # add black + interpolate smooth transition from violet to black color
    # #00000000 -> #6600FFFF
    # color choice is not optimal, by any means
    colormap <- c('black', '#330066FF', '#660099FF', '#9900CCFF', colormap0)
  } else {
    colormap <- colormap0
  }
  
  # invert color order
  if ({icolor}) colormap <- rev(colormap);
  
  colormap
}
