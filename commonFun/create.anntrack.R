create.anntrack <- function(samples=NULL, subtype=NULL, typesToPlot=NULL, text.col=NULL) {
  if (is.null(samples)) {stop("samples can not be NULL!")}
  if (length(samples)!=length(subtype)) {stop("samples and subtype do not have equal length!")}
  
  if (is.null(typesToPlot)) {
    typesToPlot <- levels(factor(unique(subtype)))
  }
  if (is.null(text.col)) {
    text.col <- rainbow(length(typesToPlot))
  }  
  names(subtype) <- samples
  names(text.col) <- levels(factor(typesToPlot))  
  col <- text.col[subtype]
  return(list(subtype=subtype, text.col=text.col, typesToPlot=typesToPlot, col=col))
}