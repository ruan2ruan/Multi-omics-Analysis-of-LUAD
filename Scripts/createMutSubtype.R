createMutSubtype <- function(indata=NULL, samples=NULL, genename=NULL) {
  if (!is.element(genename, rownames(indata))) { stop (paste(genename, "not found in indata!", sep=" ")) }
  
  comsam <- intersect(colnames(indata), samples)
  
  out <- rep("Not Available", times=length(samples))
  names(out) <- samples
  out[comsam] <- as.character(indata[genename, comsam])
  out[is.na(out)] <- "Not Available"
  out[out=="0"] <- "Normal"
  out[out=="1"] <- "Mutated"  
  
  return(out)
}
