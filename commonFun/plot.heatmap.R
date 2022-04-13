plot.heatmap=function(indata, halfwidth=NULL, outfile, trsfmFlag="none", ColCor=NULL, img.fmt="jpg", cex=1) {
## in the input indata, row:gene, column:sample
library(ClassDiscovery)
library(gplots)
indata=na.omit(indata)
hcs <- hclust(distanceMatrix(as.matrix(indata), "pearson"), "ward")
hcg <- hclust(distanceMatrix(as.matrix(t(indata)), "pearson"), "ward")
#hcs <- hclust(distanceMatrix(as.matrix(indata), "euclidean"), "ward")
#hcg <- hclust(distanceMatrix(as.matrix(t(indata)), "euclidean"), "ward")

if (toupper(trsfmFlag) == "BOTH") {
    plotdata=t(scale(t(indata), center=T, scale=T))
}
if (toupper(trsfmFlag)=="CENTER"){
  plotdata=t(scale(t(indata), center=T, scale=F))
}

if (!(toupper(trsfmFlag) %in% c("CENTER", "BOTH"))){
  plotdata=indata
}    

if (!(is.null(halfwidth))) {
    plotdata[plotdata>halfwidth]=halfwidth
    plotdata[plotdata<(-halfwidth)]= -halfwidth
}
if (img.fmt=="jpg") {
    jpeg(file=outfile, width = 900, height = 900, quality=100, pointsize=16)
}else{
    pdf(outfile)
}
if (is.null(ColCor)) {
    heatmap.2(plotdata, scale="none", Rowv=as.dendrogram(hcg), Colv=as.dendrogram(hcs), col=greenred(64), trace="none", cexRow=cex, cexCol=cex)
}else{
    heatmap.2(plotdata, scale="none", Rowv=as.dendrogram(hcg), Colv=as.dendrogram(hcs), col=greenred(64), trace="none", ColSideColors=ColCor, cexRow=cex, cexCol=cex)
}        
invisible(dev.off())
}