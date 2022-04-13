plot.boxplot =function(indata, fig.fullname, col.vector="red", xlab=NULL, ylab=NULL, type="pdf", outline=TRUE, title=NULL, las=0, cex.axis=1, xlables=NULL, pos="topright", legend=NULL, fill=NULL, y0flag=FALSE, ylim=NULL){
    if (type=="jpg") { jpeg(file=fig.fullname, width = 2000, height = 2000, quality=100, pointsize=16) }
    if (type=="pdf") { pdf(file=fig.fullname) }
    
    if(y0flag==TRUE && is.null(ylim)) {      
      tmp <- boxplot(indata, plot=FALSE)
      ylim <- c(0, max(tmp$stats[5, ], na.rm=TRUE))
    }
    
    par(las=las)
    par(cex=cex.axis)
    if(is.null(xlables)){
      boxplot(indata, col=col.vector, xlab=xlab, ylab=ylab, outline=outline, main=title, ylim=ylim)
    }else{
      boxplot(indata, col=col.vector, xlab=xlab, ylab=ylab, outline=outline, main=title, xaxt="n", ylim=ylim)
      axis(side=1,at=(1:length(xlables))[xlables!=""],labels=xlables[xlables!=""])
    }
    if (!is.null(legend)) {
      legend(pos, legend=legend, fill=fill)
    }
    if (type!="none") { invisible(dev.off()) }    
}
