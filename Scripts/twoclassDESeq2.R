twoclassDESeq2 <- function(res.path=NULL, Ginfo=NULL, countsTable=NULL, tailrows=NULL, Batchinfo=NULL, features=NULL, featType=NULL, complist=NULL, PASSFlag=NULL, overwt=FALSE) {
  
  #Groupinfo could contain "batch", which will be considered by edgeR design matrix
  sampleList <- complist[[1]]
  treatsamList <- complist[[2]]
  treatnameList <- complist[[3]]
  ctrlnameList <- complist[[4]]
  allsamples <- colnames(countsTable)
  
  options(warn=1)
  for (k in 1:length(sampleList)) {
    samples <- sampleList[[k]]
    treatsam <- treatsamList[[k]]
    treatname <- treatnameList[k]
    ctrlname <- ctrlnameList[k]
    
    compname <- paste(treatname, "_vs_", ctrlname, sep="")
    tmp = rep("others", times=length(allsamples))
    names(tmp) <- allsamples
    tmp[samples]="control"
    tmp[treatsam]="treatment"
    outfile <- file.path( res.path, paste(featType, "_deseq2_test_result.", compname, ".txt", sep="") )
    if (file.exists(outfile) & (overwt==FALSE)) {
      cat(k, ":", compname, "exists and skipped;\n")
      next
    }
    if(is.null(Batchinfo)) {
      saminfo <- data.frame("Type"=tmp[samples],"SampleID"=samples,stringsAsFactors = F)
      cts <- countsTable[setdiff(rownames(countsTable),tailrows),samples]
      coldata <- saminfo[samples,]
      dds <- DESeqDataSetFromMatrix(countData = cts,
                                    colData = coldata,
                                    design = as.formula("~ Type"))
          } else {
      saminfo <- data.frame("Type"=tmp[samples],"Batch"=Batchinfo[samples],stringsAsFactors = F)
      cts <- countsTable[setdiff(rownames(countsTable),tailrows),samples]
      coldata <- saminfo[samples,]
      dds <- DESeqDataSetFromMatrix(countData = cts,
                                    colData = coldata,
                                    design = as.formula(paste("~", paste(c("Batch","Type"), collapse=" + "))))
    }
    
    keep <- intersect(features, names(PASSFlag[PASSFlag==TRUE]) )
    dds <- dds[keep,]
    
    dds$Type <- relevel(dds$Type,ref = "control")
    
    dds <- DESeq(dds)
    res <- results(dds, contrast=c("Type","treatment","control"))
    
    resData <- as.data.frame(res[order(res$padj),])
    resData$id <- Ginfo[rownames(resData),"genename"]
    resData <- resData[!duplicated(resData$id),]
    resData <- resData[,c("id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
    
    write.table(resData, file=outfile, row.names=F, col.names=T, sep="\t", quote=F)
    cat(k, ",")
  }
  options(warn=0)
  
}