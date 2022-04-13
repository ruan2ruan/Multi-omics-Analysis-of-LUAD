#-----------------#
# MOVICS pipeline #

# set working path and creat directory
tumor.path <- "E:/LUAD/movics_pipeline"
setwd(tumor.path) #create dir
res.path    <- file.path(tumor.path, "Results")
fig.path    <- file.path(tumor.path, "Figures")
data.path  <- file.path(tumor.path, "Data")
comAnn.path <- file.path(tumor.path,"Annotation")
script.path <- file.path(tumor.path,"Scripts")
comRFun.path <- file.path(tumor.path,"commonFun")

if (!file.exists(tumor.path)) { dir.create(tumor.path) }
if (!file.exists(res.path)) { dir.create(res.path) }
if (!file.exists(fig.path)) { dir.create(fig.path) }

# set colors
blue   <- "#5bc0eb"
yellow <- "#fde74c"
green  <- "#9bc53d"
red    <- "#f25f5c"
purple <- "#531f7a"
grey   <- "#8693ab"
orange <- "#fa7921"
white  <- "#f2d7ee"
darkred   <- "#F2042C"
lightred  <- "#FF7FBF"
lightblue <- "#B2EBFF"
darkblue  <- "#1d00ff"
cherry    <- "#700353"
lightgrey <- "#dcddde"
nake <- "#F8C364"
gold <- "#ECE700"
cyan <- "#00B3D0"
sun  <- "#E53435"
peach  <- "#E43889"
violet <- "#89439B"
soil   <- "#EC7D21"
lightgreen <- "#54B642"
darkblue   <- "#21498D"
darkgreen  <- "#009047"
brown      <- "#874118"
seagreen   <- "#008B8A"
jco <- c("#2874C5","#EABF00","#868686","#C6524A","#80A7DE")
heatmap.BlWtRd <- c("#6699CC","white","#FF3C38")
heatmap.YlGnPe <- c("#440259","#345F8C","#228C8A","#78CE51","#FAE71F")
heatmap.GrWtRd <- c("#2b2d42","#8d99ae","#edf2f4","#ef233c","#d90429")
heatmap.L.BlYlRd <- c("#4281a4","#9cafb7","#ead2ac","#e6b89c","#fe938c")
heatmap.BlBkRd <- c("#54FEFF","#32ABAA","#125456","#000000","#510000","#A20000","#F30000")
heatmap.fancy <- c("#10040A", "#2A0B35", "#4D155B", "#73215B", "#9C3558", "#C34D44", "#E07038", "#F2981C", "#F2CA51", "#FAF6A3")

mutect.dataframe <- function(x){
  # delete rows of Silent
  #cut_id <- x$Variant_Classification %in% c("Silent")
  sel_id <- x$Variant_Type == "SNP" #& x$Variant_Classification %in% c("Silent","Missense_Mutation","","Nonsense_Mutation","Nonstop_Mutation")
  #x <- x[!cut_id,]
  x <- x[sel_id,]
  somatic_sum <- x %>% group_by(Tumor_Sample_Barcode) %>% summarise(TCGA_sum = n())
}

indel.dataframe <- function(x){
  # delete rows of Silent
  #cut_id <- x$Variant_Classification %in% c("Silent")
  sel_id <- x$Variant_Type %in% c("DEL","INS")
  #x <- x[!cut_id,]
  x <- x[sel_id,]
  somatic_sum <- x %>% group_by(Tumor_Sample_Barcode) %>% summarise(TCGA_sum = n())
}

# load R package
library(MOVICS)
library(ggplot2)
library(RColorBrewer)

# load MOVICS input data
load("luad.tcga.RData")

# print name of example data
names(luad.tcga)
# [1] "mRNA"      "lncRNA"    "meth"      "mut"       "count"     "tpm"       "fpkm"      "maf"       "segment"  
# [10] "surv.info"       "segment"     "clin.info"

# extract multi-omics data
mo.data   <- luad.tcga[1:4]

# extract raw count data for downstream analyses
count     <- luad.tcga$count

# extract tpm data for downstream analyses
tpm       <- luad.tcga$tpm

# extract maf for downstream analysis
maf       <- luad.tcga$maf

# extract segmented copy number for downstream analyses
segment   <- luad.tcga$segment

# extract survival information
surv.info <- luad.tcga$surv.info

#----------------------------------------------------------#
# 1. identify optimal clustering number (may take a while) #
optk.luad <- getClustNum(data        = mo.data,
                         is.binary   = c(F,F,F,T), # note: the 4th data is somatic mutation which is a binary matrix
                         try.N.clust = 2:6, # try cluster number from 2 to 6
                         fig.path    = fig.path,
                         fig.name    = "CLUSTER NUMBER OF TCGA-LUAD")
# take cluster number as 4 then

# 2. perform multi-omics integrative clustering with the 10 algorithms #
moic.res.list5 <- getMOIC(data        = mo.data,
                             N.clust     = 4,
                             methodslist = list("iClusterBayes","SNF","ConsensusClustering", "CIMLR", "MoCluster"), # specify only ONE algorithm here
                             type        = c("gaussian","gaussian","gaussian","binomial"))
save(moic.res.list5, file = file.path(res.path,"moic.res.list5.rda"))

#------------------------------------------------------#
# 3. get consensus clustering from different algorithm #
cmoic.luad <- getConsensusMOIC(moic.res.list = moic.res.list5,
                               fig.name      = "CONSENSUS HEATMAP（5methods）",
                               fig.path      = fig.path,
                               distance      = "pearson",
                               linkage       = "ward.D2")

getSilhouette(sil      = cmoic.luad$sil, # a sil object returned by getConsensusMOIC()
              fig.path = fig.path,
              fig.name = "SILHOUETTE（5methods）",
              height   = 5.5,
              width    = 5)

#------------------------------#
# 4. get comprehensive heatmap #
# convert beta value to M value for stronger signal
indata <- mo.data
indata$meth <- log2(indata$meth / (1 - indata$meth))

# set color for each omics data
mRNA.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
lncRNA.col <- c("#6699CC", "white"  , "#FF3C38")
meth.col   <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
mut.col    <- c("grey90" , "black")
col.list   <- list(mRNA.col, lncRNA.col, meth.col, mut.col)

# extract sample annotation
tmp <- cliquery[,c("gender",
                   "number_pack_years_smoked")]
colnames(tmp) <- c("gender","Smoke")
tmp1 <- cbind.data.frame(surv.info,tmp[comsam,])
tmp1 <- tmp1[,-8]
tmp1[which(tmp1$pStage %in% c("Stage IA","Stage IB")),"pStage"] <- "Stage I"
tmp1[which(tmp1$pStage %in% c("Stage IIA","Stage IIB")),"pStage"] <- "Stage II"
tmp1[which(tmp1$pStage %in% c("Stage IIIA","Stage IIIB")),"pStage"] <- "Stage III"
tmp1[tmp1 == ""] <- NA
tmp1[which(tmp1$Race %in% c("AMERICAN INDIAN OR ALASKA NATIVE","BLACK OR AFRICAN AMERICAN")),"Race"] <- "Others"
annCol    <- tmp1[,c("Age",  "Race","pStage","Gender","Smoke","TP53"), drop = FALSE]
annCol[is.na(annCol)] <- "Missing"
annCol$Age <- as.numeric(annCol$Age)
annCol$Smoke <- as.numeric(annCol$Smoke)
annColors <- list(Age    = circlize::colorRamp2(breaks = c(min(annCol$Age,na.rm = TRUE),
                                                           median(annCol$Age,na.rm = TRUE),
                                                           max(annCol$Age,na.rm = TRUE)), 
                                                colors = c("#0000AA", "#555555", "#AAAA00")),
                  pStage = c("Stage I"    = alpha(darkred,0.2),
                             "Stage II"   = alpha(darkred,0.4),
                             "Stage III"  = alpha(darkred,0.7),
                             "Stage IV"   = alpha(darkred,0.9), 
                             "Missing"    = "white"),
                  Race   = c("ASIAN" = yellow,
                             "WHITE" = "grey80",
                             "Others" = brown,
                             "Missing" = "white"),
                  Gender = c("FEMALE" = jco[1],
                             "MALE"   = jco[2]),
                  Smoke    = circlize::colorRamp2(breaks = c(min(annCol$Smoke,na.rm = TRUE),
                                                             max(annCol$Smoke,na.rm = TRUE)), 
                                                  colors = c("#2ec4b6", "#e71d36")),
                  TP53   = c("Mutated" = "black",
                             "Wild" = "grey90"))
# data normalization for heatmap
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2,2,NA), # no truncation for mutation
                     centerFlag = c(T,T,T,F), # no center for mutation
                     scaleFlag  = c(T,T,T,F)) # no scale for mutation
library(impute)
plotdata$meth <- impute.knn(plotdata$meth)$data
feat   <- moic.res.list5[["CIMLR"]][["feat.res"]]
feat1  <- feat[which(feat$dataset == "mRNA"),][1:10,"feature"] 
feat2  <- feat[which(feat$dataset == "lncRNA"),][1:10,"feature"]
feat3  <- feat[which(feat$dataset == "meth"),][1:10,"feature"]
feat4  <- feat[which(feat$dataset == "mut"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3, feat4)
# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.TPM","lncRNA.TPM","M value","Mutated"),
             clust.res     = cmoic.luad$clust.res, # consensusMOIC results
             clust.dend    = NULL, # show no dendrogram for samples
             show.rownames = c(F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             show.row.dend = c(F,F,F,F), # show dendrogram for features
             annRow        = annRow, # no selected features
             color         = col.list,
             annCol        = annCol, # annotation for samples
             annColors     = annColors, # annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.path      = fig.path,
             fig.name      = "COMPREHENSIVE HEATMAP OF CONSENSUSMOIC (5methods features+Smoke)")


surv.info[which(surv.info$pStage %in% c("Stage IA","Stage IB")),"pStage"] <- "Stage I"
surv.info[which(surv.info$pStage %in% c("Stage IIA","Stage IIB")),"pStage"] <- "Stage II"
surv.info[which(surv.info$pStage %in% c("Stage IIIA","Stage IIIB")),"pStage"] <- "Stage III"
surv.info[surv.info == ""] <- NA
surv.info[which(surv.info$Race %in% c("AMERICAN INDIAN OR ALASKA NATIVE","BLACK OR AFRICAN AMERICAN")),"Race"] <- "Others"

annCol    <- surv.info[,c("Age", "Gender", "Race","pStage","TP53"), drop = FALSE]
annCol[is.na(annCol)] <- "Missing"
annCol$Age <- as.numeric(annCol$Age)

# generate corresponding colors for sample annotation
annColors <- list(Age    = circlize::colorRamp2(breaks = c(min(annCol$Age,na.rm = TRUE),
                                                           median(annCol$Age,na.rm = TRUE),
                                                           max(annCol$Age,na.rm = TRUE)), 
                                                colors = c("#0000AA", "#555555", "#AAAA00")),
                  Gender = c("FEMALE" = jco[1],
                             "MALE"   = jco[2]),
                  pStage = c("Stage I"    = alpha(darkred,0.2),
                             "Stage II"   = alpha(darkred,0.4),
                             "Stage III"  = alpha(darkred,0.7),
                             "Stage IV"   = alpha(darkred,0.9), 
                             "Missing"    = "white"),
                  Race   = c("ASIAN" = yellow,
                             "WHITE" = "grey80",
                             "Others" = brown,
                             "Missing" = "white"),
                  TP53   = c("Mutated" = "black",
                             "Wild" = "grey90"))
# data normalization for heatmap
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2,2,NA), # no truncation for mutation
                     centerFlag = c(T,T,T,F), # no center for mutation
                     scaleFlag  = c(T,T,T,F)) # no scale for mutation

library(impute)
plotdata$meth <- impute.knn(plotdata$meth)$data

feat   <- moic.res.list5[["CIMLR"]][["feat.res"]]
feat1  <- feat[which(feat$dataset == "mRNA"),][1:10,"feature"] 
feat2  <- feat[which(feat$dataset == "lncRNA"),][1:10,"feature"]
feat3  <- feat[which(feat$dataset == "meth"),][1:10,"feature"]
feat4  <- feat[which(feat$dataset == "mut"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3, feat4)

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.TPM","lncRNA.TPM","M value","Mutated"),
             clust.res     = cmoic.luad$clust.res, # consensusMOIC results
             clust.dend    = NULL, # show no dendrogram for samples
             show.rownames = c(F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             show.row.dend = c(F,F,F,F), # show dendrogram for features
             annRow        = annRow, # no selected features
             color         = col.list,
             annCol        = annCol, # annotation for samples
             annColors     = annColors, # annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.path      = fig.path,
             fig.name      = "COMPREHENSIVE HEATMAP OF CONSENSUSMOIC (5methods features)")

#---------------------#
# 5. compare survival #
surv.luad <- compSurv(moic.res         = cmoic.luad,
                      surv.info        = surv.info,
                      convt.time       = "m", # convert day unit to month
                      surv.median.line = "h", # draw horizontal line at median survival
                      xyrs.est         = c(5,10), # estimate 5 and 10-year survival
                      surv.cut         = 120,
                      p.adjust.method  = "none",
                      fig.path         = fig.path,
                      fig.name         = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC (5methods)")

#------------------------------#
# 6. compare clinical features #
clin.luad <- compClinvar(moic.res      = cmoic.luad,
                         var2comp      = tmp1, # data.frame needs to summarize (must has row names of samples)
                         strata        = "Subtype", # stratifying variable (e.g., Subtype in this example)
                         factorVars    = c("Race","pStage","fustat","TP53","Gender"), # features that are considered categorical variables
                         nonnormalVars = c("futime","Age"), # feature(s) that are considered using nonparametric test
                         exactVars     = c("Race","pStage","fustat","TP53","Gender"), # feature(s) that are considered using exact test
                         doWord        = TRUE, # generate .docx file in local path
                         res.path      = res.path,
                         includeNA     = FALSE,
                         tab.name      = "SUMMARIZATION OF CLINICAL FEATURES (5methods)")
print(clin.luad$compTab)
write.csv(clin.luad$compTab,"clincompTab (5methods+smoke).csv")
#------------------------------------#
# 7. mutational frequency comparison #
mut.luad <- compMut(moic.res     = cmoic.luad,
                    mut.matrix   = luad.tcga$mut, # binary somatic mutation matrix
                    doWord       = TRUE, # generate table in .docx format
                    doPlot       = TRUE, # draw OncoPrint
                    freq.cutoff  = 0.18, # keep those genes that mutated in at least 3% of samples
                    #                    p.cutoff     = 0.05, # keep those genes with nominal p value < 0.25 to draw OncoPrint
                    p.adj.cutoff = 0.05, # keep those genes with adjusted p value < 1 to draw OncoPrint
                    innerclust   = TRUE, # perform clustering within each subtype
                    annCol       = annCol, # same annotation for heatmap
                    annColors    = annColors, # same annotation color for heatmap
                    width        = 8, 
                    height       = 10,
                    fig.name     = "ONCOPRINT FOR SIGNIFICANT MUTATIONS (5methods)",
                    tab.name     = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION (5methods)",
                    res.path     = res.path,
                    fig.path     = fig.path)
print(mut.luad)

#----------------------------------#
# 7. compare total mutation burden #
tmb.luad <- compTMB(moic.res     = cmoic.luad,
                    maf          = maf,
                    rmDup        = TRUE, # remove duplicated variants per sample
                    rmFLAGS      = FALSE, # keep FLAGS mutations
                    exome.size   = 38, # estimated exome size
                    test.method  = "nonparametric", # statistical testing method
                    fig.path     = fig.path,
                    fig.name     = "DISTRIBUTION OF TMB AND TITV (5methods)")

#------------------------------------#
# 8. compare fraction genome altered #
# change column names of segment data
colnames(segment) <- c("sample","chrom","start","end","value")

# compare FGA, FGG, and FGL
fga.luad <- compFGA(moic.res     = cmoic.luad,
                    segment      = segment,
                    iscopynumber = FALSE, # this is a segmented copy number file
                    cnathreshold = 0.2, # threshold to determine CNA gain or loss
                    test.method  = "nonparametric", # statistical testing method
                    fig.path     = fig.path,
                    fig.name     = "BARPLOT OF FGA (5methods)")

#--------------------------------#
# 9. drug sensitivity comparison #
drug.luad <- compDrugsen(moic.res    = cmoic.luad,
                         norm.expr   = tpm[,cmoic.luad$clust.res$samID], # double guarantee sample order
                         drugs       = c("Cisplatin","Paclitaxel", "Docetaxel", "Vinorelbine", "Gemcitabine"), # a vector of names of drug in GDSC (I get this from https://www.cancer.org/cancer/liver-cancer/treating/chemotherapy.html)
                         tissueType  = "digestive_system", # choose specific tissue type to construct ridge regression model
                         test.method = "nonparametric", # statistical testing method
                         fig.path    = fig.path,
                         prefix      = "BOXVIOLIN OF ESTIMATED IC50 (5methods)") 

drug.luad1 <- compDrugsen(moic.res    = cmoic.luad,
                         norm.expr   = tpm[,cmoic.luad$clust.res$samID], # double guarantee sample order
                         drugs       = c("Afatinib","Gefitinib","Erlotinib"), # a vector of names of drug in GDSC (I get this from https://www.cancer.org/cancer/liver-cancer/treating/chemotherapy.html)
                         tissueType  = "digestive_system", # choose specific tissue type to construct ridge regression model
                         test.method = "nonparametric", # statistical testing method
                         fig.path    = fig.path,
                         prefix      = "BOXVIOLIN OF ESTIMATED IC50 (5methods)") 

#-------------------------------------------#
# 10. compare agreement with other subtypes #
# agreement comparison (support up to 6 classifications include current subtype)
agree.luad <- compAgree(moic.res  = cmoic.luad,
                        subt2comp = annCol[,c("Gender","TP53","pStage"),drop = F],
                        doPlot    = TRUE,
                        fig.path  = fig.path,
                        box.width = 0.2,
                        fig.name  = "AGREEMENT OF CONSENSUSMOIC (5methods)")

#------------------------------------------#
# 11. run differential expression analysis #
# run DEA with edgeR
runDEA(dea.method = "edger",
       expr       = count, # raw count data
       moic.res   = cmoic.luad,
       res.path   = res.path,
       prefix     = "TCGA-LUAD") # prefix of figure name

# run DEA with DESeq2
runDEA(dea.method = "deseq2",
       expr       = count, # raw count data
       moic.res   = cmoic.luad,
       res.path   = res.path,
       prefix     = "TCGA-LUAD") # prefix of figure name

# run DEA with limma
runDEA(dea.method = "limma",
       expr       = tpm, # raw count data
       moic.res   = cmoic.luad,
       res.path   = res.path,
       prefix     = "TCGA-LUAD") # prefix of figure name

#--------------------------------------------#
# 12. run biomarker identification procedure #
# choose limma result to identify subtype-specific up-regulated biomarkers
marker.up <- runMarker(moic.res      = cmoic.luad,
                       dea.method    = "deseq2", # name of DEA method
                       prefix        = "TCGA-LUAD", # MUST be the same of argument in runDEA()
                       dat.path      = res.path, # path of DEA files
                       res.path      = res.path, # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       n.marker      = 100, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = tpm, # use normalized expression as heatmap input
                       annCol        = annCol, # sample annotation in heatmap
                       annColors     = annColors, # colors for sample annotation
                       show_rownames = FALSE, # show no rownames (biomarker name)
                       fig.path      = fig.path,
                       width         = 12,
                       height        = 12,
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP (5methods)")

# choose limma result to identify subtype-specific down-regulated biomarkers
marker.dn <- runMarker(moic.res      = cmoic.luad,
                       dea.method    = "deseq2",
                       dat.path      = res.path, # path of DEA files
                       res.path      = res.path, # path to save marker files
                       prefix        = "TCGA-LUAD",
                       dirct         = "down",
                       n.marker      = 100, # switch to 50
                       doplot        = TRUE,
                       norm.expr     = tpm,
                       annCol        = annCol,
                       annColors     = annColors,
                       fig.path      = fig.path,
                       width         = 12,
                       height        = 12,
                       fig.name      = "DOWNREGULATED BIOMARKER HEATMAP (5methods)")

#--------------------------------------#
# 13. run gene set enrichment analysis #
MSIGDB.FILE1 <- "c5.go.bp.v7.4.symbols.gmt"
MSIGDB.FILE2 <- "c6.all.v7.4.symbols.gmt"
MSIGDB.FILE3 <- "c7.immunesigdb.v7.4.symbols.gmt"

# run GSEA to identify up-regulated GO pathways using results from deseq2
gsea.up <- runGSEA(moic.res     = cmoic.luad,
                   dea.method   = "deseq2", # name of DEA method
                   prefix       = "TCGA-LUAD", # MUST be the same of argument in runDEA()
                   dat.path     = res.path, # path of DEA files
                   res.path     = res.path, # path to save GSEA files
                   msigdb.path  = MSIGDB.FILE1, # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = tpm, # use normalized expression to calculate enrichment score
                   dirct        = "up", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.25, # padj cutoff to identify significant pathways
                   gsva.method  = "ssgsea", # method to calculate single sample enrichment score
                   norm.method  = "median", # normalization method to calculate subtype-specific enrichment score
                   fig.path     = fig.path,
                   fig.name     = "UPREGULATED PATHWAY HEATMAP (GO-BP)")
# run GSEA to identify down-regulated GO pathways using results from deseq2
gsea.dn <- runGSEA(moic.res     = cmoic.luad,
                   dea.method   = "deseq2", # name of DEA method
                   prefix       = "TCGA-LUAD", # MUST be the same of argument in runDEA()
                   dat.path     = res.path, # path of DEA files
                   res.path     = res.path, # path to save GSEA files
                   msigdb.path  = MSIGDB.FILE1, # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = tpm, # use normalized expression to calculate enrichment score
                   dirct        = "down", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.25, # padj cutoff to identify significant pathways
                   gsva.method  = "ssgsea", # method to calculate single sample enrichment score
                   norm.method  = "median", # normalization method to calculate subtype-specific enrichment score
                   fig.path     = fig.path,
                   fig.name     = "DOWNREGULATED PATHWAY HEATMAP (GO-BP)")

gsea.up <- runGSEA(moic.res     = cmoic.luad,
                   dea.method   = "deseq2", # name of DEA method
                   prefix       = "TCGA-LUAD", # MUST be the same of argument in runDEA()
                   dat.path     = res.path, # path of DEA files
                   res.path     = res.path, # path to save GSEA files
                   msigdb.path  = MSIGDB.FILE2, # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = tpm, # use normalized expression to calculate enrichment score
                   dirct        = "up", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.25, # padj cutoff to identify significant pathways
                   gsva.method  = "ssgsea", # method to calculate single sample enrichment score
                   norm.method  = "median", # normalization method to calculate subtype-specific enrichment score
                   fig.path     = fig.path,
                   fig.name     = "UPREGULATED PATHWAY HEATMAP (C6)")
gsea.dn <- runGSEA(moic.res     = cmoic.luad,
                   dea.method   = "deseq2", # name of DEA method
                   prefix       = "TCGA-LUAD", # MUST be the same of argument in runDEA()
                   dat.path     = res.path, # path of DEA files
                   res.path     = res.path, # path to save GSEA files
                   msigdb.path  = MSIGDB.FILE2, # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = tpm, # use normalized expression to calculate enrichment score
                   dirct        = "down", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.25, # padj cutoff to identify significant pathways
                   gsva.method  = "ssgsea", # method to calculate single sample enrichment score
                   norm.method  = "median", # normalization method to calculate subtype-specific enrichment score
                   fig.path     = fig.path,
                   fig.name     = "DOWNREGULATED PATHWAY HEATMAP (C6)")

gsea.up <- runGSEA(moic.res     = cmoic.luad,
                   dea.method   = "deseq2", # name of DEA method
                   prefix       = "TCGA-LUAD", # MUST be the same of argument in runDEA()
                   dat.path     = res.path, # path of DEA files
                   res.path     = res.path, # path to save GSEA files
                   msigdb.path  = GSET.FILE, # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = tpm, # use normalized expression to calculate enrichment score
                   dirct        = "up", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.25, # padj cutoff to identify significant pathways
                   gsva.method  = "ssgsea", # method to calculate single sample enrichment score
                   norm.method  = "median", # normalization method to calculate subtype-specific enrichment score
                   fig.path     = fig.path,
                   fig.name     = "UPREGULATED PATHWAY HEATMAP (HALLMARK)")
gsea.dn <- runGSEA(moic.res     = cmoic.luad,
                   dea.method   = "deseq2", # name of DEA method
                   prefix       = "TCGA-LUAD", # MUST be the same of argument in runDEA()
                   dat.path     = res.path, # path of DEA files
                   res.path     = res.path, # path to save GSEA files
                   msigdb.path  = GSET.FILE, # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = tpm, # use normalized expression to calculate enrichment score
                   dirct        = "down", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.25, # padj cutoff to identify significant pathways
                   gsva.method  = "ssgsea", # method to calculate single sample enrichment score
                   norm.method  = "median", # normalization method to calculate subtype-specific enrichment score
                   fig.path     = fig.path,
                   fig.name     = "DOWNREGULATED PATHWAY HEATMAP (HALLMARK)")

GSET.FILE <- "h.all.v7.4.symbols.gmt"
gsva.res <- 
  runGSVA(moic.res      = cmoic.luad,
          norm.expr     = tpm,
          gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          annCol        = annCol,
          annColors     = annColors,
          fig.path      = fig.path,
          fig.name      = "HALLMARK OF INTEREST HEATMAP (5methods)",
          height        = 10,
          width         = 12)

#################### GSEA ####################
library(clusterProfiler)
library(enrichplot)

MSigDB=read.gmt(file.path(comAnn.path,"h.all.v7.4.symbols.gmt"))

CS23 <- read.table(file.path(res.path,"consensusMOIC_TCGA-LUAD_deseq2_test_result.CS2_vs_CS3.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
geneList <- CS23$log2fc
names(geneList) <- rownames(CS23)
geneList <- sort(geneList,decreasing = T)
GSEA.CS23 <- GSEA(geneList = geneList,TERM2GENE=MSigDB,seed = T,verbose=F,pvalueCutoff = 0.25)
pdf(file.path(fig.path,"GSEA.CS23.h_bubble.pdf"),width = 8,height = 7)
dotplot(GSEA.CS23, showCategory=15)
dev.off()

pdf(file.path(fig.path,"GSEA EMT enrichplot.pdf"),width = 8,height = 6)
gseaplot2(GSEA.CS23, geneSetID = c(3),base_size=9,
          color=c("#1b9e77","#7570b3","#66a61e","#e6ab02","#377eb8","#ec7014"),
          rel_heights = c(1.7, 0.5, 0.8))
dev.off()
pdf(file.path(fig.path,"GSEA KRAS enrichplot.pdf"),width = 8,height = 6)
gseaplot2(GSEA.CS23, geneSetID = c(13),base_size=9,
          color=c("#1b9e77","#7570b3","#66a61e","#e6ab02","#377eb8","#ec7014"),
          rel_heights = c(1.7, 0.5, 0.8))
dev.off()

tmp <- read.table(file.path(res.path,"consensusMOIC_TCGA-LUAD_deseq2_test_result.CS1_vs_CS4.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
geneList <- tmp$log2fc
names(geneList) <- rownames(tmp)
geneList <- sort(geneList,decreasing = T)
GSEA.CS14 <- GSEA(geneList = geneList,TERM2GENE=MSigDB,seed = T,verbose=F,pvalueCutoff = 0.25)
pdf(file.path(fig.path,"GSEA.CS14_bubble.pdf"),width = 9,height =8)
dotplot(GSEA.CS14, showCategory=15)
dev.off()

###############Volcano##############
library( "ggplot2" )

######## CS2 VS CS3 ########
x <- read.table(file.path(res.path,"consensusMOIC_TCGA-LUAD_deseq2_test_result.CS2_vs_CS3.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
x$label<- rownames(x)
head(x)
#colnames(x) <- c("baseMean","logFC","lfcSE","stat","pvalue","padj","Genesymbol","label")
colnames(x) <- c("baseMean","logFC","pvalue","padj","label")
#plot_mode <- "classic" 
plot_mode <- "advanced"

logFCcut <- log2(2) 
logFCcut2 <- log2(4) #for advanced mode

adjPcut <- 0.25
adjPcut2 <- 0.05 #for advanced mode

DEG_up <- x[x$log2FoldChange > log2(2) & x$padj<0.05,]
DEG_DN <- x[x$log2FoldChange < -log2(2) & x$padj<0.05,]

xmin <- -4
xmax <- 4
ymin <- 0
ymax <- 6

if (plot_mode == "classic"){
  # setting for color
  x$color_transparent <- ifelse((x$padj < adjPcut & x$logFC > logFCcut), "red", ifelse((x$padj < adjPcut & x$logFC < -logFCcut), "blue","grey30"))
  # setting for size
  size <- ifelse((x$padj < adjPcut & abs(x$logFC) > logFCcut), 4, 2)
  
} else if (plot_mode == "advanced") {
  # setting for color
  n1 <- length(x[, 1])
  cols <- rep("grey30", n1)
  names(cols)<- rownames(x)
  
  cols[x$padj < adjPcut & x$logFC >logFCcut]<- "#EBE645"
  cols[x$padj < adjPcut2 & x$logFC > logFCcut2]<- "#EBE645"
  cols[x$padj < adjPcut & x$logFC < -logFCcut]<- "#577BC1"
  cols[x$padj < adjPcut2 & x$logFC < -logFCcut2]<- "#577BC1"
  color_transparent <- adjustcolor(cols, alpha.f = 0.5)
  x$color_transparent <- color_transparent
  
  # setting for size
  n1 <- length(x[, 1])
  size <- rep(1, n1)
  
  size[x$padj < adjPcut & x$logFC > logFCcut]<- 3
  size[x$padj < adjPcut2 & x$logFC > logFCcut2]<- 5
  size[x$padj < adjPcut & x$logFC < -logFCcut]<- 3
  size[x$padj < adjPcut2 & x$logFC < -logFCcut2]<- 5
  
  
} else {
  stop("Unsupport mode")
}

# Construct the plot object
p1 <- ggplot(data=x, aes(logFC, -log10(padj), label = label, color = pathway)) +
  geom_point(alpha = 0.6, size = size, colour = x$color_transparent) +
  
  labs(x=bquote(~log[2]~"(FoldChange)"), y=bquote(~-log[10]~"FDR"), title="") + 
  scale_y_continuous(
    breaks = c(0, -round(log10(adjPcut),1),round(-log10(adjPcut2),1),2.5,5), 
    labels = c(0, -round(log10(adjPcut),1),round(-log10(adjPcut2),1),2.5,5),
    limits = c(ymin,ymax)
  ) +
  scale_x_continuous(
    breaks = c(-4,-3, -round(logFCcut2,1),-round(logFCcut,1), 0, round(logFCcut,1),round(logFCcut2,1),3, 4),
    labels = c(-4,-3, -round(logFCcut2,1),"-1", 0, "1",round(logFCcut2,1),3, 4),
    limits = c(-4, 4) 
  ) +
  geom_vline(xintercept = c(-logFCcut, logFCcut), color="grey40", 
             linetype="longdash", lwd = 0.5) + 
  geom_hline(yintercept = -log10(adjPcut), color="grey40", 
             linetype="longdash", lwd = 0.5) +
  
  theme_bw(base_size = 12#, base_family = "Times"
  ) +
  theme(panel.grid=element_blank()) +
  theme(legend.position="right",
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "Times", size=8),
        plot.title = element_text(hjust = 0.8),
        axis.text.x = element_text(face="bold", color="black", size=11),
        axis.text.y = element_text(face="bold",  color="black", size=11),
        axis.title.x = element_text(face="bold", color="black", size=11),
        axis.title.y = element_text(face="bold",color="black", size=11))

if (plot_mode == "advanced") {
  p1 <- p1 + 
    geom_vline(xintercept = c(-logFCcut2, logFCcut2), color="grey40", 
               linetype="longdash", lwd = 0.5) +
    geom_hline(yintercept = -log10(adjPcut2), color="grey40", 
               linetype="longdash", lwd = 0.5)
}

ggsave(file.path(fig.path,"volcano plot for DEGs between CS2 and CS3.pdf"),width = 8,height = 6)

######## CS1 VS CS4 ########
x <- read.table(file.path(res.path,"consensusMOIC_TCGA-LUAD_deseq2_test_result.CS1_vs_CS4.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
x$label<- rownames(x)
head(x)
#colnames(x) <- c("baseMean","logFC","lfcSE","stat","pvalue","padj","Genesymbol","label")
colnames(x) <- c("baseMean","logFC","pvalue","padj","label")
#plot_mode <- "classic" 
plot_mode <- "advanced"

logFCcut <- log2(2) 
logFCcut2 <- log2(4) #for advanced mode

adjPcut <- 0.25
adjPcut2 <- 0.05 #for advanced mode

DEG_up <- x[x$log2FoldChange > log2(2) & x$padj<0.05,]
DEG_DN <- x[x$log2FoldChange < -log2(2) & x$padj<0.05,]

xmin <- -4
xmax <- 4
ymin <- 0
ymax <- 6

if (plot_mode == "classic"){
  # setting for color
  x$color_transparent <- ifelse((x$padj < adjPcut & x$logFC > logFCcut), "red", ifelse((x$padj < adjPcut & x$logFC < -logFCcut), "blue","grey30"))
  # setting for size
  size <- ifelse((x$padj < adjPcut & abs(x$logFC) > logFCcut), 4, 2)
  
} else if (plot_mode == "advanced") {
  # setting for color
  n1 <- length(x[, 1])
  cols <- rep("grey30", n1)
  names(cols)<- rownames(x)
  
  cols[x$padj < adjPcut & x$logFC >logFCcut]<- "#EBE645"
  cols[x$padj < adjPcut2 & x$logFC > logFCcut2]<- "#EBE645"
  cols[x$padj < adjPcut & x$logFC < -logFCcut]<- "#577BC1"
  cols[x$padj < adjPcut2 & x$logFC < -logFCcut2]<- "#577BC1"
  color_transparent <- adjustcolor(cols, alpha.f = 0.5)
  x$color_transparent <- color_transparent
  
  # setting for size
  n1 <- length(x[, 1])
  size <- rep(1, n1)
  
  size[x$padj < adjPcut & x$logFC > logFCcut]<- 3
  size[x$padj < adjPcut2 & x$logFC > logFCcut2]<- 5
  size[x$padj < adjPcut & x$logFC < -logFCcut]<- 3
  size[x$padj < adjPcut2 & x$logFC < -logFCcut2]<- 5
  
  
} else {
  stop("Unsupport mode")
}

# Construct the plot object
p1 <- ggplot(data=x, aes(logFC, -log10(padj), label = label, color = pathway)) +
  geom_point(alpha = 0.6, size = size, colour = x$color_transparent) +
  
  labs(x=bquote(~log[2]~"(FoldChange)"), y=bquote(~-log[10]~"FDR"), title="") + 
  scale_y_continuous(
    breaks = c(0, -round(log10(adjPcut),1),round(-log10(adjPcut2),1),2.5,5), 
    labels = c(0, -round(log10(adjPcut),1),round(-log10(adjPcut2),1),2.5,5),
    limits = c(ymin,ymax)
  ) +
  scale_x_continuous(
    breaks = c(-4,-3, -round(logFCcut2,1),-round(logFCcut,1), 0, round(logFCcut,1),round(logFCcut2,1),3, 4),
    labels = c(-4,-3, -round(logFCcut2,1),"-1", 0, "1",round(logFCcut2,1),3, 4),
    limits = c(-4, 4) 
  ) +
  geom_vline(xintercept = c(-logFCcut, logFCcut), color="grey40", 
             linetype="longdash", lwd = 0.5) + 
  geom_hline(yintercept = -log10(adjPcut), color="grey40", 
             linetype="longdash", lwd = 0.5) +
  
  theme_bw(base_size = 12#, base_family = "Times"
  ) +
  theme(panel.grid=element_blank()) +
  theme(legend.position="right",
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "Times", size=8),
        plot.title = element_text(hjust = 0.8),
        axis.text.x = element_text(face="bold", color="black", size=11),
        axis.text.y = element_text(face="bold",  color="black", size=11),
        axis.title.x = element_text(face="bold", color="black", size=11),
        axis.title.y = element_text(face="bold",color="black", size=11))

if (plot_mode == "advanced") {
  p1 <- p1 + 
    geom_vline(xintercept = c(-logFCcut2, logFCcut2), color="grey40", 
               linetype="longdash", lwd = 0.5) +
    geom_hline(yintercept = -log10(adjPcut2), color="grey40", 
               linetype="longdash", lwd = 0.5)
}

ggsave(file.path(fig.path,"volcano plot for DEGs between CS1 and CS4.pdf"),width = 8,height = 6)

#################### GSVA ######################
library( "pheatmap" )
library(GSVA)
library("gplots")

TPM.HUGO <- read.table(file.path(res.path,"TCGA_LUAD_mRNA_TPM.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
met_sig <- read.table(file.path(res.path,"some immnue signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL,quote = "")
label <- unique(met_sig$`Metabolic Function`)
meta_sig <- list()
for (i in label) {
  meta_sig[[i]] <- met_sig[which(met_sig$`Metabolic Function` == i),"Gene Name"]  
}
meta.score <- gsva(as.matrix(TPM.HUGO),meta_sig,method="gsva")
write.table(as.data.frame(meta.score),file.path(res.path,"gsva results for some signatures by gsva.txt"),row.names = T,col.names = NA,sep = "\t")

# heatmap
Immune_Cell_Signature <- read.table(file.path(res.path,"gsva results for some signatures by gsva.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
Immune_Cell_Signature <- Immune_Cell_Signature[,rownames(annCol)]
choose_matrix = standarize.fun(indata=Immune_Cell_Signature,halfwidth = 1.5)
mycol <- colorpanel(256,low=cyan,mid = "black",high=peach) 

pheatmap(choose_matrix,color = mycol,cluster_col=F, cluster_rows = F, gaps_row=c(5,7,8,9),
         annotation_col = annCol1, annotation_colors=annColors1 , show_rownames = T,show_colnames = F, 
         annotation_legend = T, filename = file.path(fig.path,"heatmap_LUAD_gsva_some_Signatures.pdf"),
         border_color = NA,fontsize_row = 15, fontsize_col = 20,cellwidth = 1.2)

library(stringr)
TPM.HUGO <- read.table(file.path(res.path,"TCGA_LUAD_mRNA_TPM.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
met_sig <- MSigDB
met_sig$term <- str_replace_all(string = met_sig$term,pattern = "HALLMARK_",replacement = "")
met_sig$term <- str_replace_all(string = met_sig$term,pattern = "_",replacement = " ")
label <- unique(met_sig$term)
meta_sig <- list()
for (i in label) {
  meta_sig[[i]] <- met_sig[which(met_sig$term == i),"gene"]  
}
meta.score <- gsva(as.matrix(TPM.HUGO),meta_sig,method="gsva")
write.table(as.data.frame(meta.score),file.path(res.path,"gsva results for hallmark by gsva.txt"),row.names = T,col.names = NA,sep = "\t")

# heatmap
Hallmark_Signature <- read.table(file.path(res.path,"gsva results for hallmark by gsva.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
Hallmark_Signature <- Hallmark_Signature[,rownames(annCol)]
choose_matrix = standarize.fun(indata=Hallmark_Signature,halfwidth = 1.5)
mycol <- colorpanel(256,low="#07689F",mid = "black",high="#ffd424") 

pheatmap(choose_matrix,color = mycol,cluster_col=F, cluster_rows = T,
         annotation_col = annCol1, annotation_colors=annColors1 , show_rownames = T,show_colnames = F, 
         annotation_legend = T, filename = file.path(fig.path,"heatmap_LUAD_gsva_Hallmark_Signatures.pdf"),
         border_color = NA,fontsize_row = 9, fontsize_col = 5,cellwidth = 1.1,cellheight = 11)

#################### Immune genesets ####################
# MUST locate ABSOLUTE path of gene set file
GSET.FILE <- 
  system.file("extdata", "gene sets of interest.gmt", package = "MOVICS", mustWork = TRUE)
# run GSVA to estimate single sample enrichment score based on given gene set of interest
gsva.res <- 
  runGSVA(moic.res      = cmoic.luad,
          norm.expr     = tpm,
          gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          fig.path      = fig.path,
          fig.name      = "GENE SETS OF INTEREST HEATMAP",
          height        = 5,
          width         = 8)

annCol1 <- data.frame(Subtype = annCol$Subtype)
rownames(annCol1) <- rownames(annCol)
# generate corresponding colors for sample annotation
annColors1 <- list(Subtype = c("CS1" = "#2EC4B6", "CS2" = "#E71D36", 
                               "CS3" = "#FF9F1C", "CS4" = "#BDD5EA"))

library( "pheatmap" )
library(GSVA)
library("gplots")

TPM.HUGO <- read.table(file.path(res.path,"TCGA_LUAD_mRNA_TPM.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
met_sig <- read.table(file.path(comAnn.path,"CCR_Curated_Immune_Cell_Signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL,quote = "")
label <- unique(met_sig$`CellType`)
meta_sig <- list()
for (i in label) {
  meta_sig[[i]] <- met_sig[which(met_sig$`CellType` == i),"Symbol"]  
}
meta.score <- gsva(as.matrix(TPM.HUGO),meta_sig,method="gsva")
write.table(as.data.frame(meta.score),file.path(res.path,"gsva results for CCR_Curated_Immune_Cell_Signature by gsva.txt"),row.names = T,col.names = NA,sep = "\t")
# heatmap
Immune_Cell_Signature <- read.table(file.path(res.path,"gsva results for CCR_Curated_Immune_Cell_Signature by gsva.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
Immune_Cell_Signature <- Immune_Cell_Signature[,rownames(annCol)]
choose_matrix = standarize.fun(indata=Immune_Cell_Signature,halfwidth = 1.5)
mycol <- colorpanel(256,low=cyan,mid = "black",high=peach) 

pheatmap(choose_matrix,color = mycol,cluster_col=F, cluster_rows = F,
         annotation_col = annCol1, annotation_colors=annColors1 , show_rownames = T,show_colnames = F, 
         annotation_legend = T, filename = file.path(fig.path,"heatmap_LUAD_gsva_CCR_Curated_Immune_Cell_Signature.pdf"),
         border_color = NA,fontsize_row = 9, fontsize_col = 5,cellwidth = 1.1,cellheight = 12)
#dev.off()

############## Oncogenetic signature
TPM.HUGO <- read.table(file.path(res.path,"TCGA_LUAD_mRNA_TPM.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
met_sig <- read.table(file.path(comAnn.path,"Oncogenetic_signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL,quote = "")
label <- unique(met_sig$Pathway)
meta_sig <- list()
for (i in label) {
  meta_sig[[i]] <- met_sig[which(met_sig$Pathway == i),"Symbol"]  
}
meta.score <- gsva(as.matrix(TPM.HUGO),meta_sig,method="gsva")
write.table(as.data.frame(meta.score),file.path(res.path,"gsva results for Oncogenetic_signature by gsva.txt"),row.names = T,col.names = NA,sep = "\t")
# heatmap
Oncogenetic_signature <- read.table(file.path(res.path,"gsva results for Oncogenetic_signature by gsva.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
Oncogenetic_signature <- Oncogenetic_signature[,rownames(annCol)]
choose_matrix = standarize.fun(indata=Oncogenetic_signature,halfwidth = 1.5)
mycol <- colorpanel(256,low=cyan,mid = "black",high=peach) 

pheatmap(choose_matrix,color = mycol,cluster_col=F, cluster_rows = F,
         annotation_col = annCol1, annotation_colors=annColors1 , show_rownames = T,show_colnames = F, 
         annotation_legend = T, filename = file.path(fig.path,"heatmap_LUAD_gsva_Oncogenetic_signature.pdf"),
         border_color = NA,fontsize_row = 9, fontsize_col = 5,cellwidth = 1.1)

########################## C2 VS C3 Immune #############################
### CIBERSORT ###
CS1 <- rownames(annCol[which(annCol$Subtype=="CS1"),])
CS2 <- rownames(annCol[which(annCol$Subtype=="CS2"),])
CS3 <- rownames(annCol[which(annCol$Subtype=="CS3"),])
CS4 <- rownames(annCol[which(annCol$Subtype=="CS4"),])

TPM <- read.table(file.path(data.path,"TCGA_LUAD_mRNA_TPM.txt"),sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
tmp <- log2(TPM + 1)
tmp <- sweep(tmp,1,apply(tmp, 1, median))
write.table(tmp,file.path(res.path,"CIBERSOT.TPM.logTransMedCent.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

tmp <- read.table(file.path(res.path,"CIBERSORT.Output_Job2.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
tmp <- as.data.frame(t(tmp[,1:22]))
cibersort <- tmp
p <- c()
for (i in 1:22) {
  tmp1 <- as.numeric(tmp[i,CS2])
  tmp2 <- as.numeric(tmp[i,CS3])
  p <- c(p,wilcox.test(tmp1,tmp2)$p.value)
}
names(p) <- rownames(tmp)
cibersort.p<- as.data.frame(p)

### TIDE ###
tmp <- log2(TPM + 1)
tmp <- sweep(tmp,2,apply(tmp, 2, median))
tmp <- sweep(tmp,1,apply(tmp, 1, median))
write.table(tmp,file.path(res.path,"TIDE.TPM.logTransMedCentBoth.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

TIDE <- read.csv(file.path(res.path,"LUAD_TIDE_output.csv"),check.names = F, header = T,stringsAsFactors = F,row.names = 1)
TIDE <- TIDE[rownames(annCol),]
TIDE$Subtype <- rep(c("CS1","CS2","CS3","CS4"),
                    c(length(CS1),length(CS2),length(CS3),length(CS4)))
chisq.test(table(TIDE$Responder,TIDE$Subtype))  #X-squared = 66.492, df = 3, p-value = 2.405e-14
table(TIDE$Responder,TIDE$Subtype)

TIDE_result <- data.frame(table(TIDE$Responder,TIDE$Subtype)) 
TIDE_result<- TIDE_result %>% 
  group_by(Var2) %>% 
  mutate(sumVal = sum(Freq)) %>%
  ungroup() %>%
  mutate(per=Freq/sumVal) %>%
  mutate(label=paste0(round(per,2),"%"))%>%
  mutate(Var1=ifelse(TIDE_result$Var1=="False","Non-Response","Response"))
ggplot(aes(x=Var2, y=per, fill=Var1), data=TIDE_result)+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#119da4","#ffc857"))+
  xlab("")+ 
  ylab("Percentage")+
  geom_text(aes(label = label), 
            color="black", size=5,position=position_fill(0.5))+
  guides(fill=guide_legend(title=NULL))+
  theme_bw() +
  theme(axis.title.y=element_text(size=18),
        axis.text=element_text(size=13),
        legend.text=element_text(size=13))

TIDE23 <- TIDE[which(annCol$Subtype=="CS2"|annCol$Subtype=="CS3"),]
TIDE23$Subtype <- rep(c("CS2","CS3"),c(length(CS2),length(CS3)))
chisq.test(table(TIDE23$Responder,TIDE23$Subtype))  #p-value = 1.357e-13


### MCPCOUNTER and ESITIMATE ###
tmp <- log2(TPM + 1)

# library(utils)
# rforge <- "http://r-forge.r-project.org"
# install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
filterCommonGenes(input.f=file.path(res.path, "TCGA_LUAD_mRNA_TPM.txt") , output.f=file.path(res.path,"LUAD_HUGO_symbol_ESTIMATE.txt"), id="GeneSymbol")
estimateScore(file.path(res.path,"LUAD_HUGO_symbol_ESTIMATE.txt"), file.path(res.path,"LUAD_HUGO_symbol_estimate_score.txt"), platform="affymetrix")
est <- read.table(file = file.path(res.path,"LUAD_HUGO_symbol_estimate_score.txt"),header = T,row.names = NULL,check.names = F,stringsAsFactors = F,sep = "\t")
rownames(est) <- est[,2]; colnames(est) <- est[1,]; est <- est[-1,c(-1,-2)];
est <- sapply(est, as.numeric); rownames(est) <- c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity")
est.raw <- est; colnames(est.raw) <- colnames(TPM)
source(file.path(script.path,"annTrackScale.R"))
tmp <- annTrackScale(indata = est, halfwidth = 2, poolsd = F); tmp <- as.data.frame(t(tmp))
rownames(tmp) <- colnames(TPM)
annCol <- cbind.data.frame(annCol,tmp[rownames(annCol),])

library(MCPcounter)
#MCPscore <- MCPcounter.estimate(expression = tmp,featuresType = "HUGO_symbols")
#write.table(MCPscore,"LUAD_TPM_MCPscore_HUGO_symbols.txt",row.names=T, col.names=NA, sep="\t", quote=F)
MCPscore <- read.table(file.path(res.path,"LUAD_TPM_MCPscore_HUGO_symbols.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)

tmp <- as.data.frame(MCPscore[,rownames(annCol)])
MCPscore.raw <- tmp
tmp <- annTrackScale(indata = tmp, halfwidth = 2, poolsd = F)
tmp <- as.data.frame(t(tmp))
annCol <- cbind.data.frame(annCol,tmp)

# CS3 <- rownames(annCol[which(annCol$Subtype=="CS3"),])
# CS4 <- rownames(annCol[which(annCol$Subtype=="CS4"),])

p.mcp <- c()
for (i in 1:10) {
  tmp1 <- as.numeric(MCPscore.raw[i,CS2])
  tmp2 <- as.numeric(MCPscore.raw[i,CS3])
  p.mcp <- c(p.mcp,wilcox.test(tmp1,tmp2)$p.value)
}
names(p.mcp) <- rownames(MCPscore.raw)

p.est <- c()
for (i in 1:4) {
  tmp1 <- as.numeric(est.raw[i,CS2])
  tmp2 <- as.numeric(est.raw[i,CS3])
  p.est <- c(p.est,wilcox.test(tmp1,tmp2)$p.value)
}
names(p.est) <- rownames(est.raw)

# boxplot for scores
require(tidyr)
annCol_new <- annCol[which(annCol$Subtype=="CS2"|annCol$Subtype=="CS3"),]

dd <- as.data.frame(est.raw[,c(CS2,CS3)])
dd$cell = rownames(dd)
d2 <- gather(dd, sample, ES, 1:264)
tmp <- data.frame(sample=rownames(annCol_new),Subtype=rep(c("CS2","CS3"),c(length(CS2),length(CS3))),stringsAsFactors = F)
d2 <- merge(d2,tmp,by="sample",all.x=T)
pvalues <- sapply(d2$cell, function(x) {
  res <- wilcox.test(ES ~ Subtype, data = subset(d2, cell == x))$p.value
})
pv <- data.frame(cell = d2$cell, pvalue = pvalues)
d2$cell <- factor(d2$cell,levels = rownames(est.raw))
d2$Subtype <- factor(d2$Subtype,levels = c("CS2","CS3"))
pv$sigcode <- cut(pv$pvalue, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  labels=c('***', '**', '*', '.', ' '))
p1 <- ggplot(d2, aes(cell, ES, fill=Subtype)) + 
  geom_boxplot() + scale_fill_manual(values = c("#E71D36","#FF9F1C")) + 
  geom_text(aes(cell, y=max(d2$ES) * 1.1, 
                label=pv$sigcode),
            data=pv, inherit.aes=F) + 
  xlab(NULL)+ylab("ESITIMATE Score") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 8),
        axis.title = element_blank()) 
ggsave(file.path(fig.path,"Boxplot for ESITIMATE score of CS2 and CS3.pdf"),width = 8,height = 3)

dd <- as.data.frame(MCPscore.raw[,c(CS2,CS3)])
dd$cell = rownames(dd)
d2 <- gather(dd, sample, ES, 1:264)
tmp <- data.frame(sample=rownames(annCol_new),Subtype=rep(c("CS2","CS3"),c(length(CS2),length(CS3))),stringsAsFactors = F)
d2 <- merge(d2,tmp,by="sample",all.x=T)
pvalues <- sapply(d2$cell, function(x) {
  res <- wilcox.test(ES ~ Subtype, data = subset(d2, cell == x))$p.value
})
pv <- data.frame(cell = d2$cell, pvalue = pvalues)
d2$cell <- factor(d2$cell,levels = rownames(MCPscore.raw))
d2$Subtype <- factor(d2$Subtype,levels = c("CS2","CS3"))
## Signif. codes:  0 ??***?? 0.001 ??**?? 0.01 ??*?? 0.05 ??.?? 0.1 ?? ?? 1
pv$sigcode <- cut(pv$pvalue, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  labels=c('***', '**', '*', '.', ' '))
p1 <- ggplot(d2, aes(cell, ES, fill=Subtype)) + 
  geom_boxplot() + scale_fill_manual(values = c("#E71D36","#FF9F1C")) + 
  geom_text(aes(cell, y=max(d2$ES) * 1.1, 
                label=pv$sigcode),
            data=pv, inherit.aes=F) + 
  xlab(NULL)+ylab("MCPcounter Score") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 8),
        axis.title = element_blank())
ggsave(file.path(fig.path,"Boxplot for MCPcounter score of CS2 and CS3.pdf"),width = 8,height = 3)

dd <- as.data.frame(cibersort[,c(CS2,CS3)])
dd$cell = rownames(dd)
d2 <- gather(dd, sample, ES, 1:264)
tmp <- data.frame(sample=rownames(annCol_new),Subtype=rep(c("CS2","CS3"),
                                                          c(length(CS2),length(CS3))),stringsAsFactors = F)
d2 <- merge(d2,tmp,by="sample",all.x=T)
pvalues <- sapply(d2$cell, function(x) {
  res <- wilcox.test(ES ~ Subtype, data = subset(d2, cell == x))$p.value
})
pv <- data.frame(cell = d2$cell, pvalue = pvalues)
d2$cell <- factor(d2$cell,levels = rownames(cibersort))
d2$Subtype <- factor(d2$Subtype,levels = c("CS2","CS3"))
## Signif. codes:  0 ??***?? 0.001 ??**?? 0.01 ??*?? 0.05 ??.?? 0.1 ?? ?? 1
pv$sigcode <- cut(pv$pvalue, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  labels=c('***', '**', '*', '.', ' '))
library(ggplot2)
p2 <- ggplot(d2, aes(cell, ES, fill=Subtype)) + 
  geom_boxplot() + scale_fill_manual(values = c("#E71D36","#FF9F1C")) + 
  geom_text(aes(cell, y=max(d2$ES) * 1.1, 
                label=pv$sigcode),
            data=pv, inherit.aes=F) + 
  xlab(NULL)+ylab("CIBERSORT Score") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 8),axis.title = element_blank())
ggsave(file.path(fig.path,"Boxplot for CIBERSORT of CS2 and CS3.pdf"),width = 8,height = 3)


################ Immune_Cell_Signature ################
# boxplot
Immune <- read.table(file = file.path(res.path,"gsva results for Immune_Cell_Signature by gsva.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
#Immune <- Immune[which(rownames(Immune) %in% c("CD8 T cells","T cells","T helper cells","Th2 cells","Th17 cells")),]
#Immune <- sapply(Immune, as.numeric); rownames(Immune) <- c("T cells","T helper cells","Th2 cells","Th17 cells","CD8 T cells")
Immune.raw <- Immune
#tmp <- annTrackScale(indata = Immune, halfwidth = 2, poolsd = F)
tmp <- as.data.frame(t(Immune[,rownames(annCol)]))
tmp$Cluster <- annCol$Subtype
p.Immune <- c()
for (i in 1:28) {
  tmp1 <- as.numeric(Immune.raw[i,CS2])
  tmp2 <- as.numeric(Immune.raw[i,CS3])
  p.Immune <- c(p.Immune,wilcox.test(tmp1,tmp2)$p.value)
}
names(p.Immune) <- rownames(Immune.raw)

dd <- as.data.frame(Immune.raw[,c(CS2,CS3)])
dd$cell = rownames(dd)
d2 <- gather(dd, sample, ES, 1:236)
tmp <- data.frame(sample=rownames(annCol_new),Subtype=rep(c("CS2","CS3"),c(length(CS2),length(CS3))),stringsAsFactors = F)
d2 <- merge(d2,tmp,by="sample",all.x=T)
pvalues <- sapply(d2$cell, function(x) {
  res <- wilcox.test(ES ~ Subtype, data = subset(d2, cell == x))$p.value
})
pv <- data.frame(cell = d2$cell, pvalue = pvalues)
d2$cell <- factor(d2$cell,levels = rownames(Immune.raw))
d2$Subtype <- factor(d2$Subtype,levels = c("CS2","CS3"))
pv$sigcode <- cut(pv$pvalue, c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  labels=c('***', '**', '*', '.', ' '))
p1 <- ggplot(d2, aes(cell, ES, fill=Subtype)) + 
  geom_boxplot() + scale_fill_manual(values = c("#E71D36","#FF9F1C")) + 
  geom_text(aes(cell, y=max(d2$ES) * 1.1, 
                label=pv$sigcode),
            data=pv, inherit.aes=F) + 
  xlab(NULL)+ylab("GSVA Immune Cell Signature") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 8),
        axis.title = element_blank())
ggsave(file.path(fig.path,"Boxplot for gsva_Immune_Cell_Signature of CS2 and CS3.pdf"),width = 10,height = 3)

########################## GSVA #############################
### EMT(NTP) ###
library(CMScaller)
EMT_signature <- read.table(file.path(comAnn.path,"EMT_signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
emat <- read.table(file.path(res.path,"TCGA_LUAD_mRNA_TPM.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
ntp.matrix <- ntp(emat=emat, templates=EMT_signature, nPerm=1000,distance = "pearson", doPlot =  FALSE)
ntp <- ntp.matrix[rownames(annCol),]
ntp$Cluster <- annCol$Subtype
#write.table(ntp,file.path(res.path,"LUAD_EMT_result.txt"),row.names = T,col.names = NA,sep = "\t")
table(ntp$prediction,ntp$Cluster)
fisher.test(table(ntp$prediction,ntp$Cluster))   #p-value = 0.2137
ntp23 <- ntp[which(ntp$Cluster=="CS2"|ntp$Cluster=="CS3"),]
fisher.test(table(ntp23$prediction,ntp23$Cluster))   #p-value =0.2642
ntp14 <- ntp[which(ntp$Cluster=="CS1"|ntp$Cluster=="CS4"),]
fisher.test(table(ntp14$prediction,ntp14$Cluster))   #p-value = 0.08861

# ssGSEA--Metabolism_signature
library(GSVA)
TPM.HUGO <- read.table(file.path(res.path,"TCGA_LUAD_mRNA_TPM.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
met_sig <- read.table(file.path(comAnn.path,"Metabolism_signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL,quote = "")
label <- unique(met_sig$`Metabolic Function`)
meta_sig <- list()
for (i in label) {
  meta_sig[[i]] <- met_sig[which(met_sig$`Metabolic Function` == i),"Gene Name"]  
}
meta.score <- gsva(as.matrix(TPM.HUGO),meta_sig,method="gsva")
write.table(as.data.frame(meta.score),file.path(res.path,"gsva results for Metabolism_signature by gsva.txt"),row.names = T,col.names = NA,sep = "\t")

# ssGSEA--Metabolism_signature2
met_sig2 <- read.table(file.path(comAnn.path,"Metabolism_signature2.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL,quote = "")
label <- unique(met_sig2$`Metabolic Function`)
meta_sig2 <- list()
for (i in label) {
  meta_sig2[[i]] <- met_sig2[which(met_sig2$`Metabolic Function` == i),"Gene Name"]  
}
meta.score2 <- gsva(as.matrix(TPM.HUGO),meta_sig2,method="gsva")
write.table(as.data.frame(meta.score2),file.path(res.path,"gsva results for Metabolism_signature2(NEW) by gsva.txt"),row.names = T,col.names = NA,sep = "\t")

# heatmap--Metabolism_signature
Metabolism_signature1 <-read.table(file.path(res.path,"gsva results for Metabolism_signature by gsva.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
Metabolism_signature1 <- Metabolism_signature1[,rownames(annCol1)]
Metabolism_signature2 <-read.table(file.path(res.path,"gsva results for Metabolism_signature2(NEW) by gsva.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
Metabolism_signature2 <- Metabolism_signature2[,rownames(annCol1)]
Metabolism_signature <- rbind(Metabolism_signature1,Metabolism_signature2)

library( "pheatmap" )
library("gplots")
choose_matrix = standarize.fun(indata=Metabolism_signature,halfwidth = 1.5)
#mycol <- colorpanel(256,low=cyan,mid = "black",high=peach) 
mycol <- colorpanel(256,low="#07689F",mid = "black",high="#ffd424") 
#mycol <- colorpanel(256,low=blue,mid = "black",high=gold)

pheatmap(choose_matrix,color = mycol,cluster_col=F, cluster_rows = F,
         annotation_col = annCol1, annotation_colors=annColors1 , show_rownames = T,show_colnames = F, 
         annotation_legend = T, filename = file.path(fig.path,"heatmap_LUAD_gsva_Metabolism_signature(NEW).pdf"),
         border_color = NA,fontsize_row = 9, fontsize_col = 5,cellwidth = 1.1,cellheight = 12)

# ssGSEA--Metabolism_signature（7sig）
library(GSVA)
TPM.HUGO <- read.table(file.path(res.path,"TCGA_LUAD_mRNA_TPM.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
met_sig <- read.table(file.path(comAnn.path,"Metabolism_signature7.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL,quote = "")
label <- unique(met_sig$`Metabolic Function`)
meta_sig <- list()
for (i in label) {
  meta_sig[[i]] <- met_sig[which(met_sig$`Metabolic Function` == i),"Gene Name"]  
}
meta.score <- gsva(as.matrix(TPM.HUGO),meta_sig,method="gsva")
write.table(as.data.frame(meta.score),file.path(res.path,"gsva results for Metabolism_signature7 by gsva.txt"),row.names = T,col.names = NA,sep = "\t")

# heatmap--Metabolism_signature7
Metabolism_signature7 <-read.table(file.path(res.path,"gsva results for Metabolism_signature7 by gsva.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
Metabolism_signature7 <- Metabolism_signature7[,rownames(annCol1)]

library( "pheatmap" )
library("gplots")
choose_matrix = standarize.fun(indata=Metabolism_signature7,halfwidth = 1.5)
#mycol <- colorpanel(256,low=cyan,mid = "black",high=peach) 
mycol <- colorpanel(256,low="#07689F",mid = "black",high="#ffd424") 
#mycol <- colorpanel(256,low=blue,mid = "black",high=gold)

pheatmap(choose_matrix,color = mycol,cluster_col=F, cluster_rows = F,
         annotation_col = annCol1, annotation_colors=annColors1 , show_rownames = T,show_colnames = F, 
         annotation_legend = T, filename = file.path(fig.path,"heatmap_LUAD_gsva_Metabolism_signature7.pdf"),
         border_color = NA,fontsize_row = 9, fontsize_col = 5,cellwidth = 1.1,cellheight = 12)

# ssGSEA--Immune_Cell_Signature
met_sig <- read.table(file.path(comAnn.path,"Immune_Cell_Signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL,quote = "")
label <- unique(met_sig$`CellType`)
meta_sig <- list()
for (i in label) {
  meta_sig[[i]] <- met_sig[which(met_sig$`CellType` == i),"Symbol"]  
}
meta.score <- gsva(as.matrix(TPM.HUGO),meta_sig,method="gsva")
write.table(as.data.frame(meta.score),file.path(res.path,"gsva results for Immune_Cell_Signature by gsva.txt"),row.names = T,col.names = NA,sep = "\t")

# heatmap
Immune_Cell_Signature <- read.table(file.path(res.path,"gsva results for Immune_Cell_Signature by gsva.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
Immune_Cell_Signature <- Immune_Cell_Signature[,rownames(annCol1)]
choose_matrix = standarize.fun(indata=Immune_Cell_Signature,halfwidth = 1.5)
mycol <- colorpanel(256,low="#07689F",mid = "black",high="#ffd424") 

pheatmap(choose_matrix,color = mycol,cluster_col=F, cluster_rows = F,
         annotation_col = annCol1, annotation_colors=annColors1 , show_rownames = T,show_colnames = F, 
         annotation_legend = T, filename = file.path(fig.path,"heatmap_LUAD_gsva_Immune_Cell_Signature.pdf"),
         border_color = NA,fontsize_row = 9, fontsize_col = 5,cellwidth = 1.1,cellheight = 14)

################### Ferroptosis ################### 
library(GSVA)
TPM.HUGO <- read.table(file.path(res.path,"TCGA_LUAD_mRNA_TPM.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
met_sig <- read.table(file.path(comAnn.path,"Ferroptosis.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL,quote = "")
label <- unique(met_sig$`Metabolic Function`)
meta_sig <- list()
for (i in label) {
  meta_sig[[i]] <- met_sig[which(met_sig$`Metabolic Function` == i),"Gene Name"]  
}
meta.score <- gsva(as.matrix(TPM.HUGO),meta_sig,method="gsva")
write.table(as.data.frame(meta.score),file.path(res.path,"gsva results for Ferroptosis by gsva.txt"),row.names = T,col.names = NA,sep = "\t")

# heatmap--Ferroptosis
Ferroptosis_signature <-read.table(file.path(res.path,"gsva results for Ferroptosis by gsva.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
Ferroptosis_signature <- Ferroptosis_signature[,rownames(annCol1)]

library( "pheatmap" )
library("gplots")
choose_matrix = standarize.fun(indata=Ferroptosis_signature,halfwidth = 1.5)
#mycol <- colorpanel(256,low=cyan,mid = "black",high=peach) 
mycol <- colorpanel(256,low="#07689F",mid = "black",high="#ffd424") 
#mycol <- colorpanel(256,low=blue,mid = "black",high=gold)

pheatmap(choose_matrix,color = mycol,cluster_col=F, cluster_rows = F,
         annotation_col = annCol1, annotation_colors=annColors1 , show_rownames = T,show_colnames = F, 
         annotation_legend = T, filename = file.path(fig.path,"heatmap_LUAD_gsva_Ferroptosis_signature.pdf"),
         border_color = NA,fontsize_row = 9, fontsize_col = 5,cellheight = 15)

#############################
### Copy Number Variation ###
# compare copy number variation
seg <- read.table(file.path(data.path,"LUAD.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt"),sep = "\t",row.names = NULL,header = T,check.names = F,stringsAsFactors = F)
seg <- seg[grepl("-01A-",seg$Sample),]
seg$Sample <- substr(seg$Sample,1,15)

CS1.cnv <- intersect(seg$Sample,rownames(annCol[which(annCol$Subtype == "CS1"),]))
CS2.cnv <- intersect(seg$Sample,rownames(annCol[which(annCol$Subtype == "CS2"),]))
CS3.cnv <- intersect(seg$Sample,rownames(annCol[which(annCol$Subtype == "CS3"),]))
CS4.cnv <- intersect(seg$Sample,rownames(annCol[which(annCol$Subtype == "CS4"),]))

seg.CS1 <- seg[which(seg$Sample %in% CS1.cnv),]
seg.CS2 <- seg[which(seg$Sample %in% CS2.cnv),]
seg.CS3 <- seg[which(seg$Sample %in% CS3.cnv),]
seg.CS4 <- seg[which(seg$Sample %in% CS4.cnv),]

seg <- seg[which(seg$Sample %in% c(CS1.cnv,CS2.cnv,CS3.cnv,CS4.cnv)),]
seg$Sample <- as.character(seg$Sample)
rownames(seg) <- 1:nrow(seg)
write.table(as.data.frame(seg),file.path(res.path,"LUAD_all_segment_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)
write.table(as.data.frame(seg.CS1),file.path(res.path,"LUAD_CS1_segment_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)
write.table(as.data.frame(seg.CS2),file.path(res.path,"LUAD_CS2_segment_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)
write.table(as.data.frame(seg.CS3),file.path(res.path,"LUAD_CS3_segment_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)
write.table(as.data.frame(seg.CS4),file.path(res.path,"LUAD_CS4_segment_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)

marker <- seg[,1:4]
a <- b <- c <- c()
for(i in 1:nrow(marker)) {
  a <- c(a,rep(marker[i,"Sample"],2))
  b <- c(b,rep(marker[i,"Chromosome"],2))
  c <- c(c,marker[i,"Start"],marker[i,"End"])
}

marker <- data.frame(Marker.Name = a,Chromosome = b,Marker.Position=c,stringsAsFactors = F)
write.table(as.data.frame(marker),file.path(res.path,"luad_all_marker_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)
rm(a); rm(b); rm(c)

marker <- seg.CS1[,1:4]
a <- b <- c <- c()
for(i in 1:nrow(marker)) {
  a <- c(a,rep(marker[i,"Sample"],2))
  b <- c(b,rep(marker[i,"Chromosome"],2))
  c <- c(c,marker[i,"Start"],marker[i,"End"])
}

marker <- data.frame(Marker.Name = a,Chromosome = b,Marker.Position=c,stringsAsFactors = F)
write.table(as.data.frame(marker),file.path(res.path,"luad_CS1_marker_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)
rm(a); rm(b); rm(c)

marker <- seg.CS2[,1:4]
a <- b <- c <- c()
for(i in 1:nrow(marker)) {
  a <- c(a,rep(marker[i,"Sample"],2))
  b <- c(b,rep(marker[i,"Chromosome"],2))
  c <- c(c,marker[i,"Start"],marker[i,"End"])
}

marker <- data.frame(Marker.Name = a,Chromosome = b,Marker.Position=c,stringsAsFactors = F)
write.table(as.data.frame(marker),file.path(res.path,"luad_CS2_marker_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)
rm(a); rm(b); rm(c)

marker <- seg.CS3[,1:4]
a <- b <- c <- c()
for(i in 1:nrow(marker)) {
  a <- c(a,rep(marker[i,"Sample"],2))
  b <- c(b,rep(marker[i,"Chromosome"],2))
  c <- c(c,marker[i,"Start"],marker[i,"End"])
}

marker <- data.frame(Marker.Name = a,Chromosome = b,Marker.Position=c,stringsAsFactors = F)
write.table(as.data.frame(marker),file.path(res.path,"luad_CS3_marker_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)
rm(a); rm(b); rm(c)

marker <- seg.CS4[,1:4]
a <- b <- c <- c()
for(i in 1:nrow(marker)) {
  a <- c(a,rep(marker[i,"Sample"],2))
  b <- c(b,rep(marker[i,"Chromosome"],2))
  c <- c(c,marker[i,"Start"],marker[i,"End"])
}

marker <- data.frame(Marker.Name = a,Chromosome = b,Marker.Position=c,stringsAsFactors = F)
write.table(as.data.frame(marker),file.path(res.path,"luad_CS4_marker_forGISTIC2.0.txt"),sep = "\t",row.names = F,quote = F)
rm(a); rm(b); rm(c)

# generate gistic plot
# Create a chromosomes reference objects function
chrom_extract <- function(BSgenome.hg  = NULL) {
  if (is.null(BSgenome.hg )) stop("NULL object !", call. = FALSE)
  obj <- list(species = GenomeInfoDb::organism(BSgenome.hg), genomebuild = BSgenome::providerVersion(BSgenome.hg))
  df <- data.frame(chrom = BSgenome::seqnames(BSgenome.hg), chrN = seq_along(BSgenome::seqnames(BSgenome.hg)), chr.length = GenomeInfoDb::seqlengths(BSgenome.hg), stringsAsFactors = FALSE)
  df <- df[1:24,]
  df$chr.length.sum <- cumsum(as.numeric(df$chr.length))
  df$chr.length.cumsum <- c(0, df$chr.length.sum[-nrow(df)])
  df$middle.chr <- round(diff(c(0, df$chr.length.sum)) /2)
  df$middle.chr.genome <- df$middle.chr + df$chr.length.cumsum
  obj$chromosomes <- df
  obj$chrom2chr <- sapply(obj$chromosomes$chrom, function(k) { obj$chromosomes$chrN[obj$chromosomes$chrom == k]}, simplify = FALSE)
  obj$chr2chrom <- sapply(obj$chromosomes$chrN, function(k) { obj$chromosomes$chrom[obj$chromosomes$chrN == k]}, simplify = FALSE)
  names(obj$chr2chrom) <- obj$chromosomes$chrN
  obj$genome.length <- sum(as.numeric(obj$chromosomes$chr.length), na.rm = TRUE)
  return(obj)
}

# Extract a chromosomes reference loci
BSgenome.hg = "BSgenome.Hsapiens.UCSC.hg19"
BSg.obj <- getExportedValue(BSgenome.hg, BSgenome.hg)
genome.version <- BSgenome::providerVersion(BSg.obj)
chrom <- chrom_extract(BSg.obj)
#str(chrom)
col1 <- adjustcolor("#800f2f", alpha.f = .7)
col2 <- adjustcolor("#013a63", alpha.f = .7)

### CS1 ###
pdf(file.path(fig.path,"CS1.cnv.scores.gistic.pdf"),14,5)

scores <- read.table(file.path(data.path,"CS1_.scores.gistic"), sep="\t",header=T,stringsAsFactors = F)

# Important step for accurate length to match back to continual chrom loci
scores[scores$Chromosome==23,"Chromosome"]="X"
scores[scores$Chromosome==24,"Chromosome"]="Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$G.score <- scores.amp$G.score * 1
scores.del <- scores[scores$Type=="Del",]
scores.del$G.score <- scores.del$G.score * -1
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$G.score)-0.1,max(scores$G.score)+0.1)
title=paste0("LUAD CS1 copy number"," ","n=",length(CS1.cnv))

plot(scores.amp$Start.geno, scores.amp$G.score,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", 
     xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 2, ylab = "gistic score", xlab = NA,
     cex.lab = 1.5, col = adjustcolor("#800f2f", alpha.f = .7),
     xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$G.score, type='h', lwd = 2, 
      col = adjustcolor("#013a63", alpha.f = .7))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+0.08,ylim[2]-0.08)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], 
         labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("#800f2f", alpha.f = .7)
col2 <- adjustcolor("#0466c8", alpha.f = .7)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.7, bty="n", fill=c(col1,col2))
dev.off()

### CS2 ###
pdf(file.path(fig.path,"CS2.cnv.scores.gistic.pdf"),14,5)

scores <- read.table(file.path(data.path,"CS2_.scores.gistic"), sep="\t",header=T,stringsAsFactors = F)

# Important step for accurate length to match back to continual chrom loci
scores[scores$Chromosome==23,"Chromosome"]="X"
scores[scores$Chromosome==24,"Chromosome"]="Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$G.score <- scores.amp$G.score * 1
scores.del <- scores[scores$Type=="Del",]
scores.del$G.score <- scores.del$G.score * -1
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$G.score)-0.1,max(scores$G.score)+0.1)
title=paste0("LUAD CS2 copy number"," ","n=",length(CS2.cnv))

plot(scores.amp$Start.geno, scores.amp$G.score,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", 
     xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 2, ylab = "gistic score", xlab = NA,
     cex.lab = 1.5, 
     col = adjustcolor("#800f2f", alpha.f = .7), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$G.score, type='h', lwd = 2, 
      col = adjustcolor("#013a63", alpha.f = .7))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+0.05,ylim[2]-0.05)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], 
         labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("#800f2f", alpha.f = .7)
col2 <- adjustcolor("#013a63", alpha.f = .7)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.7, bty="n", fill=c(col1,col2))

dev.off()

### CS3 ###
pdf(file.path(fig.path,"CS3.cnv.scores.gistic.pdf"),14,5)

scores <- read.table(file.path(data.path,"CS3_.scores.gistic"), sep="\t",header=T,stringsAsFactors = F)

# Important step for accurate length to match back to continual chrom loci
scores[scores$Chromosome==23,"Chromosome"]="X"
scores[scores$Chromosome==24,"Chromosome"]="Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$G.score <- scores.amp$G.score * 1
scores.del <- scores[scores$Type=="Del",]
scores.del$G.score <- scores.del$G.score * -1
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$G.score)-0.1,max(scores$G.score)+0.1)
title=paste0("LUAD CS3 copy number"," ","n=",length(CS3.cnv))

plot(scores.amp$Start.geno, scores.amp$G.score,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", 
     xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 2, ylab = "gistic score", xlab = NA,
     cex.lab = 1.5, col = adjustcolor("#800f2f", alpha.f = .7), 
     xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$G.score, type='h', lwd = 2, 
      col = adjustcolor("#013a63", alpha.f = .7))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+0.05,ylim[2]-0.05)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], 
         labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("#800f2f", alpha.f = .7)
col2 <- adjustcolor("#0466c8", alpha.f = .7)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.7, bty="n", fill=c(col1,col2))
dev.off()

### CS4 ###
pdf(file.path(fig.path,"CS4.cnv.scores.gistic.pdf"),14,5)

scores <- read.table(file.path(data.path,"CS4_.scores.gistic"), sep="\t",header=T,stringsAsFactors = F)

# Important step for accurate length to match back to continual chrom loci
scores[scores$Chromosome==23,"Chromosome"]="X"
scores[scores$Chromosome==24,"Chromosome"]="Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$G.score <- scores.amp$G.score * 1
scores.del <- scores[scores$Type=="Del",]
scores.del$G.score <- scores.del$G.score * -1
scores <- rbind.data.frame(scores.amp,scores.del)

# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$G.score)-0.1,max(scores$G.score)+0.1)
title=paste0("LUAD CS4 copy number"," ","n=",length(CS4.cnv))

plot(scores.amp$Start.geno, scores.amp$G.score,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", 
     xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 2, ylab = "gistic score", xlab = NA,
     cex.lab = 1.5, col = adjustcolor("#800f2f", alpha.f = .7), 
     xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$G.score, type='h', lwd = 2, 
      col = adjustcolor("#013a63", alpha.f = .7))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+0.1,ylim[2]-0.1)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) + 2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], 
         labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("#800f2f", alpha.f = .7)
col2 <- adjustcolor("#0466c8", alpha.f = .7)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.7, bty="n", fill=c(col1,col2))
dev.off()

########### Mutation ###########
################################
library(tidyverse)
library(magrittr)
library(readxl)
library(stringr)
library(forcats)

maf <- read_tsv(file.path(data.path,"mc3.v0.2.8.PUBLIC.maf"), comment = "#")
label <- c("Tumor_Sample_Barcode","Hugo_Symbol","NCBI_Build","Chromosome","Start_Position","End_Position","Strand",
           "Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2",
           "HGVSc","HGVSp","HGVSp_Short","BIOTYPE")

maf$bcr_patient_barcode <- substr(maf$Tumor_Sample_Barcode,start = 1,stop = 12)

tmp <- read.table(file.path(comAnn.path,"pancancerSurvivalData_XLu.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
tmp <- tmp[which(tmp$type %in% c("LUAD")),"bcr_patient_barcode"]
luad.mut.maf <- maf[which(maf$bcr_patient_barcode %in% tmp),label]
luad.mut.maf$Tumor_Sample_Barcode <- substr(luad.mut.maf$Tumor_Sample_Barcode,start = 1,stop = 15)

mut.maf <- as.data.frame(luad.mut.maf)
mut.maf <- mut.maf[which(mut.maf$Tumor_Sample_Barcode %in% rownames(annCol)),] 
write.table(as.data.frame(mut.maf),file.path(res.path,"LUAD_MAF.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

# TMB with only SNP
mutload.luad <- as.data.frame(mutect.dataframe(mut.maf)); rownames(mutload.luad) <- mutload.luad$Tumor_Sample_Barcode
indel.luad <- as.data.frame(indel.dataframe(mut.maf)); rownames(indel.luad) <- indel.luad$Tumor_Sample_Barcode

dim(mutload.luad)
head(mutload.luad)

CS1.mut <- intersect(rownames(mutload.luad),CS1) #81
CS2.mut <- intersect(rownames(mutload.luad),CS2) #120
CS3.mut <- intersect(rownames(mutload.luad),CS3) #144
CS4.mut <- intersect(rownames(mutload.luad),CS4) #92

# Independent test
tmp <- as.data.frame(mut.maf[which(mut.maf$Tumor_Sample_Barcode %in% c(CS1.mut,CS2.mut,CS3.mut,CS4.mut)),
                             c("Variant_Classification","Tumor_Sample_Barcode","Hugo_Symbol","Variant_Type")])
mut.binary <- matrix(0,nrow = length(unique(tmp$Hugo_Symbol)),ncol = length(c(CS1.mut,CS2.mut,CS3.mut,CS4.mut)),
                     dimnames = list(c(unique(tmp$Hugo_Symbol)),c(CS1.mut,CS2.mut,CS3.mut,CS4.mut)))

for (i in colnames(mut.binary)) {
  tmp1 <- tmp[which(tmp$Tumor_Sample_Barcode == i),]
  tmp1 <- tmp1[which(tmp1$Variant_Type == "SNP"),]
  # if(is.element("Silent",tmp$Variant_Classification)) {
  #   tmp1 <- tmp1[-which(tmp1$Variant_Classification %in% "Silent"),]
  # }
  for (j in tmp1$Hugo_Symbol)
    mut.binary[j,i] <- 1
}
mut.binary <- as.data.frame(mut.binary); #rownames(mut.binary) <- toupper(rownames(mut.binary))
write.table(mut.binary,file.path(res.path,"mut binary SNP.txt"),sep = "\t",row.names = T,col.names = NA)

source(file.path(script.path,"create.anntrack.R")) 
source(file.path(script.path,"createMutSubtype.R")) 
ans <- rep(c("CS1","CS2","CS3","CS4"),c(81,120,144,92))
names(ans) <- c(CS1.mut,CS2.mut,CS3.mut,CS4.mut)
genelist <- rownames(mut.binary[rowSums(mut.binary) > 0.05*ncol(mut.binary),])
binarymut <- as.data.frame(matrix(0, nrow=437, ncol=length(genelist)))
rownames(binarymut) <- names(ans)
colnames(binarymut) <- genelist
for (k in 1:length(genelist)) {
  res <- create.anntrack(samples=names(ans), subtype=createMutSubtype(mut.binary, names(ans), genelist[k]))
  binarymut[, genelist[k]] <- res$subtype
}

out <- matrix(0, nrow=length(genelist), ncol=5)
colnames(out) <- c(levels(factor(ans)), "pvalue")
rownames(out) <- paste(genelist, "Mutated", sep="_")
for (k in 1:length(genelist)) {
  genek <- genelist[k]
  x <- ans
  y <- binarymut[names(x), genek, drop=T]; y <- as.character(y); names(y) <- names(x)
  tmp <- setdiff(names(x), names(y)[y=="Not Available"])
  x <- x[tmp]
  y <- y[tmp]  
  res <- table(y, x)
  if (!all(colnames(res)==colnames(out)[1:(ncol(out)-1)])) {stop(paste("colnames mismatch for ", k, sep=""))}
  out[k, 1:(ncol(out)-1)] <- res["Mutated", ]
  
  out[k, "pvalue"] <- (fisher.test(x, y))$p.value  
}
out <- as.data.frame(out)
#out$FDR <- p.adjust(out$pvalue)
write.table(out, file.path(res.path, "Independence test between 5% cut gene mutation and CS 1234 with SNP.txt"), row.names=T, col.names=NA, sep="\t", quote=F)

# MutSigCV
CS1.maf <- mut.maf[which(mut.maf$Tumor_Sample_Barcode %in% CS1.mut),]
CS1.maf$Tumor_Seq_Allele1 <- CS1.maf$Tumor_Seq_Allele2
all(CS1.maf$Tumor_Seq_Allele1 == CS1.maf$Tumor_Seq_Allele2)
all(CS1.maf$Tumor_Seq_Allele1 == CS1.maf$Reference_Allele)
outTable <- CS1.maf %>%
  select(Hugo_Symbol,Tumor_Sample_Barcode,Variant_Classification,Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2)
write.table(outTable,file.path(res.path,"CS1.maf mutation.full_for_MutSig.txt"),sep = "\t",row.names = F,quote = F)

CS2.maf <- mut.maf[which(mut.maf$Tumor_Sample_Barcode %in% CS2.mut),]
CS2.maf$Tumor_Seq_Allele1 <- CS2.maf$Tumor_Seq_Allele2
all(CS2.maf$Tumor_Seq_Allele1 == CS2.maf$Tumor_Seq_Allele2)
all(CS2.maf$Tumor_Seq_Allele1 == CS2.maf$Reference_Allele)
outTable <- CS2.maf %>%
  select(Hugo_Symbol,Tumor_Sample_Barcode,Variant_Classification,Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2)
write.table(outTable,file.path(res.path,"CS2.maf mutation.full_for_MutSig.txt"),sep = "\t",row.names = F,quote = F)

CS3.maf <- mut.maf[which(mut.maf$Tumor_Sample_Barcode %in% CS3.mut),]
CS3.maf$Tumor_Seq_Allele1 <- CS3.maf$Tumor_Seq_Allele2
all(CS3.maf$Tumor_Seq_Allele1 == CS3.maf$Tumor_Seq_Allele2)
all(CS3.maf$Tumor_Seq_Allele1 == CS3.maf$Reference_Allele)
outTable <- CS3.maf %>%
  select(Hugo_Symbol,Tumor_Sample_Barcode,Variant_Classification,Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2)
write.table(outTable,file.path(res.path,"CS3.maf mutation.full_for_MutSig.txt"),sep = "\t",row.names = F,quote = F)

CS4.maf <- mut.maf[which(mut.maf$Tumor_Sample_Barcode %in% CS4.mut),]
CS4.maf$Tumor_Seq_Allele1 <- CS4.maf$Tumor_Seq_Allele2
all(CS4.maf$Tumor_Seq_Allele1 == CS4.maf$Tumor_Seq_Allele2)
all(CS4.maf$Tumor_Seq_Allele1 == CS4.maf$Reference_Allele)
outTable <- CS4.maf %>%
  select(Hugo_Symbol,Tumor_Sample_Barcode,Variant_Classification,Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2)
write.table(outTable,file.path(res.path,"CS4.maf mutation.full_for_MutSig.txt"),sep = "\t",row.names = F,quote = F)

CS.maf <- mut.maf[which(mut.maf$Tumor_Sample_Barcode %in% c(CS1.mut,CS2.mut,CS3.mut,CS4.mut)),]
CS.maf$Tumor_Seq_Allele1 <- CS.maf$Tumor_Seq_Allele2
all(CS.maf$Tumor_Seq_Allele1 == CS.maf$Tumor_Seq_Allele2)
all(CS.maf$Tumor_Seq_Allele1 == CS.maf$Reference_Allele)
outTable <- CS.maf %>%
  select(Hugo_Symbol,Tumor_Sample_Barcode,Variant_Classification,Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2)
write.table(outTable,file.path(res.path,"CS1234 mutation.full_for_MutSig.txt"),sep = "\t",row.names = F,quote = F)

# mutation signature
library(deconstructSigs)
maf <- as.data.frame(mut.maf[which(mut.maf$Variant_Classification %in% c("Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Silent") & mut.maf$Variant_Type == "SNP"),])
maf$Chromosome <- paste0("chr",maf$Chromosome)
sigs.input <- mut.to.sigs.input(mut.ref = maf, 
                                sample.id = "Tumor_Sample_Barcode", 
                                chr = "Chromosome", 
                                pos = "Start_Position", 
                                ref = "Reference_Allele", 
                                alt = "Tumor_Seq_Allele2")
write.table(sigs.input,file.path(res.path,"mutation.sig.input.snp.bydeconstructSigs.txt"),sep = "\t",row.names = T,col.names = NA)

cut.off <- 0.05
mut.wt <- data.frame()
sigs.out.list <- list()
for (sample in rownames(sigs.input)) {
  tmp <- whichSignatures(tumor.ref = sigs.input, 
                         signatures.ref = signatures.cosmic, 
                         sample.id = sample, 
                         contexts.needed = TRUE,
                         tri.counts.method = 'exome2genome',
                         signature.cutoff = cut.off)
  
  #Plot output
  # pdf(file.path(fig.path,paste0(sample,"_plotSignatures.pdf")))
  # plotSignatures(tmp)
  # invisible(dev.off())
  # 
  # pdf(file.path(fig.path,paste0(sample,"_weightPie.pdf")))
  # makePie(tmp)
  # invisible(dev.off())
  # 
  sigs.out.list[[sample]] <- tmp
  tmp <- data.frame(c(tmp$weights,unknown=tmp$unknown),row.names = sample)
  mut.wt <- rbind.data.frame(mut.wt,tmp)
}
write.table(mut.wt,file.path(res.path,"mutation.snp.signature.weightMatrix.bydeconstructSigs.txt"),sep = "\t",row.names = T,col.names = NA)

p.mut <- c()
for (i in 1:30) {
  tmp1 <- mut.wt[CS3.mut,i]
  tmp2 <- mut.wt[CS4.mut,i]
  p.mut <- c(p.mut,wilcox.test(tmp1,tmp2)$p.value)
}
names(p.mut) <- colnames(mut.wt)[1:30]
p.mut

#####################  maftools #####################
require(maftools)

flags <- read.table(file.path(comAnn.path,"Mutation Flags 100.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
tmp <- as.data.frame(mut.maf)
tmp1 <- tmp[-which(tmp$Hugo_Symbol %in% flags$FLAGS),]

CS1.mut.maf <- tmp1[which(tmp1$Tumor_Sample_Barcode %in% CS1),]
CS2.mut.maf <- tmp1[which(tmp1$Tumor_Sample_Barcode %in% CS2),]
CS3.mut.maf <- tmp1[which(tmp1$Tumor_Sample_Barcode %in% CS3),]
CS4.mut.maf <- tmp1[which(tmp1$Tumor_Sample_Barcode %in% CS4),]

clin <- read.table(file.path(comAnn.path,"pancancerSurvivalData_XLu.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
clin <- as.data.frame(na.omit(clin[which(clin$type == "LUAD"),c("OS","OS.time")]))
clin <- clin[which(clin$OS.time > 0),] # get survival time greater than 0
colnames(clin) <- c("fustat","futime")
rownames(clin) <- paste0(rownames(clin),"-01")
clin$Tumor_Sample_Barcode <- rownames(clin)
clin <- clin[rownames(annCol),]
clin$Subtype <- annCol$Subtype

clin1 <- clin[which(rownames(clin) %in% tmp1$Tumor_Sample_Barcode),]
clin1$Tumor_Sample_Barcode <- rownames(clin1)
laml <- read.maf(maf = tmp1, clinicalData= clin1)

library("RColorBrewer")
col = c("#1F78B4","#B2DF8A", "#33A02C" ,"#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")
names(col) = c('Missense_Mutation', 'Nonsense_Mutation', "Frame_Shift_Del" ,'Splice_Site','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Translation_Start_Site','Nonstop_Mutation')
titvCol = RColorBrewer::brewer.pal(n = 6, name = 'Set2')
names(titvCol) = c('T>G','T>A','T>C','C>T','C>G','C>A')
pdf(file.path(fig.path,"LUAD_mutation_plotmafSummary.pdf"),width=8, height=6)
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE , color = col, titvColor = titvCol, fs=0.9,top=8)
invisible(dev.off())
#We will draw oncoplots for top ten mutated genes.
#oncoplot(maf = laml, top = 10, fontSize = 12)

#Changing colors for variant classifications
#brewer.pal(10,"Paired")
col =  RColorBrewer::brewer.pal(n = 8, name = 'Paired')
col =  c("#A6CEE3","#1599ba","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00")
names(col) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
#Color coding for prediction,myasthenia_gravis
cincolors = list(Subtype = c(CS1 = "#2EC4B6", CS2 = "#E71D36", 
              CS3 = "#FF9F1C", CS4 = "#BDD5EA"))
pdf(file.path(fig.path,"LUAD_mutation_oncoplot.pdf"),width=10, height=8)
oncoplot(maf = laml,top = 20, removeNonMutated=F,clinicalFeatures = c("Subtype"), 
         colors = col, sortByAnnotation = TRUE, groupAnnotationBySize=FALSE,annotationColor = cincolors,
         legendFontSize = 1.5 ,titleFontSize = 1.5)
invisible(dev.off())

pdf(file.path(fig.path,"LUAD_mutation_exclusive and Co-occurrence.pdf"),width=6.5, height=6)
somaticInteractions(maf = laml, top = 30, fontSize=0.7,colPal="PiYG",
                    showSum=FALSE, pvalue = c(0.05, 0.1))
invisible(dev.off())

################### CS2 VS CS3 ###################
clinCS2 <- clin[which(clin$Subtype == "CS2"),]
clinCS3 <- clin[which(clin$Subtype == "CS3"),]

mafCS2 = read.maf(maf = CS2.mut.maf, clinicalData= clinCS2 )
mafCS3 = read.maf(maf = CS3.mut.maf, clinicalData= clinCS3 )
pt.vs.rt <- mafCompare(m1 = mafCS3, m2 = mafCS2, m1Name = 'CS3', m2Name = 'CS2', minMut = 25)
pdf(file.path(fig.path,"forestPlot for CS2 vs CS3 mutation.pdf"),width=8, height=8)
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05, color = c('royalblue', 'maroon'))
invisible(dev.off())

col =  c("#1F78B4","#B2DF8A", "#FDBF6F","#33A02C")
names(col) = c('Missense_Mutation', 'Nonsense_Mutation', "In_Frame_Ins" ,'Frame_Shift_Del')
pdf(file.path(fig.path,"LUAD_mutation_oncoplot for CS2 vs CS3.pdf"),width=6, height=5)
coBarplot(m1 = mafCS2, m2 = mafCS3, m1Name = 'CS2', m2Name = 'CS3')
invisible(dev.off())

################### CS2 VS CS3 ###################
clinCS1 <- clin[which(clin$Subtype == "CS1"),]
clinCS4 <- clin[which(clin$Subtype == "CS4"),]

mafCS1 = read.maf(maf = CS1.mut.maf, clinicalData= clinCS1 )
mafCS4 = read.maf(maf = CS4.mut.maf, clinicalData= clinCS4 )
pt.vs.rt <- mafCompare(m1 = mafCS4, m2 = mafCS1, m1Name = 'CS4', m2Name = 'CS1', minMut = 17)
pdf(file.path(fig.path,"forestPlot for CS1 vs CS4 mutation.pdf"),width=8, height=8)
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05, color = c('royalblue', 'maroon'))
invisible(dev.off())

#col =  c("#1F78B4","#B2DF8A", "#FDBF6F","#33A02C")
#names(col) = c('Missense_Mutation', 'Nonsense_Mutation', "In_Frame_Ins" ,'Frame_Shift_Del')
pdf(file.path(fig.path,"LUAD_mutation_oncoplot for CS1 vs CS4.pdf"),width=6, height=5)
coBarplot(m1 = mafCS1, m2 = mafCS4, m1Name = 'CS1', m2Name = 'CS4')
invisible(dev.off())

################### ITH ################### 
library("mclust")
maf <- read_tsv(file.path(data.path,"mc3.v0.2.8.PUBLIC.maf"), comment = "#")
label <- c("Tumor_Sample_Barcode","Hugo_Symbol","NCBI_Build","Chromosome","Start_Position","End_Position","Strand",
           "Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2",
           "HGVSc","HGVSp","HGVSp_Short","BIOTYPE","t_ref_count","t_alt_count","n_alt_count","n_ref_count")
maf$bcr_patient_barcode <- substr(maf$Tumor_Sample_Barcode,start = 1,stop = 12)
tmp <- read.table(file.path(comAnn.path,"pancancerSurvivalData_XLu.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
tmp <- tmp[which(tmp$type %in% c("LUAD")),"bcr_patient_barcode"]
luad.mut.maf <- maf[which(maf$bcr_patient_barcode %in% tmp),label]
luad.mut.maf$Tumor_Sample_Barcode <- substr(luad.mut.maf$Tumor_Sample_Barcode,start = 1,stop = 15)
mut.maf <- as.data.frame(luad.mut.maf)
mut.maf <- mut.maf[which(mut.maf$Tumor_Sample_Barcode %in% rownames(annCol)),]

mut.maf$vaf <- (mut.maf$t_alt_count+mut.maf$n_alt_count)/(mut.maf$t_alt_count+mut.maf$n_alt_count+mut.maf$t_ref_count+mut.maf$n_ref_count)
flags <- read.table(file.path(comAnn.path,"Mutation Flags 100.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
tmp <- as.data.frame(mut.maf)
#tmp1 <- tmp[-which(tmp$Hugo_Symbol %in% flags$FLAGS),]
tmp1 <- tmp
laml <- read.maf(maf = tmp1, clinicalData= clin )

ITH = inferHeterogeneity(maf = laml, tsb = rownames(annCol) , vafCol = 'vaf')
ITH1 <- unique(ITH$clusterData[,c("Tumor_Sample_Barcode","MATH","MedianAbsoluteDeviation")])
ITH.CS1 = ITH1[which(ITH1$Tumor_Sample_Barcode %in% CS1),]
ITH.CS1$Subtype <- "CS1"
ITH.CS2 = ITH1[which(ITH1$Tumor_Sample_Barcode %in% CS2),]
ITH.CS2$Subtype <- "CS2"
ITH.CS3 = ITH1[which(ITH1$Tumor_Sample_Barcode %in% CS3),]
ITH.CS3$Subtype <- "CS3"
ITH.CS4 = ITH1[which(ITH1$Tumor_Sample_Barcode %in% CS4),]
ITH.CS4$Subtype <- "CS4"
ITH2 <- rbind(ITH.CS2,ITH.CS3)

kruskal.test(MATH~Subtype, ITH2)
t.test(MATH~Subtype, ITH2)

library(doBy)  
library(ggplot2)
dat <- summaryBy(MATH~Subtype, ITH2, FUN = c(mean, sd))
p <- ggplot(dat, aes(Subtype, MATH.mean, fill = Subtype)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(aes(ymin = MATH.mean - MATH.sd, ymax = MATH.mean + MATH.sd), width = 0.15, size = 0.5) +
  scale_fill_manual(values=c("CS2" = "#E71D36", 
                             "CS3" = "#FF9F1C")) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5)) +
  labs(x = 'Cluster', y = 'MATH', title = 't-test: p-value = 0.058')
ggsave(file.path(fig.path,"plot for MATH between 2Subtype.pdf"),width = 5,height = 4)

################### SubMap ################### 
source(file.path(comRFun.path,"generateInputFileForSubMap.R"))

gct_file=file.path(res.path,'Submap_gct.gct')
cls_file=file.path(res.path,'Submap_cls.cls')

sam_info=data.frame('rank'=annCol_new$Subtype)
rownames(sam_info)=rownames(annCol_new)
sam_info$type=sam_info$rank
data <- read.table(file.path(res.path,"TCGA_LUAD_mRNA_TPM.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
in_gct=log2(data+1)
in_gct <- in_gct[,rownames(sam_info)]
generateInputFileForSubMap(in_gct, gct_file, cls_file, sam_info, type_name = "type")

submap <- read.table(file.path(res.path,"SubMapResult.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F)
rownames(submap) <- submap$Subtype
library( "pheatmap" )
choose_matrix <- submap[,2:5]
annotation_row = data.frame(Pvalue = factor(rep(c("Nominal.p", "Bonferroni"), c(2,2))))  
rownames(annotation_row) = rownames(choose_matrix)
#mycol <- colorpanel(256,low=blue,mid = "black",high=gold)
mycol <- c("#e5361f","#d16c69","#5587ec","#7968cc","#423293")
library("RColorBrewer")
#brewer.pal(8,"Set2")
ann_colors = list( Pvalue = c(Nominal.p = "#61bb30", Bonferroni = gold))
pheatmap(choose_matrix,color = mycol,cluster_col=F, cluster_rows = F, annotation_row=annotation_row,
         annotation_colors=ann_colors , show_rownames = T,show_colnames = T, 
         gaps_row = 2, annotation_legend = T, legend_breaks=c(0,0.25,0.5,0.75,1),
         filename = file.path(fig.path,"heatmap_submap.pdf"),
         border_color = "white", height=4)
dev.off()

####################### Drug ###########################
library(oncoPredict)
library(pRRophetic)

library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)

set.seed(12345)
trainingExprData = readRDS(file =paste0 (data.path,"/CTRP2_Expr (TPM, log2(x+1) Transformed).rds"))
dim(trainingExprData) #51847 829
trainingPtype = readRDS(file = paste0 (data.path,"/CTRP2_Res.rds"))
dim(trainingPtype) #829 545 
trainingPtype <- trainingPtype[,c("ML162","ML210","erastin")]

testExpr<- read.table(file.path(res.path,"TCGA_LUAD_mRNA_TPM.txt"),sep = "\t",check.names = F, header = T,stringsAsFactors = F,row.names = 1)
testExpr<- as.matrix(testExpr)
dim(testExpr)  

calcPhenotype(trainingExprData = trainingExprData,
              trainingPtype = trainingPtype,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'homogenizeData' )


prefix      = "BOXVIOLIN OF ESTIMATED IC50 (5methods)"
DrugPredictions <- read.csv(file.path("E:/LUAD/movics_pipeline/calcPhenotype_Output/DrugPredictions.csv"),check.names = F, header = T,stringsAsFactors = F,row.names = 1)
DrugPredictions <- DrugPredictions[rownames(annCol),]
Drug <- cbind(annCol,DrugPredictions)
statistic = "kruskal.test"

ic50.test  <- kruskal.test(Drug$ML162 ~ Drug$Subtype)$p.value
pairwise.ic50.test <- pairwise.wilcox.test(Drug$ML162,Drug$Subtype,p.adjust.method = "BH")
cat(paste0("ML162",": Kruskal-Wallis rank sum test p value = ", formatC(ic50.test, format = "e", digits = 2),"\npost-hoc pairwise wilcoxon rank sum test with Benjamini-Hochberg adjustment presents below:\n"))
print(formatC(pairwise.ic50.test$p.value, format = "e", digits = 2))

colvec <- clust.col[1:length(unique(moic.res$clust.res$clust))]
names(colvec) <- paste0("CS",unique(moic.res$clust.res$clust))
p <- ggplot(data = Drug,
            aes(x = Subtype, y = ML162, fill = Subtype)) +
  scale_fill_manual(values = colvec) +
  geom_violin(alpha = 0.4, position = position_dodge(width = .75),
              size = 0.8, color = "black") +
  geom_boxplot(notch = TRUE, outlier.size = -1,
               color = "black", lwd = 0.8, alpha = 0.7) +
  geom_point(shape = 21, size = 2,
             position = position_jitterdodge(),
             color = "black", alpha = 1) +
  theme_classic() +
  ylab(bquote("Estimated IC"[50]~"of"~"ML162")) + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.ticks = element_line(size = 0.2, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10)) +
  # add statistical inference
  stat_compare_means(method = statistic,
                     hjust = ifelse(n.moic %% 2 == 0, 0.5, 0),
                     label.x = ifelse(n.moic %% 2 == 0, n.moic / 2 + 0.5, n.moic / 2),
                     label.y = min(Drug$ML162))
p
outFig <- paste0(prefix, " for ML162", ".pdf")
ggsave(file.path(fig.path, outFig), width = 5, height = 5)


ic50.test  <- kruskal.test(Drug$ML210 ~ Drug$Subtype)$p.value
pairwise.ic50.test <- pairwise.wilcox.test(Drug$ML210,Drug$Subtype,p.adjust.method = "BH")
cat(paste0("ML210",": Kruskal-Wallis rank sum test p value = ", formatC(ic50.test, format = "e", digits = 2),"\npost-hoc pairwise wilcoxon rank sum test with Benjamini-Hochberg adjustment presents below:\n"))
print(formatC(pairwise.ic50.test$p.value, format = "e", digits = 2))

colvec <- clust.col[1:length(unique(moic.res$clust.res$clust))]
names(colvec) <- paste0("CS",unique(moic.res$clust.res$clust))
p <- ggplot(data = Drug,
            aes(x = Subtype, y = ML210, fill = Subtype)) +
  scale_fill_manual(values = colvec) +
  geom_violin(alpha = 0.4, position = position_dodge(width = .75),
              size = 0.8, color = "black") +
  geom_boxplot(notch = TRUE, outlier.size = -1,
               color = "black", lwd = 0.8, alpha = 0.7) +
  geom_point(shape = 21, size = 2,
             position = position_jitterdodge(),
             color = "black", alpha = 1) +
  theme_classic() +
  ylab(bquote("Estimated IC"[50]~"of"~"ML210")) + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.ticks = element_line(size = 0.2, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10)) +
  # add statistical inference
  stat_compare_means(method = statistic,
                     hjust = ifelse(n.moic %% 2 == 0, 0.5, 0),
                     label.x = ifelse(n.moic %% 2 == 0, n.moic / 2 + 0.5, n.moic / 2),
                     label.y = 1.2)
p
outFig <- paste0(prefix, " for ML210", ".pdf")
ggsave(file.path(fig.path, outFig), width = 5, height = 5)


ic50.test  <- kruskal.test(Drug$erastin ~ Drug$Subtype)$p.value
pairwise.ic50.test <- pairwise.wilcox.test(Drug$erastin,Drug$Subtype,p.adjust.method = "BH")
cat(paste0("erastin",": Kruskal-Wallis rank sum test p value = ", formatC(ic50.test, format = "e", digits = 2),"\npost-hoc pairwise wilcoxon rank sum test with Benjamini-Hochberg adjustment presents below:\n"))
print(formatC(pairwise.ic50.test$p.value, format = "e", digits = 2))

colvec <- clust.col[1:length(unique(moic.res$clust.res$clust))]
names(colvec) <- paste0("CS",unique(moic.res$clust.res$clust))
p <- ggplot(data = Drug,
            aes(x = Subtype, y = erastin, fill = Subtype)) +
  scale_fill_manual(values = colvec) +
  geom_violin(alpha = 0.4, position = position_dodge(width = .75),
              size = 0.8, color = "black") +
  geom_boxplot(notch = TRUE, outlier.size = -1,
               color = "black", lwd = 0.8, alpha = 0.7) +
  geom_point(shape = 21, size = 2,
             position = position_jitterdodge(),
             color = "black", alpha = 1) +
  theme_classic() +
  ylab(bquote("Estimated IC"[50]~"of"~"erastin")) + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.ticks = element_line(size = 0.2, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10)) +
  # add statistical inference
  stat_compare_means(method = statistic,
                     hjust = ifelse(n.moic %% 2 == 0, 0.5, 0),
                     label.x = ifelse(n.moic %% 2 == 0, n.moic / 2 + 0.5, n.moic / 2),
                     label.y = 6)
p
outFig <- paste0(prefix, " for erastin", ".pdf")
ggsave(file.path(fig.path, outFig), width = 5, height = 5)