rm(list = ls())
options(stringsAsFactors = F)
a=read.table('../RNAseq/all.counts.id.txt',header = T)
dim(a)
dat=a[,7:16]
rownames(dat)=a[,1]
dat[1:4,1:4]
library(stringr)
group_list=str_split(colnames(dat),'_',simplify = T)[,1]
table(group_list)

######################################################################
###################      Firstly for DEseq2      #####################
######################################################################
if(T){
  exprSet=dat
  suppressMessages(library(DESeq2)) 
  (colData <- data.frame(row.names=colnames(exprSet), 
                         group_list=group_list) )
  dds <- DESeqDataSetFromMatrix(countData = exprSet,
                                colData = colData,
                                design = ~ group_list)
  dds <- DESeq(dds)
  png("qc_dispersions.png", 1000, 1000, pointsize=20)
  plotDispEsts(dds, main="Dispersion plot")
  dev.off()
  
  
  rld <- rlogTransformation(dds)
  exprMatrix_rlog=assay(rld) 
  #write.csv(exprMatrix_rlog,'exprMatrix.rlog.csv' )
  
  normalizedCounts1 <- t( t(counts(dds)) / sizeFactors(dds) )
  # normalizedCounts2 <- counts(dds, normalized=T) # it's the same for the tpm value
  # we also can try cpm or rpkm from edgeR pacage
  exprMatrix_rpm=as.data.frame(normalizedCounts1) 
  head(exprMatrix_rpm)
  #write.csv(exprMatrix_rpm,'exprMatrix.rpm.csv' )
  
  png("DEseq_RAWvsNORM.png",height = 800,width = 800)
  par(cex = 0.7)
  n.sample=ncol(exprSet)
  if(n.sample>40) par(cex = 0.5)
  cols <- rainbow(n.sample*1.2)
  par(mfrow=c(2,2))
  boxplot(exprSet, col = cols,main="expression value",las=2)
  boxplot(exprMatrix_rlog, col = cols,main="expression value",las=2)
  hist(as.matrix(exprSet))
  hist(exprMatrix_rlog)
  dev.off()
  
  library(RColorBrewer)
  (mycols <- brewer.pal(8, "Dark2")[1:length(unique(group_list))])
  cor(as.matrix(exprSet))
  # Sample distance heatmap
  sampleDists <- as.matrix(dist(t(exprMatrix_rlog)))
  #install.packages("gplots",repos = "http://cran.us.r-project.org")
  library(gplots)
  png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
  heatmap.2(as.matrix(sampleDists), key=F, trace="none",
            col=colorpanel(100, "black", "white"),
            ColSideColors=mycols[group_list], RowSideColors=mycols[group_list],
            margin=c(10, 10), main="Sample Distance Matrix")
  dev.off()
  
  cor(exprMatrix_rlog) 
  
  table(group_list)
  res <- results(dds, 
                 contrast=c("group_list","SppsKO","WT"))
  resOrdered <- res[order(res$padj),]
  head(resOrdered)
  DEG_SppsKO=as.data.frame(resOrdered)
  DEG_SppsKO=na.omit(DEG_SppsKO)
 
  table(group_list)
  res <- results(dds, 
                 contrast=c("group_list","PhoKO","WT"))
  resOrdered <- res[order(res$padj),]
  head(resOrdered)
  DEG_PhoKO=as.data.frame(resOrdered)
  DEG_PhoKO=na.omit(DEG_PhoKO)
  save(DEG_PhoKO,DEG_SppsKO,file = 'deg_output.Rdata')
}


ac=data.frame(group=group_list)
rownames(ac)=colnames(dat) 


