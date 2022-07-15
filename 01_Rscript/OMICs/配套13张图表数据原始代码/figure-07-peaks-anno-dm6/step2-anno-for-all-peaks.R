require(ChIPseeker) 
require(TxDb.Dmelanogaster.UCSC.dm6.ensGene  )
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
require(clusterProfiler)  
anno_bed <- function(bedPeaksFile){
  peak <- readPeakFile( bedPeaksFile )  
  seqlevels(peak)
  keepChr= seqlevels(peak) %in% c('2L','2R','3L','3R','4','X','Y','M')
  seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
  seqlevels(peak, pruning.mode="coarse") <- paste0('chr',seqlevels(peak))
  seqlevels(peak)
  cat(paste0('there are ',length(peak),' peaks for this data' ))
  peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), 
                           TxDb=txdb, annoDb="org.Dm.eg.db") 
  peakAnno_df <- as.data.frame(peakAnno)
  sampleName=basename(strsplit(bedPeaksFile,'\\.')[[1]][1])
  return(peakAnno_df)
}

f=list.files(path = './',pattern = 'WT',full.names = T)
f
tmp = lapply(f, anno_bed)
head(tmp[[1]])
df=do.call(rbind,lapply(tmp, function(x){
  #table(x$annotation)
  num1=length(grep('Promoter',x$annotation))
  num2=length(grep("5' UTR",x$annotation))
  num3=length(grep('Exon',x$annotation))
  num4=length(grep('Intron',x$annotation))
  num5=length(grep("3' UTR",x$annotation))
  num6=length(grep('Intergenic',x$annotation))
  return(c(num1,num2,num3,num4,num5,num6 ))
}))
colnames(df)=c('Promoter',"5' UTR",'Exon','Intron',"3' UTR",'Intergenic')
rownames(df)=unlist(lapply(f, function(x){ 
  (strsplit(basename(x),'\\.')[[1]][1])
}))
library(reshape2)
df2=melt(apply(df,1,function(x) x/sum(x)))
colnames(df2)=c('dis','sample','fraction')
library(ggpubr)
ggbarplot(df2, "sample", "fraction",
          fill = "dis", color = "dis", palette = "jco"  )



