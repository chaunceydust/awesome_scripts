rm(list=ls())
require(ChIPseeker) 
require(TxDb.Dmelanogaster.UCSC.dm3.ensGene  )
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
require(clusterProfiler) 
(bedFiles=list.files(pattern = '*dm6'))

library(ChIPpeakAnno)
fs=lapply(bedFiles,function(x){
  peak <- readPeakFile( x )  
  keepChr= !grepl('Het',seqlevels(peak)) 
  seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
  peak
})
library(stringr)
ol <- findOverlapsOfPeaks(fs[[1]],fs[[2]],fs[[3]] )
png('3_factors_overlapVenn.png')
makeVennDiagram(ol,
                NameOfPeaks=str_split(bedFiles,'_',simplify = T)[,1],
                TxDb=txdb)
dev.off()
