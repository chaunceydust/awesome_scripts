rm(list=ls())
require(ChIPseeker) 
require(TxDb.Dmelanogaster.UCSC.dm3.ensGene  )
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
require(clusterProfiler) 
(bedFiles=list.files(pattern = '*dm6'))
bedFiles=bedFiles[-2]
library(ChIPpeakAnno)
fs=lapply(bedFiles,function(x){
  peak <- readPeakFile( x )  
  keepChr= !grepl('Het',seqlevels(peak)) 
  seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
  peak
})
library(stringr)
str_split(bedFiles,'_',simplify = T)[,1]
H3K27 <- readPeakFile('H3K27_peaks.narrowPeak')  
keepChr= seqlevels(H3K27) %in% c('2L','2R','3L','3R','4','X','Y','M')
seqlevels(H3K27, pruning.mode="coarse") <- seqlevels(H3K27)[keepChr]
seqlevels(H3K27, pruning.mode="coarse") <- paste0('chr',seqlevels(H3K27))
seqlevels(H3K27)
tmp1=findOverlapsOfPeaks(fs[[1]], H3K27)
tmp2=findOverlapsOfPeaks(fs[[2]], H3K27)
tmp3=findOverlapsOfPeaks(fs[[3]], H3K27)

ol <- findOverlapsOfPeaks(tmp1$peaklist$`fs..1..///H3K27`,
                          tmp2$peaklist$`fs..2..///H3K27`,
                          tmp3$peaklist$`fs..3..///H3K27`)
png('3_factors_overlapVenn_in_domain.png')
makeVennDiagram(ol,
                NameOfPeaks=str_split(bedFiles,'_',simplify = T)[,1],
                TxDb=txdb)
dev.off()

 
lapply(1:7,function(i){
  gr=ol$peaklist[[i]]
  dat=as.data.frame(gr)[,1:3]
  dat[,1]=gsub('chr','',dat[,1])
  write.table(dat,sep = '\t',
              col.names = F,file = paste0('G',i,'_in.bed')
               ,quote = F,row.names = F)
})

ol <- findOverlapsOfPeaks(tmp1$peaklist$fs..1..,
                          tmp2$peaklist$fs..2..,
                          tmp3$peaklist$fs..3..)
png('3_factors_overlapVenn_not_domain.png')
makeVennDiagram(ol,
                NameOfPeaks=str_split(bedFiles,'_',simplify = T)[,1],
                TxDb=txdb)
dev.off()

lapply(1:7,function(i){
  gr=ol$peaklist[[i]]
  dat=as.data.frame(gr)[,1:3]
  dat[,1]=gsub('chr','',dat[,1])
  write.table(dat,sep = '\t',
              col.names = F,file = paste0('G',i,'_out.bed')
              ,quote = F,row.names = F)
})


