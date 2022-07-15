bedPeaksFile 		= 'Ez_SppsKO_1_raw_peaks.narrowPeak'; 
bedPeaksFile
## loading packages
require(ChIPseeker)
# BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
# BiocManager::install("org.Dm.eg.db")
require(TxDb.Dmelanogaster.UCSC.dm6.ensGene  )
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
require(clusterProfiler) 
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
print(sampleName)
plotAnnoPie(peakAnno) 
plotAnnoBar(peakAnno)








