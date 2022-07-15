bedPeaksFile 		= 'oldBedFiles/Cg_WT.narrowPeak.bed'; 
bedPeaksFile
## loading packages
require(ChIPseeker)
# BiocManager::install("TxDb.Dmelanogaster.UCSC.dm3.ensGene")
# BiocManager::install("org.Dm.eg.db")
require(TxDb.Dmelanogaster.UCSC.dm3.ensGene  )
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
require(clusterProfiler) 
peak <- readPeakFile( bedPeaksFile )  
keepChr= !grepl('Het',seqlevels(peak)) 
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
cat(paste0('there are ',length(peak),' peaks for this data' ))
peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), 
                         TxDb=txdb, annoDb="org.Dm.eg.db") 
peakAnno_df <- as.data.frame(peakAnno)
sampleName=basename(strsplit(bedPeaksFile,'\\.')[[1]][1])
print(sampleName)
plotAnnoPie(peakAnno) 
plotAnnoBar(peakAnno)








