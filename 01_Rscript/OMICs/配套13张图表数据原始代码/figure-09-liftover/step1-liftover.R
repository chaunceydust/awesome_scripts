# https://genome.ucsc.edu/cgi-bin/hgLiftOver
library(rtracklayer)
path='dm3ToDm6.over.chain'
ch = import.chain(path)
ch
 
require(ChIPseeker) 
require(TxDb.Dmelanogaster.UCSC.dm3.ensGene  )
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
require(clusterProfiler) 
library(ChIPpeakAnno)
peak_dm3 <- readPeakFile('Pho_WT.narrowPeak.bed')  

seqlevelsStyle(peak_dm3) = "UCSC"  # necessary
peak_dm6 = liftOver(peak_dm3, ch)
class(peak_dm6)
peak_dm3
peak_dm6



