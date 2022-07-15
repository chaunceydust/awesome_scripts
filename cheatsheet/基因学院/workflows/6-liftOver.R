## ----doit,echo=FALSE,results="hide"--------------------------------------
library(gwascat)
library(GenomicRanges)
library(rtracklayer)
library(Homo.sapiens)
library(BiocGenerics)
library(liftOver)

## ----lkOne,eval=FALSE----------------------------------------------------
## library(gwascat)
## cur = makeCurrentGwascat()  # result varies by day

## ----lkcur---------------------------------------------------------------
data(cur)
cur

## ----getch---------------------------------------------------------------
library(rtracklayer)
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)
ch
str(ch[[1]])

## ----dolift--------------------------------------------------------------
seqlevelsStyle(cur) = "UCSC"  # necessary
cur19 = liftOver(cur, ch)
class(cur19)

## ----ul------------------------------------------------------------------
cur19 = unlist(cur19)
genome(cur19) = "hg19"
cur19 = new("gwaswloc", cur19)
cur19

## ----lkloss--------------------------------------------------------------
length(cur)-length(cur19)
setdiff(cur$SNPs, cur19$SNPs)

