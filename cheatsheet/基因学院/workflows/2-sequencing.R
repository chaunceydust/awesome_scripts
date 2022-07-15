## ----style, echo=FALSE, message=FALSE, warning=FALSE, results="asis"-----
options(width=100)
knitr::opts_chunk$set(message = FALSE, error = FALSE, warning = FALSE,
fig.width=6, fig.height=4)
BiocStyle::markdown()

## ---- echo=FALSE, results="hide", warning=FALSE--------------------------
suppressPackageStartupMessages({
   library(GenomicRanges)
   library(GenomicAlignments) 
   library(Biostrings)
   library(Rsamtools)
   library(ShortRead)
   library(BiocParallel)
   library(rtracklayer)
   library(VariantAnnotation)
   library(AnnotationHub)
   library(BSgenome.Hsapiens.UCSC.hg19)
   library(RNAseqData.HNRNPC.bam.chr14)
})
ah = AnnotationHub()


## ----eval=FALSE----------------------------------------------------------
## help(package="GenomicRanges")
## vignette(package="GenomicRanges")

## ----message=FALSE-------------------------------------------------------
library(Biostrings)
d <- DNAString("TTGAAAA-CTC-N")
length(d)  #no of letters in the DNAString

## ----eval=FALSE----------------------------------------------------------
## library(AnnotationHub)
## ah <- AnnotationHub()

## ------------------------------------------------------------------------
ah2 <- query(ah, c("fasta", "homo sapiens", "Ensembl"))
fa <- ah2[["AH18522"]]
fa

## ------------------------------------------------------------------------
readDNAStringSet(path(fa))

## ------------------------------------------------------------------------
library(Rsamtools)
idx <- scanFaIndex(fa)
idx

## ------------------------------------------------------------------------
long <- idx[width(idx) > 82000]
getSeq(fa, param=long)

## ------------------------------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg19)

chr14_range = GRanges("chr14", IRanges(1, seqlengths(Hsapiens)["chr14"]))
chr14_dna <- getSeq(Hsapiens, chr14_range)
letterFrequency(chr14_dna, "GC", as.prob=TRUE)

## ----eval=FALSE----------------------------------------------------------
## ## 1. attach ShortRead and BiocParallel
## library(ShortRead)
## library(BiocParallel)
## 
## ## 2. create a vector of file paths
## fls <- dir("~/fastq", pattern="*fastq", full=TRUE)
## 
## ## 3. collect statistics
## stats0 <- qa(fls)
## 
## ## 4. generate and browse the report
## if (interactive())
##     browseURL(report(stats0))

## ------------------------------------------------------------------------
## 1. load software packages
library(GenomicRanges)
library(GenomicAlignments)

## 2. load sample data
library('RNAseqData.HNRNPC.bam.chr14')
bf <- BamFile(RNAseqData.HNRNPC.bam.chr14_BAMFILES[[1]], asMates=TRUE)

## 3. define our 'region of interest'
roi <- GRanges("chr14", IRanges(19653773, width=1)) 

## 4. alignments, junctions, overlapping our roi
paln <- readGAlignmentsList(bf)
j <- summarizeJunctions(paln, with.revmap=TRUE)
j_overlap <- j[j %over% roi]

## 5. supporting reads
paln[j_overlap$revmap[[1]]]

## ------------------------------------------------------------------------
library(VariantAnnotation)
fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
vcf <- readVcf(fl, "hg19")

## ------------------------------------------------------------------------
library(rtracklayer)
test_path <- system.file("tests", package = "rtracklayer")
test_bed <- file.path(test_path, "test.bed")
  
test <- import(test_bed, format = "bed")
test

## ------------------------------------------------------------------------
sessionInfo()

