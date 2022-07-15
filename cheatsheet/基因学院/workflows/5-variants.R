## ---- echo=FALSE---------------------------------------------------------
library('variants')

## ---- eval=FALSE---------------------------------------------------------
## library(VariantAnnotation)
## library(cgdv17)
## library(org.Hs.eg.db)
## library(TxDb.Hsapiens.UCSC.hg19.knownGene)
## library(BSgenome.Hsapiens.UCSC.hg19)
## library(PolyPhen.Hsapiens.dbSNP131)

## ---- eval=FALSE---------------------------------------------------------
## ## try http:// if https:// URLs are not supported
## source("https://bioconductor.org/biocLite.R")
## biocLite("mypackage")

## ---- eval=FALSE---------------------------------------------------------
## browseVignettes("cgdv17")

## ------------------------------------------------------------------------
file <- system.file("vcf", "NA06985_17.vcf.gz", package = "cgdv17")

## ------------------------------------------------------------------------
hdr <- scanVcfHeader(file)
 
info(hdr) 
 
geno(hdr) 

## ------------------------------------------------------------------------
meta(hdr)$META

## ------------------------------------------------------------------------
## get entrez ids from gene symbols
genesym <- c("TRPV1", "TRPV2", "TRPV3")
geneid <- select(org.Hs.eg.db, keys=genesym, keytype="SYMBOL",
                 columns="ENTREZID")
geneid

## ------------------------------------------------------------------------
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb

## ------------------------------------------------------------------------
txdb <- renameSeqlevels(txdb, gsub("chr", "", seqlevels(txdb)))
txdb <- keepSeqlevels(txdb, "17")

## ------------------------------------------------------------------------
txbygene = transcriptsBy(txdb, "gene")

## ------------------------------------------------------------------------
gnrng <- unlist(range(txbygene[geneid$ENTREZID]), use.names=FALSE)
names(gnrng) <- geneid$SYMBOL

## ------------------------------------------------------------------------
param <- ScanVcfParam(which = gnrng, info = "DP", geno = c("GT", "cPd"))
param
 
## Extract the TRPV ranges from the VCF file 
vcf <- readVcf(file, "hg19", param)
## Inspect the VCF object with the 'fixed', 'info' and 'geno' accessors
vcf
 
head(fixed(vcf))

geno(vcf)

## ---- eval=FALSE---------------------------------------------------------
## ## Use the 'region' argument to define the region
## ## of interest. See ?locateVariants for details.
## cds <- locateVariants(vcf, txdb, CodingVariants())
## five <- locateVariants(vcf, txdb, FiveUTRVariants())
## splice <- locateVariants(vcf, txdb, SpliceSiteVariants())
## intron <- locateVariants(vcf, txdb, IntronVariants())

## ------------------------------------------------------------------------
all <- locateVariants(vcf, txdb, AllVariants())

## ------------------------------------------------------------------------
## Did any variants match more than one gene?
table(sapply(split(mcols(all)$GENEID, mcols(all)$QUERYID), 
      function(x) length(unique(x)) > 1))
 
## Summarize the number of variants by gene:
idx <- sapply(split(mcols(all)$QUERYID, mcols(all)$GENEID), unique)
sapply(idx, length)
 
## Summarize variant location by gene:
sapply(names(idx), 
    function(nm) {
        d <- all[mcols(all)$GENEID %in% nm, c("QUERYID", "LOCATION")]
        table(mcols(d)$LOCATION[duplicated(d) == FALSE])
    })

## ------------------------------------------------------------------------
seqlevelsStyle(vcf) <- "UCSC"
seqlevelsStyle(txdb) <- "UCSC"
aa <- predictCoding(vcf, txdb, Hsapiens)

## ------------------------------------------------------------------------
## Did any variants match more than one gene?
table(sapply(split(mcols(aa)$GENEID, mcols(aa)$QUERYID), 
        function(x) length(unique(x)) > 1))

## Summarize the number of variants by gene:
idx <- sapply(split(mcols(aa)$QUERYID, mcols(aa)$GENEID, drop=TRUE), unique)
sapply(idx, length)

## Summarize variant consequence by gene:
sapply(names(idx), 
       function(nm) {
           d <- aa[mcols(aa)$GENEID %in% nm, c("QUERYID","CONSEQUENCE")]
           table(mcols(d)$CONSEQUENCE[duplicated(d) == FALSE])
       })

## ----eval=FALSE----------------------------------------------------------
## browseVignettes(package="VariantAnnotation")

## ----eval=FALSE----------------------------------------------------------
## help.start()

## ------------------------------------------------------------------------
sessionInfo()

