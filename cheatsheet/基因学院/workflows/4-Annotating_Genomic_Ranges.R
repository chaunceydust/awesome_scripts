## ---- echo=FALSE---------------------------------------------------------
suppressPackageStartupMessages(library(annotation))

## ------------------------------------------------------------------------
hub <- AnnotationHub()

## ------------------------------------------------------------------------
mcols(query(hub, "clinvar.vcf", "GRCh37"))[,"sourceurl", drop=FALSE]

## ------------------------------------------------------------------------
fl <- query(hub, "clinvar.vcf", "GRCh37")[[1]]

## ------------------------------------------------------------------------
vcf <- readVcf(fl, "hg19")
dim(vcf)

## ------------------------------------------------------------------------
txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
head(seqlevels(txdb_hg19))
seqlevels(vcf)
seqlevels(vcf) <- paste0("chr", seqlevels(vcf))

## ------------------------------------------------------------------------
seqlevels(vcf, pruning.mode="coarse") <- c("chr3", "chr18")
seqlevels(txdb_hg19) <- c("chr3", "chr18")

## ------------------------------------------------------------------------
intersect(seqlevels(txdb_hg19), seqlevels(vcf))

## ------------------------------------------------------------------------
unique(genome(txdb_hg19))
unique(genome(vcf))

## ------------------------------------------------------------------------
gr_hg19 <- rowRanges(vcf)

## ------------------------------------------------------------------------
txdb_mm10 <- keepStandardChromosomes(TxDb.Mmusculus.UCSC.mm10.ensGene)

## ------------------------------------------------------------------------
head(seqlevels(txdb_mm10))
gr_mm10 <- GRanges("chr4", IRanges(c(4000000, 107889000), width=1000))

## ------------------------------------------------------------------------
unique(genome(txdb_mm10))
genome(gr_mm10) <- "mm10"

## ------------------------------------------------------------------------
loc_hg19 <- locateVariants(gr_hg19, txdb_hg19, AllVariants())
table(loc_hg19$LOCATION)
loc_mm10 <- locateVariants(gr_mm10, txdb_mm10, AllVariants()) 
table(loc_mm10$LOCATION)

## ------------------------------------------------------------------------
cols <- c("UNIPROT", "PFAM")
keys <- na.omit(unique(loc_hg19$GENEID))
head(select(org.Hs.eg.db, keys, cols, keytype="ENTREZID"))

## ------------------------------------------------------------------------
keys <- unique(loc_mm10$GENEID)
head(select(org.Mm.eg.db, keys, cols, keytype="ENSEMBL"))

## ------------------------------------------------------------------------
hub <- AnnotationHub()
hub_hg19 <- subset(hub, 
                  (hub$species == "Homo sapiens") & (hub$genome == "hg19"))
length(hub_hg19)

## ---- echo=FALSE---------------------------------------------------------
ov_hg19 <- lapply(1:3, function(i) subsetByOverlaps(hub_hg19[[i]], gr_hg19))

## ------------------------------------------------------------------------
ov_hg19 <- lapply(1:3, function(i) subsetByOverlaps(hub_hg19[[i]], gr_hg19))

## ------------------------------------------------------------------------
names(ov_hg19) <- names(hub_hg19)[1:3]
lapply(ov_hg19, head, n=3)

## ------------------------------------------------------------------------
head(predictCoding(vcf, txdb_hg19, Hsapiens), 3)

## ----sess----------------------------------------------------------------
sessionInfo()

