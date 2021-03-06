## ---- eval=FALSE---------------------------------------------------------
## ## try http:// if https:// URLs are not supported
## source("https://bioconductor.org/biocLite.R")
## biocLite(c("AnnotationHub", "Homo.sapiens",
##            "TxDb.Hsapiens.UCSC.hg19.knownGene",
##            "BSgenome.Hsapiens.UCSC.hg19", "biomaRt",
##            "TxDb.Athaliana.BioMart.plantsmart22"))

## ---- echo=FALSE---------------------------------------------------------
library(annotation)
ah <- AnnotationHub()

## ---- eval=FALSE---------------------------------------------------------
## ah <- AnnotationHub()

## ------------------------------------------------------------------------
ah

## ------------------------------------------------------------------------
unique(ah$dataprovider)

## ------------------------------------------------------------------------
unique(ah$rdataclass)

## ------------------------------------------------------------------------
grs <- query(ah, "GRanges")
grs

## ------------------------------------------------------------------------
grs <- ah[ah$rdataclass == "GRanges",]

## ------------------------------------------------------------------------
orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgs

## ------------------------------------------------------------------------
meta <- mcols(ah)

## ---- eval=FALSE---------------------------------------------------------
## sah <- display(ah)

## ------------------------------------------------------------------------
res <- grs[[1]]
head(res, n=3)

## ------------------------------------------------------------------------
dog <- query(orgs, "Canis familiaris")[[1]]
dog

## ------------------------------------------------------------------------
columns(dog)

## ------------------------------------------------------------------------
keytypes(dog)

## ------------------------------------------------------------------------
head(keys(dog, keytype="ENTREZID"))

## ------------------------------------------------------------------------
keys(dog, keytype="SYMBOL", pattern="COX")

## ------------------------------------------------------------------------
keys(dog, keytype="ENTREZID", pattern="COX", column="SYMBOL")

## ------------------------------------------------------------------------
select(dog, keys="804478", columns=c("SYMBOL","REFSEQ"), keytype="ENTREZID")

## ------------------------------------------------------------------------
select(dog, keys="804478", columns="GO", keytype="ENTREZID")

## ------------------------------------------------------------------------
mapIds(dog, keys="804478", column="GO", keytype="ENTREZID")

## ------------------------------------------------------------------------
mapIds(dog, keys="804478", column="GO", keytype="ENTREZID", multiVals="list")

## ------------------------------------------------------------------------
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb

## ------------------------------------------------------------------------
txs <- transcripts(txdb)
txs

## ------------------------------------------------------------------------
txby <- transcriptsBy(txdb, by="gene")
txby

## ------------------------------------------------------------------------
si <- seqinfo(txdb)
si

## ------------------------------------------------------------------------
txby <- transcriptsBy(txdb, by="gene")
si <- seqinfo(txby)

## ------------------------------------------------------------------------
head(seqlevels(txdb))
seqlevelsStyle(txdb)
seqlevelsStyle(txdb) <- "NCBI"
head(seqlevels(txdb))

## then change it back
seqlevelsStyle(txdb) <- "UCSC"
head(seqlevels(txdb))

## ------------------------------------------------------------------------
head(isActiveSeq(txdb), n=30)

## ----eval=FALSE----------------------------------------------------------
## isActiveSeq(txdb)["chrY"] <- FALSE
## head(isActiveSeq(txdb), n=26)

## ------------------------------------------------------------------------
Homo.sapiens

## ------------------------------------------------------------------------
select(Homo.sapiens, keys="4488", columns=c("SYMBOL","TXNAME"), keytype="ENTREZID")

## ------------------------------------------------------------------------
txs <- transcripts(Homo.sapiens, columns=c("TXNAME","SYMBOL"))
txs

## ------------------------------------------------------------------------
head(available.genomes())

## ------------------------------------------------------------------------
ls(2)
Hsapiens

## ------------------------------------------------------------------------
seqNms <- seqnames(Hsapiens)
head(seqNms)
getSeq(Hsapiens, seqNms[1:2])

## ----eval=FALSE----------------------------------------------------------
## txby <- transcriptsBy(txdb, by="gene")
## geneOfInterest <- txby[["4488"]]
## res <- getSeq(Hsapiens, geneOfInterest)
## res

## ------------------------------------------------------------------------
head(listMarts())
ensembl <- useMart("ensembl")
ensembl

## ------------------------------------------------------------------------
head(listDatasets(ensembl))
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
ensembl

## ------------------------------------------------------------------------
head(listAttributes(ensembl))

## ------------------------------------------------------------------------
head(getBM(attributes="chromosome_name", mart=ensembl))

## ------------------------------------------------------------------------
head(listFilters(ensembl))

## ------------------------------------------------------------------------
res <- getBM(attributes=c("hgnc_symbol", "entrezgene"),
                    filters = "chromosome_name",
                    values = "1", mart = ensembl)
head(res)

## ------------------------------------------------------------------------
head(columns(ensembl))

## ----eval=FALSE----------------------------------------------------------
## help("makeTxDbPackageFromUCSC")

## ------------------------------------------------------------------------
sessionInfo()

## ------------------------------------------------------------------------
ahs <- query(ah, "UCSC")

## ------------------------------------------------------------------------
ahs <- subset(ahs, ahs$genome=='hg19')
length(ahs)
ahs <- subset(ahs, ahs$species=='Homo sapiens')
length(ahs)

## ------------------------------------------------------------------------
ahs <- query(ah, 'oreganno')
ahs
ahs[1]
oreg <- ahs[['AH5087']]
oreg

## ------------------------------------------------------------------------
keys <- "MSX2"
columns <- c("ENTREZID", "CHR")
select(org.Hs.eg.db, keys, columns, keytype="SYMBOL")

## ------------------------------------------------------------------------
## 1st get all the gene symbols
orgSymbols <- keys(org.Hs.eg.db, keytype="SYMBOL")
## and then use that to get all gene symbols matched to all entrez gene IDs
egr <- select(org.Hs.eg.db, keys=orgSymbols, "ENTREZID", "SYMBOL")
length(egr$ENTREZID)
length(unique(egr$ENTREZID))
## VS:
length(egr$SYMBOL)
length(unique(egr$SYMBOL))
## So lets trap these symbols that are redundant and look more closely...
redund <- egr$SYMBOL
badSymbols <- redund[duplicated(redund)]
select(org.Hs.eg.db, badSymbols, "ENTREZID", "SYMBOL")

## ------------------------------------------------------------------------
res1 <- select(TxDb.Hsapiens.UCSC.hg19.knownGene,
               keys(TxDb.Hsapiens.UCSC.hg19.knownGene, keytype="TXID"),
       	       columns=c("GENEID","TXNAME","TXCHROM"), keytype="TXID")

head(res1)

## ------------------------------------------------------------------------
res2 <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene,
                    columns = c("gene_id","tx_name"))
head(res2)

## ------------------------------------------------------------------------
res <- transcripts(TxDb.Athaliana.BioMart.plantsmart22, columns = c("gene_id"))

## ----eval=FALSE----------------------------------------------------------
## keys <- keys(Homo.sapiens, keytype="TXID")
## res1 <- select(Homo.sapiens,
##                keys= keys,
##        	       columns=c("SYMBOL","TXSTART","TXCHROM"), keytype="TXID")
## 
## head(res1)

## ------------------------------------------------------------------------
res2 <- transcripts(Homo.sapiens, columns="SYMBOL")
head(res2)

## ------------------------------------------------------------------------
columns(Homo.sapiens)
columns(org.Hs.eg.db)
columns(TxDb.Hsapiens.UCSC.hg19.knownGene)
## You might also want to look at this:
transcripts(Homo.sapiens, columns=c("SYMBOL","CHRLOC"))

## ------------------------------------------------------------------------
xk = head(keys(Homo.sapiens, keytype="ENTREZID", pattern="X", column="SYMBOL"))
xk

## ------------------------------------------------------------------------
select(Homo.sapiens, xk, "SYMBOL", "ENTREZID")

## ----eval=FALSE----------------------------------------------------------
## ## Get the transcript ranges grouped by gene
## txby <- transcriptsBy(Homo.sapiens, by="gene")
## ## look up the entrez ID for the gene symbol 'PTEN'
## select(Homo.sapiens, keys='PTEN', columns='ENTREZID', keytype='SYMBOL')
## ## subset that genes transcripts
## geneOfInterest <- txby[["5728"]]
## ## extract the sequence
## res <- getSeq(Hsapiens, geneOfInterest)
## res

## ------------------------------------------------------------------------
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
ids=c("1")
getBM(attributes=c('go_id', 'entrezgene'),
		    filters = 'entrezgene',
                    values = ids, mart = ensembl)


## ------------------------------------------------------------------------
ids=c("1")
select(org.Hs.eg.db, keys=ids, columns="GO", keytype="ENTREZID")

