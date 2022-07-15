## ----style, echo=FALSE, results='hide', message=FALSE--------------------
library(BiocStyle)
library(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)
opts_chunk$set(fig.width=7, fig.height=7)
opts_chunk$set(dpi=300, dev="png", dev.args=list(pointsize=15))

remap <- FALSE # should we remap everything?
redownload <- TRUE # should we re-download from BioC's AWS servers?
keep.files <- FALSE # should we keep BAM files afterwards?
clear.memory <- TRUE # should we clear objects from memory?

## ---- eval=clear.memory, results='hide', echo=FALSE----------------------
protected <- c("blacklist", "gax", "greg",
    "remap", "redownload", "keep.files", "clear.memory",
    "core.loc", "protected", "infile")

## ---- eval=redownload, results='asis', echo=FALSE------------------------
cat('***Note:*** *to cite this article, please refer to https://f1000research.com/articles/4-1080/v2 for instructions.*')

## ------------------------------------------------------------------------
sra.numbers <- c("SRR499718", "SRR499719", "SRR499720", "SRR499721",
    "SRR499734", "SRR499735", "SRR499736", "SRR499737", "SRR499738")
grouping <- c("proB-8113", "proB-8113", "proB-8108", "proB-8108",
    "matureB-8059", "matureB-8059", "matureB-8059", "matureB-8059", "matureB-8086")
all.sra <- paste0(sra.numbers, ".lite.sra")
data.frame(SRA=all.sra, Condition=grouping)

## ---- eval=remap---------------------------------------------------------
## for (sra in all.sra) {
##     code <- system(paste("fastq-dump", sra))
##     stopifnot(code==0L)
## }
## all.fastq <- paste0(sra.numbers, ".fastq")

## ---- eval=remap---------------------------------------------------------
## by.group <- split(all.fastq, grouping)
## for (group in names(by.group)) {
##     code <- system(paste(c("cat", by.group[[group]], ">",
##         paste0(group, ".fastq")), collapse=" "))
##     stopifnot(code==0L)
## }
## group.fastq <- paste0(names(by.group), ".fastq")

## ---- results="hide", eval=remap, message=FALSE--------------------------
## library(Rsubread)
## bam.files <- paste0(names(by.group), ".bam")
## align(index="index/mm10", readfile1=group.fastq, TH1=2, type=1,
##     input_format="FASTQ", output_file=bam.files)

## ---- message=FALSE, results="hide", eval=remap--------------------------
## library(Rsamtools)
## for (bam in bam.files) {
##     out <- suppressWarnings(sortBam(bam, "h3k9ac_temp"))
##     file.rename(out, bam)
## }

## ----killindex, echo=FALSE, results="hide", message=FALSE, eval=remap----
## # For some reason, MarkDuplicates uses BAM index files if they're available; but we don't
## # want it using old indices, so we delete them beforehand if any are present.
## indices <- paste0(bam.files, ".bai")
## exist.indices <- file.exists(indices)
## if (any(exist.indices)) { unlink(indices[exist.indices]) }

## ---- eval=remap---------------------------------------------------------
## temp.bam <- "h3k9ac_temp.bam"
## temp.file <- "h3k9ac_metric.txt"
## temp.dir <- "h3k9ac_working"
## dir.create(temp.dir)
## for (bam in bam.files) {
##     code <- system(sprintf("MarkDuplicates I=%s O=%s M=%s \\
##         TMP_DIR=%s AS=true REMOVE_DUPLICATES=false \\
##         VALIDATION_STRINGENCY=SILENT", bam, temp.bam,
##         temp.file, temp.dir))
##     stopifnot(code==0L)
##     file.rename(temp.bam, bam)
## }

## ---- eval=!remap, results="hide", echo=FALSE, message=FALSE-------------
library(Rsubread)
library(Rsamtools)
by.group <- split(all.sra, grouping)
bam.files <- paste0(names(by.group), ".bam")

## ---- eval=redownload, results="hide", echo=FALSE, message=FALSE---------
core.loc <- "http://s3.amazonaws.com/chipseqdb-bamfiles/"
for (bam in bam.files) { # Downloading all files.
    bam.url <- paste0(core.loc, bam)
    download.file(bam.url, bam)
    download.file(paste0(bam.url, ".bai"), paste0(bam, ".bai"))
}

## ----mapstat-------------------------------------------------------------
diagnostics <- list()
for (bam in bam.files) {
    total <- countBam(bam)$records
    mapped <- countBam(bam, param=ScanBamParam(
        flag=scanBamFlag(isUnmapped=FALSE)))$records
    marked <- countBam(bam, param=ScanBamParam(
        flag=scanBamFlag(isUnmapped=FALSE, isDuplicate=TRUE)))$records
    diagnostics[[bam]] <- c(Total=total, Mapped=mapped, Marked=marked)
}
diag.stats <- data.frame(do.call(rbind, diagnostics))
diag.stats$Prop.mapped <- diag.stats$Mapped/diag.stats$Total*100
diag.stats$Prop.marked <- diag.stats$Marked/diag.stats$Mapped*100
diag.stats

## ---- eval=remap, results = 'hide'---------------------------------------
## indexBam(bam.files)

## ---- echo = FALSE, results = 'hide', eval=remap-------------------------
## # Cleaning up
## unlink(all.fastq)
## unlink(group.fastq)
## unlink(temp.dir, recursive=TRUE)
## unlink(temp.file)

## ---- eval=redownload, message=FALSE, echo=FALSE, results="hide"---------
bed.file <- "mm10.blacklist.bed.gz"
download.file(paste0(core.loc, bed.file), bed.file)

## ---- message=FALSE------------------------------------------------------
library(rtracklayer)
blacklist <- import("mm10.blacklist.bed.gz")
blacklist

## ------------------------------------------------------------------------
celltype <- sub("-.*", "", bam.files)
data.frame(BAM=bam.files, CellType=celltype)

## ---- message=FALSE------------------------------------------------------
library(csaw)
standard.chr <- paste0("chr", c(1:19, "X", "Y"))
param <- readParam(minq=50, discard=blacklist, restrict=standard.chr)

## ------------------------------------------------------------------------
x <- correlateReads(bam.files, param=reform(param, dedup=TRUE))
frag.len <- maximizeCcf(x)
frag.len

## ----ccfplot, fig.cap="**Figure 1:** Cross-correlation function (CCF) against delay distance for the H3K9ac data set. The delay with the maximum correlation is shown as the red line."----
plot(1:length(x)-1, x, xlab="Delay (bp)", ylab="CCF", type="l")
abline(v=frag.len, col="red")
text(x=frag.len, y=min(x), paste(frag.len, "bp"), pos=4, col="red")

## ------------------------------------------------------------------------
win.data <- windowCounts(bam.files, param=param, width=150, ext=frag.len)
win.data

## ------------------------------------------------------------------------
bins <- windowCounts(bam.files, bin=TRUE, width=2000, param=param)
filter.stat <- filterWindows(win.data, bins, type="global")
min.fc <- 3
keep <- filter.stat$filter > log2(min.fc)
summary(keep)

## ----bghistplot, fig.cap="**Figure 2:** Histogram of average abundances across all 2 kbp genomic bins. The filter threshold is shown as the red line."----
hist(filter.stat$back.abundances, main="", breaks=50,
    xlab="Background abundance (log2-CPM)")
threshold <- filter.stat$abundances[1] - filter.stat$filter[1] + log2(min.fc)
abline(v=threshold, col="red")

## ------------------------------------------------------------------------
filtered.data <- win.data[keep,]

## ----trendplot, fig.cap="**Figure 3:** Abundance-dependent trend in the log-fold change between two H3K9ac libraries (mature B over pro-B), across all windows retained after filtering."----
library(edgeR)
win.ab <- aveLogCPM(asDGEList(filtered.data))
adjc <- log2(assay(filtered.data)+0.5)
logfc <- adjc[,1] - adjc[,4]
smoothScatter(win.ab, logfc, ylim=c(-6, 6), xlim=c(0, 5),
    xlab="Average abundance", ylab="Log-fold change")

## ------------------------------------------------------------------------
offsets <- normOffsets(filtered.data, type="loess")
head(offsets)

## ----normplot, fig.cap="**Figure 4:** Effect of non-linear normalization on the trended bias between two H3K9ac libraries. Normalized log-fold changes are shown for all windows retained after filtering."----
norm.adjc <- adjc - offsets/log(2)
norm.fc <- norm.adjc[,1]-norm.adjc[,4]
smoothScatter(win.ab, norm.fc, ylim=c(-6, 6), xlim=c(0, 5),
    xlab="Average abundance", ylab="Log-fold change")

## ------------------------------------------------------------------------
celltype <- factor(celltype)
design <- model.matrix(~0+celltype)
colnames(design) <- levels(celltype)
design

## ---- message=FALSE, echo=FALSE, results="hide"--------------------------
library(statmod)
library(locfit)

## ---- message=FALSE------------------------------------------------------
library(edgeR)
y <- asDGEList(filtered.data)
y <- scaleOffset(y, offsets)
y <- estimateDisp(y, design)
summary(y$trended.dispersion)

## ----bcvplot, fig.cap="**Figure 5:** Abundance-dependent trend in the BCV for each window, represented by the blue line. Common (red) and tagwise estimates (black) are also shown."----
plotBCV(y)

## ----qlplot, fig.cap="**Figure 6:** Effect of EB shrinkage on the raw QL dispersion estimate for each window (black) towards the abundance-dependent trend (blue) to obtain squeezed estimates (red)."----
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)

## ------------------------------------------------------------------------
summary(fit$df.prior)

## ----mdsplot, fig.cap="**Figure 7:** MDS plot with two dimensions for all libraries in the H3K9ac data set. Libraries are labelled and coloured according to the cell type."----
plotMDS(norm.adjc, labels=celltype,
    col=c("red", "blue")[as.integer(celltype)])

## ------------------------------------------------------------------------
contrast <- makeContrasts(proB-matureB, levels=design)
res <- glmQLFTest(fit, contrast=contrast)
head(res$table)

## ------------------------------------------------------------------------
merged <- mergeWindows(rowRanges(filtered.data), tol=100, max.width=5000)

## ------------------------------------------------------------------------
tabcom <- combineTests(merged$id, res$table)
head(tabcom)

## ------------------------------------------------------------------------
is.sig <- tabcom$FDR <= 0.05
summary(is.sig)

## ------------------------------------------------------------------------
table(tabcom$direction[is.sig])

## ------------------------------------------------------------------------
tabbest <- getBestTest(merged$id, res$table)
head(tabbest)

## ------------------------------------------------------------------------
is.sig.pos <- (tabbest$logFC > 0)[is.sig]
summary(is.sig.pos)

## ------------------------------------------------------------------------
out.ranges <- merged$region
elementMetadata(out.ranges) <- data.frame(tabcom,
    best.pos=mid(ranges(rowRanges(filtered.data[tabbest$best]))),
    best.logFC=tabbest$logFC)
saveRDS(file="h3k9ac_results.rds", out.ranges)

## ------------------------------------------------------------------------
simplified <- out.ranges[is.sig]
simplified$score <- -10*log10(simplified$FDR)
export(con="h3k9ac_results.bed", object=simplified)

## ------------------------------------------------------------------------
save(file="h3k9ac_objects.Rda", win.data, bins, y)

## ---- message=FALSE------------------------------------------------------
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
anno <- detailRanges(out.ranges, orgdb=org.Mm.eg.db,
    txdb=TxDb.Mmusculus.UCSC.mm10.knownGene)
head(anno$overlap)

## ------------------------------------------------------------------------
head(anno$left)
head(anno$right)

## ------------------------------------------------------------------------
meta <- elementMetadata(out.ranges)
elementMetadata(out.ranges) <- data.frame(meta, anno)

## ---- message=FALSE------------------------------------------------------
library(ChIPpeakAnno)
data(TSS.mouse.GRCm38)
minimal <- out.ranges
elementMetadata(minimal) <- NULL
anno.regions <- annotatePeakInBatch(minimal, AnnotationData=TSS.mouse.GRCm38)
colnames(elementMetadata(anno.regions))

## ------------------------------------------------------------------------
prom <- suppressWarnings(promoters(TxDb.Mmusculus.UCSC.mm10.knownGene,
    upstream=3000, downstream=1000, columns=c("tx_name", "gene_id")))
entrez.ids <- sapply(prom$gene_id, FUN=function(x) x[1]) # Using the first Entrez ID.
gene.name <- select(org.Mm.eg.db, keys=entrez.ids, keytype="ENTREZID", column="SYMBOL")
prom$gene_name <- gene.name$SYMBOL[match(entrez.ids, gene.name$ENTREZID)]
head(prom)

## ------------------------------------------------------------------------
olap <- findOverlaps(prom, rowRanges(filtered.data))
tabprom <- combineOverlaps(olap, res$table)
head(data.frame(ID=prom$tx_name, Gene=prom$gene_name, tabprom)[!is.na(tabprom$PValue),])

## ---- message=FALSE------------------------------------------------------
library(Gviz)
gax <- GenomeAxisTrack(col="black", fontsize=15, size=2)
greg <- GeneRegionTrack(TxDb.Mmusculus.UCSC.mm10.knownGene, showId=TRUE,
    geneSymbol=TRUE, name="", background.title="transparent")
symbols <- unlist(mapIds(org.Mm.eg.db, gene(greg), "SYMBOL",
    "ENTREZID", multiVals = "first"))
symbol(greg) <- symbols[gene(greg)]

## ------------------------------------------------------------------------
o <- order(out.ranges$PValue)
cur.region <- out.ranges[o[1]]
cur.region

## ---- echo=FALSE, results="hide"-----------------------------------------
if (!overlapsAny(cur.region, GRanges("chr17", IRanges(34285101, 34289950)), type="equal")) {
    warning("first region does not match expectations")
}

## ----simplebroadplot, fig.width=8, fig.height=6, fig.cap="**Figure 8:** Coverage tracks for a simple DB event between pro-B and mature B cells, across a broad region in the H3K9ac data set. Read coverage for each library is shown as a per-million value at each base."----
collected <- list()
lib.sizes <- filtered.data$totals/1e6
for (i in 1:length(bam.files)) {
    reads <- extractReads(bam.file=bam.files[i], cur.region, param=param)
    cov <- as(coverage(reads)/lib.sizes[i], "GRanges")
    collected[[i]] <- DataTrack(cov, type="histogram", lwd=0, ylim=c(0,10),
        name=bam.files[i], col.axis="black", col.title="black",
        fill="darkgray", col.histogram=NA)
}
plotTracks(c(gax, collected, greg), chromosome=as.character(seqnames(cur.region)),
    from=start(cur.region), to=end(cur.region))

## ------------------------------------------------------------------------
complex <- out.ranges$logFC.up > 0 & out.ranges$logFC.down > 0
cur.region <- out.ranges[o[complex[o]][2]]
cur.region

## ---- echo=FALSE, results="hide"-----------------------------------------
if (!overlapsAny(cur.region, GRanges("chr5", IRanges(122987201, 122991450)), type="equal")) {
    warning("second region does not match expectations")
}

## ----complexplot, fig.width=8, fig.height=6, fig.cap="**Figure 9:** Coverage tracks for a complex DB event in the H3K9ac data set, shown as per-million values."----
collected <- list()
for (i in 1:length(bam.files)) {
    reads <- extractReads(bam.file=bam.files[i], cur.region, param=param)
    cov <- as(coverage(reads)/lib.sizes[i], "GRanges")
    collected[[i]] <- DataTrack(cov, type="histogram", lwd=0, ylim=c(0,3),
        name=bam.files[i], col.axis="black", col.title="black",
        fill="darkgray", col.histogram=NA)
}
plotTracks(c(gax, collected, greg), chromosome=as.character(seqnames(cur.region)),
    from=start(cur.region), to=end(cur.region))

## ------------------------------------------------------------------------
sharp <- out.ranges$nWindows < 20
cur.region <- out.ranges[o[sharp[o]][1]]
cur.region

## ---- echo=FALSE, results="hide"-----------------------------------------
if (!overlapsAny(cur.region, GRanges("chr16", IRanges(36665551, 36666200)), type="equal")) {
    warning("second region does not match expectations")
}

## ----simplesharpplot, fig.width=8, fig.height=6, fig.cap="**Figure 10:** Coverage tracks for a sharp and simple DB event in the H3K9ac data set, shown as per-million values."----
collected <- list()
for (i in 1:length(bam.files)) {
    reads <- extractReads(bam.file=bam.files[i], cur.region, param=param)
    cov <- as(coverage(reads)/lib.sizes[i], "GRanges")
    collected[[i]] <- DataTrack(cov, type="histogram", lwd=0, ylim=c(0,3),
        name=bam.files[i], col.axis="black", col.title="black",
        fill="darkgray", col.histogram=NA)
}
plotTracks(c(gax, collected, greg), chromosome=as.character(seqnames(cur.region)),
    from=start(cur.region), to=end(cur.region))

## ---- eval=!keep.files, echo=FALSE, results="hide"-----------------------
unlink(bam.files)
unlink(paste0(bam.files, ".bai"))

## ---- eval=clear.memory, echo=FALSE, results="hide"----------------------
rm(list=setdiff(ls(), protected))
gc()

## ------------------------------------------------------------------------
sra.numbers <- c("SRR1145787", "SRR1145788", "SRR1145789", "SRR1145790")
genotype <- c("wt", "wt", "ko", "ko")
all.sra <- paste0(sra.numbers, ".sra")
data.frame(SRA=all.sra, Condition=genotype)

## ---- eval=remap---------------------------------------------------------
## for (sra in all.sra) {
##     code <- system(paste("fastq-dump", sra))
##     stopifnot(code==0L)
## }
## all.fastq <- paste0(sra.numbers, ".fastq")

## ---- eval=remap, message=FALSE, results="hide"--------------------------
## bam.files <- paste0(sra.numbers, ".bam")
## align(index="index/mm10", readfile1=all.fastq, type=1, phredOffset=64,
##     input_format="FASTQ", output_file=bam.files)

## ---- ref.label="killindex", echo=FALSE, results="hide", message=FALSE, eval=remap----
## NA

## ---- eval=remap, results="hide"-----------------------------------------
## temp.bam <- "cbp_temp.bam"
## temp.file <- "cbp_metric.txt"
## temp.dir <- "cbp_working"
## dir.create(temp.dir)
## for (bam in bam.files) {
##     out <- suppressWarnings(sortBam(bam, "cbp_temp"))
##     file.rename(out, bam)
##     code <- system(sprintf("MarkDuplicates I=%s O=%s M=%s \\
##         TMP_DIR=%s AS=true REMOVE_DUPLICATES=false \\
##         VALIDATION_STRINGENCY=SILENT",
##         bam, temp.bam, temp.file, temp.dir))
##     stopifnot(code==0L)
##     file.rename(temp.bam, bam)
## }
## indexBam(bam.files)

## ---- eval=remap, echo = FALSE, results = 'hide'-------------------------
## # Cleaning up
## unlink(all.fastq)
## unlink(temp.dir, recursive=TRUE)
## unlink(temp.file)

## ---- echo=FALSE, results="hide", eval=!remap, message=FALSE-------------
bam.files <- paste0(sra.numbers, ".bam")

## ---- echo=FALSE, results="hide", eval=redownload, message=FALSE---------
for (bam in bam.files) {
    bam.url <- paste0(core.loc, bam)
    download.file(bam.url, bam)
    download.file(paste0(bam.url, ".bai"), paste0(bam, ".bai"))
}

## ---- ref.label="mapstat", echo=FALSE------------------------------------

## ------------------------------------------------------------------------
param <- readParam(minq=50, discard=blacklist)

## ------------------------------------------------------------------------
x <- correlateReads(bam.files, param=reform(param, dedup=TRUE))
frag.len <- maximizeCcf(x)
frag.len

## ------------------------------------------------------------------------
win.data <- windowCounts(bam.files, param=param, width=10, ext=frag.len)
win.data

## ------------------------------------------------------------------------
bins <- windowCounts(bam.files, bin=TRUE, width=10000, param=param)
normfacs <- normOffsets(bins)
normfacs

## ----compoplot, fig.width=12, fig.height=6, fig.cap="**Figure 11:** Mean-difference plots for the bin counts, comparing library 4 to all other libraries. The red line represents the log-ratio of the normalization factors between libraries."----
y.bin <- asDGEList(bins)
bin.ab <- aveLogCPM(y.bin)
adjc <- cpm(y.bin, log=TRUE)
par(cex.lab=1.5, mfrow=c(1,3))
smoothScatter(bin.ab, adjc[,1]-adjc[,4], ylim=c(-6, 6),
    xlab="Average abundance", ylab="Log-ratio (1 vs 4)")
abline(h=log2(normfacs[1]/normfacs[4]), col="red")
smoothScatter(bin.ab, adjc[,2]-adjc[,4], ylim=c(-6, 6),
    xlab="Average abundance", ylab="Log-ratio (2 vs 4)")
abline(h=log2(normfacs[2]/normfacs[4]), col="red")
smoothScatter(bin.ab, adjc[,3]-adjc[,4], ylim=c(-6, 6),
    xlab="Average abundance", ylab="Log-ratio (3 vs 4)")
abline(h=log2(normfacs[3]/normfacs[4]), col="red")

## ------------------------------------------------------------------------
filter.stat <- filterWindows(win.data, bins, type="global")
min.fc <- 3
keep <- filter.stat$filter > log2(min.fc)
summary(keep)
filtered.data <- win.data[keep,]

## ------------------------------------------------------------------------
genotype <- factor(genotype)
design <- model.matrix(~0+genotype)
colnames(design) <- levels(genotype)
design

## ------------------------------------------------------------------------
y <- asDGEList(filtered.data, norm.factors=normfacs)
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$df.prior)

## ----mdsplot2, fig.cap="**Figure 12:** MDS plot with two dimensions for all libraries in the CBP data set. Libraries are labelled and coloured according to the genotype. A larger top set of windows was used to improve the visualization of the genome-wide differences between the WT libraries."----
plotMDS(cpm(y, log=TRUE), top=10000, labels=genotype,
    col=c("red", "blue")[as.integer(genotype)])

## ------------------------------------------------------------------------
contrast <- makeContrasts(wt-ko, levels=design)
res <- glmQLFTest(fit, contrast=contrast)
merged <- mergeWindows(rowRanges(filtered.data), tol=100, max.width=5000)
tabcom <- combineTests(merged$id, res$table)
tabbest <- getBestTest(merged$id, res$table)
is.sig <- tabcom$FDR <= 0.05
summary(is.sig)
table(tabcom$direction[is.sig])
is.sig.pos <- (tabbest$logFC > 0)[is.sig]
summary(is.sig.pos)

## ------------------------------------------------------------------------
out.ranges <- merged$region
elementMetadata(out.ranges) <- data.frame(tabcom,
    best.pos=mid(ranges(rowRanges(filtered.data[tabbest$best]))),
    best.logFC=tabbest$logFC)
saveRDS(file="cbp_results.rds", out.ranges)
save(file="cbp_objects.Rda", win.data, bins, y)

## ------------------------------------------------------------------------
anno <- detailRanges(out.ranges, orgdb=org.Mm.eg.db,
            txdb=TxDb.Mmusculus.UCSC.mm10.knownGene)
meta <- elementMetadata(out.ranges)
elementMetadata(out.ranges) <- data.frame(meta, anno)

## ------------------------------------------------------------------------
o <- order(out.ranges$PValue)    
cur.region <- out.ranges[o[2]]
cur.region

## ---- results="hide", echo=FALSE-----------------------------------------
if (!overlapsAny(cur.region, GRanges("chr16", IRanges(70313851, 70314860)), type="equal")) {
        warning("first region does not match expectations")
}

## ----tfplot, fig.height=6, fig.width=8, fig.cap="**Figure 13:** Coverage tracks for TF binding sites that are differentially bound in the WT (top two tracks) against the KO (last two tracks). Blue and red tracks represent forward- and reverse-strand coverage, respectively, on a per-million scale (capped at 5 in SRR1145788, for visibility)."----
collected <- list()
lib.sizes <- filtered.data$totals/1e6
for (i in 1:length(bam.files)) {
    reads <- extractReads(bam.file=bam.files[i], cur.region, param=param)
    pcov <- as(coverage(reads[strand(reads)=="+"])/lib.sizes[i], "GRanges")
    ncov <- as(coverage(reads[strand(reads)=="-"])/-lib.sizes[i], "GRanges")
    ptrack <- DataTrack(pcov, type="histogram", lwd=0, ylim=c(-5, 5),
        name=bam.files[i], col.axis="black", col.title="black",
        fill="blue", col.histogram=NA)
    ntrack <- DataTrack(ncov, type="histogram", lwd=0, ylim=c(-5, 5),
        fill="red", col.histogram=NA)
    collected[[i]] <- OverlayTrack(trackList=list(ptrack, ntrack))
}
plotTracks(c(gax, collected, greg), chromosome=as.character(seqnames(cur.region)),
    from=start(cur.region), to=end(cur.region))

## ---- eval=!keep.files, echo=FALSE, results="hide"-----------------------
unlink(bam.files)
unlink(paste0(bam.files, ".bai"))

## ---- eval=clear.memory, echo=FALSE, results="hide"----------------------
rm(list=setdiff(ls(), protected))
gc()

## ---- eval=!redownload, results='asis', echo=FALSE-----------------------
## cat('Installation of required packages and execution of the workflow can be performed by following the instructions at https://www.bioconductor.org/help/workflows/chipseqDB.')

## ------------------------------------------------------------------------
sessionInfo()

