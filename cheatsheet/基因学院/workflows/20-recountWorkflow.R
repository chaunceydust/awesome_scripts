## ----bioc_setup, include = FALSE-----------------------------------------
on.bioc <- TRUE
knitr::opts_chunk$set(fig.path = "")

## ----Figure1, out.width="100%", fig.align="center", fig.cap = "Overview of the data available in recount2. Reads (pink boxes) aligned to the reference genome can be used to compute a base-pair coverage curve and identify exon-exon junctions (split reads). Gene and exon count matrices are generated using annotation information providing the gene (green boxes) and exon (blue boxes) coordinates together with the base-level coverage curve. The reads spanning exon-exon junctions (jx) are used to compute a third count matrix that might include unannotated junctions (jx 3 and 4). Without using annotation information, expressed regions (orange box) can be determined from the base-level coverage curve to then construct data-driven count matrices.", echo = FALSE----
knitr::include_graphics("Figure1.png")

## ----Figure2, out.width="100%", fig.align="center", fig.cap = "recount2 provides coverage count matrices in RangedSummarizedExperiment (rse) objects. Once the rse object has been downloaded and loaded into R, the feature information is accessed with rowRanges(rse) (blue box), the counts with assays(rse)\\$counts (pink box) and the sample metadata with colData(rse) (green box). The sample metadata can be expanded using add\\_predictions(rse) (orange box) or with custom code (brown box) matching by a unique sample identifier such as the SRA Run ID. The rse object is inside the purple box and matching data is highlighted in each box.", echo = FALSE----
knitr::include_graphics("Figure2.png")

## ----"install", eval = FALSE---------------------------------------------
## ## Install packages from Bioconductor
## source("https://bioconductor.org/biocLite.R")
## biocLite(c("recount", "GenomicRanges", "limma", "edgeR", "DESeq2",
##     "regionReport", "clusterProfiler", "org.Hs.eg.db", "gplots",
##     "derfinder", "rtracklayer", "GenomicFeatures", "bumphunter",
##     "derfinderPlot", "devtools"))

## ----"load libraries", message = FALSE, warning = FALSE------------------
library("recount")
library("GenomicRanges")
library("limma")
library("edgeR")
library("DESeq2")
library("regionReport")
library("clusterProfiler")
library("org.Hs.eg.db")
library("gplots")
library("derfinder")
library("rtracklayer")
library("GenomicFeatures")
library("bumphunter")
library("derfinderPlot")
library("devtools")

## ----Figure3, out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "RNA-seq starting data. 16 RNA-seq un-aligned RNA-seq reads 3 base-pairs long are shown (pink boxes) alongside a reference genome that is 16 base-pairs long (white box).", echo = FALSE----
knitr::include_graphics("Figure3.png")

## ----Figure4, out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Aligned RNA-seq reads. Spice-aware RNA-seq aligners such as Rail-RNA are able to find the coordinates to which the reads map, even if they span exon-exon junctions (connected boxes). Rail-RNA soft clips some reads (purple boxes with rough edges) such that a portion of these reads align to the reference genome.", echo = FALSE----
knitr::include_graphics("Figure4.png")

## ----Figure5, out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Gene annotation. A single gene with two isoforms composed by three distinct exons (blue boxes) is illustrated. Exons 1 and 3 share the first five base-pairs while exon 2 is common to both isoforms.", echo = FALSE----
knitr::include_graphics("Figure5.png")

## ----"disjoin", message = FALSE------------------------------------------
library("GenomicRanges")
exons <- GRanges("seq", IRanges(start = c(1, 1, 13), end = c(5, 8, 15)))
exons
disjoin(exons)

## ----Figure6, out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Disjoint exons. Windows of distinct exonic sequence for the example gene. Disjoint exons 1 and 2 form exon 1.", echo = FALSE----
knitr::include_graphics("Figure6.png")

## ----Figure7, out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Base-pair coverage counting for exonic base-pairs. At each exonic base-pair we compute the number of reads overlapping that given base-pair. The first base (orange arrow) has 3 reads overlapping that base-pair. Base-pair 11 has a coverage of 3 but does not overlap known exonic sequence, so that information is not used for the gene and exon count matrices (grey arrow). If a read partially overlaps exonic sequence, only the portion that overlaps is used in the computation (see right most read).", echo = FALSE----
knitr::include_graphics("Figure7.png")

## ----Figure8, out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Exon and gene coverage counts. The coverage counts for each disjoint exon are the sum of the base-pair coverage. The gene coverage count is the sum of the disjoint exons coverage counts.", echo = FALSE----
knitr::include_graphics("Figure8.png")

## ----"coverage"----------------------------------------------------------
## Take the example and translate it to R code
library("GenomicRanges")
reads <- GRanges("seq", IRanges(
    start = rep(
        c(1, 2, 3, 4, 5, 7, 8, 9, 10, 13, 14), 
        c(3, 1, 2, 1, 2, 1, 2, 1, 2, 4, 1)
    ), width = rep(
        c(1, 3, 2, 3, 1, 2, 1, 3, 2, 3, 2, 1, 3),
        c(1, 4, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 1)
    )
))
## Get the base-level genome coverage curve
cov <- as.integer(coverage(reads)$seq)

## AUC
sum(cov)

## ----"coverage-reproduce", eval = FALSE----------------------------------
## ## Code for reproducing the bottom portion of Figure 8.
## pdf("base_pair_coverage.pdf", width = 20)
## par(mar = c(5, 6, 4, 2) + 0.1)
## plot(cov, type = "o", col = "violetred1", lwd = 10, ylim = c(0, 5),
##     xlab = "Genome", ylab = "Coverage", cex.axis = 2, cex.lab = 3,
##     bty = "n")
## polygon(c(1, seq_len(length(cov)), length(cov)), c(0, cov, 0),
##     border = NA, density = -1, col = "light blue")
## points(seq_len(length(cov)), cov, col = "violetred1", type = "o",
##     lwd = 10)
## dev.off()

## ----Figure9, out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Area under coverage (AUC). The area under coverage is the sum of the base-pair coverage for all positions in the genome regardless of the annotation. It is the area under the base-level coverage curve shown as the light blue area under the pink curve.", echo = FALSE----
knitr::include_graphics("Figure9.png")

## ----"example_scaled", message = FALSE-----------------------------------
## Check that the number of reads is less than or equal to 40 million
## after scaling.
library("recount")
rse_scaled <- scale_counts(rse_gene_SRP009615, round = FALSE)
summary(colSums(assays(rse_scaled)$counts)) / 1e6

## ----Figure10, out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Exon-exon junctions go beyond the annotation. Reads spanning exon-exon junctions are highlighted and compared against the annotation. Three of them match the annotated junctions, but one (blue and orange read) spans an unannotated exon-exon junction with the left end matching the annotation and the right end hinting at a possible new isoform for this gene (blue and orange isoform).", echo = FALSE----
knitr::include_graphics("Figure10.png")

## ----Figure11, out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Intron retention events. Some reads might align with known intronic segments of the genome and provide information for exploring intron retention events (pink read). Some might support an intron retention event or a new isoform when coupled with exon-exon junction data (orange read).", echo = FALSE----
knitr::include_graphics("Figure11.png")

## ----Figure12, out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Exon boundaries. Reads that go beyond the known exon boundaries can inform us of whether the annotated boundaries are correct or if there was a run-off transcription event.", echo = FALSE----
knitr::include_graphics("Figure12.png")

## ----"download gene"-----------------------------------------------------
library("recount")

## Find the project ID by searching abstracts of studies
abstract_search("human brain development by age")

## Download the data if it is not there
if(!file.exists(file.path("SRP045638", "rse_gene.Rdata"))) {
    download_study("SRP045638", type = "rse-gene")
}

## Check that the file was downloaded
file.exists(file.path("SRP045638", "rse_gene.Rdata"))

## Load the data
load(file.path("SRP045638", "rse_gene.Rdata"))

## ----colData-------------------------------------------------------------
## One row per sample, one column per phenotype variable
dim(colData(rse_gene))

## Mostly technical variables are included
colnames(colData(rse_gene))

## ----"explore colData"---------------------------------------------------
## Input reads: number reported by SRA might be larger than number
## of reads Rail-RNA downloaded
colData(rse_gene)[,
    c("read_count_as_reported_by_sra", "reads_downloaded")]
summary(
    colData(rse_gene)$proportion_of_reads_reported_by_sra_downloaded
)

## AUC information used by scale_counts() by default
head(colData(rse_gene)$auc)

## Alternatively, scale_scounts() can use the number of mapped reads
## and other information
colData(rse_gene)[, c("mapped_read_count", "paired_end",
    "avg_read_length")]

## ----sharq---------------------------------------------------------------
## SHARQ tissue predictions: not present for all studies
head(colData(rse_gene)$sharq_beta_tissue)
head(colData(rse_gene_SRP009615)$sharq_beta_tissue)

## ----characteristics-----------------------------------------------------
## GEO information was absent for the SRP045638 data set
colData(rse_gene)[, c("geo_accession", "title", "characteristics")]

## GEO information for the SRP009615 data set
head(colData(rse_gene_SRP009615)$geo_accession)
head(colData(rse_gene_SRP009615)$title, 2)
head(colData(rse_gene_SRP009615)$characteristics, 2)

## Similar but not exactly the same wording used for two different samples
colData(rse_gene_SRP009615)$characteristics[[1]]
colData(rse_gene_SRP009615)$characteristics[[11]]

## Extract the target information
target <- sapply(colData(rse_gene_SRP009615)$characteristics, "[", 2)
target

## Build a useful factor vector, set the reference level and append the result 
## to the colData() slot
target_factor <- sapply(strsplit(target, "targeting "), "[", 2)
target_factor[is.na(target_factor)] <- "none"
target_factor <- factor(target_factor)
target_factor <- relevel(target_factor, "none")
target_factor
colData(rse_gene_SRP009615)$target_factor <- target_factor

## ----"add_predictions"---------------------------------------------------
## Before adding predictions
dim(colData(rse_gene))

## Add the predictions
rse_gene <- add_predictions(rse_gene)

## After adding the predictions
dim(colData(rse_gene))

## Explore the variables
colData(rse_gene)[, 22:ncol(colData(rse_gene))]

## ----"sra_run_table"-----------------------------------------------------
## Save the information from 
## https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP045638
## to a table. We saved the file as SRP045638/SraRunTable.txt.
file.exists(file.path("SRP045638", "SraRunTable.txt"))

## Read the table
sra <- read.table(file.path("SRP045638", "SraRunTable.txt"),
    header = TRUE, sep = "\t")

## Explore it
head(sra)

## We will remove some trailing '_s' from the variable names
colnames(sra) <- gsub("_s$", "", colnames(sra))

## Choose some variables we want to add
sra_vars <- c("sex", "race", "RIN", "age", "isolate", "disease",
    "tissue")

## Re-organize the SRA table based on the SRA Run IDs we have
sra <- sra[match(colData(rse_gene)$run, sra$Run), ]

## Double check the order
identical(colData(rse_gene)$run, as.character(sra$Run))

## Append the variables of interest
colData(rse_gene) <- cbind(colData(rse_gene), sra[, sra_vars])

## Final dimensions
dim(colData(rse_gene))

## Explore result
colData(rse_gene)[, 34:ncol(colData(rse_gene))]

## ----"sex preds"---------------------------------------------------------
table("Predicted" = colData(rse_gene)$predicted_sex,
    "Observed" = colData(rse_gene)$sex)

## ----"age_groups"--------------------------------------------------------
## Create the original 6 age groups
age_bins <- cut( colData(rse_gene)$age, c(-1, 0, 1, 10, 20, 50, Inf),
    include.lowest=TRUE )
levels( age_bins ) <- c("prenatal", "infant", "child", "teen", "adult",
    "late life")
colData(rse_gene)$age_group <- age_bins

## ----"prenatal_factor"---------------------------------------------------
## Create prenatal factor
colData(rse_gene)$prenatal <- factor(
    ifelse(colData(rse_gene)$age_group == "prenatal", "prenatal",
        "postnatal"),
    levels = c("prenatal", "postnatal"))

## ----"scale_counts"------------------------------------------------------
## Scale counts
rse_gene_scaled <- scale_counts(rse_gene)

## To highlight that we scaled the counts
rm(rse_gene)

## ----"filter_low"--------------------------------------------------------
## Extract counts and filter out lowly expressed geens
counts <- assays(rse_gene_scaled)$counts
filter <- rowMeans(counts) > 0.5

## ----"limmade1", out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Multi-dimensional scaling plot of the gene-level data by age group."----
library("limma")
library("edgeR")

## Build DGEList object
dge <- DGEList(counts = counts[filter, ])

## Calculate normalization factors
dge <- calcNormFactors(dge)

## Explore the data
plotMDS(dge, labels = substr(colData(rse_gene_scaled)$prenatal, 1, 2) )

## ----"limmade2", out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Multi-dimensional scaling plot of the gene-level data by sex."----
plotMDS(dge, labels = substr(colData(rse_gene_scaled)$sex, 1, 1) )
tapply(colData(rse_gene_scaled)$RIN, colData(rse_gene_scaled)$prenatal,
    summary)

## Specify our design matrix
design <- with(colData(rse_gene_scaled),
    model.matrix(~ sex + RIN + prenatal))

## ----"limmade3", out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "voom mean-variance plot of the gene-level data."----
## Run voom
v <- voom(dge, design, plot = TRUE)

## Run remaining parts of the DE analysis
fit <- lmFit(v, design)
fit <- eBayes(fit)

## ----"limmaplots1", out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "MA plot of the gene-level data. Testing for prenatal and postnatal DE adjusting for sex and RIN."----
## Visually explore DE results
limma::plotMA(fit, coef = 4)

## ----"limmaplots2", out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Volcano plot of the gene-level data. Testing for prenatal and postnatal DE adjusting for sex and RIN."----
limma::volcanoplot(fit, coef = 4)

## ----"report_setup"------------------------------------------------------
## Extract data from limma-voom results
top <- topTable(fit, number = Inf, sort.by = "none",
    coef = "prenatalpostnatal")

## Build a DESeqDataSet with the count data and model we used
library("DESeq2")
dds <- DESeqDataSet(rse_gene_scaled[filter, ], ~ sex + RIN + prenatal)

## Add gene names keeping only the Ensembl part of the Gencode IDs
rownames(dds) <- gsub("\\..*", "", rownames(dds))

## Build a DESeqResults object with the relevant information
## Note that we are transforming the baseMean so it will look ok
## with DESeq2's plotting functions.
limma_res <- DESeqResults(DataFrame(pvalue = top[, "P.Value"], 
    log2FoldChange = top[, "logFC"], 
    baseMean = exp(top[, "AveExpr"]),
    padj = top[, "adj.P.Val"]))
rownames(limma_res) <- rownames(dds)

## Specify FDR cutoff to use
metadata(limma_res)[["alpha"]] <- 0.001

## Add gene symbols so they will be displayed in the report
limma_res$symbol <- rowRanges(rse_gene_scaled)$symbol[filter]

## Some extra information used by the report function
mcols(dds) <- limma_res
mcols(mcols(dds)) <- DataFrame(type = "results",
    description = "manual incomplete conversion from limma-voom to DESeq2")

## ----"create_report", cache = !on.bioc, eval = !on.bioc, message = FALSE, warning = FALSE, results = "hide"----
## library("regionReport")
## ## This takes about 20 minutes to run
## report <- DESeq2Report(dds,
##     project = "SRP045638 gene results with limma-voom",
##     output = "gene_report", outdir = "SRP045638",
##     intgroup = c("prenatal", "sex"), res = limma_res,
##     software = "limma")

## ----'browse_report', eval = FALSE---------------------------------------
## browseURL(file.path("SRP045638", "gene_report.html"))

## ----"goanalysis", out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Biological processes enriched in the DE genes."----
library("clusterProfiler")
library("org.Hs.eg.db")

## Remember that limma_res had ENSEMBL IDs for the genes
head(rownames(limma_res))

## Perform enrichment analysis for Biological Process (BP)
## Note that the argument is keytype instead of keyType in Bioconductor 3.5
enrich_go <- enrichGO(
    gene = rownames(limma_res)[limma_res$padj < 0.001],
    OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP",
    pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05,
    universe = rownames(limma_res))

## Visualize enrichment results
dotplot(enrich_go, font.size = 7)

## ----"exondeanalysis1", out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "voom mean-variance plot of the exon-level data."----
## Download the data if it is not there
if(!file.exists(file.path("SRP045638", "rse_exon.Rdata"))) {
    download_study("SRP045638", type = "rse-exon")
}

## Load the data
load(file.path("SRP045638", "rse_exon.Rdata"))

## Scale and add the metadata (it is in the same order)
identical(colData(rse_exon)$run, colData(rse_gene_scaled)$run)
colData(rse_exon) <- colData(rse_gene_scaled)
rse_exon_scaled <- scale_counts(rse_exon)
## To highlight that we scaled the counts
rm(rse_exon)

## Filter lowly expressed exons
filter_exon <- rowMeans(assays(rse_exon_scaled)$counts) > 0.5
round(table(filter_exon) / length(filter_exon) * 100, 2)

## Build DGEList object
dge_exon <- DGEList(
    counts = assays(rse_exon_scaled)$counts[filter_exon, ])

## Calculate normalization factors
dge_exon <- calcNormFactors(dge_exon)

## Run voom
v_exon <- voom(dge_exon, design, plot = TRUE)

## ----"exondeanalysis2", out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Volcano plot of the exon-level data. Testing for prenatal and postnatal DE adjusting for sex and RIN."----
## Run remaining parts of the DE analysis
fit_exon <- lmFit(v_exon, design)
fit_exon <- eBayes(fit_exon)

## Visualize inspect results
limma::volcanoplot(fit_exon, coef = 4)

## Get p-values and other statistics
top_exon <- topTable(fit_exon, number = Inf, sort.by = "none",
    coef = "prenatalpostnatal")
table(top_exon$adj.P.Val < 0.001)

## ----"geneexon", out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Venn diagram of the overlap between DE genes and genes with at least one exon DE."----
## Get the gene IDs for genes that are DE at the gene-level or that have at
## least one exon with DE signal.
genes_w_de_exon <- unique(
    rownames(rse_exon_scaled)[top_exon$adj.P.Val < 0.001])
genes_de <- rownames(rse_gene_scaled)[
    which(filter)[top$adj.P.Val < 0.001]]

## Make a venn diagram
library("gplots")
vinfo <- venn(list("genes" = genes_de, "exons" = genes_w_de_exon),
    names = c("genes", "exons"), show.plot = FALSE) 
plot(vinfo) +
    title("Genes/exons with DE signal")

## ----"geneexonmatch", out.width="100%", fig.align="center", fig.cap = "Log fold change (FC) for DE genes compared against the most extreme exon log FC among exons that are DE for the given gene."----
## Keep only the DE exons that are from a gene that is also DE
top_exon_de <- top_exon[top_exon$adj.P.Val < 0.001 & 
    top_exon$ID %in% attr(vinfo, "intersections")[["genes:exons"]], ]
    
## Find the fold change that is the most extreme among the DE exons of a gene
exon_max_fc <- tapply(top_exon_de$logFC, top_exon_de$ID, function(x) { 
    x[which.max(abs(x))] })

## Keep only the DE genes that match the previous selection
top_gene_de <- top[match(names(exon_max_fc), rownames(top)), ]

## Make the plot
plot(top_gene_de$logFC, exon_max_fc, pch = 20, col = adjustcolor("black", 1/5),
    ylab = "Most extreme exon log FC",
    xlab = "Gene log FC",
    main = "DE genes with at least one DE exon")
abline(a = 0, b = 1, col = "red")
abline(h = 0, col = "grey80")
abline(v = 0, col = "grey80")

## ----"identify regions", eval = .Platform$OS.type != "windows"-----------
## Define expressed regions for study SRP045638, only for chromosome 21
regions <- expressed_regions("SRP045638", "chr21", cutoff = 5L,
    maxClusterGap = 3000L)
    
## Explore the resulting expressed regions
regions
summary(width(regions))
table(width(regions) >= 100)

## Keep only the ones that are at least 100 bp long
regions <- regions[width(regions) >= 100]
length(regions)

## ----"build_rse_ER", eval = .Platform$OS.type != "windows"---------------
## Compute coverage matrix for study SRP045638, only for chromosome 21
## Takes about 4 minutes
rse_er <- coverage_matrix("SRP045638", "chr21", regions,
    chunksize = 2000, verboseLoad = FALSE, scale = FALSE)

## Use the expanded metadata we built for the gene model
colData(rse_er) <- colData(rse_gene_scaled)

## Scale the coverage matrix
rse_er_scaled <- scale_counts(rse_er)

## To highlight that we scaled the counts
rm(rse_er)

## ----"erdeanalysis1", eval = .Platform$OS.type != "windows", out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Multi-dimensional scaling plot of the expressed regions level data by age group."----
## Build DGEList object
dge_er <- DGEList(counts = assays(rse_er_scaled)$counts)

## Calculate normalization factors
dge_er <- calcNormFactors(dge_er)

## Explore the data
plotMDS(dge_er, labels = substr(colData(rse_er_scaled)$prenatal, 1, 2) )

## ----"erdeanalysis2", eval = .Platform$OS.type != "windows", out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Multi-dimensional scaling plot of the expressed regions level data by sex."----
plotMDS(dge_er, labels = substr(colData(rse_er_scaled)$sex, 1, 1) )

## ----"erdeanalysis3", eval = .Platform$OS.type != "windows", out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "voom mean-variance plot of the expressed regions level data."----
## Run voom
v_er <- voom(dge_er, design, plot = TRUE)

## Run remaining parts of the DE analysis
fit_er <- lmFit(v_er, design)
fit_er <- eBayes(fit_er)

## ----"erdeanalysis4", eval = .Platform$OS.type != "windows", out.width=ifelse(on.bioc, "100%", "70%"), fig.align="center", fig.cap = "Volcano plot of the expressed regions level data. Testing for prenatal and postnatal DE adjusting for sex and RIN."----
## Visually explore the results
limma::volcanoplot(fit_er, coef = 4)

## Number of DERs
top_er <- topTable(fit_er, number = Inf, sort.by = "none",
    coef = "prenatalpostnatal")
table(top_er$adj.P.Val < 0.001)

## ----"sort_qvalue", eval = .Platform$OS.type != "windows"----------------
## Sort regions by q-value
regions_by_padj <- regions[order(top_er$adj.P.Val, decreasing = FALSE)]

## Look at the top 10
regions_by_padj[1:10]
width(regions_by_padj[1:10])

## ----"find_bws", eval = .Platform$OS.type != "windows"-------------------
## Construct the list of bigWig URLs
## They have the following form:
## http://duffel.rail.bio/recount/
## project id
## /bw/
## sample run id
## .bw
bws <- paste0("http://duffel.rail.bio/recount/SRP045638/bw/",
    colData(rse_er_scaled)$bigwig_file)

## Note that they are also present in the recount_url data.frame
bws <- recount_url$url[match(colData(rse_er_scaled)$bigwig_file,
    recount_url$file_name)]

## Use the sample run IDs as the sample names
names(bws) <- colData(rse_er_scaled)$run

## ----"add_padding", eval = .Platform$OS.type != "windows"----------------
## Add 100 bp padding on each side
regions_resized <- resize(regions_by_padj[1:10],
    width(regions_by_padj[1:10]) + 200, fix = "center")

## ----"regionCov", eval = .Platform$OS.type != "windows"------------------
## Get the bp coverage data for the plots
library("derfinder")
regionCov <- getRegionCoverage(regions = regions_resized, files = bws,
    targetSize = 40 * 1e6 * 100,
    totalMapped = colData(rse_er_scaled)$auc,
    verbose = FALSE)

## ----"gencode_txdb", eval = .Platform$OS.type != "windows"---------------
## Import the Gencode v25 hg38 gene annotation
library("rtracklayer")
gencode_v25_hg38 <- import(paste0(
    "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/",
    "gencode.v25.annotation.gtf.gz"))
            
## Keep only the chr21 info
gencode_v25_hg38 <- keepSeqlevels(gencode_v25_hg38, "chr21",
    pruning.mode = "coarse")

## Get the chromosome information for hg38
library("GenomicFeatures")
chrInfo <- getChromInfoFromUCSC("hg38")
chrInfo$chrom <- as.character(chrInfo$chrom)
chrInfo <- chrInfo[chrInfo$chrom %in% seqlevels(regions), ]
chrInfo$isCircular <- FALSE

## Assign the chromosome information to the object we will use to
## create the txdb object
si <- with(chrInfo, Seqinfo(as.character(chrom), length, isCircular,
    genome = "hg38"))
seqinfo(gencode_v25_hg38) <- si

## Switch from Gencode gene IDs to Ensembl gene IDs
gencode_v25_hg38$gene_id <- gsub("\\..*", "", gencode_v25_hg38$gene_id)

## Create the TxDb object
gencode_v25_hg38_txdb <- makeTxDbFromGRanges(gencode_v25_hg38)

## Explore the TxDb object
gencode_v25_hg38_txdb

## ----"bump_ann", eval = .Platform$OS.type != "windows"-------------------
library("bumphunter")
## Annotate all transcripts for gencode v25 based on the TxDb object
## we built previously.
ann_gencode_v25_hg38 <- annotateTranscripts(gencode_v25_hg38_txdb,
    annotationPackage = "org.Hs.eg.db",
    mappingInfo = list("column" = "ENTREZID", "keytype" = "ENSEMBL",
    "multiVals" = "first"))
    
## Annotate the regions of interest
## Note that we are using the original regions, not the resized ones
nearest_ann <- matchGenes(regions_by_padj[1:10], ann_gencode_v25_hg38)

## ----"make_gs", eval = .Platform$OS.type != "windows"--------------------
## Create the genomic state object using the gencode TxDb object
gs_gencode_v25_hg38 <- makeGenomicState(gencode_v25_hg38_txdb,
    chrs = seqlevels(regions))
    
## Annotate the original regions
regions_ann <- annotateRegions(regions_resized,
    gs_gencode_v25_hg38$fullGenome)

## ----"regionplots", eval = .Platform$OS.type != "windows", out.width=ifelse(on.bioc, "100%", "75%"), fig.align="center", fig.cap = 'Base-pair resolution plot of differentially expressed region 2.'----
library("derfinderPlot")
pdf('region_plots.pdf')
plotRegionCoverage(regions = regions_resized, regionCoverage = regionCov, 
   groupInfo = colData(rse_er_scaled)$prenatal,
   nearestAnnotation = nearest_ann, 
   annotatedRegions = regions_ann,
   txdb = gencode_v25_hg38_txdb,
   scalefac = 1, ylab = "Coverage (RP40M, 100bp)",
   ask = FALSE, verbose = FALSE)
dev.off()

## Visualize DER #2
plotRegionCoverage(regions = regions_resized, regionCoverage = regionCov, 
   groupInfo = colData(rse_er_scaled)$prenatal,
   nearestAnnotation = nearest_ann, 
   annotatedRegions = regions_ann,
   txdb = gencode_v25_hg38_txdb,
   scalefac = 1, ylab = "Coverage (RP40M, 100bp)",
   ask = FALSE, verbose = FALSE, whichRegions = 2)

## ----sessionInfo---------------------------------------------------------
## Final list of files created
dir("SRP045638")

## Pandoc information
library("rmarkdown")
pandoc_version()

## Time for reproducing this workflow, in minutes
round(proc.time()[3] / 60, 1)

options(width = 100)
library("devtools")
session_info()

