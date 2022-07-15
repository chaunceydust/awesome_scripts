## ----style, echo=FALSE, results='hide', message=FALSE--------------------
library(BiocStyle)
library(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)
opts_chunk$set(fig.width=7, fig.height=7)
opts_chunk$set(dpi=300, dev="png", dev.args=list(pointsize=15))
options(bitmapType="cairo", width=100)

# Setting single-core unless explicitly specified otherwise.
library(BiocParallel)
register(SerialParam())

# Deciding whether we want to re-download everything or not.
on.bioc <- TRUE

## ---- message=FALSE, echo=FALSE, results='hide'--------------------------
library(Rtsne)
library(mvoutlier)
library(destiny)

## ---- eval=on.bioc, echo=FALSE, results='hide'---------------------------
all.urls <- c("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE61533&format=file&file=GSE61533%5FHTSEQ%5Fcount%5Fresults%2Exls%2Egz", 
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE29087&format=file&file=GSE29087%5FL139%5Fexpression%5Ftab%2Etxt%2Egz",
"https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt",
"https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mito_17-Aug-2014.txt",
"https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_spikes_17-Aug-2014.txt",
"http://www.ebi.ac.uk/teichmann-srv/espresso/static/counttable_es.csv", 
"http://www.nature.com/nbt/journal/v33/n2/extref/nbt.3102-S7.xlsx")
all.basenames <- basename(all.urls)
all.basenames[1] <- "GSE61533_HTSEQ_count_results.xls.gz"
all.basenames[2] <- "GSE29087_L139_expression_tab.txt.gz"
all.modes <- rep("w", length(all.urls))
all.modes[!grepl("(txt|csv)$", all.basenames)] <- "wb"
for (x in seq_along(all.urls)) { 
    download.file(all.urls[x], all.basenames[x], mode=all.modes[x])
}

## ---- results='asis', eval=on.bioc, echo=FALSE---------------------------
cat("***Note:*** *to cite this article, please refer to http://f1000research.com/articles/5-2122/v1 for instructions.*")

## ------------------------------------------------------------------------
library(R.utils)
gunzip("GSE61533_HTSEQ_count_results.xls.gz", remove=FALSE, overwrite=TRUE)
library(readxl)
all.counts <- as.data.frame(read_excel('GSE61533_HTSEQ_count_results.xls', sheet=1))
rownames(all.counts) <- all.counts$ID
all.counts <- as.matrix(all.counts[,-1])

## ------------------------------------------------------------------------
library(SingleCellExperiment)
sce <- SingleCellExperiment(list(counts=all.counts))
dim(sce)

## ------------------------------------------------------------------------
is.spike <- grepl("^ERCC", rownames(sce))
isSpike(sce, "ERCC") <- is.spike
summary(is.spike)

## ------------------------------------------------------------------------
is.mito <- grepl("^mt-", rownames(sce))
summary(is.mito)

## ------------------------------------------------------------------------
library(scater)
sce <- calculateQCMetrics(sce, feature_controls=list(ERCC=is.spike, Mt=is.mito))
head(colnames(colData(sce)))

## ----qcplothsc, fig.height=12, fig.width=12, fig.cap="Histograms of various QC metrics for all cells in the HSC data set. This includes the library sizes, number of expressed genes, and proportion of reads mapped to spike-in transcripts or mitochondrial genes."----
par(mfrow=c(2,2))
hist(sce$total_counts/1e6, xlab="Library sizes (millions)", main="", 
    breaks=20, col="grey80", ylab="Number of cells")
hist(sce$total_features, xlab="Number of expressed genes", main="", 
    breaks=20, col="grey80", ylab="Number of cells")
hist(sce$pct_counts_ERCC, xlab="ERCC proportion (%)", 
    ylab="Number of cells", breaks=20, main="", col="grey80")
hist(sce$pct_counts_Mt, xlab="Mitochondrial proportion (%)", 
    ylab="Number of cells", breaks=20, main="", col="grey80")

## ------------------------------------------------------------------------
libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE)

## ------------------------------------------------------------------------
spike.drop <- isOutlier(sce$pct_counts_ERCC, nmads=3, type="higher")

## ------------------------------------------------------------------------
sce <- sce[,!(libsize.drop | feature.drop | spike.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
    BySpike=sum(spike.drop), Remaining=ncol(sce))

## ----pcaqualplothsc, fig.cap="PCA plot for cells in the HSC dataset, constructed using quality metrics. The first and second components are shown on each axis, along with the percentage of total variance explained by each component. Bars represent the coordinates of the cells on each axis."----
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotPCA(sce, pca_data_input="pdata") + fontsize

## ------------------------------------------------------------------------
set.seed(100)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
library(org.Mm.eg.db)
ensembl <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="SYMBOL", column="ENSEMBL")
library(scran)
assignments <- cyclone(sce, mm.pairs, gene.names=ensembl)

## ----phaseplothsc, message=FALSE, fig.cap="Cell cycle phase scores from applying the pair-based classifier on the HSC dataset, where each point represents a cell."----
plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)

## ------------------------------------------------------------------------
sce$phases <- assignments$phases
table(sce$phases)

## ----topgenehsc, fig.height=9, fig.width=6, fig.cap="Percentage of total counts assigned to the top 50 most highly-abundant features in the HSC dataset. For each feature, each bar represents the percentage assigned to that feature for a single cell, while the circle represents the average across all cells. Bars are coloured by the total number of expressed features in each cell, while circles are coloured according to whether the feature is labelled as a control feature."----
plotQC(sce, type = "highest-expression", n=50) + fontsize

## ----abhisthsc, fig.cap="Histogram of log-average counts for all genes in the HSC dataset."----
ave.counts <- calcAverage(sce)
hist(log10(ave.counts), breaks=100, main="", col="grey80", 
    xlab=expression(Log[10]~"average count"))

## ------------------------------------------------------------------------
demo.keep <- ave.counts >= 1
filtered.sce <- sce[demo.keep,]
summary(demo.keep)

## ----nexprshisthsc, fig.cap="The number of cells expressing each gene in the HSC data set, plotted against the log-average count. Intensity of colour corresponds to the number of genes at any given location."----
num.cells <- nexprs(sce, byrow=TRUE)
smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells", 
    xlab=expression(Log[10]~"average count"))

## ------------------------------------------------------------------------
to.keep <- num.cells > 0
sce <- sce[to.keep,]
summary(to.keep)

## ---- warning=FALSE------------------------------------------------------
sce <- computeSumFactors(sce)
summary(sizeFactors(sce))

## ----normplothsc, fig.cap="Size factors from deconvolution, plotted against library sizes for all cells in the HSC dataset. Axes are shown on a log-scale."----
plot(sizeFactors(sce), sce$total_counts/1e6, log="xy",
    ylab="Library size (millions)", xlab="Size factor")

## ------------------------------------------------------------------------
sce <- computeSpikeFactors(sce, type="ERCC", general.use=FALSE)

## ------------------------------------------------------------------------
sce <- normalize(sce)

## ---- echo=FALSE, results="hide"-----------------------------------------
gc()

## ----explvarplothsc, fig.cap="Density plot of the percentage of variance explained by the (log-transformed) total spike-in counts across all genes in the HSC dataset. For each gene, the percentage of the variance of the normalized log-expression values across cells that is explained by each factor is calculated. Each curve corresponds to one factor and represents the distribution of percentages across all genes."----
plotExplanatoryVariables(sce, variables=c("total_counts_ERCC", 
    "log10_total_counts_ERCC")) + fontsize

## ------------------------------------------------------------------------
var.fit <- trendVar(sce, parametric=TRUE, span=0.2)

## ------------------------------------------------------------------------
var.out <- decomposeVar(sce, var.fit)
head(var.out)

## ----hvgplothsc, fig.cap="Variance of normalized log-expression values for each gene in the HSC dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of the spike-in transcripts (red)."----
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
    ylab="Variance of log-expression")
curve(var.fit$trend(x), col="dodgerblue", lwd=2, add=TRUE)
cur.spike <- isSpike(sce)
points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)

## ----hvgvioplothsc, fig.cap="Violin plots of normalized log-expression values for the top 10 genes with the largest biological components in the HSC dataset. Each point represents the log-expression value in a single cell."----
chosen.genes <- order(var.out$bio, decreasing=TRUE)[1:10]
plotExpression(sce, features=rownames(var.out)[chosen.genes]) + fontsize

## ------------------------------------------------------------------------
var.fit.nospike <- trendVar(sce, parametric=TRUE, use.spikes=FALSE, span=0.2)
var.out.nospike <- decomposeVar(sce, var.fit.nospike)

## ----hvgplothsc2, fig.cap="Variance of normalized log-expression values for each gene in the HSC dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of the endogenous genes (black), with spike-in transcripts shown in red."----
plot(var.out.nospike$mean, var.out.nospike$total, pch=16, cex=0.6, 
    xlab="Mean log-expression", ylab="Variance of log-expression")
curve(var.fit.nospike$trend(x), col="dodgerblue", lwd=2, add=TRUE)
points(var.out.nospike$mean[cur.spike], var.out.nospike$total[cur.spike], col="red", pch=16)

## ------------------------------------------------------------------------
sce <- denoisePCA(sce, technical=var.fit$trend) 
dim(reducedDim(sce, "PCA")) 

## ------------------------------------------------------------------------
sce2 <- denoisePCA(sce, technical=var.fit$trend, value="lowrank") 
assayNames(sce2)

## ---- echo=FALSE, results="hide"-----------------------------------------
rm(sce2)
gc()

## ----pcaplothsc, fig.cap="Pairwise PCA plots of the first three PCs in the HSC data set, constructed from normalized log-expression values of genes with positive biological components. Each point represents a cell, coloured according to its total number of expressed features. Bars represent the coordinates of the cells on each axis.", fig.height=9, fig.width=9----
plotReducedDim(sce, use_dimred="PCA", ncomponents=3, colour_by="total_features") + fontsize

## ----tsneplothsc, fig.cap="_t_-SNE plots constructed from the denoised PCs in the HSC data set, using a range of perplexity values. Each point represents a cell, coloured according to its total number of expressed features. Bars represent the coordinates of the cells on each axis.", fig.width=12, fig.height=6----
out5 <- plotTSNE(sce, use_dimred="PCA", perplexity=5, colour_by="total_features", 
    rand_seed=100) + fontsize + ggtitle("Perplexity = 5")
out10 <- plotTSNE(sce, use_dimred="PCA", perplexity=10, colour_by="total_features",
    rand_seed=100) + fontsize + ggtitle("Perplexity = 10")
out20 <- plotTSNE(sce, use_dimred="PCA", perplexity=20, colour_by="total_features",
    rand_seed=100) + fontsize + ggtitle("Perplexity = 20")
multiplot(out5, out10, out20, cols=3)

## ------------------------------------------------------------------------
pc1 <- reducedDim(sce, "PCA")[,1]
design <- model.matrix(~pc1)
library(limma)
fit <- lmFit(logcounts(sce), design)
fit <- eBayes(fit, trend=TRUE, robust=TRUE)
topTable(fit)

## ----heatmaphsc, fig.height=10, fig.width=6, fig.cap="Heatmap of the top 50 DE genes along the first PC in the HSC data set. The colour for each cell (column) represents the log-fold change from the average log-expression for each gene (row), bounded to [-2, 2] for visualization purposes. Cells are ordered by their location on the first PC."----
de.genes <- rownames(topTable(fit, coef=2, n=50))
heat.vals <- logcounts(sce)[de.genes,]
heat.vals <- heat.vals - rowMeans(heat.vals)
heat.vals[heat.vals > 2] <- 2
heat.vals[heat.vals < -2] <- -2
library(pheatmap)
pheatmap(heat.vals[,order(pc1)], cluster_cols=FALSE)

## ------------------------------------------------------------------------
saveRDS(file="hsc_data.rds", sce)

## ---- echo=FALSE, results='hide'-----------------------------------------
rm(sce, all.counts)
gc()

## ------------------------------------------------------------------------
readFormat <- function(infile) { 
    # First column is empty.
    metadata <- read.delim(infile, stringsAsFactors=FALSE, header=FALSE, nrow=10)[,-1] 
    rownames(metadata) <- metadata[,1]
    metadata <- metadata[,-1]
    metadata <- as.data.frame(t(metadata))
    # First column after row names is some useless filler.
    counts <- read.delim(infile, stringsAsFactors=FALSE, header=FALSE, row.names=1, skip=11)[,-1] 
    counts <- as.matrix(counts)
    return(list(metadata=metadata, counts=counts))
}

## ------------------------------------------------------------------------
endo.data <- readFormat("expression_mRNA_17-Aug-2014.txt")
spike.data <- readFormat("expression_spikes_17-Aug-2014.txt")
mito.data <- readFormat("expression_mito_17-Aug-2014.txt")

## ------------------------------------------------------------------------
m <- match(endo.data$metadata$cell_id, mito.data$metadata$cell_id)
mito.data$metadata <- mito.data$metadata[m,]
mito.data$counts <- mito.data$counts[,m]

## ---- echo=FALSE---------------------------------------------------------
stopifnot(identical(endo.data$metadata$cell_id, spike.data$metadata$cell_id)) # should be the same.
stopifnot(all(endo.data$metadata$cell_id==mito.data$metadata$cell_id)) # should now be the same.

## ------------------------------------------------------------------------
raw.names <- sub("_loc[0-9]+$", "", rownames(endo.data$counts))
new.counts <- rowsum(endo.data$counts, group=raw.names, reorder=FALSE)
endo.data$counts <- new.counts

## ------------------------------------------------------------------------
all.counts <- rbind(endo.data$counts, mito.data$counts, spike.data$counts)
sce <- SingleCellExperiment(list(counts=all.counts), colData=endo.data$metadata)
dim(sce)

## ------------------------------------------------------------------------
nrows <- c(nrow(endo.data$counts), nrow(mito.data$counts), nrow(spike.data$counts))
is.spike <- rep(c(FALSE, FALSE, TRUE), nrows)
is.mito <- rep(c(FALSE, TRUE, FALSE), nrows)
isSpike(sce, "Spike") <- is.spike
sce

## ---- echo=FALSE, results='hide'-----------------------------------------
# Save some memory.
rm(mito.data, endo.data, spike.data, new.counts)
gc()

## ------------------------------------------------------------------------
sce <- calculateQCMetrics(sce, feature_controls=list(Spike=is.spike, Mt=is.mito)) 

## ----libplotbrain, fig.width=12, fig.height=12, fig.cap="Histograms of QC metrics including the library sizes, number of expressed genes and proportion of UMIs assigned to spike-in transcripts or mitochondrial genes for all cells in the brain dataset."----
par(mfrow=c(2,2))
hist(sce$total_counts/1e3, xlab="Library sizes (thousands)", main="", 
    breaks=20, col="grey80", ylab="Number of cells")
hist(sce$total_features, xlab="Number of expressed genes", main="", 
    breaks=20, col="grey80", ylab="Number of cells")
hist(sce$pct_counts_Spike, xlab="ERCC proportion (%)",
    ylab="Number of cells", breaks=20, main="", col="grey80")
hist(sce$pct_counts_Mt, xlab="Mitochondrial proportion (%)", 
    ylab="Number of cells", breaks=20, main="", col="grey80")

## ------------------------------------------------------------------------
libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE)
spike.drop <- isOutlier(sce$pct_counts_Spike, nmads=3, type="higher")

## ------------------------------------------------------------------------
sce <- sce[,!(libsize.drop | feature.drop | spike.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), 
    BySpike=sum(spike.drop), Remaining=ncol(sce))

## ----echo=FALSE, results='hide'------------------------------------------
gc()

## ---- echo=FALSE, results='hide', message=FALSE--------------------------
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
library(org.Mm.eg.db)

## ----phaseplotbrain, message=FALSE, fig.cap="Cell cycle phase scores from applying the pair-based classifier on the brain dataset, where each point represents a cell."----
ensembl <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="SYMBOL", column="ENSEMBL")
assignments <- cyclone(sce, mm.pairs, gene.names=ensembl)
plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)

## ----echo=FALSE, results='hide'------------------------------------------
gc()

## ----topgenebrain, fig.height=9, fig.width=6, fig.cap="Percentage of total counts assigned to the top 50 most highly-abundant features in the brain dataset. For each feature, each bar represents the percentage assigned to that feature for a single cell, while the circle represents the average across all cells. Bars are coloured by the total number of expressed features in each cell, while circles are coloured according to whether the feature is labelled as a control feature."----
plotQC(sce, type = "highest-expression", n=50) + fontsize

## ----abhistbrain, fig.cap="Histogram of log-average counts for all genes in the brain dataset. The filter threshold is represented by the blue line."----
ave.counts <- calcAverage(sce)
hist(log10(ave.counts), breaks=100, main="", col="grey",
    xlab=expression(Log[10]~"average count"))
abline(v=log10(0.1), col="blue", lwd=2, lty=2)

## ------------------------------------------------------------------------
rowData(sce)$ave.count <- ave.counts
to.keep <- ave.counts > 0
sce <- sce[to.keep,]
summary(to.keep)

## ----echo=FALSE, results='hide'------------------------------------------
gc()

## ------------------------------------------------------------------------
high.ave <- rowData(sce)$ave.count >= 0.1
clusters <- quickCluster(sce, subset.row=high.ave, method="igraph")
sce <- computeSumFactors(sce, cluster=clusters, 
    subset.row=high.ave, min.mean=NULL)
summary(sizeFactors(sce))

## ----echo=FALSE, results='hide'------------------------------------------
gc()

## ----normplotbrain, fig.cap="Size factors from deconvolution, plotted against library sizes for all cells in the brain dataset. Axes are shown on a log-scale."----
plot(sizeFactors(sce), sce$total_counts/1e3, log="xy",
    ylab="Library size (thousands)", xlab="Size factor")

## ------------------------------------------------------------------------
sce <- computeSpikeFactors(sce, type="Spike", general.use=FALSE)

## ------------------------------------------------------------------------
sce <- normalize(sce)

## ----echo=FALSE, results='hide'------------------------------------------
gc()

## ---- echo=FALSE, results='hide', message=FALSE--------------------------
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))

## ----explvarplotbrain, fig.cap="Density plot of the percentage of variance explained by each factor across all genes in the brain dataset. For each gene, the percentage of the variance of the normalized log-expression values that is explained by the (log-transformed) total spike-in counts, the sex or age of the mouse, or the tissue of origin is calculated. Each curve corresponds to one factor and represents the distribution of percentages across all genes."----
plotExplanatoryVariables(sce, variables=c("log10_total_counts_Spike", 
    "log10_total_counts_Spike", "sex", "tissue", "age")) + fontsize

## ------------------------------------------------------------------------
design <- model.matrix(~sce$sex)

## ------------------------------------------------------------------------
var.fit <- trendVar(sce, parametric=TRUE, span=0.4, design=design)
var.out <- decomposeVar(sce, var.fit)

## ----hvgplotbrain, fig.cap="Variance of normalized log-expression values against the mean for each gene, calculated across all cells in the brain data set after blocking on the sex effect. The blue line represents the mean-dependent trend in the technical variance of the spike-in transcripts (also highlighted as red points)."----
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
    ylab="Variance of log-expression")
points(var.out$mean[isSpike(sce)], var.out$total[isSpike(sce)], col="red", pch=16)
curve(var.fit$trend(x), col="dodgerblue", add=TRUE, lwd=2)

## ----hvgvioplotbrain, fig.cap="Violin plots of normalized log-expression values for the top 10 HVGs in the brain dataset. For each gene, each point represents the log-expression value for an individual cell."----
chosen.genes <- order(var.out$bio, decreasing=TRUE)[1:10]
plotExpression(sce, rownames(var.out)[chosen.genes], 
    alpha=0.05, jitter="jitter") + fontsize

## ------------------------------------------------------------------------
sce <- denoisePCA(sce, technical=var.fit$trend, design=design, approximate=TRUE)
ncol(reducedDim(sce, "PCA"))

## ---- echo=FALSE, results='hide', message=FALSE--------------------------
gc()

## ------------------------------------------------------------------------
collected <- list()
for (block in levels(sce$sex)) {
    cur.sce <- sce[,sce$sex==block]
    cur.sce <- normalize(cur.sce) 
    var.fit <- trendVar(cur.sce, parametric=TRUE, span=0.4)
    collected[[block]] <- decomposeVar(cur.sce, var.fit)
}
var.out <- do.call(combineVar, collected)

## ---- echo=FALSE, results='hide', message=FALSE--------------------------
rm(cur.sce)
gc()

## ------------------------------------------------------------------------
library(limma)
adj.exprs <- logcounts(sce)
adj.exprs <- removeBatchEffect(adj.exprs, batch=sce$sex)
norm_exprs(sce) <- adj.exprs 

## ---- echo=FALSE, results='hide', message=FALSE--------------------------
rm(adj.exprs)
gc()

## ----tsneplotbrain, fig.cap="_t_-SNE plots constructed from the denoised PCs of the brain dataset. Each point represents a cell and is coloured according to its expression of _Neurod6_ (left) or _Mog_ (right).", fig.width=15, fig.height=6----
tsne1 <- plotTSNE(sce, use_dimred="PCA", colour_by="Neurod6",
    perplexity=10, rand_seed=100) + fontsize
tsne2 <- plotTSNE(sce, use_dimred="PCA", colour_by="Mog",
    perplexity=10, rand_seed=100) + fontsize
multiplot(tsne1, tsne2, cols=2)

## ----pcaplotbrain, fig.cap="PCA plots constructed from the denoised PCs of the brain dataset. Each point represents a cell and is coloured according to its expression of the _Neurod6_ (left) or _Mog_ (right).", fig.width=15, fig.height=6----
pca1 <- plotReducedDim(sce, use_dimred="PCA", colour_by="Neurod6") + fontsize
pca2 <- plotReducedDim(sce, use_dimred="PCA", colour_by="Mog") + fontsize
multiplot(pca1, pca2, cols=2)

## ----echo=FALSE, results='hide'------------------------------------------
rm(tsne1, tsne2, pca1, pca2)
gc()

## ------------------------------------------------------------------------
pcs <- reducedDim(sce, "PCA")
my.dist <- dist(pcs)
my.tree <- hclust(my.dist, method="ward.D2")

## ------------------------------------------------------------------------
library(dynamicTreeCut)
my.clusters <- unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist), verbose=0))

## ----tsneclusterbrain, message=FALSE, fig.width=7, fig.height=10, fig.cap="_t_-SNE plot of the denoised PCs of the brain data set. Each point represents a cell and is coloured according to the cluster identity to which it was assigned."----
sce$cluster <- factor(my.clusters)
plotTSNE(sce, use_dimred="PCA", colour_by="cluster",
    perplexity=10, rand_seed=100) + fontsize

## ----silhouettebrain, message=FALSE, fig.cap="Barplot of silhouette widths for cells in each cluster. Each cluster is assigned a colour and cells with positive widths are coloured according to the colour of its assigned cluster. Any cell with a negative width is coloured according to the colour of the cluster that it is closest to. The average width for all cells in each cluster is shown, along with the average width for all cells in the data set."----
library(cluster)
clust.col <- scater:::.get_palette("tableau10medium") # hidden scater colours
sil <- silhouette(my.clusters, dist = my.dist)
sil.cols <- clust.col[ifelse(sil[,3] > 0, sil[,1], sil[,2])]
sil.cols <- sil.cols[order(-sil[,1], sil[,3])]
plot(sil, main = paste(length(unique(my.clusters)), "clusters"), 
    border=sil.cols, col=sil.cols, do.col.sort=FALSE) 

## ----echo=FALSE, results='hide'------------------------------------------
gc()

## ------------------------------------------------------------------------
snn.gr <- buildSNNGraph(sce, use.dimred="PCA")
gr.clusters <- igraph::cluster_fast_greedy(snn.gr)
table(gr.clusters$membership)

## ----echo=FALSE, results='hide'------------------------------------------
rm(snn.gr, gr.clusters)
gc()

## ------------------------------------------------------------------------
markers <- findMarkers(sce, my.clusters, design=design)

## ---- echo=FALSE, results="hide"-----------------------------------------
old.digits <- options()$digits
options(digits=3)

## ------------------------------------------------------------------------
marker.set <- markers[["1"]]
head(marker.set, 10)

## ---- echo=FALSE, results="hide"-----------------------------------------
options(digits=old.digits)

## ------------------------------------------------------------------------
write.table(marker.set, file="brain_marker_1.tsv", sep="\t", quote=FALSE, col.names=NA)

## ----heatmapmarkerbrain, fig.cap="Heatmap of mean-centred normalized and corrected log-expression values for the top set of markers for cluster 1 in the brain dataset. Column colours represent the cluster to which each cell is assigned, as indicated by the legend."----
top.markers <- marker.set$Gene[marker.set$Top <= 10]
top.exprs <- norm_exprs(sce)[top.markers,,drop=FALSE]
heat.vals <- top.exprs - rowMeans(top.exprs)
pheatmap(heat.vals, cluster_cols=my.tree,
    annotation_col=data.frame(Cluster=factor(my.clusters), row.names=colnames(sce)),
    annotation_colors=list(Cluster=setNames(clust.col, seq_along(unique(my.clusters)))))

## ------------------------------------------------------------------------
library(edgeR)
y <- convertTo(sce, type="edgeR")

## ------------------------------------------------------------------------
saveRDS(file="brain_data.rds", sce)

## ---- echo=FALSE, results='hide'-----------------------------------------
gc()

## ------------------------------------------------------------------------
counts <- read.table("GSE29087_L139_expression_tab.txt.gz", colClasses=c(list("character", 
    NULL, NULL, NULL, NULL, NULL, NULL), rep("integer", 96)), skip=6, sep='\t', row.names=1)
is.spike <- grep("SPIKE", rownames(counts)) 
sce <- SingleCellExperiment(list(counts=as.matrix(counts)))
isSpike(sce, "spike") <- is.spike
sce$grouping <- rep(c("mESC", "MEF", "Neg"), c(48, 44, 4))
sce <- sce[,sce$grouping!="Neg"] # Removing negative control wells.
sce <- calculateQCMetrics(sce, feature_controls=list(spike=is.spike))
sce

## ------------------------------------------------------------------------
sce <- computeSpikeFactors(sce, general.use=TRUE)

## ------------------------------------------------------------------------
sce <- normalize(sce)

## ----normplotspikemef, fig.cap="Size factors from spike-in normalization, plotted against the size factors from deconvolution for all cells in the mESC/MEF dataset. Axes are shown on a log-scale, and cells are coloured according to their identity. Deconvolution size factors were computed with small pool sizes owing to the low number of cells of each type."----
colours <- c(mESC="red", MEF="grey")
deconv.sf <- computeSumFactors(sce, sf.out=TRUE, cluster=sce$grouping)
plot(sizeFactors(sce), deconv.sf, col=colours[sce$grouping], pch=16, log="xy", 
    xlab="Size factor (spike-in)", ylab="Size factor (deconvolution)")
legend("bottomleft", col=colours, legend=names(colours), pch=16)

## ------------------------------------------------------------------------
sce <- readRDS("hsc_data.rds")
var.fit <- trendVar(sce, parametric=TRUE, span=0.2)
var.out <- decomposeVar(sce, var.fit)
hvg.out <- var.out[which(var.out$FDR <= 0.05),]
nrow(hvg.out)

## ------------------------------------------------------------------------
hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),] 
write.table(file="hsc_hvg.tsv", hvg.out, sep="\t", quote=FALSE, col.names=NA)
head(hvg.out)

## ------------------------------------------------------------------------
set.seed(100)
var.cor <- correlatePairs(sce, subset.row=grep("^H2-", rownames(sce)))
head(var.cor)

## ------------------------------------------------------------------------
sig.cor <- var.cor$FDR <= 0.05
summary(sig.cor)

## ------------------------------------------------------------------------
correlatePairs(sce, subset.row=cbind("Fos", "Jun"))

## ----fosjuncorplot, fig.cap="Expression of _Fos_ plotted against the expression of _Jun_ for all cells in the HSC data set."----
plotExpression(sce, features="Fos", x="Jun")

## ------------------------------------------------------------------------
incoming <- as.data.frame(read_excel("nbt.3102-S7.xlsx", sheet=1))
rownames(incoming) <- incoming[,1]
incoming <- incoming[,-1]
incoming <- incoming[,!duplicated(colnames(incoming))] # Remove duplicated genes.
sce <- SingleCellExperiment(list(logcounts=t(incoming)))

## ---- echo=FALSE, results='hide', message=FALSE--------------------------
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
library(org.Mm.eg.db)
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))

## ----phaseplotth2, message=FALSE, fig.cap="Cell cycle phase scores from applying the pair-based classifier on the T~H~2 dataset, where each point represents a cell."----
anno <- select(org.Mm.eg.db, keys=rownames(sce), keytype="SYMBOL", column="ENSEMBL")
ensembl <- anno$ENSEMBL[match(rownames(sce), anno$SYMBOL)]
set.seed(100)
assignments <- cyclone(sce, mm.pairs, gene.names=ensembl, assay.type="logcounts")
plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)

## ------------------------------------------------------------------------
design <- model.matrix(~ G1 + G2M, assignments$score)
fit.block <- trendVar(sce, design=design, parametric=TRUE, use.spikes=NA)
sce.block <- denoisePCA(sce, technical=fit.block$trend, design=design) 

## ----pcaplotth2, fig.width=12, fig.height=6, fig.cap="PCA plots before (left) and after (right) removal of the cell cycle effect in the T~H~2 dataset. Each cell is represented by a point with colour and size determined by the G1 and G2/M scores, respectively."----
sce$G1score <- sce.block$G1score <- assignments$score$G1
sce$G2Mscore <- sce.block$G2Mscore <- assignments$score$G2M

# Without blocking on phase score.
fit <- trendVar(sce, parametric=TRUE, use.spikes=NA) 
sce <- denoisePCA(sce, technical=fit$trend)
out <- plotReducedDim(sce, use_dimred="PCA", ncomponents=2, colour_by="G1score", 
    size_by="G2Mscore") + fontsize + ggtitle("Before removal")

# After blocking on the phase score.
out2 <- plotReducedDim(sce.block, use_dimred="PCA", ncomponents=2, colour_by="G1score", 
    size_by="G2Mscore") + fontsize + ggtitle("After removal")
multiplot(out, out2, cols=2)

## ----diffusionth2, fig.cap="A diffusion map for the T~H~2 dataset, where each cell is coloured by its expression of _Gata3_."----
sce.block <- denoisePCA(sce, technical=fit.block$trend, design=design, 
    min.rank=10) # avoid odd-looking plot due to too few PCs.
plotDiffusionMap(sce.block, use_dimred="PCA", colour_by="Gata3") + fontsize

## ---- echo=FALSE, results='hide'-----------------------------------------
saveRDS(file="th2_data.rds", sce)
gc()

## ---- eval=!on.bioc, echo=FALSE, results="hide"--------------------------
## # Cleaning out DLLs, with some protection to get it to work.
## repeat {
##     present <- setdiff(loadedNamespaces(),
##         c(rownames(installed.packages(.Library, priority="base")),
##         "AnnotationDbi", "GenomeInfoDb", "scran"))
##     discarded <- 0L
##     for (pkg in present) {
##        try({
##            unloadNamespace(pkg)
##            discarded <- discarded + 1L
##        }, silent=TRUE)
##     }
##     if (discarded==0) break
## }
## library(R.utils)
## library(org.Mm.eg.db)
## library(BiocStyle)
## gcDLLs()

## ------------------------------------------------------------------------
incoming <- read.table("counttable_es.csv", header=TRUE, row.names=1)
my.ids <- rownames(incoming)
symb <- mapIds(org.Mm.eg.db, keys=my.ids, keytype="ENSEMBL", column="SYMBOL")
head(symb)

## ------------------------------------------------------------------------
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
location <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=my.ids, 
    column="CDSCHROM", keytype="GENEID")
is.mito <- location == "chrM" & !is.na(location)
sum(is.mito)

## ------------------------------------------------------------------------
is.spike <- grepl("^ERCC", my.ids)
sum(is.spike)

## ------------------------------------------------------------------------
sce <- SingleCellExperiment(list(counts=as.matrix(incoming)), 
    rowData=DataFrame(Symbol=symb, Chr=location))
isSpike(sce, "ERCC") <- is.spike
sce <- calculateQCMetrics(sce, feature_controls=list(ERCC=is.spike, Mito=is.mito))
sce

## ------------------------------------------------------------------------
sce <- sce[grepl("ENSMUS", rownames(sce)) | isSpike(sce),]
dim(sce)

## ---- echo=FALSE, results='hide'-----------------------------------------
saveRDS(file="mesc_data.rds", sce)
gc()

## ---- echo=FALSE, results='asis'-----------------------------------------
if (!on.bioc) { 
    cat("Users can install all required packages and execute the workflow by following the instructions at https://www.bioconductor.org/help/workflows/simpleSingleCell.\n")
}
cat("The workflow takes less than an hour to run on a desktop computer with 8 GB of memory.\n")

## ------------------------------------------------------------------------
sessionInfo()

## ---- eval=on.bioc, echo=FALSE, results='hide'---------------------------
unlink(all.basenames)
unlink("GSE61533_HTSEQ_count_results.xls")

