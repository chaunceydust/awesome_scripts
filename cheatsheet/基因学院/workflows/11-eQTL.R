## ----echo=FALSE,results="hide"-------------------------------------------
date()

## ----echo=FALSE,results="hide"-------------------------------------------
scoresCis = function(...){NULL}  # does commented unevald code get checked?
suppressMessages({  # include some warnings on symbol replacements
if (!("GGdata" %in% installed.packages()[,1])) {
 source("http://www.bioconductor.org/biocLite.R")
 biocLite("GGdata")
}
if (!("SNPlocs.Hsapiens.dbSNP144.GRCh37" %in% installed.packages()[,1])) {
 source("http://www.bioconductor.org/biocLite.R")
 biocLite("SNPlocs.Hsapiens.dbSNP144.GRCh37")
}
library(knitcitations)
library(bibtex)
allbib = read.bibtex("allbib.bib")
library(GenomeInfoDb)
library(S4Vectors)
library(GGtools)
library(GGdata)
library(yri1kgv)
library(snpStats)
library(scatterplot3d)
library(lumi)
library(parallel)
library(foreach)
library(doParallel)
library(biglm)
library(lumiHumanAll.db)
library(rmeta)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(eQTL)
})

## ----dolo,echo=FALSE-----------------------------------------------------
   cc = new("CisConfig")
   chrnames(cc) = "21"
   radius(cc) = 25000L
   lkp = try(library(parallel))
   if (!inherits(lkp, "try-error")) {
      nc = 2 # max(c(1,min(c(10, detectCores()-1)))) # attempt to trim for workflow builder
      options(mc.cores=nc)
      geneApply(cc) = mclapply
      if (.Platform$OS.type == "windows") geneApply(cc) = lapply
      }
   estimates(cc) = FALSE
   set.seed(1234)
   f1 <- All.cis( cc )  # devel: cisScores

## ----lkex,eval=FALSE-----------------------------------------------------
## cc = new("CisConfig") # take a default configuration
## chrnames(cc) = "21"   # confine to chr21
## estimates(cc) = FALSE # no point estimates neede
## f1 <- All.cis( cc )   # compute the tests; can be slow without attendance
##                       # to parallelization

## ----lookf1,echo=TRUE----------------------------------------------------
length(f1)
f1[1:3]
metadata(f1)

## ----demoy,fig=TRUE,fig.width=7,fig.height=4,echo=FALSE,results="hide"----
suppressMessages({
library(yri1kgv)
if (!exists("c20")) c20 = getSS("yri1kgv", "chr20")
par(mfrow=c(1,2))
plot_EvG(probeId("o67h4JQSuEa02CJJIQ"), rsid("rs2259928"), c20,
  main="observed expr.")
if (!exists("c20f")) c20f = clipPCs(c20, 1:10)
plot_EvG(probeId("o67h4JQSuEa02CJJIQ"), rsid("rs2259928"), c20f,
  main="10 expr. PC removed")
})

## ----shoit,eval=FALSE----------------------------------------------------
## plot_EvG(probeId("o67h4JQSuEa02CJJIQ"), rsid("rs2259928"), c20f,
##   main="10 expr. PC removed")

## ----bag,fig=TRUE,fig.width=4,fig.height=4,echo=FALSE--------------------
library(snpStats)
library(scatterplot3d)
tmp = as.raw(1:253)
yy = g2post(tmp)
EB = yy %*% c(0,1,2) 
scatterplot3d(yy[,1], yy[,3], EB, xlab="Pr(A/A)", ylab="Pr(B/B)", zlab="mean num. B")

## ----bag2----------------------------------------------------------------
library(GGtools)
library(yri1kgv)
library(lumiHumanAll.db)
if (!exists("y22")) y22 = getSS("yri1kgv", "chr22")
y22
dim(exprs(y22))
fn = featureNames(y22)[1:5]

## ----getseq--------------------------------------------------------------
library(lumi)
id2seq(fn) # get the 50mer for each probe
# and some annotation

## ----getann--------------------------------------------------------------
select( lumiHumanAll.db, keys=fn, keytype="PROBEID", columns=c("SYMBOL", "CHR", "ENTREZID"))

## ----getgen--------------------------------------------------------------
gt22 <- smList(y22)[[1]]  # access to genotypes
as( gt22[1:5,1:5], "character" )
cs22 = col.summary(gt22)  # some information on genotypes
cs22[1:10,]

## ----showscript,eval=FALSE-----------------------------------------------
## library(parallel)
## newcl = makePSOCKcluster(c("master", paste0("node00", 1:3)))
## library(foreach)
## library(doParallel)
## registerDoParallel(cores=8)  # may want to keep at 5
## 
## library(GGtools)
## ceuDemoRecov = try(ciseqByCluster( newcl,
##    chromsToRun=19:22, finaltag="partceu100k",
##    outprefix="ceurun",
##    ncoresPerNode=8, targetfolder="/freshdata/CEU_DEMO"  ))
## save(ceuDemoRecov, file="ceuDemoRecov.rda")
## stopCluster(newcl)
## stopImplicitCluster()
## sessionInfo()

## ----lkqq,eval=FALSE-----------------------------------------------------
## binnedQQ(partceu100k_dt, ylim=c(0,30), xlim=c(-2,15), end45=12)

## ----co,echo=FALSE-------------------------------------------------------
update_fdr_filt = function (tab, filt = function(x) x, by = c("pairs", "snps", 
    "probes")[1]) {
    require(GGtools, quietly = TRUE)
    tab = filt(tab)
    psinds = grep("permScore", colnames(tab), value = TRUE)
    nr = nrow(tab)
    pscores = vector("numeric", nr * length(psinds))
    for (np in 1:length(psinds)) pscores[(((np - 1) * nr) + 1):(np * 
        nr)] = tab[[psinds[np]]]
    if (by == "pairs") {
        newfd = pifdr(tab$score, pscores)
    }
    else {
        if (by == "snps") 
            byvar = "snp"
        else if (by == "probes") 
            byvar = "probeid"
        base = tab[, max(score), by = byvar]
        maxbysnp = base$V1
        ol = list()
        pnames = grep("permScore", names(tab))
        for (i in 1:length(pnames)) {
            tab$score = tab[, pnames[i], with = FALSE]
            ol[[i]] = tab[, max(score), by = byvar]$V1
        }
        newfd = pifdr(maxbysnp, as.numeric(unlist(ol)))
        tab = base
    }
    tab$fdr = newfd
    tab
}

#
# defining here so that release version can be used
#
filtgen.maf.dist = function (maf.dist, 
     validate.tab = function(tab) all(c("mindist", 
    "MAF") %in% colnames(tab))) 
{
    stopifnot(is.atomic(maf.dist))
    stopifnot(length(maf.dist) == 2)
    maf = maf.dist[1]
    dist = maf.dist[2]
    function(tab) {
        stopifnot(isTRUE(validate.tab(tab)))
        tab[tab$mindist <= dist & tab$MAF >= maf, ]
    }
}

eqsens_dt = function (dtab, filtgen = filtgen.maf.dist, by = c("pairs", "snps", 
    "probes")[1], targfdrs = c(0.05, 0.01, 0.005), parmslist = list(mafs = c(0.025, 
    0.05, 0.075, 0.1, 0.125), dists = c(1000, 5000, 10000, 25000, 
    50000, 1e+05))) 
{
    parmset = data.matrix(do.call(expand.grid, parmslist))
    ntune = nrow(parmset)
    ans = foreach(curp = 1:ntune) %dopar% {
        tmp = update_fdr_filt(tab = dtab, filt = filtgen(parmset[curp, 
            ]), by = by)
        sapply(targfdrs, function(x) sum(tmp$fdr <= x))
    }
    hold = t(sapply(ans, force))
    colnames(hold) = paste0("at_", targfdrs)
    cbind(parmset, hold)
}
filtgen.maf.dist = function (maf.dist, validate.tab = function(tab) all(c("mindist", 
    "MAF") %in% colnames(tab))) 
{
    stopifnot(is.atomic(maf.dist))
    stopifnot(length(maf.dist) == 2)
    maf = maf.dist[1]
    dist = maf.dist[2]
    function(tab) {
        stopifnot(isTRUE(validate.tab(tab)))
        tab[tab$mindist <= dist & tab$MAF >= maf, ]
    }
}
plotsens = function (eqsout, ylab = "count of eQTL at given FDR", 
    title = "cis radius in bp") {
    require(reshape2)
    require(ggplot2)
    mdf = melt(data.frame(eqsout), id.vars = c("mafs", "dists"))
    vind = which(names(mdf) == "variable")
    names(mdf)[vind] = "FDR"
    mdf[, vind] = gsub("at_", "", mdf[, vind])
    ggplot(data = mdf) + geom_point(aes(x = mafs, y = value, 
        colour = FDR)) + facet_grid(~dists) + theme(axis.text.x = element_text(angle = 90)) + 
        ylab(ylab) + labs(title = title)
}

## ----dosens,fig=TRUE-----------------------------------------------------
data(partceu100k_dt)
library(foreach)  # basic function includes a %dopar%
library(doParallel)
library(parallel)
# dcdown = max(c(detectCores()-1,1)) ## use with lots of RAM
dcdown = 1
registerDoSEQ()
if (.Platform$OS.type != "windows") {
 #cl = makeForkCluster(dcdown)
 registerDoParallel(cores=dcdown) # nesting?
 }

## ----doeq,eval=FALSE-----------------------------------------------------
## eq1 = eqsens_dt(partceu100k_dt)

## ----nnn,echo=FALSE------------------------------------------------------
data(sensSave)
eq1 = sensSave

## ----lkar----------------------------------------------------------------
args(eqsens_dt)

## ----coded,eval=FALSE----------------------------------------------------
## library(data.table)
## data("partceu100k_dt.rda")
## scoresCis("CPNE1", partceu100k_dt)

## ----disc----------------------------------------------------------------
distcat = cut(partceu100k_dt$mindist,c(-1, 1, 1000, 5000, 10000, 50000, 100001))
fdrcat = cut(partceu100k_dt$fdr,c(-.01,.005, .05, .1, .2, 1.01))
fdrcat = relevel(fdrcat, "(0.2,1.01]")
mafcat = cut(partceu100k_dt$MAF,c(0,.05, .1, .2, .3, .51))
approm = 1*partceu100k_dt$chromcat878 %in% c("1_Active_Promoter", "3_Poised_Promoter")

## ----fit-----------------------------------------------------------------
partceu100k_dt = cbind(partceu100k_dt, distcat, fdrcat, mafcat, approm)
set.seed(1234)
train = sample(1:nrow(partceu100k_dt), 
   size=floor(nrow(partceu100k_dt)/2), replace=FALSE)
library(biglm)
b1 = bigglm(isgwashit~distcat+fdrcat+mafcat+approm, fam=binomial(),
 data=partceu100k_dt[train,], maxit=30)

## ----cali----------------------------------------------------------------
pp = predict(b1, newdata=partceu100k_dt[-train,], type="response")
summary(pp)
cpp = cut(pp, c(0,.025, .05, .12, .21))
table(cpp)
sapply(split(partceu100k_dt$isgwashit[-train], cpp), mean)

## ----demomodco,fig=TRUE,fig.width=7,fig.height=4-------------------------
tmat = matrix(rownames(summary(b1)$mat),nc=1)
est = summary(b1)$mat[,1]
library(rmeta)
forestplot(tmat, est, est-.01, est+.01, xlog=TRUE,
  boxsize=.35, graphwidth=unit(3, "inches"),
  xticks=exp(seq(-4,2,2)))

## ----cleanup, echo=FALSE-------------------------------------------------
rm(list=ls())
gc()

## ----results='asis',echo=FALSE-------------------------------------------
bibliography() #style="markdown")

## ----sess----------------------------------------------------------------
sessionInfo()

