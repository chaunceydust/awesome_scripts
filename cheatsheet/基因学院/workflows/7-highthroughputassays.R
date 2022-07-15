## ----echo=FALSE, results="hide", warning=FALSE---------------------------
suppressPackageStartupMessages(library(flowCore))
suppressPackageStartupMessages(library(flowStats))
suppressPackageStartupMessages(library(flowViz))

## ------------------------------------------------------------------------
## Load packages
library(flowCore)
library(flowStats)
library(flowViz) # for flow data visualization

## Load data
data(ITN)
ITN

## ------------------------------------------------------------------------
## Create a workflow instance and transform data using asinh
wf <- workFlow(ITN)
asinh <- arcsinhTransform()
tl <- transformList(colnames(ITN)[3:7], asinh, 
                      transformationId = "asinh")
add(wf, tl)

## ------------------------------------------------------------------------
## Identify T-cells population
lg <- lymphGate(Data(wf[["asinh"]]), channels=c("SSC", "CD3"),
         preselection="CD4", filterId="TCells", eval=FALSE,
         scale=2.5)
add(wf, lg$n2gate, parent="asinh")
print(xyplot(SSC ~ CD3| PatientID, wf[["TCells+"]],
             par.settings=list(gate=list(col="red", 
             fill="red", alpha=0.3))))

## ----eval=FALSE----------------------------------------------------------
## ## try http:// if https:// URLs are not supported
## source("https://bioconductor.org/biocLite.R")
## biocLite("flowCore")

## ----eval=FALSE----------------------------------------------------------
## library("flowCore")

## ----eval=FALSE----------------------------------------------------------
## browseVignettes(package="flowCore")

## ----eval=FALSE----------------------------------------------------------
## help.start()

## ------------------------------------------------------------------------
sessionInfo()

