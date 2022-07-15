#=====================================packages basic=====================================

# display your library path
.libPaths()
# display your packages
library()
# display your loaded packages
search()

#=====================================packages install====================================
## From CRAN
# Install the package named "ggplot2"
install.packages("tidyverse",repos="http://mirror.bjtu.edu.cn/ ")
# install Rtools first
# Install the package named "XML"
install.packages("Your/Path/to/XML_3.98-1.20.tar.gz", repos = NULL, type = "source")

## From Bioconductor
# Install the package named "Iranges"
if (!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("biomaRt")
BiocManager::install("clusterProfiler")

source("http://bioconductor.org/biocLite.R")
#### 有安全隐患
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
biocLite("RGalaxy")##这样就用中科大的镜像来下载包啦
## bioconductor还有很多其它镜像：https://www.bioconductor.org/about/mirrors/

## From Github
# Install the package named "flipPlots"
install.packages("devtools")
library(devtools)
install_github("Displayr/flipPlots")


#=====================================packages library=====================================
# library ggplt2
library(ggplot2)


#=====================================packages update======================================
# list packagescould be updated
old.packages()
# update all packages
update.packages(ask = FALSE)
# update one packages
update.packages("ggplot2")



#=====================================packages remove=======================================
# remove ggplot2
remove.packages("ggplot2")




#=====================================packages backup=======================================

# list all installed pacckages
installed.packages()
oldip <- installed.packages()[,1]
save(oldip, file = "Your/Path/to/InstalledPackages.Rdata")
# start your new R
load("Your/Path/to/InstalledPackages.Rdata")
newip <- installed.packages()[,1]
for (i in setdiff(oldip, newip)) {
  install.packages(i)
}

