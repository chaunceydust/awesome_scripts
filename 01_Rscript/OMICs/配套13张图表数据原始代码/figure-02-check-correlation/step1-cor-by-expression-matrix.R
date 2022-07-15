rm(list = ls())
options(stringsAsFactors = F)
a=read.table('../RNAseq/all.counts.id.txt',header = T)
dim(a)
dat=a[,7:16]
library(stringr)
ac=data.frame(group=str_split(colnames(dat),'_',simplify = T)[,1])
rownames(ac)=colnames(dat)
M=cor(log(dat+1))
pheatmap::pheatmap(M,
                   annotation_col = ac) 
pheatmap::pheatmap(scale(M),
                   annotation_col = ac) 

cg=dat[,colnames(dat)[grepl('_1',colnames(dat))]]
library(ggpubr)
cg=log(cg+1)
pairs(cg)






