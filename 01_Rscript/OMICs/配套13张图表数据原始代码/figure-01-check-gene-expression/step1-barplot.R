rm(list = ls())
options(stringsAsFactors = F)
a=read.table('all.counts.id.txt',header = T)
dim(a)


cg=a[a[,1]=='pho',7:16]
library(ggpubr)
library(stringr)
dat=data.frame(gene=as.numeric(cg),
               sample=names(cg),
               group=str_split(names(cg),'_',simplify = T)[,1]
               )
ggbarplot(dat,x='sample',y='gene',color = 'group')




cg=a[a[,1]=='Spps',7:16]
library(ggpubr)
library(stringr)
dat=data.frame(gene=as.numeric(cg),
               sample=names(cg),
               group=str_split(names(cg),'_',simplify = T)[,1]
)
ggbarplot(dat,x='sample',y='gene',color = 'group')






