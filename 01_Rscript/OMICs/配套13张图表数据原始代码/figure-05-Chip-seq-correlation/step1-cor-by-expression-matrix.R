rm(list = ls())
options(stringsAsFactors = F)
a=read.table('SpearmanCorr_readCounts.tab')
pheatmap::pheatmap(a)
library(stringr)
ac=data.frame(group=str_split(rownames(a),'_',simplify = T)[,1])
rownames(ac)=colnames(a)
M=a
pheatmap::pheatmap(M,
                   annotation_col = ac) 
 

a=read.table('merge_SpearmanCorr_readCounts.tab')
pheatmap::pheatmap(a)


