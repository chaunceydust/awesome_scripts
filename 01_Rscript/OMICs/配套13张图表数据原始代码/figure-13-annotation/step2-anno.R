rm(list = ls())
load(file = '6-genesets.Rdata')

library(ggplot2)
library(clusterProfiler)
library(org.Dm.eg.db)

kegg_anno <- function(gs,pro){
  #pro='SppsKO_up_uniq'
  #gs=SppsKO_up_uniq
  gene_up= bitr(gs, fromType = "SYMBOL",
                toType = c( "ENTREZID"),
                OrgDb = org.Dm.eg.db)[,2] 
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'dme', 
                      keyType =  'ncbi-geneid',
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
  head(kk.up)[,1:6]
  kk=kk.up
  dotplot(kk)
  ggsave(paste0(pro,'_kk.png'))
  kk=DOSE::setReadable(kk, OrgDb='org.Dm.eg.db',keyType='ENTREZID')
  write.csv(kk@result,paste0(pro,'_kk.up.csv'))
  
}
l=list(SppsKO_up_uniq,SppsKO_down_uniq,
       PhoKO_up_uniq, PhoKO_down_uniq,
       SppsKO_up_PhoKO_up,SppsKO_down_PhoKO_down)
names(l)=c('SppsKO_up_uniq','SppsKO_down_uniq',
           'PhoKO_up_uniq', 'PhoKO_down_uniq',
           'SppsKO_up_PhoKO_up','SppsKO_down_PhoKO_down')
lapply(1:length(l), function(i){
  kegg_anno(l[[i]],names(l)[i])
})

library(ggplot2)
library(clusterProfiler)
library(org.Dm.eg.db)
gcSample=lapply(l, function(x){
  return( bitr(x, fromType = "SYMBOL",
               toType = c( "ENTREZID"),
               OrgDb = org.Dm.eg.db)[,2] )
})
xx <- compareCluster(gcSample, fun="enrichKEGG",
                     keyType =  'ncbi-geneid',
                     organism="dme", pvalueCutoff=0.05)
dotplot(xx)



