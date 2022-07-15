load(file = 'deg_output.Rdata')
library(ggpubr)
colnames(DEG_PhoKO)
DEG_PhoKO$log=log(DEG_PhoKO$baseMean+1)
DEG_PhoKO$change=ifelse(DEG_PhoKO$padj>0.05,'stable',
                        ifelse(DEG_PhoKO$log2FoldChange > 0,'up','down'))
table(DEG_PhoKO$change)
ggscatter(DEG_PhoKO,x="log" ,y="log2FoldChange",color = 'change')

DEG_SppsKO$log=log(DEG_SppsKO$baseMean+1)
DEG_SppsKO$change=ifelse(DEG_SppsKO$padj>0.05,'stable',
                         ifelse(DEG_SppsKO$log2FoldChange > 0,'up','down'))
table(DEG_SppsKO$change)
ggscatter(DEG_SppsKO,x="log" ,y="log2FoldChange",color = 'change')

library(UpSetR)

SppsKO_up=rownames(DEG_SppsKO[DEG_SppsKO$change=='up',])
SppsKO_down=rownames(DEG_SppsKO[DEG_SppsKO$change=='down',])
PhoKO_up=rownames(DEG_PhoKO[DEG_PhoKO$change=='up',])
PhoKO_down=rownames(DEG_PhoKO[DEG_PhoKO$change=='down',])

allG=unique(c(SppsKO_up,SppsKO_down,PhoKO_up,PhoKO_down))

df=data.frame(allG=allG,
              SppsKO_up=as.numeric(allG %in% SppsKO_up),
              SppsKO_down=as.numeric(allG %in% SppsKO_down),
              PhoKO_up=as.numeric(allG %in% PhoKO_up),
              PhoKO_down=as.numeric(allG %in% PhoKO_down))

upset(df)

colnames(df)
# "SppsKO_up"   "SppsKO_down" "PhoKO_up"    "PhoKO_down" 
SppsKO_up_uniq=df[df[,2]==1 & df[,3]==0 & df[,4]==0 & df[,5]==0 ,1] 
SppsKO_down_uniq=df[df[,2]==0 & df[,3]==1 & df[,4]==0 & df[,5]==0 ,1]
PhoKO_up_uniq=df[df[,2]==0 & df[,3]==0 & df[,4]==1 & df[,5]==0 ,1]
PhoKO_down_uniq=df[df[,2]==0 & df[,3]==0 & df[,4]==0 & df[,5]==1 ,1]
SppsKO_up_PhoKO_up=df[df[,2]==1 & df[,3]==0 & df[,4]==1 & df[,5]==0 ,1]
SppsKO_down_PhoKO_down=df[df[,2]==0 & df[,3]==1 & df[,4]==0 & df[,5]==1 ,1]
  
library(org.Dm.eg.db)
t1=toTable(org.Dm.egCHRLOC)
t2=toTable(org.Dm.egCHRLOCEND)
t3=toTable(org.Dm.egSYMBOL)
 

#　作业 


require(TxDb.Dmelanogaster.UCSC.dm6.ensGene  )
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
gs=as.data.frame(genes(txdb))
library(org.Dm.eg.db)
t100=toTable(org.Dm.egFLYBASE)
t101=toTable(org.Dm.egSYMBOL)
t200=merge(t100,t101,by='gene_id')
colnames(gs)[6]="flybase_id"
t300=merge(t200,gs,by="flybase_id")

SppsKO_up_uniq
SppsKO_bed=t300[match(SppsKO_up_uniq[SppsKO_up_uniq %in% t300$symbol ],
                      t300$symbol),4:6]

