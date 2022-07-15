rm(list = ls())
options(stringsAsFactors = F)
a=read.table('../RNAseq/all.counts.id.txt',header = T)
dim(a)
dat=a[,7:16]
rownames(dat)=a[,1]
dat[1:4,1:4]
library(stringr)
group_list=str_split(colnames(dat),'_',simplify = T)[,1]
table(group_list)

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

load(file = 'deg_output.Rdata')
library(ggpubr)
DEG_PhoKO=DEG_PhoKO[rownames(DEG_SppsKO),]
po=data.frame(SppsKO=DEG_SppsKO$log2FoldChange,
              PhoKO=DEG_PhoKO$log2FoldChange)
ggscatter(po,'SppsKO','PhoKO')
sp <- ggscatter(po,'SppsKO','PhoKO',
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 3, label.y = 30)



