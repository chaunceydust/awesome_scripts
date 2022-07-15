data=read.table('Parse_List.txt',header=T)
dd=unique(data)
write.table(dd,'Parse_List_new.txt',sep='\t',quote=F,row.names=F,col.names=T)

group=unique(as.character(data[,3]))

for(i in group){
	i1=as.character(data[,3]) %in% i
	d1=data[i1,]
	write.table(d1,paste0(i,'.txt'),quote=F,row.names=F,col.names=T,sep='\t')
}

i1=as.character(data[,3]) %in% group[6]
dd=data[i1,]
gg=unique(as.character(dd[,1]))
for(i in gg){
	i1=as.character(dd[,1]) %in% i
	d1=dd[i1,]
	write.table(d1,paste0(i,'.txt'),quote=F,row.names=F,col.names=T,sep='\t')
}

