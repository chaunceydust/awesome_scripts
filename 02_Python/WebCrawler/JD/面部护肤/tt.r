data=read.table('Parse_List.txt',header=T)
dd=unique(data)

new=read.table('Parse_List_new.txt',header=T)

gg=c('面部护肤','香水彩妆')

for(i in gg){
	g1=grep(i,as.character(data[,3]))
	d1=data[g1,]
	d2=unique(d1)
	write.table(d2,paste0(i,'.txt'),sep='\t',quote=F,row.names=F,col.names=T)
}


