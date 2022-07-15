barplot(VADeaths, 
        col = NULL, #col: 设置条形底纹或者填充颜色。
        border =par("fg"),#border：设置条形边框颜色。如果设置为NA，则消除了边缘
        main = NULL, sub = NULL, 
        xlab = NULL,ylab = NULL, #xlab和ylab：设置x轴与y轴的lable
        xlim = NULL, ylim = NULL,  #xlim和ylim:设置图形x轴与y轴的范围。
        
        beside = FALSE,#beside:逻辑参数。如果FALSE，那么将绘画堆叠式的条形；如果是TRUE，将绘画并列式条形。
        horiz = FALSE, #horiz：逻辑参数。设置图形是水平或是垂直
        
        width = 1, #width：设置条形的宽度
        space = NULL, #space：设置各个条形间的间隔宽度。相当于各个条形宽度的一部分。
        names.arg = NULL,  #names.arg:设置条形标签（barlabels）。
        
        
        #density 和 angle : 设置柱子用线条填充，density 控制线条的密度， angel 控制线条的角度
        density = NULL,  #density:底纹的密度。默认值为NULL
        angle =45,  #angle：设置底纹的斜率
        
        
        axes = TRUE, #axes:逻辑参数。设置图形是否显示坐标轴。
        las=1,            #设置刻度值的方向, 0表示总是平行于坐标轴；1表示总是水平方向；2表示总是垂直于坐标轴；3表示总是垂直方向。
        
        yaxt= "s", #是否绘制Y坐标轴，s 绘制，n不绘制
        
        axisnames = TRUE, #axisnames：逻辑参数。设置是否显示x轴条形标签
        cex.axis=par("cex.axis"),#cex.axis:设置坐标轴数值的大小。
        cex.names=par("cex.axis"), #cex.names: 设置条形标签（barlabels）的大小。
        
        add = FALSE #add = “TRUE”将barplot加在目前已经有的图上
)


#载入调色板包
require("RColorBrewer")


#读取数据
x<-read.table("new_tax_abundance.txt",sep="\t",header = TRUE,check.names=FALSE,comment.char="@")

#获取组合颜色
mycol<-c(brewer.pal(9, "Set1"),brewer.pal(12, "Paired"),brewer.pal(12, "Set3"))

#数据处理一下，以适合用barplot绘制柱状图
x<-x[order(x[,2],decreasing =T),]
a<-x[,2:ncol(x)]
rownames(a)<-x[,1]
t<-as.matrix(a)


#修改绘图参数，调整绘制满意的柱状图

par(mar=c(18, 5, 4, 1.1))
p<-barplot(t,main="abundance",xaxs="i",ylim=c(0,1),axisnames=F,border = NA ,space=0.05,col=mycol[1:nrow(t)],ylab="Relative abundance ")
#自定义坐标轴
axis(side=1,p,labels=F)
labs <- colnames(t)
text(cex=1, x=p-0.25, y=-0.05, labs, xpd=TRUE, srt=45,adj=c(1,1))
#调整legend位置
box(bty="l")
legend("bottom",ncol=3,xpd=T,bty = "n",rownames(t),fill=mycol[1:nrow(t)],inset=c(0,-0.75),cex=1)