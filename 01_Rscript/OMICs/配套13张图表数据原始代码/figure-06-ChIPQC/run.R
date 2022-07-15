library(ChIPQC)
samples = read.csv(file.path(system.file("extdata", package="ChIPQC"),
                             "example_QCexperiment.csv"))
samples
exampleExp = ChIPQC(samples,annotaiton="hg19")
QCmetrics(exampleExp)  #shows a summary of the main QC metrics
ChIPQCreport(exampleExp)


# SampleID: 样本ID
# Tissue, Factor, Condition: 不同的实验设计对照信息
# 三列信息必须包含在sampleSheet里，如果没有某一列的信息设为NA。
# Replicate : 重复样本的编号
# bamReads : 实验组BAM 文件的路径（data/bams）
# ControlID : 对照组样本ID
# bamControl :对照组样本的bam文件路径
# Peaks :样本peaks文件的路径
# PeakCaller ：peak类型的字符串，可以是raw,bed,narrow,macs等。

(bams=list.files(path = '~/fly/CHIP-SEQ/align/',
                 pattern = '*.raw.bam$', full.names = T))
bams=bams[grepl('WT',bams)]
peaks=list.files(path = '~/fly/CHIP-SEQ/peaks-single/',pattern = '*_raw_peaks.narrowPeak', full.names = T)
peaks=peaks[grepl('WT',peaks)]


library(stringr)
SampleID=sub('.raw.bam','',basename(bams))
Replicate=str_split(basename(bams),'_',simplify = T)[,3]
Factor=str_split(basename(bams),'_',simplify = T)[,1]


samples=data.frame(SampleID=SampleID,
                   Tissue='WT', 
                   Factor=Factor, 
                   Replicate=1,            
                   bamReads=bams,                           
                   Peaks=peaks) 

exampleExp = ChIPQC(samples,
                    chromosomes='2L',
                    annotaiton="dm6")
QCmetrics(exampleExp)  #shows a summary of the main QC metrics
ChIPQCreport(exampleExp)
