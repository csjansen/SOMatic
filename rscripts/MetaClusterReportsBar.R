library("optparse")
option_list = list(
  make_option(c("--MetaClusterFile"),type="character",default=NULL,help="Meta Cluster File from MetaClusters in SOMatic",metavar="character"),
  make_option(c("--SOMFile"),type="character",default=NULL,help="SOM file from SOMatic",metavar="character"),
  make_option(c("--ClusterNum"),type="integer",default=1,help="Number of Meta Clusters in Meta CLuster File",metavar="character"),
  make_option(c("--SampleList"),type="character",default=NULL,help="Sample List used in SOM run"),
  make_option(c("--TrainingMatrix"),type="character",default=NULL,help="Training Matrix file used in SOM run"),
  make_option(c("--GeneFilePrefix"),type="character",default=NULL,help="Folder that contains output fron MetaClusterGene"),
  make_option(c("--TypeFile"),type="character",default=NULL,help="Type File"),
  make_option(c("--OutputPrefix"),type="character",default="Traits.pdf",help="Outfile File Name")
);
opt_parser = OptionParser(option_list=option_list);
opt=parse_args(opt_parser);
#sample <- read.delim("/bio/csjansen/TraitFiles/sample3.Xeno.list")
print(1)
#clusters2<- read.delim(opt$MetaClusterFile,header=F)
#typeBar <- read.delim("/bio/csjansen/SOMatic/FusionATAC.CosDist.40x60.v11/data/types", header=F)
typeBar <- read.delim(opt$TypeFile, header=F)
#typeBar <- read.delim("/bio/csjansen/SOMatic/xenoATAC.40x60.v2/data/xenoTypes", header=F)

#clusters2<- read.delim("/bio/csjansen/SOMatic/xenoATAC.40x60.v2/data/MetaClusters", header=F)
clusters2<- read.delim(opt$MetaClusterFile, header=F,stringsAsFactors=F)
clusters2<-clusters2[-1,]
#print(as.numeric(clusters2))
print(1)
Atac <- read.delim(opt$SOMFile, header=F, comment.char="#")
#Atac <- read.delim("/bio/csjansen/SOMatic/xenoATAC.40x60.v2.som", header=F, comment.char="#")
#Atac <- read.delim("/bio/csjansen/SOMatic/encode.40x60.v6.som", header=F, comment.char="#")
Atac2=Atac

Atac2=Atac2[,-1]
Atac2=Atac2[,-1]
#Atac=cbind(Atac$V1,Atac$V2,Atac2+.0001)
clusternum <- opt$ClusterNum-1
#clusternum <- 115-1
print(1)
#samples <- read.delim(opt$SampleList, header=F, comment.char="#")
print(1)
#samples <- read.delim("/bio/csjansen/SOMatic/xenoATAC.40x60.v2/data/sample.list", header=F, comment.char="#")
samples <- read.delim(opt$SampleList, header=F, comment.char="#")
trainingMatrix <- read.delim(opt$TrainingMatrix, header=F, comment.char="#")
Signal = trainingMatrix[,-1]>1
colSignal = colSums(Signal)
#trainingMatrix <- read.delim("/bio/csjansen/SOMatic/scripts/ATAC.training.matrix",header=F, comment.char="#")
#trainingMatrix <- read.delim("/bio/csjansen/SOMatic/scripts/xenoAtacTrainingMatrix2",header=F, comment.char="#")
#GeneFilePrefix = "/bio/csjansen/SOMatic/xenoATAC.40x60.v2/data/Genes_"
#GeneFilePrefix = paste0(opt$GeneFilePrefix,"/Genes_")
#colnames(trainingMatrix)=c("ProbeID",as.matrix(samples))
#clusternum <- 101
#Atac <- read.delim("/bio/zengw/SOM/SOMatic/Bcl11b_SOM_combat_Gata3KD_removed_log_scale_40by60.som", header=F, comment.char="#")
#clusters2<- read.delim("/bio/csjansen/SOM_Meta_Clusters/BennyRNAAIC40606.cluster", header=F)
#clusters2<-clusters2[-1,]
#samples <- read.delim("/bio/zengw/SOM/SOMatic/Bcl11b_SOM_combat_Gata3KD_removed_log_scale_40by60/data/sample.list", header=F, comment.char="#")
#trainingMatrix <- read.delim("/bio/zengw/SOM/SOMatic/scripts/Combat_adjusted_matrix_Gata3KD_removed_modified_log_scale.txt", header=T, comment.char="#")
#GeneFilePrefix = "/bio/csjansen/SOM_Meta_Clusters_GO/Cluster-Final/Genes_"
colsAtac<-c()

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggdendro))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(proxy))
if (!require("plyr")) {
  install.packages("plyr", dependencies = TRUE)
  library(plyr)
}

clusters2$V1=as.numeric(clusters2$V1)
rows = max(clusters2$V1)
print(rows)
print(clusters2$V1)
cols = max(clusters2$V2)
print(cols)
print(clusters2$V2)
rowcol = c()
colcol = c()
for(i in 0:rows) {
        for(j in 0:cols) {
                rowcol = c(rowcol, i)
                colcol = c(colcol, j)
        }
}
#print(rowcol)
Atac = cbind(rowcol,colcol,Atac)

colnames(Atac)=c("V1","V2",t(samples))
blueColor=rgb(75,75,200,maxColorValue=255,alpha=255)
blackColor=rgb(0,0,70,maxColorValue=255,alpha=255)
yellowColor=rgb(200,200,20,maxColorValue=255,alpha=255)
matrix_palette <- colorRampPalette(c("blue3", "cyan", "yellow2","red3"))(n = 299)
matrix_palette3 <- colorRampPalette(c("lightblue","cyan", "yellow2","orange3"))(n = 299)
matrix_palette2 <- colorRampPalette(c("white","green"))(n = 299)

HeatAtac3=c()
Maxsig=0
for(i in 0:clusternum) {
  cluster = clusters2[clusters2$V3==i,]
  HeatAtac = join(cluster,Atac,by=c("V1","V2"))
  HeatAtac = data.matrix(HeatAtac[,4:ncol(HeatAtac)])
 Maxsig=max(c(Maxsig,max(HeatAtac))) 
#  Maxsig=max(c(Maxsig,max(Atac2$V3+.0001)))
  HeatAtac2 = colSums(HeatAtac)/nrow(HeatAtac)
  HeatAtac3 = rbind(HeatAtac3,HeatAtac2)
}
#trainingMatrix[trainingMatrix>Maxsig]=Maxsig
base_size = 8
Minsig = 0
for(i in 0:clusternum) {
#i=33
#geneListFile = read.delim(paste0(GeneFilePrefix,toString(i)), header=F, comment.char="#")
#colnames(geneListFile)=c("ProbeID","row","col")
#geneListNames = c()
#for(j in 1:nrow(geneListFile)) {
#  geneListNames=c(geneListNames,paste0(geneListFile[j],paste0(" ",paste0(paste0(toString(HeatAtac[j,1]),','),toString(HeatAtac[j,2])))))
#}
#geneList = join(geneListFile,trainingMatrix)
#geneList2 = geneList[,-1]
#geneList2 = geneList2[,-1]
#geneList2 = geneList2[,-1]
#dfgenes = melt(as.matrix(geneList2))
#for(j in 1:nrow(dfgenes)) {
#  dfgenes[j,1]=toString(geneList[dfgenes[j,1],1])
#}
cluster = clusters2[clusters2$V3==i,]
print(i)
print(nrow(cluster))
HeatAtac = join(cluster,Atac,by=c("V1","V2"))
rownames1 = c()
for(j in 1:nrow(HeatAtac)) {
  rownames1=c(rownames1,paste0('(',toString(HeatAtac[j,1]),',',toString(HeatAtac[j,2]),')'))
}
HeatAtac = data.matrix(HeatAtac[,4:ncol(HeatAtac)])
rownames(HeatAtac)=rownames1
HeatAtac2 = colSums(HeatAtac)/nrow(HeatAtac)
df = melt(as.matrix(HeatAtac))
df2 = melt(as.matrix(HeatAtac3))
dfbox = melt(as.matrix(HeatAtac2))

colSide_by = NULL#}
col_metadata = NULL
col_labels = NULL
row_labels = NULL
colSide_by = NULL
merge_col_mdata_on = NULL
row_metadata = NULL
merge_row_mdata_on = NULL
col_dendro = TRUE
row_dendro = TRUE
dist = "euclidean"
hclust = "complete"
height = NULL
width = NULL
rowSide_by=NULL

matrix_legend_title="value"
if (!is.null(rowSide_by)) {rowSide_by = strsplit(rowSide_by, ",")[[1]]} else {rowSide_by = NULL}
if (!is.null(col_labels) && col_labels != "none") {col_label_fields = strsplit(col_labels,",")[[1]]} else {col_label_fields=NULL}
if (!is.null(row_labels) && row_labels != "none") {row_label_fields = strsplit(row_labels,",")[[1]]} else {row_label_fields=NULL}



col_mdata_header = unique(c(merge_col_mdata_on, colSide_by, col_label_fields))
row_mdata_header = unique(c(merge_row_mdata_on, rowSide_by, row_label_fields))
row_mdata_header = c(merge_row_mdata_on, setdiff(row_mdata_header, intersect(col_mdata_header, row_mdata_header)))


same_mdata = FALSE
#colDist = dist(t(HeatAtac3), method="cosine")
colDist = dist(t(HeatAtac3))
colHC = hclust(colDist, method=hclust)
colHC_data = dendro_data(as.dendrogram(colHC))
col_ggdendro = ggplot(segment(colHC_data))
col_ggdendro = col_ggdendro + geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
#col_ggdendro = col_ggdendro + coord_flip()
col_ggdendro = col_ggdendro + scale_x_continuous(expand=c(0.0, 0.5), labels=NULL) 
col_ggdendro = col_ggdendro + scale_y_continuous(expand=c(0.0, 0.0), labels=NULL)
col_ggdendro = col_ggdendro + theme(plot.margin=unit(c(0.00, 0.00, 0.00, 0.00), "inch"))
col_ggdendro = col_ggdendro + theme_dendro()
col_ggdendro = col_ggdendro + labs(x=NULL, y=NULL)
rowDist = dist(HeatAtac, method="cosine")
dorows = 1

if(nrow(cluster)<2) {
  dorows=0
}
rowHC=c()
rowHC_data=c()
row_ggdendro=c()
if(dorows) {
  rowHC = hclust(rowDist, method=hclust)
  rowHC_data = dendro_data(as.dendrogram(rowHC))
  row_ggdendro = ggplot(segment(rowHC_data))
  row_ggdendro = row_ggdendro + geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
  row_ggdendro = row_ggdendro + coord_flip()
  row_ggdendro = row_ggdendro + scale_x_continuous(expand=c(0.0, 0.5), labels=NULL) 
  row_ggdendro = row_ggdendro + scale_y_continuous(expand=c(0.0, 0.0), labels=NULL)
  row_ggdendro = row_ggdendro + theme(plot.margin=unit(c(0.00, 0.10, 0.00, 0.00), "inch"))
  row_ggdendro = row_ggdendro + theme_dendro()
  row_ggdendro = row_ggdendro + labs(x=NULL, y=NULL)
}


col_labels = colnames(HeatAtac)
col_limits = colnames(HeatAtac)
col_limits = colnames(HeatAtac)[colHC$order]
col_labels = col_labels[colHC$order]
row_labels = rownames1
row_limits = rownames1
if(dorows) {
  row_limits = rownames1[rowHC$order]
  row_labels = row_labels[rowHC$order]
}
#clust = hclust(dist(t(Atac2)))
#colCut = cutreeDynamic(colHC,cutHeight=max(colHC$height)*.5)#,cutHeight=max(colHC$height)*.35)
#colCut = cutreeDynamicTree(colHC, minModuleSize=2)#,maxTreeHeight=1)
n=length(colHC$height)
colCut = cutree(colHC,h=max(colHC$height)*.4)
typeBar$V1=max(typeBar$V1)-typeBar$V1
colMat = t(rbind(seq(from=1,to=length(colCut)),t(typeBar$V1[colHC$order]),colCut[colHC$order],rep(1.5,length(colCut)),colSignal[colHC$order]))
colnames(colMat)=c("Number","Type","Cut","Y","Signal")
p2 = ggplot(data.frame(colMat), aes(x=Number,y=Y))
p2 = p2 + geom_tile(aes(fill=Cut))
p2 = p2 + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
p2 = p2 + scale_x_discrete(expand=c(0,0), limits=col_limits, labels=NULL)
p2 = p2 + scale_y_discrete(expand=c(0,0), limits=c(1,1), labels=NULL)
p2 = p2 + scale_fill_gradientn(colours=matrix_palette3)
#p1 = p1 + scale_fill_gradientn(colours=matrix_palette)
p2 = p2 + theme(plot.margin=unit(c(0.00, 0.00, 0.00, 0.00),"inch"),axis.ticks=element_blank())
p2 = p2 + labs(x=NULL, y=NULL)
p2 = p2 + guides(fill=guide_colourbar(
  title.position="top", 
  direction="horizontal", 
  title.hjust=0, 
  title=paste0("Meta Cluster ",toString(i))
  #  draw.ulim=FALSE,
  #  draw.llim=FALSE
))
#p1_legend = g_legend(p1)
p2 = p2 + theme(legend.position = "none",panel.border=element_rect(colour="black",fill=NA, size=1))

p3 = ggplot(data.frame(colMat), aes(x=Number,y=Y))
p3 = p3 + geom_tile(aes(fill=Type))
p3 = p3 + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
p3 = p3 + scale_x_discrete(expand=c(0,0), limits=col_limits, labels=NULL)
p3 = p3 + scale_y_discrete(expand=c(0,0), limits=c(1,1), labels=NULL)
p3 = p3 + scale_fill_gradientn(colours=matrix_palette)
#p1 = p1 + scale_fill_gradientn(colours=matrix_palette)
p3 = p3 + theme(plot.margin=unit(c(0.00, 0.00, 0.00, 0.00),"inch"),axis.ticks=element_blank())
p3 = p3 + labs(x=NULL, y=NULL)
p3 = p3 + guides(fill=guide_colourbar(
  title.position="top", 
  direction="horizontal", 
  title.hjust=0 
  #  draw.ulim=FALSE,
  #  draw.llim=FALSE
))
#p1_legend = g_legend(p1)
p3 = p3 + theme(legend.position = "none",panel.border=element_rect(colour="black",fill=NA, size=1))

p4 = ggplot(data.frame(colMat), aes(x=Number,y=Y))
p4 = p4 + geom_tile(aes(fill=Signal))
p4 = p4 + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
p4 = p4 + scale_x_discrete(expand=c(0,0), limits=col_limits, labels=NULL)
p4 = p4 + scale_y_discrete(expand=c(0,0), limits=c(1,1), labels=NULL)
p4 = p4 + scale_fill_gradientn(colours=matrix_palette)
#p1 = p1 + scale_fill_gradientn(colours=matrix_palette)
p4 = p4 + theme(plot.margin=unit(c(0.00, 0.00, 0.00, 0.00),"inch"),axis.ticks=element_blank())
p4 = p4 + labs(x=NULL, y=NULL)
p4 = p4 + guides(fill=guide_colourbar(
  title.position="top",
  direction="horizontal",
  title.hjust=0
  #  draw.ulim=FALSE,
  #  draw.llim=FALSE
))
#p1_legend = g_legend(p1)
p4 = p4 + theme(legend.position = "none",panel.border=element_rect(colour="black",fill=NA, size=1))


row_labels_inches = 1.5*max(strwidth(row_labels, units="in", cex=base_size*(as.numeric(theme_get()$axis.text$size))*par()$cex/par()$ps))
col_labels_inches = 1.50*max(strwidth(col_labels, units="in", cex=base_size*(as.numeric(theme_get()$axis.text$size))*par()$cex/par()$ps))
#row_labels_inches2 = strwidth("timepoint", units="in", cex=base_size*(as.numeric(theme_get()$axis.text$size))*par()$cex/par()$ps)
matrix_fill_limits = NULL  
df$value[df$value>100]=100
p1 = ggplot(df, aes(x=Var2, y=Var1))
p1 = p1 + geom_tile(aes(fill=value))
p1 = p1 + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
p1 = p1 + scale_x_discrete(expand=c(0,0), limits=col_limits, labels=col_labels)
p1 = p1 + scale_y_discrete(expand=c(0,0), limits=row_limits, labels=row_labels)
print(Minsig)
print(Maxsig)
p1 = p1 + scale_fill_gradientn(colours=matrix_palette, limits=c(0,Maxsig))
#p1 = p1 + scale_fill_gradientn(colours=matrix_palette)
p1 = p1 + theme(plot.margin=unit(c(0.00, 0.00, 0.00, 0.00),"inch"))
p1 = p1 + labs(x=NULL, y=NULL)
p1 = p1 + guides(fill=guide_colourbar(
  title.position="top", 
  direction="horizontal", 
  title.hjust=0, 
  title=paste0("Meta Cluster ",toString(i))
  #	draw.ulim=FALSE,
  #	draw.llim=FALSE
))
p1_legend = g_legend(p1)
p1 = p1 + theme(legend.position = "none")
gcol_labels = colnames(HeatAtac)
gcol_limits = colnames(HeatAtac)
gcol_limits = colnames(HeatAtac)[colHC$order]
gcol_labels = col_labels[colHC$order]
#grow_labels = make.unique(as.character(paste0(paste0(paste0(paste0(paste0(geneListFile$ProbeID," ("),geneListFile$row),","),geneListFile$col),")")))
#grow_limits = grow_labels
grow_labels_inches = 1.5#1.5*max(strwidth(grow_labels, units="in", cex=base_size*(as.numeric(theme_get()$axis.text$size))*par()$cex/par()$ps))

#genebox=ggplot(dfgenes, aes(x=Var2, y=Var1))
#genebox=genebox+geom_tile(aes(fill=value))
#genebox=genebox+theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
#genebox = genebox + scale_x_discrete(expand=c(0,0), limits=gcol_limits, labels=col_labels)
#genebox = genebox + scale_y_discrete(expand=c(0,0), limits=grow_limits, labels=grow_labels)
#genebox = genebox + scale_fill_gradientn(colours=matrix_palette, limits=c(Minsig,Maxsig))
#genebox = genebox + theme(plot.margin=unit(c(0.00, 0.00, 0.00, 0.00),"inch"))
#genebox = genebox + labs(x=NULL, y=NULL)
#genebox = genebox + guides(fill=guide_colourbar(
#  title.position="top", 
#  direction="horizontal", 
#  title.hjust=0
  #title=opt$matrix_legend_title
  #  draw.ulim=FALSE,
  #	draw.llim=FALSE
#))
#genebox_legend = g_legend(p1)
#genebox = genebox + theme(legend.position = "none")
#genebox = genebox + theme()
library(scales)
dfbox$value=dfbox$value-mean(dfbox$value)
box=ggplot(dfbox, aes(x=Var1, y=value))
box = box + geom_bar(stat="identity")
box=box+theme(axis.text.x = element_blank())
box = box + scale_x_discrete(expand=c(0,0), limits=col_limits, labels=col_labels)
box = box + scale_y_continuous(limits=c(-.75,.75),oob=rescale_none)
#box = box + scale_fill_gradientn(colours=matrix_palette, limits=matrix_fill_limits)
box = box + theme(plot.margin=unit(c(0.00, 0.00, 0.00, 0.00),"inch"))
box = box + labs(x=NULL, y=NULL)
#box = box + guides(fill=guide_colourbar(
#  title.position="top", 
##  direction="horizontal", 
#  title.hjust=0, 
#  title=opt$matrix_legend_title))
#box_legend = g_legend(box)
#box = box + theme(legend.position = "none")
#print(box)
matrix_genes_x = .1
matrix_genes_y = 0
matrix_genes_h = base_size/72.27 #* nrow(geneList2) + col_labels_inches
matrix_genes_w = .15* ncol(HeatAtac) + row_labels_inches
matrix_gene = viewport(
  y = matrix_genes_y,
  x = matrix_genes_x, 
  h = matrix_genes_h, 
  w = matrix_genes_w, 
  default.units="inch", 
  just=c("left","bottom")
)
matrix_box_x = matrix_genes_x+grow_labels_inches-.1
matrix_box_y = matrix_genes_h
matrix_box_h = 3
matrix_box_w = matrix_genes_w-grow_labels_inches+.1
matrix_box = viewport(
  y = matrix_box_y,
  x = matrix_box_x, 
  h = matrix_box_h, 
  w = matrix_box_w, 
  default.units="inch", 
  just=c("left","bottom")
)
ColSides = list(); ColSide_legends = list()
RowSides = list(); RowSide_legends = list()
matrix_vp_y = max(0, as.numeric(strwidth(rowSide_by, "in")) - col_labels_inches)+matrix_box_h+matrix_box_y+2
#matrix_vp_x = max(0, as.numeric(strwidth(colSide_by, "in")) - row_labels_inches)
matrix_vp_x = matrix_genes_x+grow_labels_inches-row_labels_inches
#matrix_vp_h =0
#matrix_vp_w =0
matrix_vp_h = base_size/72.27 * nrow(HeatAtac)*2 + col_labels_inches
matrix_vp_w = matrix_genes_w-grow_labels_inches+row_labels_inches
matrix_vp = viewport(
  y = matrix_vp_y,
  x = matrix_vp_x, 
  h = matrix_vp_h, 
  w = matrix_vp_w, 
  default.units="inch", 
  just=c("left","bottom")
)
#matrix_vp_y2 = max(0, as.numeric(strwidth(rowSide_by, "in")) - col_labels_inches)+matrix_box_h+matrix_box_y
#matrix_vp_x2 = max(0, as.numeric(strwidth(colSide_by, "in")) - row_labels_inches)
#if (is.null(height)) {matrix_vp_h2 = base_size/72.27 * 1 + col_labels_inches} else {matrix_vp_h = height}
#if (is.null(width)) {matrix_vp_w2 = .05 * ncol(HeatAtac) + row_labels_inches+matrix_vp_x} else {matrix_vp_w = width}
#matrix_vp2 = viewport(
#  y = matrix_vp_y2+matrix_vp_y+matrix_vp_h-col_labels_inches-.05,
#  x = matrix_vp_x2, 
###  h = matrix_vp_h2, 
#  w = matrix_vp_w2, 
#  default.units="inch", 
#  just=c("left","bottom")
#)
matrix_scale_h = sum(sapply(p1_legend$heights, convertUnit, "in"))
matrix_scale_w = sum(sapply(p1_legend$widths, convertUnit, "in"))
matrix_scale_vp = viewport(
  y = matrix_vp_y + matrix_vp_h + 0.01,
  x = matrix_vp_x + matrix_vp_w + 0.05,
  h = matrix_scale_h,
  w = matrix_scale_w,
  default.units = "inch",
  just = c("left", "bottom")
)
matrix_signal_vp = viewport(
  y = matrix_vp_y + matrix_vp_h+2*base_size/72.27,
  x = matrix_vp_x+ row_labels_inches,
  h = base_size/72.27*2,
  w = matrix_vp_w - row_labels_inches-.01,
  default.units = "inch",
  just = c("left", "bottom")
  )
matrix_type_vp = viewport(
  y = matrix_vp_y + matrix_vp_h+base_size/72.27,
  x = matrix_vp_x+ row_labels_inches,
  h = base_size/72.27*2,
  w = matrix_vp_w - row_labels_inches-.01,
  default.units = "inch",
  just = c("left", "bottom")
  )
matrix_cut_vp = viewport(
  y = matrix_vp_y + matrix_vp_h,
  x = matrix_vp_x+ row_labels_inches,
  h = base_size/72.27*2,
  w = matrix_vp_w - row_labels_inches-.01,
  default.units = "inch",
  just = c("left", "bottom")
  )
col_dendro_h = .5
colDendro_vp = viewport(
  y = matrix_vp_y + matrix_vp_h-0.25*length(colSide_by)+4*base_size/72.27,#-1.5*base_size/72.27
  x = matrix_vp_x + row_labels_inches+.01,
  h = col_dendro_h,
  w = matrix_vp_w - row_labels_inches-.01,
  default.units = "inch",
  just = c("left", "bottom")
)
row_dendro_w = 2.0
rowDendro_vp = viewport(
  y = matrix_vp_y + col_labels_inches,
  x = matrix_vp_x + matrix_vp_w + 0.25*length(rowSide_by), 
  h = matrix_vp_h - col_labels_inches,
  w = row_dendro_w,
  default.units = "inch",
  just = c("left", "bottom")
)
total_h = matrix_vp_y + matrix_vp_h + 0.25*length(colSide_by) + max(col_dendro_h, matrix_scale_h)
total_w = matrix_vp_x + matrix_vp_w + 0.25*length(rowSide_by) + max(row_dendro_w, matrix_scale_w)
legend_width_inch = 0
total_w = total_w + legend_width_inch
pdf(paste0(paste0(opt$OutputPrefix,toString(i)),".pdf"), h = total_h, w=total_w)
#pdf(paste0(paste0("/pub/public-www/csjansen/EncodeReportsBar-",toString(i)),".pdf"), h = total_h, w=total_w)
#print(genebox, matrix_gene, newpage=FALSE)
print(box, matrix_box, newpage=FALSE)
print(col_ggdendro, vp=colDendro_vp, newpage=FALSE)
print(p4, matrix_signal_vp, newpage=FALSE)
print(p3, matrix_type_vp, newpage=FALSE)
print(p2, matrix_cut_vp, newpage=FALSE)
print(p1, matrix_vp, newpage=FALSE)
pushViewport(matrix_scale_vp); grid.draw(p1_legend); upViewport()




# Print row dendrogram
if(dorows) {
  print(row_ggdendro, vp=rowDendro_vp, newpage=FALSE)
}
dev.off()

}
