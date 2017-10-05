library("optparse")
option_list = list(
  make_option(c("--TraitFile"),type="character",default=NULL,help="Trait File Name",metavar="character"),
  make_option(c("--MetaClusterFile"),type="character",default=NULL,help="Meta Cluster File from MetaClusters in SOMatic",metavar="character"),
  make_option(c("--SOMFile"),type="character",default=NULL,help="SOM file from SOMatic",metavar="character"),
  make_option(c("--ClusterNum"),type="integer",default=1,help="Number of Meta Clusters in Meta CLuster File",metavar="character"),
  make_option(c("--OutputName"),type="character",default="Traits.pdf",help="Outfile File Name")
  );
opt_parser = OptionParser(option_list=option_list);
opt=parse_args(opt_parser);
samples<-read.delim(opt$TraitFile);
#sample <- read.delim("/bio/csjansen/TraitFiles/sample3.Xeno.list")
#sample=sample[,-2]#sample <- read.delim("/bio/csjansen/TraitFiles/Trait2.txt")

#sample <- read.delim("/bio/csjansen/TraitFiles/sample3.Xeno.list")
clusters2<- read.delim(opt$MetaClusterFile, header=F)
#clusters2<- read.delim("/bio/csjansen/SOMatic/XenoRNAFusion.20x30.v10/data/MetaClusters", header=F)
#clusters2<- read.delim("/bio/csjansen/SOM_Meta_Clusters/Xeno.v9.AIC.75.cluster", header=F)
clusters2<-clusters2[-1,]
Atac <- read.delim(opt$SOMFile, header=F, comment.char="#")
#Atac <- read.delim("/bio/csjansen/SOMatic/XenoRNAFusion.20x30.v10.som", header=F, comment.char="#")
#samples <- read.delim("/bio/zengw/SOM/SOMatic/Bcl11b_SOM_combat_Gata3KD_removed_log_scale_40by60/data/sample.list", header=F, comment.char="#")

#Atac <- read.delim("/samlab/csjansen/SOMatic/XenoRNAFusion.20x30.v9.som", header=F, comment.char="#")
clusternum = opt$ClusterNum-1
#clusternum <- 34
#clusternum <- 299
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggdendro))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(RColorBrewer))
if (!require("plyr")) {
  install.packages("plyr", dependencies = TRUE)
  library(plyr)
}
print(samples)
sample = samples[,-1]
#colnames(sample) = samples[1,-1]
rownames(sample) = samples[,1]
print(sample)
#sample=sample[,-1]
traitnames = colnames(sample)

cols = ncol(sample)
#i=2
HeatAtac3=c()
for(i in 0:clusternum) {
  cluster = clusters2[clusters2$V3==i,]
  HeatAtac = join(cluster,Atac,by=c("V1","V2"))
  HeatAtac = data.matrix(HeatAtac[,4:ncol(HeatAtac)])
  HeatAtac2 = colSums(HeatAtac)/nrow(HeatAtac)
  HeatAtac3 = rbind(HeatAtac3,HeatAtac2)
}
cors=cor(t(HeatAtac3),sample)
print(cors)
pvals = 2*pt(-abs(cors/(sqrt((1-cors^2)/(clusternum-1)))),df=clusternum+1)
#for each trait
#traitpvals=c()
#traiteigens = c()
#traitnegs=c()
#for(i in 1:cols) {
#  traitvect=sample[,i]
#  traitsamples=c()
  #nottraitvect = 1-sample[,i]
#  nottraitvect=c()
  #for each sample, check if trait
#  for(j in 1:nrow(sample)) {
#    if(traitvect[j]==1) {
#      traitsamples=rbind(traitsamples,HeatAtac3[,j])
#    } else {
##      nottraitvect=rbind(nottraitvect,HeatAtac3[,j])
#    }
#  }
#  cor(colMeans(traitsamples),t(nottraitvect))
  #traiteigenvector = colMeans(traitsamples)
  #traiteigens=rbind(traiteigens,traiteigenvector)
#}
  # for each metacluster,check mean and std dev
  #meantrait=c()
  #meannottrait=c()
  #stddevtrait=c()
  #stddevnottrait=c()
  #pvalrow=c()
  #negsrow=c()
 # for(j in 1:ncol(traitsamples)) {
#    meantrait=mean(traitsamples[,j])
    #stddevtrait=c(stddevtrait,sd(traitsamples[,j]))
#    meannottrait=mean(nottraitvect[,j])
    #stddevnottrait=c(stddevnottrait,sd(nottraitvect[,j]))
    #ttest = t.test(traitsamples[,j],nottraitvect[,j],var.equal=TRUE)
#    utest=wilcox.test(traitsamples[,j],nottraitvect[,j])
    #if(utest$p.value<.01) {
    #  print(i)
    #  print(j)
    ##  print(utest)
#      boxplot(traitsamples[,j],nottraitvect2[,j],names=c("With Trait", "Without Trait"),ylab="Expression",main=paste0(paste0(paste0("Meta Cluster: ",toString(j-1)),paste0(" Trait: ",traitnames[i])),paste0(" pval: ",toString(utest$p.value))))
    #  #title(main=paste0(paste0(paste0("Meta Cluster: ",toString(j)),paste0(" Trait: ",traitnames[i+1])),paste0(" pval: ",toString(ttest$p.value))))
      #cat("Press [enter] to continue")
      #line=readline()
    #}
 #   pvalrow=c(pvalrow,utest$p.value)
  #  if(meantrait>meannottrait) {
  #    negsrow=c(negsrow,1)
  #  } else {
  #    negsrow=c(negsrow,-1)
  #  }
  #}
  #traitpvals=rbind(traitpvals,pvalrow)
  #traitnegs=rbind(traitnegs,negsrow)
#}
#traitpvals=t(traitpvals)
#traitnegs=t(traitnegs)
pvals[pvals>.05/((clusternum+1))]=1
pvals=-1*log10(pvals)
maxval = min(10,max(pvals))
pvals=pvals/maxval
pvals[pvals>1]=1

cors = cors/abs(cors)
pvals=pvals*cors
#cors[pvals==1]=0
rownames(pvals)=0:clusternum
colnames(pvals)=traitnames
rownames(pvals) <- make.names(rownames(pvals))
pvals=pvals[rowSums(abs(pvals))>0,]
print(pvals)
matrix_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
matrix_fill_limits = NULL  

rowDist = dist(pvals, method="euclidean")
rowHC = hclust(rowDist, method="complete")
rowHC_data = dendro_data(as.dendrogram(rowHC))


row_labels=row.names(pvals)
row_limits=row_labels
row_limits = row_limits[rowHC$order]
row_labels = row_labels[rowHC$order]

col_labels=colnames(pvals)
col_limits=col_labels

df = melt(as.matrix(pvals))
p1 = ggplot(df, aes(x=Var2, y=Var1))
p1 = p1 + geom_tile(aes(fill=value))
p1 = p1 + scale_fill_gradientn(colours=matrix_palette, limits=c(-1,1))
p1 = p1 + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
p1 = p1 + scale_y_discrete(expand=c(0,0), limits=row_limits, labels=row_labels)
p1 = p1 + scale_x_discrete(expand=c(0,0), limits=col_limits, labels=col_labels)
p1 = p1 + theme(plot.margin=unit(c(0.00, 0.00, 0.00, 0.00),"inch"))
p1 = p1 + labs(x="Trait", y="Cluster #")

matrix_vp_y = 0
#matrix_vp_x = max(0, as.numeric(strwidth(colSide_by, "in")) - row_labels_inches)
matrix_vp_x = 0
#matrix_vp_h =0
#matrix_vp_w =0
matrix_vp_h = 8/72.27 * nrow(cors) +2
matrix_vp_w = .15* ncol(cors) +2
matrix_vp = viewport(
  y = matrix_vp_y,
  x = matrix_vp_x, 
  h = matrix_vp_h, 
  w = matrix_vp_w, 
  default.units="inch", 
  just=c("left","bottom")
)
pdf(opt$OutputName, h = matrix_vp_h+.5, w=matrix_vp_w)
#pdf("/pub/public-www/csjansen/XenoRNATrait.pdf", h = matrix_vp_h+.5, w=matrix_vp_w)
print(p1, matrix_vp,newpage=FALSE)
dev.off()
