#library("optparse")
#option_list = list(
#  make_option(c("--HypoFile"),type="character",default=NULL,help="Hypo File Name",metavar="character"),
#  make_option(c("--MetaClusterFile"),type="character",default=NULL,help="Meta Cluster File from MetaClusters in SOMatic",metavar="character"),
#  make_option(c("--SOMFile"),type="character",default=NULL,help="SOM file from SOMatic",metavar="character"),
#  make_option(c("--ClusterNum"),type="integer",default=1,help="Number of Meta Clusters in Meta CLuster File",metavar="character"),
#  make_option(c("--SampleList"),type="character",default=NULL,help="SOM Sample List"),
#  make_option(c("--LinkedReplicates"),type="logical",default=FALSE,help="If Controls and Experiment reps are linked"),
#  make_option(c("--OutputName"),type="character",default="Hypo.pdf",help="Outfile File Name")
#);
#opt_parser = OptionParser(option_list=option_list);
#opt=parse_args(opt_parser);
#HypoFile = read.delim(opt$HypoFile,header=T,comment.char="#")
#clusters2<- read.delim("/bio/csjansen/SOM_Meta_Clusters/Xeno.v10.AIC.75.cluster", header=F)
#clusters2<- read.delim("/bio/csjansen/SOM_Meta_Clusters/BennyRNAAIC4060.300.cluster", header=F)
#clusters2<- read.delim("/bio/csjansen/SOMatic/xenoRNA.40x60.v1/data/MetaClusters", header=F)
#clusters2<- read.delim("/bio/csjansen/SOMatic/FusionATAC.CosDist.40x60.v11/data/MetaClusters", header=F)
#Atac <- read.delim(opt$SOMFile,header=F,comment.char="#")
#clusters2 <- read.delim(opt$MetaClusterFile,header=F,comment.char="#")
#Atac <- read.delim("/samlab/csjansen/SOMatic/statATAC.20x30.v3.som", header=F, comment.char="#")
#Atac <- read.delim("/bio/zengw/SOM/SOMatic/Bcl11b_SOM_combat_Gata3KD_removed_log_scale_40by60.som", header=F, comment.char="#")
#Atac <- read.delim("/bio/csjansen/SOMatic/FusionATAC.CosDist.40x60.v11.som", header=F, comment.char="#")
#Atac <- read.delim("/bio/csjansen/SOMatic/xenoRNA.40x60.v1.som", header=F, comment.char="#")
#Atac <- read.delim("/bio/csjansen/SOMatic/XenoRNAFusion.20x30.v10.som", header=F, comment.char="#")

#samples <- read.delim("/bio/zengw/SOM/SOMatic/Bcl11b_SOM_combat_Gata3KD_removed_log_scale_40by60/data/sample.list", header=F, comment.char="#")
#samples <- read.delim("/samlab/csjansen/SOMatic/FusionATAC.CosDist.40x60.v11/data/sample.list", header=F, comment.char="#")
#samples <- read.delim(opt$SampleList, header=F, comment.char="#")
samples <- read.table("/bio/csjansen/SOMatic/statATAC.20x30.v3/data/sample.list", header=F, comment.char="#",stringsAsFactors=F)
#asamples <- read.delim("/bio/csjansen/SOMatic/XenoRNAFusion.20x30.v10/data/Beta-catenin-Morphs7", header=F, comment.char="#")
#bsamples <- read.delim("/bio/csjansen/SOMatic/XenoRNAFusion.20x30.v10/data/Beta-catenin-Control7", header=F, comment.char="#")
#OutputName <- "/pub/public-www/csjansen/xenoRNA.40x60.v1-Hypo.pdf"
#clusternum <- 17-1
Atac <- read.table("/share/samdata/csjansen/SOMatic/statATAC.20x30.v3.som", header=F,stringsAsFactors=F)
#Atac=Atac[-1,]
MetaClusters <- read.delim("/share/samdata/csjansen/SOMatic/statATAC.20x30.v3/data/MetaClusters", header=F)
MetaClusters = MetaClusters[-1,]
clusters2=MetaClusters
#OutputName <- opt$OutputName
OutputName <- "/pub/public-www/csjansen/statATAC.20x30.v3-Hypo.pdf"
#clusternum <- opt$ClusterNum-1
clusternum =  102
SOMrows = 19
SOMcols = 29
rowcol=c()
colcol=c()
Atac=data.matrix(Atac)
for(i in 0:SOMrows) {
  for(j in 0:SOMcols) {
    rowcol=c(rowcol,i)
    colcol=c(colcol,j)
  }
}
Atac = cbind(rowcol,colcol,Atac)
colnames(Atac)=c('V1','V2',samples$V1)
names0h = rownames(traits)[traits$X0h==1]
names24h = rownames(traits)[traits$X24h==1]
top200h=names(head(sort(colSums(Atac)[names0h],decreasing=TRUE), n = 20))
top2024h=names(head(sort(colSums(Atac)[names24h],decreasing=TRUE), n = 20))
#LinkedReplicates = opt$LinkedReplicates
LinkedReplicates = F
#rows <- 40
#cols <- 60
#clusters2<-clusters2[-1,]
#asamples=HypoFile[HypoFile[,2*(1)]==1,1]
#bsamples=HypoFile[HypoFile[,2*(1)+1]==1,1]
#for(i in 2:((ncol(HypoFile)-1)/2)) {
# asamples = cbind.data.frame(asamples,HypoFile[HypoFile[,2*(i)]==1,1])
#  bsamples = cbind.data.frame(bsamples,HypoFile[HypoFile[,2*(i)+1]==1,1])
#}
asamples=top200h
bsamples=top2024h

#print(ncol(HypoFile))
if (!require("reshape2")) {
  install.packages("reshape2", dependencies = TRUE)
  suppressPackageStartupMessages(library(reshape2))
}
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  suppressPackageStartupMessages(library(ggplot2))
}
if (!require("ggdendro")) {
  install.packages("ggdendro", dependencies = TRUE)
  suppressPackageStartupMessages(library(ggdendro))
}
if (!require("grid")) {
  install.packages("grid", dependencies = TRUE)
  suppressPackageStartupMessages(library(grid))
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  suppressPackageStartupMessages(library(RColorBrewer))
}
if (!require("plyr")) {
  install.packages("plyr", dependencies = TRUE)
  suppressPackageStartupMessages(library(plyr))
}
#asamples=list(as.matrix(asamples))
#bsamples=list(as.matrix(bsamples))
alist=c()
names=c()
#print(bsamples)
#for(i in 1:10) {
#i=121
#print(ncol(asamples))
  anums=match(asamples,samples$V1)
  bnums=match(bsamples,samples$V1)
#print(asamples)
  traitpvals=c()
  traitnegs=c()
  for(j in 0:clusternum) {
    cluster = clusters2[clusters2$V3==j,]
    HeatAtac = join(cluster,as.data.frame(Atac),by=c("V1","V2"))
  
    HeatAtac2 = HeatAtac[,4:ncol(HeatAtac)]

    #as=exp(as.matrix(HeatAtac2[,anums]))-1
    as=as.matrix(HeatAtac2[,anums])
    #bs=exp(as.matrix(HeatAtac2[,bnums]))-1
    bs=as.matrix(HeatAtac2[,bnums])
    asubs=c()
    bsubs=c()
    csubs=c()
    dsubs=c()
    if(!LinkedReplicates) {
      for(n in 1:(ncol(as)-1)) {
          if(n<ncol(as))
          for(k in (n+1):ncol(as)) {
           
            asubs=c(asubs,abs(as[,n]-as[,k]))
            #asubs=abs(as[,n]-as[,k])
          }
        }
        for(n in 1:(ncol(bs)-1)) {
          if(n<ncol(bs))
          for(k in (n+1):ncol(bs)) {
            bsubs=c(bsubs,abs(bs[,n]-bs[,k]))
          }
        }
        for(n in 1:ncol(as)) {
          for(k in n:ncol(bs)) {
            csubs=c(csubs,abs(as[,n]-bs[,k]))
          }
        }
      for(n in 1:ncol(as)) {
        for(k in n:ncol(bs)) {
          dsubs=c(dsubs,as[,n]-bs[,k])
        }
      }
    } else {
      for(n in 1:(ncol(as)-1)) {
        if(n<ncol(as))
          for(k in (n+1):ncol(as)) {
            asubs=c(asubs,abs(as[,n]-as[,k]))
          }
      }
      for(n in 1:(ncol(bs)-1)) {
        if(n<ncol(bs))
          for(k in (n+1):ncol(bs)) {
            bsubs=c(bsubs,abs(bs[,n]-bs[,k]))
          }
      }
      for(n in 1:(ncol(as))) {
        csubs=c(csubs,abs(as[,n]-bs[,n]))
      }
      for(n in 1:ncol(as)) {
        dsubs=c(dsubs,as[,n]-bs[,n])
      }
    }
    #print(i)
    #print(j)
    if(length(csubs)<3){
      print(i)
      print(j)
      csubs=c(0,0)
    }
    if(length(asubs)<3){
      print(i)
      print(j)
      asubs=c(0,0)
    }
    utest=wilcox.test(csubs,c(asubs,bsubs),alternative="greater")
    #boxplot(csubs,c(asubs,bsubs),names=c("With Trait", "Without Trait"),ylab="Expression",ylim=c(0,100))
    
    if(abs(mean(asubs))>0) {
      traitpvals=c(traitpvals,utest$p.value)
    } else {
      traitpvals=c(traitpvals,1)
    }
    
    if(mean(dsubs)>0) {
      traitnegs=c(traitnegs,1)
    } else {
      traitnegs=c(traitnegs,-1)
    }
  }
  traitpvals[traitpvals>.05/(clusternum)]=1
  
  traitpvals=-1*log10(traitpvals)
  maxval = 20
  traitpvals[traitpvals>maxval]=maxval
  
  
  traitpvals=traitpvals*traitnegs
 #print(m)
#write.table(traitpvals,file=paste0(),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE) 
  alist=cbind(alist,traitpvals)
#}
names=c("0hr vs 24hr")
#for(i in seq(2,ncol(HypoFile),2)) {
#  names=c(names,colnames(HypoFile)[i])
#}
colnames(alist)=names
rownames(alist)=make.names(0:clusternum)
#rownames(alist)=names
ind <- rowSums(alist == 0) != ncol(alist)
if(sum(ind)==0) {
  print("no sig metaclusters")
  return
} else {
  alist2=data.frame(alist[ind,])
  colnames(alist2) = names
  matrix_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
  matrix_fill_limits = NULL  

  row_labels=row.names(alist2)
  row_limits=row_labels
  col_labels=colnames(alist2)
  col_limits=col_labels
  if(ncol(alist)>1) {
    colDist = dist(t(alist2), method="euclidean")
    colHC = hclust(colDist, method="complete")
    colHC_data = dendro_data(as.dendrogram(colHC))
    
    
    col_limits = col_limits[colHC$order]
    col_labels = col_labels[colHC$order]
  }
  if(nrow(alist)>1) {
    rowDist = dist(alist2, method="euclidean")
    rowHC = hclust(rowDist, method="complete")
    rowHC_data = dendro_data(as.dendrogram(rowHC))


    row_limits = row_limits[rowHC$order]
    row_labels = row_labels[rowHC$order]
  }
  df = melt(as.matrix(alist2))
  p1 = ggplot(df, aes(x=Var2, y=Var1))
  p1 = p1 + geom_tile(aes(fill=value))
  p1 = p1 + scale_fill_gradientn(colours=matrix_palette,limits=c(-20,20))
  p1 = p1 + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
  p1 = p1 + scale_y_discrete(expand=c(0,0), limits=row_limits, labels=row_labels)
  p1 = p1 + scale_x_discrete(expand=c(0,0), limits=col_limits, labels=col_labels)
  p1 = p1 + theme(plot.margin=unit(c(0.00, 0.00, 0.00, 0.00),"inch"),axis.ticks.x=element_blank())
  p1 = p1 + labs(x="Hypothesis", y="Cluster #")
  row_labels_inches = 1.5*max(strwidth(row_labels, units="in", cex=8*(as.numeric(theme_get()$axis.text$size))*par()$cex/par()$ps))
  col_labels_inches = 1.5*max(strwidth(col_labels, units="in", cex=8*(as.numeric(theme_get()$axis.text$size))*par()$cex/par()$ps))#

  matrix_vp_y = 0
  matrix_vp_x = 0

  matrix_vp_h = 12/72.27 * nrow(alist2)+ col_labels_inches
  matrix_vp_w = .15* ncol(alist2) + row_labels_inches+2
  matrix_vp = viewport(
    y = matrix_vp_y,
    x = matrix_vp_x, 
    h = matrix_vp_h, 
    w = matrix_vp_w, 
    default.units="inch", 
    just=c("left","bottom")
  )

  pdf(OutputName, h = matrix_vp_h+.25, w=matrix_vp_w)
  print(p1, matrix_vp, newpage=FALSE)

  dev.off()
  #colHC_data = dendro_data(as.dendrogram(colHC))
  #col_ggdendro = ggplot(segment(colHC_data))
  #col_ggdendro = col_ggdendro + geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
  #p1 = p1 + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
  #col_ggdendro = col_ggdendro + scale_x_continuous(expand=c(0,0), labels=col_labels, limits=col_limits) 
  #col_ggdendro = col_ggdendro + scale_y_continuous(expand=c(0,0), labels=NULL)
  #col_ggdendro = col_ggdendro + theme(plot.margin=unit(c(0.10, 0.00, 0.00, 0.01), "inch"))
  #col_ggdendro = col_ggdendro + theme_dendro()
  #col_ggdendro = col_ggdendro + labs(x=NULL, y=NULL)
  #pdf("ColDendro.pdf",h=5,w=8)
  
}
