library("optparse")
option_list = list(
  make_option(c("--HypoFile"),type="character",default=NULL,help="Hypo File Name",metavar="character"),
  make_option(c("--MetaClusterFile"),type="character",default=NULL,help="Meta Cluster File from MetaClusters in SOMatic",metavar="character"),
  make_option(c("--SOMFile"),type="character",default=NULL,help="SOM file from SOMatic",metavar="character"),
  make_option(c("--ClusterNum"),type="integer",default=1,help="Number of Meta Clusters in Meta CLuster File",metavar="character"),
  make_option(c("--SampleList"),type="character",default=NULL,help="SOM Sample List"),
  make_option(c("--LinkedReplicates"),type="logical",default=FALSE,help="If Controls and Experiment reps are linked"),
  make_option(c("--OutputName"),type="character",default="Hypo.pdf",help="Outfile File Name")
);
opt_parser = OptionParser(option_list=option_list);
opt=parse_args(opt_parser);
HypoFile = read.delim(opt$HypoFile,header=T,comment.char="#")
#clusters2<- read.delim("/bio/csjansen/SOM_Meta_Clusters/Xeno.v10.AIC.75.cluster", header=F)
#clusters2<- read.delim("/bio/csjansen/SOM_Meta_Clusters/BennyRNAAIC4060.300.cluster", header=F)
#clusters2<- read.delim("/bio/csjansen/SOMatic/xenoRNA.40x60.v1/data/MetaClusters", header=F)
#clusters2<- read.delim("/bio/csjansen/SOMatic/FusionATAC.CosDist.40x60.v11/data/MetaClusters", header=F)
Atac <- read.delim(opt$SOMFile,header=F,comment.char="#")
clusters2 <- read.delim(opt$MetaClusterFile,header=F,comment.char="#")
#Atac <- read.delim("/samlab/csjansen/SOMatic/XenoRNAFusion.20x30.v10.som", header=F, comment.char="#")
#Atac <- read.delim("/bio/zengw/SOM/SOMatic/Bcl11b_SOM_combat_Gata3KD_removed_log_scale_40by60.som", header=F, comment.char="#")
#Atac <- read.delim("/bio/csjansen/SOMatic/FusionATAC.CosDist.40x60.v11.som", header=F, comment.char="#")
#Atac <- read.delim("/bio/csjansen/SOMatic/xenoRNA.40x60.v1.som", header=F, comment.char="#")
#Atac <- read.delim("/bio/csjansen/SOMatic/XenoRNAFusion.20x30.v10.som", header=F, comment.char="#")

#samples <- read.delim("/bio/zengw/SOM/SOMatic/Bcl11b_SOM_combat_Gata3KD_removed_log_scale_40by60/data/sample.list", header=F, comment.char="#")
#samples <- read.delim("/samlab/csjansen/SOMatic/FusionATAC.CosDist.40x60.v11/data/sample.list", header=F, comment.char="#")
samples <- read.delim(opt$SampleList, header=F, comment.char="#")
#samples <- read.delim("/bio/csjansen/SOMatic/xenoRNA.40x60.v1/data/sample.list", header=F, comment.char="#")
#asamples <- read.delim("/bio/csjansen/SOMatic/XenoRNAFusion.20x30.v10/data/Beta-catenin-Morphs7", header=F, comment.char="#")
#bsamples <- read.delim("/bio/csjansen/SOMatic/XenoRNAFusion.20x30.v10/data/Beta-catenin-Control7", header=F, comment.char="#")
#OutputName <- "/pub/public-www/csjansen/xenoRNA.40x60.v1-Hypo.pdf"
#clusternum <- 17-1
OutputName <- opt$OutputName
clusternum <- opt$ClusterNum-1

LinkedReplicates = opt$LinkedReplicates
#rows <- 40
#cols <- 60
#clusters2<-clusters2[-1,]
asamples=c()
bsamples=c()
for(i in 1:((ncol(HypoFile)-1)/2)) {
  asamples = cbind(asamples,HypoFile[HypoFile[,2*(i)]==1,1])
  bsamples = cbind(bsamples,HypoFile[HypoFile[,2*(i)+1]==1,1])
}
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
for(i in 1:ncol(asamples)) {
  anums=match(HypoFile[asamples[,i],1],samples$V1)
  bnums=match(HypoFile[bsamples[,i],1],samples$V1)
  traitpvals=c()
  traitnegs=c()
  for(j in 0:clusternum) {
    cluster = clusters2[clusters2$V3==j,]
    HeatAtac = join(cluster,Atac,by=c("V1","V2"))
  
    HeatAtac2 = data.matrix(HeatAtac[,4:ncol(HeatAtac)])

    as=exp(as.matrix(HeatAtac2[,anums]))-1
    bs=exp(as.matrix(HeatAtac2[,bnums]))-1
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
}
names=c()
for(i in seq(2,ncol(HypoFile),2)) {
  names=c(names,colnames(HypoFile)[i])
}
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
    just=c("left","bottom") 
  )

  pdf(OutputName, h = matrix_vp_h+.25, w=matrix_vp_w)
  print(p1, matrix_vp, newpage=FALSE)

  dev.off()
}
