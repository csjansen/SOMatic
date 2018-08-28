#!/bin/bash
#$ -N Som
#$ -q sam,bio
#$ -pe openmp 1


./Som_Meta_Clusters -Rows 20 -Cols 30 -SOMFile ../../ -Metaclusters 100 -Trials 20 -Outfile FusionATAC.AIC.v11.100.cluster

# ./Som_density_kmeans rows cols somfile mincluster maxcluster numberoftrials outfile1 AIC_output BIC_output BIC_centroids_output BIC_centroids_after_reclustering_output_file


#./Som_density_kmeans ../Somatic_Splice/Fusion.HL60.txt 20 30 20 30 HL60Clusters.txt /pub/public-www/csjansen/HL60RNA.v2/data/som/units/unit /samlab/csjansen/SOMatic/HL60ATAC.v2.som /samlab/csjansen/SOMatic/HL60RNA.v2.som GOHL60/GO_
#for i in `seq 1 20`; do
#	for j in `seq 15 20`; do
#./Som_density_kmeans 40 60 /samlab/csjansen/SOMatic/FusionATAC.CosDist.40x60.v6.som 14 14 40 outputradii3.txt ATAC3.cluster ATAC4.cluster Atac_centroids Atac_PostCentroids
#./Som_density_kmeans 40 60 /samlab/csjansen/SOMatic/FusionATAC.CosDist.40x60.v6.som 20 20 1 outputradii3.txt ATAC3.cluster ATAC4.cluster Atac_centroids Atac_PostCentroids
#./Som_density_kmeans 40 60 /bio/zengw/SOM/SOMatic/Bcl11b_SOM_combat_Gata3KD_removed_log_scale_40by60.som 200 200 1 BennyRNAAIC4060.200.cluster BennyRNABIC4060.200.cluster RNA_centroids RNA_PostCentroids
#./Som_density_kmeans 40 60 /bio/zengw/SOM/SOMatic/Bcl11b_SOM_combat_Gata3KD_removed_log_scale_40by60.som 40 40 5 BennyRNAAIC4060.40.cluster BennyRNABIC4060.40.cluster RNA_centroids RNA_PostCentroids
#./Som_density_kmeans 40 60 /bio/zengw/SOM/SOMatic/Bcl11b_SOM_combat_Gata3KD_removed_log_scale_40by60.som 60 60 5 BennyRNAAIC4060.60.cluster BennyRNABIC4060.60.cluster RNA_centroids RNA_PostCentroids
#./Som_density_kmeans 20 30 /samlab/csjansen/SOMatic/XenoRNAFusion.20x30.v10.som 75 75 20 Xeno.v10.AIC.75.cluster Xeno.v10.BIC.75.cluster RNA_centroids RNA_PostCentroids
#./Som_density_kmeans 20 30 /samlab/csjansen/SOMatic/XenoRNAFusion.20x30.v10.som 2 40 20 Xeno.v9.AIC.cluster Xeno.v9.BIC.cluster RNA_centroids RNA_PostCentroids
#./Som_density_kmeans 20 30 /samlab/csjansen/SOMatic/XenoRNAFusion.20x30.v9.som 75 75 5 Xeno.v9.AIC.75.cluster Xeno.v9.BIC.75.cluster RNA_centroids RNA_PostCentroids
#./Som_Meta_Clusters 20 35 /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v7.som 70 70 20 FusionRNA.AIC.v7.70.cluster FusionRNA.BIC.v7.60.cluster RNA_centroids RNA_PostCentroids
#	done
#done
#./Som_density_kmeans ../Somatic_Splice/Fusion.CosDist.txt 40 60 20 35 FusionClusters5.txt /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4/data/som/units/unit /samlab/csjansen/SOMatic/FusionATAC.CosDist.40x60.v4.som /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4.som GOHL60/GO_ 5
#./Som_density_kmeans ../Somatic_Splice/Fusion.CosDist.txt 40 60 20 35 FusionClusters6.txt /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4/data/som/units/unit /samlab/csjansen/SOMatic/FusionATAC.CosDist.40x60.v4.som /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4.som GOHL60/GO_ 6
#./Som_density_kmeans ../Somatic_Splice/Fusion.CosDist.txt 40 60 20 35 FusionClusters7.txt /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4/data/som/units/unit /samlab/csjansen/SOMatic/FusionATAC.CosDist.40x60.v4.som /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4.som GOHL60/GO_ 7
#./Som_density_kmeans ../Somatic_Splice/Fusion.CosDist.txt 40 60 20 35 FusionClusters8.txt /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4/data/som/units/unit /samlab/csjansen/SOMatic/FusionATAC.CosDist.40x60.v4.som /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4.som GOHL60/GO_ 8
#./Som_density_kmeans ../Somatic_Splice/Fusion.CosDist.txt 40 60 20 35 FusionClusters9.txt /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4/data/som/units/unit /samlab/csjansen/SOMatic/FusionATAC.CosDist.40x60.v4.som /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4.som GOHL60/GO_ 9
#./Som_density_kmeans ../Somatic_Splice/Fusion.CosDist.txt 40 60 20 35 FusionClusters10.txt /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4/data/som/units/unit /samlab/csjansen/SOMatic/FusionATAC.CosDist.40x60.v4.som /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4.som GOHL60/GO_ 10
#./Som_density_kmeans ../Somatic_Splice/Fusion.CosDist.txt 40 60 20 35 FusionClusters11.txt /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4/data/som/units/unit /samlab/csjansen/SOMatic/FusionATAC.CosDist.40x60.v4.som /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4.som GOHL60/GO_ 11
#./Som_density_kmeans ../Somatic_Splice/Fusion.CosDist.txt 40 60 20 35 FusionClusters12.txt /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4/data/som/units/unit /samlab/csjansen/SOMatic/FusionATAC.CosDist.40x60.v4.som /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4.som GOHL60/GO_ 12
#./Som_density_kmeans ../Somatic_Splice/Fusion.CosDist.txt 40 60 20 35 FusionClusters13.txt /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4/data/som/units/unit /samlab/csjansen/SOMatic/FusionATAC.CosDist.40x60.v4.som /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4.som GOHL60/GO_ 13
#./Som_density_kmeans ../Somatic_Splice/Fusion.CosDist.txt 40 60 20 35 FusionClusters14.txt /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4/data/som/units/unit /samlab/csjansen/SOMatic/FusionATAC.CosDist.40x60.v4.som /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4.som GOHL60/GO_ 14
#./Som_density_kmeans ../Somatic_Splice/Fusion.CosDist.txt 40 60 20 35 FusionClusters15.txt /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4/data/som/units/unit /samlab/csjansen/SOMatic/FusionATAC.CosDist.40x60.v4.som /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4.som GOHL60/GO_ 15
#./Som_density_kmeans ../Somatic_Splice/Fusion.CosDist.txt 40 60 20 35 FusionClusters16.txt /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4/data/som/units/unit /samlab/csjansen/SOMatic/FusionATAC.CosDist.40x60.v4.som /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4.som GOHL60/GO_ 16
#./Som_density_kmeans ../Somatic_Splice/Fusion.CosDist.txt 40 60 20 35 FusionClusters17.txt /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4/data/som/units/unit /samlab/csjansen/SOMatic/FusionATAC.CosDist.40x60.v4.som /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4.som GOHL60/GO_ 17
#./Som_density_kmeans ../Somatic_Splice/Fusion.CosDist.txt 40 60 20 35 FusionClusters18.txt /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4/data/som/units/unit /samlab/csjansen/SOMatic/FusionATAC.CosDist.40x60.v4.som /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4.som GOHL60/GO_ 18
#./Som_density_kmeans ../Somatic_Splice/Fusion.CosDist.txt 40 60 20 35 FusionClusters19.txt /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4/data/som/units/unit /samlab/csjansen/SOMatic/FusionATAC.CosDist.40x60.v4.som /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4.som GOHL60/GO_ 19
#./Som_density_kmeans ../Somatic_Splice/Fusion.CosDist.txt 40 60 20 35 FusionClusters20.txt /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4/data/som/units/unit /samlab/csjansen/SOMatic/FusionATAC.CosDist.40x60.v4.som /samlab/csjansen/SOMatic/FusionRNA.CosDist.20x35.v4.som GOHL60/GO_ 20
