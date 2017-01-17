#!/bin/bash

if [ $# -eq 0 ]    # Script invoked with no command-line args?
then
  echo "Usage: `basename $0` [required options] "
  echo "Required Options:<Description>[default]{options}"
  echo "-SOMName1 <SOM name for first SOM>"
  echo "-SOMName2 <SOM name for second SOM>"
  echo "-Row1 <Rows for first SOM> "
  echo "-Row2 <Rows for second SOM> "
  echo "-Col1 <Cols for first SOM> "
  echo "-Col2 <Cols for second SOM> "
  echo "-Type <Type of Data Used.  The first type is before the x.>[ATACxRNA]{ATACxRNA}"
  echo "-Algorithm <Type of Fusion Function Algorithm Used.  Will change depending on the Type chosen>"
  echo "	ATACxRNA-GREAT Algorithm-[OneClosest]{OneClosest,TwoClosest}"
  echo "-GTFFile <Location of GTF File for your organism>"
  echo "-GeneIDType <Type of ID used in RNA SOMs>[gene_id]{gene_id,gene_name}"
  echo "-ClusterNum1 <Number of Meta-Clusters in first SOM>"
  echo "-ClusterNum2 <Number of Meta-Clusters in second SOM>"
  echo "-OutputFolder <Output Folder for Pairwise Meta-Clusters>"
  exit          # Exit and explain usage.
                            # Usage: scriptname -options
                            # Note: dash (-) necessary
fi
Type="ATACxRNA"
Algorithm="OneClosest"
geneIDType="gene_id"
while (( "$#" )); 
do
  case "$1" in
    -SOMName1) SOMName1=$2;;
	-SOMName2) SOMName2=$2;;
	-Row1) Row1=$2;;
	-Row2) Row2=$2;;
	-Col1) Col1=$2;;
	-Col2) Col2=$2;;
    -Type) Type=$2;;
    -Algorithm) Algorithm=$2;;
    -GTFFile) gtfFile=$2;;
	-GeneIDType) geneIDType=$2;;
	-ClusterNum1) ClusterNum1=$2;;
	-ClusterNum2) ClusterNum2=$2;;
	-OutputFolder) OutputFolder=$2;;
  esac

  shift
done

../bin/Fusion -UnitPrefix1 ../$SOMName1/data/som/units/unit -Row1 $Row1 -Col1 $Col1 -UnitPrefix2 ../$SOMName2/data/som/units/unit -Row2 $Row2 -Col2 $Col2 -Type $Type -Algorithm $Algorithm -GTFFile $gtfFile -Output ../$SOMName1/data/Fusion.txt -GeneIDType $geneIDType

../bin/FusionCluster -FusionFile ../$SOMName1/Fusion.txt -ClusterFile1 ../$SOMName1/data/MetaClusters -ClusterFile2 ../$SOMName2/data/MetaClusters -Output ../$SOMName1/data/FusionCluster.txt

mkdir $OutputFolder

../bin/FusionBreakup -FusionClusterFile ../$SOMName1/data/FusionCluster.txt -FusionFile ../$SOMName1/data/Fusion.txt -ClusterNum1 $ClusterNum1 -ClusterNum2 $ClusterNum2 -OutputPrefix $OutputFolder  -Row1 $Row1 -Col1 $Col1 -Row2 $Row2 -Col2 $Col2

