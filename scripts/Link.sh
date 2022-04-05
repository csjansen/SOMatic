#!/bin/bash

if [ $# -eq 0 ]    # Script invoked with no command-line args?
then
  echo "Usage: `basename $0` [required options] "
  echo "Required Options:"
  echo "-SOMName1 <SOM name for DNA SOM>"
  echo "-Row1 <Number of Rows for SOM1> "
  echo "-Col1 <Number of Cols for SOM1> "
  echo "-SOMName2 <SOM name for RNA SOM> "
  echo "-Row2 <Number of Rows for SOM2> "
  echo "-Col2 <Number of Cols for SOM2> "
  echo "-GTFFile <Location of gene annotations> "
  echo "-Output <Linking Output Folder> "
  exit          # Exit and explain usage.
                            # Usage: scriptname -options
                            # Note: dash (-) necessary
fi
while (( "$#" ));
do
  case "$1" in
    -SOMName1) SOMName1=$2;;
    -SOMName2) SOMName2=$2;;
    -Row1) Row1=$2;;
    -Col1) Col1=$2;;
    -Row2) Row2=$2;;
    -Col2) Col2=$2;;
    -GTFFile) GTFFile=$2;;
    -Output) Output=$2;;
  esac

  shift
done


mkdir $Output
../bin/LinkUnits/LinkUnits -UnitPrefix1 ../$SOMName1/data/som/units/unit -Row1 $Row1 -Col1 $Col1 -UnitPrefix2 ../$SOMName2/data/som/units/unit -Row2 $Row2 -Col2 $Col2 -Output $Output/LinkUnitsFile -GTFFile $GTFFile -GeneIDType gene_id -AddChr -Underscore
../bin/LinkMeta/LinkMeta -FusionFile $Output/LinkUnitsFile -ClusterFile1 ../$SOMName1/data/MetaClusters -ClusterFile2 ../$SOMName2/data/MetaClusters -Output $Output/LinkMetaFile
Meta1=`head -n 1 ../$SOMName1/data/MetaClusters | sed "s/.*: \(.*\)/\1/g"`
Meta2=`head -n 1 ../$SOMName2/data/MetaClusters | sed "s/.*: \(.*\)/\1/g"`
../bin/LinkAnalyze/LinkAnalyze -FusionFile $Output/LinkUnitsFile -FusionClusterFile $Output/LinkMetaFile -ClusterNum1 $Meta1 -ClusterNum2 $Meta2 -OutputPrefix $Output -Row1 $Row1 -Col1 $Col1 -Row2 $Row2 -Col2 $Col2

