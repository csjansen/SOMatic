#!/bin/bash

if [ $# -eq 0 ]    # Script invoked with no command-line args?
then
  echo "Usage: `basename $0` [required options] "
  echo "Required Options:"
  echo "-SOMName <SOM name>"
  echo "-Rows <Number of Rows in your SOM> "
  echo "-Cols <Number of Columns in your SOM> "
  echo "-Metaclusters <Number of metaclusters you would like to start with> "
  echo "-MetaclustersEnd <Number of metaclusters you would like to end with> "
  echo "-Trials <Number of trials you'd like to run.>"
  echo "-Sparse"
  echo "-Dimensionality <Dimensionality>"
  exit          # Exit and explain usage.
                            # Usage: scriptname -options
                            # Note: dash (-) necessary
fi
Options=""
Dimensionality=-1
while (( "$#" ));
do
  case "$1" in
    -SOMName) SOMName=$2;;
    -Rows) Rows=$2;;
    -Cols) Cols=$2;;
    -Metaclusters) Metaclusters=$2;;
    -MetaclustersEnd) MetaclustersEnd=$2;;
    -Trials) Trials=$2;;
	-Dimensionality) Dimensionality=$2;;
    -DistanceMetric) Options="$Options -DistanceMetric $2";;
  esac

  shift
done

echo "../bin/metasomThread/metasom -Rows $Rows -Cols $Cols -SOMFile ../$SOMName.som -Metaclusters $Metaclusters -MetaclustersEnd $MetaclustersEnd -Trials $Trials -Outfile ../$SOMName/data/MetaClusters -genePrefix ../$SOMName/data/som/units/unit -Dimensionality $Dimensionality $Options"
../bin/metasomThread/metasom -Rows $Rows -Cols $Cols -SOMFile ../$SOMName.som -Metaclusters $Metaclusters -MetaclustersEnd $MetaclustersEnd -Trials $Trials -Outfile ../$SOMName/data/MetaClusters -genePrefix ../$SOMName/data/som/units/unit -Dimensionality $Dimensionality $Options

#../bin/metasomgene/metasomgene -UnitPrefix ../$SOMName/data/som/units/unit -MetaclusterFile ../$SOMName/data/MetaClusters -OutPrefix ../$SOMName/data/ -Rows $Rows -Cols $Cols -Metaclusters $Metaclusters
Meta=`head -n 1 ../$SOMName/data/MetaClusters | sed "s/.*: \(.*\)/\1/g"`
sed -i -e "s/var MetaOn = 0/var MetaOn = 1/g" ../$SOMName/options.js
sed -i -e "s/var MetaClusterNumber = .*/var MetaClusterNumber = $Meta/g" ../$SOMName/options.js

../bin/metasomgene/metasomgene -UnitPrefix ../$SOMName/data/som/units/unit -MetaclusterFile ../$SOMName/data/MetaClusters -OutPrefix ../$SOMName/data/ -Rows $Rows -Cols $Cols -Metaclusters $Meta
