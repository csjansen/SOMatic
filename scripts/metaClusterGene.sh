#!/bin/bash

if [ $# -eq 0 ]    # Script invoked with no command-line args?
then
  echo "Usage: `basename $0` [required options] "
  echo "Required Options:"
  echo "-SOMName <SOM name>"
  echo "-Rows <Number of Rows in your SOM> "
  echo "-Cols <Number of Columns in your SOM> "
  echo "-Metaclusters <Number of metaclusters returned by metaClusterSOM.sh> "
  exit          # Exit and explain usage.
                            # Usage: scriptname -options
                            # Note: dash (-) necessary
fi
Sparse=0
Dimensionality=-1
while (( "$#" ));
do
  case "$1" in
    -SOMName) SOMName=$2;;
    -Rows) Rows=$2;;
    -Cols) Cols=$2;;
    -Metaclusters) Metaclusters=$2;;
  esac

  shift
done
echo $Metaclusters
sed -i -e "s/var MetaOn = 0/var MetaOn = 1/g" ../$SOMName/options.js
sed -i -e "s/var MetaClusterNumber = 75/var MetaClusterNumber = $Metaclusters/g" ../$SOMName/options.js


../bin/metasomgene/metasomgene -UnitPrefix ../$SOMName/data/som/units/unit -MetaclusterFile ../$SOMName/data/MetaClusters -OutPrefix ../$SOMName/data/ -Rows $Rows -Cols $Cols -Metaclusters $Metaclusters
