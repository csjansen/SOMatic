#!/bin/bash

if [ $# -eq 0 ]    # Script invoked with no command-line args?
then
  echo "Usage: `basename $0` [required options] "
  echo "Required Options:"
  echo "-SOMName <SOM name>"
  echo "-Rows <Number of Rows in your SOM> "
  echo "-Cols <Number of Columns in your SOM> "
  echo "-Metaclusters <Number of metaclusters you would like to try> "
  echo "-Trials <Number of trials you'd like to run.>"
  echo "-Sparse"
  exit          # Exit and explain usage.
                            # Usage: scriptname -options
                            # Note: dash (-) necessary
fi
Sparse=0
while (( "$#" ));
do
  case "$1" in
    -SOMName) SOMName=$2;;
    -Rows) Rows=$2;;
    -Cols) Cols=$2;;
    -Metaclusters) Metaclusters=$2;;
    -Trials) Trials=$2;;
    -Sparse) Sparse=1;;
  esac

  shift
done

sed -i -e "s/var MetaOn = 0/var MetaOn = 1/g" ../$SOMName/options.js
sed -i -e "s/var MetaClusterNumber = 75/var MetaClusterNumber = $Metaclusters/g" ../$SOMName/options.js

if [ "$Sparse" = 0 ]
then
../bin/metasom/metasom -Rows $Rows -Cols $Cols -SOMFile ../$SOMName.som -Metaclusters $Metaclusters -Trials $Trials -Outfile ../$SOMName/data/MetaClusters
fi

if [ "$Sparse" = 1 ]
then
../bin/metasom/metasom -Rows $Rows -Cols $Cols -SOMFile ../$SOMName.som -Metaclusters $Metaclusters -Trials $Trials -Outfile ../$SOMName/data/MetaClusters -Sparse
fi

../bin/metasomgene/metasomgene -UnitPrefix ../$SOMName/data/som/units/unit -MetaclusterFile ../$SOMName/data/MetaClusters -OutPrefix ../$SOMName/data/ -Rows $Rows -Cols $Cols -Metaclusters $Metaclusters
