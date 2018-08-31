#!/bin/bash

if [ $# -eq 0 ]    # Script invoked with no command-line args?
then
  echo "Usage: `basename $0` [required options]"
  echo "Required Options:"
  echo "-SOMName: SOM name"
  echo "-CutLevel: % of max height to make cut (suggest 40).  This tool can be biased by having a sample far away from the others. Choose a higher value in that case."
  echo "getclusters.sh needs to be run first."
  exit          # Exit and explain usage.
                            # Usage: scriptname -options
                            # Note: dash (-) necessary
fi
while (( "$#" )); 
do
  case "$1" in
    -SOMName) SOMName=$2;;
	-CutLevel) CutLevel=$2;;
  esac
  shift
done
	echo $CutLevel
	../bin/cutree/cutree -ClusterFile ../$SOMName/data/Cluster2.txt -CutLevel $CutLevel -SampleList ../$SOMName/data/sample.list -Output ../$SOMName/data/cutree.txt

