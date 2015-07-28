#!/bin/bash

if [ $# -eq 0 ]    # Script invoked with no command-line args?
then
  echo "Usage: `basename $0` [required options] [options]"
  echo "Required Options:"
  echo "-SOMName: SOM name"

  exit          # Exit and explain usage.
                            # Usage: scriptname -options
                            # Note: dash (-) necessary
fi
Method="TwoClosest"
AddToChrom=""
while (( "$#" )); 
do
  case "$1" in
    -SOMName) SOMName=$2;;
  esac

  shift
done

../bin/cluster/cluster -SOMFile ../$SOMName.som -Clusters1 ../$SOMName/data/Cluster1.txt -Clusters2 ../$SOMName/data/Cluster2.txt

sed -i -e "s/var ClusterUnits = 0;/var ClusterUnits = 1;/g" ../$SOMName/options.js
sed -i -e "s/var ClusterProfiles = 0;/var ClusterProfiles = 1;/g" ../$SOMName/options.js
sed -i -e "s/var ClusterBoth = 0;/var ClusterBoth = 1;/g" ../$SOMName/options.js
