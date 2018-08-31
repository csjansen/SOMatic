#!/bin/bash

if [ $# -eq 0 ]    # Script invoked with no command-line args?
then
	echo "Usage: ./regionCounts.sh [options] [options] -RawDataFile <raw sam file list location> -Partitions <partition file> -Output <output file location>"
	echo "Options: <default>"
	echo "-Log2: Log2(x+1) scale RPKM"

  	exit          # Exit and explain usage.
                            # Usage: scriptname -options
                            # Note: dash (-) necessary
fi
LogScale=0
while (( "$#" )); 
do
  case "$1" in
    -RawDataFile) RawDataFile=$2;;
    -Partitions) Partitions=$2;;
    -Output) Outfile=$2;;
    -Log2) LogScale=1
  esac

  shift
done

if [ "$LogScale" = 0 ]
then
../bin/regionCounts/regionCounts -RawDataFile $RawDataFile -Partitions $Partitions -Output $Outfile
fi
if [ "$LogScale" = 1 ]
then
../bin/regionCounts/regionCounts -RawDataFile $RawDataFile -Partitions $Partitions -Output $Outfile -Log2
fi
