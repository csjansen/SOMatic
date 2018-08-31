#!/bin/bash

if [ $# -eq 0 ]    # Script invoked with no command-line args?
then
	echo "Usage: ./partition.sh [options] -PeakDataFile <peak file list location> -Output <output file location>"
	echo "Options: <default>"
	echo "-MinFeature: Size of smallest partition. <200>"

  	exit          # Exit and explain usage.
                            # Usage: scriptname -options
                            # Note: dash (-) necessary
fi
MinFeature="200"
while (( "$#" )); 
do
  case "$1" in
    -PeakDataFile) PeakDataFile=$2;;
    -MinFeature) MinFeature=$2;;
    -Output) Outfile=$2;;
  esac

  shift
done

echo "../bin/partition/partition -PeakDataFile $PeakDataFile -MinFeature $MinFeature -Output $Outfile"
../bin/partition/partition -PeakDataFile $PeakDataFile -MinFeature $MinFeature -Output $Outfile
