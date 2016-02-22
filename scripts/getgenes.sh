#!/bin/bash

if [ $# -eq 0 ]    # Script invoked with no command-line args?
then
  echo "Usage: `basename $0` [required options] [options]"
  echo "Required Options:"
  echo "-SOMName: SOM name"
  echo "-Rows: Number of rows you'd like in your SOM "
  echo "-Cols: Number of columns you'd like in your SOM "
  echo "-GTFFile: Gene annotations file.  See README.txt for file format "
  echo "Options: [choices] <default>"
  echo "-Method: GREAT algorithm of choice. [TwoClosest] <TwoClosest>"
  echo "--AddToChrom: If your gtf file uses a different format for it's chromosomes, this option allows you to add text to all the chromosomes in the gtf file. <>"

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
    -Rows) Rows=$2;;
    -Cols) Cols=$2;;
	-GTFFile) GTFFile=$2;;
	-Method) Method=$2;;
	-AddToChrom) AddToChrom=$2;;
  esac

  shift
done

../bin/genes/getgenes -GTFFile $GTFFile -InputPrefix ../$SOMName/data/som/units/unit -OutputPrefix ../$SOMName/data/som/genes/gene -Rows $Rows -Cols $Cols -Method $Method #-AddToChrom $AddToChrom 
sed -i -e "s/var GenesOn = 0;/var GenesOn = 1;/g" ../$SOMName/options.js
