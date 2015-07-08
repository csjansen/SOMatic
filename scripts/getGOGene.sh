#!/bin/bash

if [ $# -eq 0 ]    # Script invoked with no command-line args?
then
  echo "Usage: `basename $0` [required options] [options]"
  echo "Required Options:"
  echo "-SOMName: SOM name"
  echo "-Rows: Number of rows you'd like in your SOM "
  echo "-Cols: Number of columns you'd like in your SOM "
  echo "-Gene2GO: Gene2GO file.  See README.txt for file format "
  echo "-GeneInfo: Gene Info file.  See README.txt for file format"
  echo "-GOFile: GO obo file.  Download from http://geneontology.org/ontology/go.obo"
  echo "Options: [choices] <default>"
  echo "-Sanity: If set to true, only GO terms with 5 genes in the unit will be reported. [true, false] <true>"
  exit          # Exit and explain usage.
                            # Usage: scriptname -options
                            # Note: dash (-) necessary
fi
Sanity="true"
while (( "$#" )); 
do
  case "$1" in
    -SOMName) SOMName=$2;;
    -Rows) Rows=$2;;
    -Cols) Cols=$2;;
	-Gene2GO) Gene2GO=$2;;
	-GeneInfo) GeneInfo=$2;;
	-GOFile) GOFile=$2;;
	-Sanity) Sanity=$2;;
  esac

  shift
done

../bin/GO/getGO -GenePrefix ../$SOMName/data/som/units/unit -OutputPrefix ../$SOMName/data/som/GO/GO -Rows $Rows -Cols $Cols -Sanity $Sanity -Gene2GO $Gene2GO -GeneInfo $GeneInfo -GOFile $GOFile
sed -i -e "s/var GOTermsOn = 0;/var GOTermsOn = 1;/g" ../$SOMName/options.js
