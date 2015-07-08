#!/bin/bash

if [ $# -eq 0 ]    # Script invoked with no command-line args?
then
  echo "Usage: `basename $0` [required options] "
  echo "Required Options:"
  echo "-SOMName <SOM name>"
  echo "-Matrix <Training Matrix File Location> "
  echo "-Rows <Number of rows you'd like in your SOM> "
  echo "-Cols <Number of Columns you'd like in your SOM> "
  echo "-SampleList <File with list of samples> "
  echo "-Timesteps <Number of timesteps for your SOM> "
  echo "-Trials <Number of trials you'd like to run.  The best SOM will be chosen.>"
  exit          # Exit and explain usage.
                            # Usage: scriptname -options
                            # Note: dash (-) necessary
fi
while (( "$#" )); 
do
  case "$1" in
    -SOMName) SOMName=$2;;
    -Matrix) Matrix=$2;;
    -Rows) Rows=$2;;
    -Cols) Cols=$2;;
	-SampleList) SampleList=$2;;
	-Timesteps) Timesteps=$2;;
	-Trials) Trials=$2;;
  esac

  shift
done

rm -rf ../Som_Package
rm -rf ../$SOMName
tar -C ../ -zxf ../Som_Package.tgz
mv ../Som_Package ../$SOMName
sed -i -e "s/= 30/= $Rows/g" ../$SOMName/options.js
sed -i -e "s/= 50/= $Cols/g" ../$SOMName/options.js
../bin/train/trainsom -Rows $Rows -Cols $Cols -TrainingMatrix $Matrix -SOMFile ../$SOMName.som -Trials $Trials -Timesteps $Timesteps -Topology toroid
../bin/score/scoresom -SOMFile ../$SOMName.som -TrainingMatrix $Matrix -ScoreFile ../$SOMName.score
../bin/map/mapsom -SOMFile ../$SOMName.som -SampleList $SampleList -Prefix ../$SOMFile/data/som/map_
mv ../$SOMName/data/som/map_summery.map ../$SOMName/data/map_summery.map
../bin/units/getunits -ScoreFile ../$SOMName.score -Rows $Rows -Cols $Cols -Prefix ../$SOMName/data/som/units/unit

