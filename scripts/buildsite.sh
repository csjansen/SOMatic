#!/bin/bash

if [ $# -eq 0 ]    # Script invoked with no command-line args?
then
  echo "Usage: `basename $0` [required options] "
  echo "Required Options:"
  echo "-SOMName <SOM name>"
  echo "-Matrix <Training Matrix File Location> "
  echo "-Rows <Number of Rows you'd like in your SOM - choose auto if you'd like the size to be chosen for you> "
  echo "-Cols <Number of Columns you'd like in your SOM> "
  echo "-SampleList <File with list of samples> "
  echo "-Timesteps <Number of timesteps for your SOM> "
  echo "-Trials <Number of trials you'd like to run.  The best SOM will be chosen.>"
  echo "-Sparse"
  echo "-LearningRate <Learning Rate for your program. default .2>"
  exit          # Exit and explain usage.
                            # Usage: scriptname -options
                            # Note: dash (-) necessary
fi
Sparse=0
LearningRate=.2
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
	-Sparse) Sparse=1;;
	-LearningRate) LearningRate=$2;;
  esac

  shift
done

sed -i -e 's/\r/\n/g' $Matrix
mac2unix $SampleList
awk '{$1=$1;print}' $SampleList > temp.samples
SampleList="temp.samples"
#sed -i -e 's/\r/\n/g' $SampleList

if [ "$Rows" = "auto" ]
then
	../bin/estimateSize/estimateSize -TrainingMatrix $Matrix -Output ../suggestedSize.txt
	Rows=1
	Cols=1
	while read -r ro co; do
		Rows=$ro
		Cols=$co
	done < ../suggestedSize.txt
fi
		echo $Rows
		echo $Cols
rm -rf ../Som_Package
rm -rf ../$SOMName
tar -C ../ -zxf ../Som_Package.tgz
mv ../Som_Package ../$SOMName
cp $SampleList ../$SOMName/data/sample.list
sed -i -e "s/= 20/= $Rows/g" ../$SOMName/options.js
sed -i -e "s/= 50/= $Cols/g" ../$SOMName/options.js
if [ "$Sparse" = 0 ]
then
echo ../bin/train/trainsom -Rows $Rows -Cols $Cols -TrainingMatrix $Matrix -SOMFile ../$SOMName.som -Trials $Trials -Timesteps $Timesteps -Topology toroid
	../bin/train/trainsom -Rows $Rows -Cols $Cols -TrainingMatrix $Matrix -SOMFile ../$SOMName.som -Trials $Trials -Timesteps $Timesteps -Topology toroid -LearningRate $LearningRate
../bin/score/scoresom -SOMFile ../$SOMName.som -TrainingMatrix $Matrix -ScoreFile ../$SOMName.score
fi
if [ "$Sparse" = 1 ]
then
echo ../bin/train/trainsom -Rows $Rows -Cols $Cols -TrainingMatrix $Matrix -SOMFile ../$SOMName.som -Trials $Trials -Timesteps $Timesteps -Topology toroid -Sparse 
	../bin/train/trainsom -Rows $Rows -Cols $Cols -TrainingMatrix $Matrix -SOMFile ../$SOMName.som -Trials $Trials -Timesteps $Timesteps -Topology toroid -Sparse -LearningRate $LearningRate
../bin/score/scoresom -SOMFile ../$SOMName.som -TrainingMatrix $Matrix -ScoreFile ../$SOMName.score -Sparse
fi
../bin/map/mapsom -SOMFile ../$SOMName.som -SampleList $SampleList -Prefix ../$SOMName/data/som/
mv ../$SOMName/data/som/summery.map ../$SOMName/data/map_summery.map
../bin/units/getunits -ScoreFile ../$SOMName.score -Rows $Rows -Cols $Cols -Prefix ../$SOMName/data/som/units/unit
cp ../$SOMName.som ../$SOMName/data/out.som
