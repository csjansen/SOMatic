if [ $# -eq 0 ]    # Script invoked with no command-line args?
then
  echo "Usage: `basename $0` [required options] "
  echo "Required Options:"
  echo "-SOMName <SOM name>"
  echo "-Rows <Number of Rows in your SOM> "
  echo "-Cols <Number of Columns in your SOM> "
  echo "-Metaclusters <Number of metaclusters you would like to try.> "
  echo "-TrainingMatrix <Original Training Matrix used in SOM creation.> "
  echo "-ShowSegments <Turn off if metaclusters are huge> [0,1] {1}"
  echo "-OutputPrefix <Report output file location.> "
  exit          # Exit and explain usage.
                            # Usage: scriptname -options
                            # Note: dash (-) necessary
fi

ShowSegments=1
while (( "$#" ));
do
  case "$1" in
    -SOMName) SOMName=$2;;
    -Rows) Rows=$2;;
    -Cols) Cols=$2;;
    -Metaclusters) Metaclusters=$2;;
    -TrainingMatrix) TrainingMatrix=$2;;
    -OutputPrefix) OutputPrefix=$2;;
	-ShowSegments) ShowSegments=$2;;
  esac

  shift
done


Rscript ../rscripts/MetaClusterReportsBar.R --MetaClusterFile ../$SOMName/data/MetaClusters --SOMFile ../$SOMName.som --ClusterNum $Metaclusters --SampleList ../$SOMName/data/sample.list --TrainingMatrix $TrainingMatrix --GeneFilePrefix ../$SOMName/data --OutputPrefix $OutputPrefix --TypeFile ../$SOMName/data/types
