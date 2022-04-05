if [ $# -eq 0 ]    # Script invoked with no command-line args?
then
  echo "Usage: `basename $0` [required options] "
  echo "Required Options:"
  echo "-SOMName <SOM name>"
  echo "-TraitFile <Trait file location> "
  echo "-Metaclusters <Number of metaclusters> "
  echo "-Output <Trait report output file location.> "
  exit          # Exit and explain usage.
                            # Usage: scriptname -options
                            # Note: dash (-) necessary
fi

while (( "$#" ));
do
  case "$1" in
    -SOMName) SOMName=$2;;
    -TraitFile) TraitFile=$2;;
    -Metaclusters) Metaclusters=$2;;
    -Output) Output=$2;;
  esac

  shift
done

echo Rscript ../rscripts/Traits.R --TraitFile $TraitFile --MetaClusterFile ../$SOMName/data/MetaClusters --SOMFile ../$SOMName.som --ClusterNum $Metaclusters --OutputName $Output
Rscript ../rscripts/Traits.R --TraitFile $TraitFile --MetaClusterFile ../$SOMName/data/MetaClusters --SOMFile ../$SOMName.som --ClusterNum $Metaclusters --OutputName $Output
