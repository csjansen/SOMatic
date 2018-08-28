# SOMatic
##Installing
Clone github folder
`github commit http://github.com/csjansen/SOMatic`

Enter SOMatic bin folder
`cd SOMatic/bin`

Compile SOMatic
`make`

If this doesn't work, be sure your g++ version is >2.8.2
`gcc -version`

##Tutorial: RNA-seq data after RSEM quantification
###Prior Requirements
This tutoral assumes you have a number of RSEM outputs that are locatable by a regular expression i.e. in a folder together.

###Build Training matrix and Sample List
`cd scripts`
`./rsemToTrainingMatrix (regular expression for output files in quotes) (Sample List output location) (Training Matrix output location)`
`./rsemToTrainingMatrix *.rsem.genes.results sample.list trainingMatrix`

###Decide which buildsite script to use
There are 2 buildSite scripts.  If you have a machine capable of making use of multithreading, such as having multiple cores, you should use buildSiteMT.sh, otherwise use buildSite.sh.  They make use of the same options.


###Run buildSite.sh
Usage: buildSite.sh [required options]
Required Options:
-SOMName <SOM name>
-Matrix <Training Matrix File Location>
-Rows <Number of rows you'd like in your SOM>
-Cols <Number of Columns you'd like in your SOM>
-SampleList <File with list of samples>
-Timesteps <Number of timesteps for your SOM>
-Trials <Number of trials you'd like to run.  The best SOM will be chosen>

From base SOMatic folder:
`cd examples`
`tar -zxf *.tgz`
`cd ../scripts`
`./buildSite.sh -SOMName Example -Matrix ../examples/example.matrix -Rows 20 -Cols 30 -SampleList ../examples/sample.list -Timesteps 1000000 -Trials 5`

###Copy Example SOM to your webserver
`cp ../Example 
