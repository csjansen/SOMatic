# SOMatic
## Installing

Clone github folder

`git clone http://github.com/csjansen/SOMatic`

Enter SOMatic bin folder

`cd SOMatic/bin`

Compile SOMatic

`make`

If this doesn't work, be sure your g++ version is >2.8.2

`gcc -version`

## Tutorial: Using example data

### Requirements

A websever needs to be available for you to use.  Install and setup Apache if a webserver isn't available.  Know where you need to copy your website folder.  By default, it is /var/www/html.

### Unzip example data

From base SOMatic folder:

`cd examples`

`tar -zxf *.tgz`

### Run buildSite

`cd ../scripts`

`./buildSite.sh -SOMName Example -Matrix ../examples/example.matrix -Rows 20 -Cols 30 -SampleList ../examples/sample.list -Timesteps 4000000 -Trials 1`

`cp ../Example (webserver location)`

## Tutorial: RNA-seq data after RSEM quantification

### Prior Requirements

This tutoral assumes you have a number of RSEM outputs that are locatable by a regular expression i.e. in a folder together.  

### Build Training matrix and Sample List

`cd scripts`

`./rsemToTrainingMatrix_TPM.sh (regular expression for output files in quotes) (Sample List output location) (Training Matrix output location)`

`./rsemToTrainingMatrix_TPM.sh *.rsem.genes.results sample.list trainingMatrix`

This sample list file should be edited to give proper titles to all of your SOM maps.  Also, be sure that the sample names have no special characters as they can mess up the website. (. or - are fine)

### Decide which buildsite script to use

There are 2 buildSite scripts.  If you have a machine capable of making use of multithreading, such as having multiple cores, you should use buildSiteMT.sh, otherwise use buildSite.sh.  They make use of the same options.


### Run buildSite.sh
```
Usage: buildSite.sh [required options]
Required Options:
-SOMName <SOM name>
-Matrix <Training Matrix File Location>
-Rows <Number of rows you'd like in your SOM>
-Cols <Number of Columns you'd like in your SOM>
-SampleList <File with list of samples>
-Epochs <Number of Epochs for your SOM (Number of times that the trainer will be shown each segment; recommend 100 on RNA data and 10 on DNA data)>
-Trials <Number of trials you'd like to run.  The best SOM will be chosen (recommend 100, but fewer is fine for initial analysis)>
-Log2 <Log2 correct data>
```
From base SOMatic folder:

`cd ../scripts`

`./buildSite.sh -SOMName RNAdata -Matrix trainingMatrix -Rows 20 -Cols 30 -SampleList sample.list -Epochs 100 -Trials 5 -Log2`

### Copy Example SOM to your webserver and check for smoothness of Summary

`cp ../RNAdata (webserver location)`

Take your browser to your newly created website.  The first map that comes up should be your summary map.  Be sure that the units look smooth with their neighbors.  If they don't, be sure that you don't have regions or genes that are extemely out-of-scale or train again with more epochs.

The genes in each unit will be in the files in the RNAdata/data/som/units folder.

### Add GO terms overlay

#### Requirement files:

##### NCBI gene2go File:

If your organism is supported:

`wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz`

Format:

tax_id  GeneID  GO_ID   Evidence        Qualifier       GO_term PubMed  Category

(tab is used as a separator, pound sign - start of a comment)

If you're organism is not supported, you can create your own gene2go file.  The following fields are used by SOMatic:

tax_id: Species ID.  Must match the ID in the GeneInfo file.

GeneID: Gene ID.  Must match the IDs in the GeneInfo file.

GO_ID: GO ID.  These must match the IDs in the GO file.

GO_term: GO term description.  These must match the IDs in the GO File.

##### go obo file

`wget http://purl.obolibrary.org/obo/go.obo`

##### GeneInfo file
If your organism is supported, find your organism here:

ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/

And download the proper file.

Format:

tax_id  GeneID  Symbol  LocusTag        Synonyms        dbXrefs chromosome      map_location    description     type_of_gene    Symbol_from_nomenclature_authority      Full_name_from_nomenclature_authority   Nomenclature_status     Other_designations      Modification_date

(tab is used as a separator, pound sign - start of a comment)

If you're organism is not supported, you can create your own gene_info file.

The following fields are used by SOMatic:

tax_id: Species ID.  Must match the ID in the Gene2GO file

GeneID: Gene ID.  Must match the IDs in the Gene2GO file

Symbol: Gene Name.  Must match the gene names in the GTF file.

#### Process
```
Usage: getGORNA.sh [required options] [options]
Required Options:
-SOMName: SOM name
-Rows: Number of rows you'd like in your SOM
-Cols: Number of columns you'd like in your SOM
-Gene2GO: Gene2GO file.  See below for file format
-GeneInfo: Gene Info file.  See below for file format
-GOFile: OBO file from geneontology.org. http://geneontology.org/ontology/go.obo
Options: [choices] <default>
Sanity: If set to true, only GO terms with 5 genes in the unit will be
reported. [true, false] <true>
```

From the base SOMatic folder:

`cd scripts`

`./getGORNA.sh -SOMName RNAdata -Rows 20 -Cols 30 -Gene2GO gene2go -GeneInfo (your organism's gene info file) -GOFile go.obo`

Recopy your website to see go terms.

`cp ../RNAdata (your webserver)`

The GO tab should now be visible and GO terms will appear in units.

The files in RNAdata/data/som/GO will contain all of the GO terms.

## Tutorial: DNA data (such as ChIP or ATAC)

### Requirements

Mapped DNA sequencing experiments need to be in sam format and peaks need to be called for them in bed format.  The locations of these files need to be placed in 2 text files (1 experiment per line).  We'll name them bedFiles and samFiles in this tutorial.  A sample.list file also needs to be made that provides titles to each of the experiments in samFiles.  Be sure that these sample names have no special characters as it could ruin the website (. and - are fine).

### Partition the genome
```
Usage: ./partition.sh [options] -PeakDataFile <peak file list location> -Output <output file location>
Options: <default>
-MinFeature: Size of smallest partition. <200>
```

From base SOMatic folder:

`cd scripts`

`./partition.sh -PeakDataFile bedFiles -Output partition.list`

### Generate Training Matrix from RPKMs on the partitions
```
Usage: ./regionCounts.sh [options] -RawDataFile <raw sam file list location> -Partitions <partition file> -Output <output file location>
       Options: <default>"
	-LogScale: Log2(x+1) scale RPKM
```

`./regionCounts.sh -RawDataFile samFiles -Partitions partition.lish -Output DNAMatrix -LogScale`

This program makes a number of temporary files in case it crashes to restart where it left off.  You can remove them after this program finishes.

### Run buildSite
```
Usage: ./buildSite.sh [required options]
       Required Options:
       -SOMName <SOM name>
       -Matrix <Training Matrix File Location>
       -Rows <Number of rows you'd like in your SOM>
       -Cols <Number of Columns you'd like in your SOM>
       -SampleList <File with list of samples>
       -Epochs <Number of Epochs for your SOM (Number of times that the trainer will be shown each segment; recommend 100 on RNA data and 10 on DNA data)>
       -Trials <Number of trials you'd like to run.  The best SOM will be chosen (recommend 100, but fewer is fine for initial analysis)>
       -Log2 <Log2(x+1) correct data>
```

`./buildSite.sh -SOMName DNAdata -Matrix DNAMatrix -Rows 20 -Cols 30 -SampleList sample.list -Epochs 10 -Trials 5`

#### Copy your website to your web server

`cp ../DNAdata (your webserver)`

Take your browser to your newly created website.  The first map that comes up should be your summary map.  Be sure that the units look smooth with their neighbors.  If they don't, be sure that you don't have regions or genes that are extemely out-of-scale or train again with more epochs.

The partitions in each unit will be in the files in the RNAdata/data/som/units folder.


### Add Gene Overlay

#### Required Files

#### Process
```
Usage: ./getGenes.sh [required options] [options]
       Required Options:
       -SOMName: SOM name
       -Rows: Number of rows you'd like in your SOM
       -Cols: Number of columns you'd like in your SOM
       -GTFFile: Gene annotations file.  See README.txt for file format
       Options: [choices] <default>
       -Method: GREAT algorithm of choice. [TwoClosest,OneClosest] <OneClosest>
       -AddToChrom: If your gtf file uses a different format for it's chromosomes than your reference genome, this option allows you to add text to all the chromosomes in the gtf file. <>
```

`./getGenes.sh -SOMName DNAdata -Rows 20 -Cols 30 -GTFFile (Your GTF file)`

Recopy your website to see go terms.

`cp ../DNAdata (your webserver)`

The Genes tab should now be visible and Genes will appear in units.

The files in RNAdata/data/som/genes will contain all of the genes in each unit.

### Add GO Overlay

#### Requirement files:

##### NCBI gene2go File:

If your organism is supported:

`wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz`

Format:

tax_id  GeneID  GO_ID   Evidence        Qualifier       GO_term PubMed  Category

(tab is used as a separator, pound sign - start of a comment)

If you're organism is not supported, you can create your own gene2go file.  The following fields are used by SOMatic:

tax_id: Species ID.  Must match the ID in the GeneInfo file.

GeneID: Gene ID.  Must match the IDs in the GeneInfo file.

GO_ID: GO ID.  These must match the IDs in the GO file.

GO_term: GO term description.  These must match the IDs in the GO File.

##### go obo file

`wget http://purl.obolibrary.org/obo/go.obo`

##### GeneInfo file
If your organism is supported, find your organism here:

ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/

And download the proper file.

Format:

tax_id  GeneID  Symbol  LocusTag        Synonyms        dbXrefs chromosome      map_location    description     type_of_gene    Symbol_from_nomenclature_authority      Full_name_from_nomenclature_authority   Nomenclature_status     Other_designations      Modification_date

(tab is used as a separator, pound sign - start of a comment)

If you're organism is not supported, you can create your own gene_info file.

The following fields are used by SOMatic:

tax_id: Species ID.  Must match the ID in the Gene2GO file

GeneID: Gene ID.  Must match the IDs in the Gene2GO file

Symbol: Gene Name.  Must match the gene names in the GTF file.

#### Process
```
Usage: getGODNA.sh [required options] [options]
Required Options:
-SOMName: SOM name
-Rows: Number of rows you'd like in your SOM
-Cols: Number of columns you'd like in your SOM
-Gene2GO: Gene2GO file.  See below for file format
-GeneInfo: Gene Info file.  See below for file format
-GOFile: OBO file from geneontology.org. http://geneontology.org/ontology/go.obo
Options: [choices] <default>
Sanity: If set to true, only GO terms with 5 genes in the unit will be
reported. [true, false] <true>
```

From the base SOMatic folder:

`cd scripts`

`./getGODNA.sh -SOMName DNAdata -Rows 20 -Cols 30 -Gene2GO gene2go -GeneInfo (your organism's gene info file) -GOFile go.obo`

## Tutorial: Metaclustering (and GO/motif analysis)

### Requirements

The base buildSite.sh has to have been run.  It is assumed for this tutorial that the SOM name is RNAdata and that it was a 20x30 SOM.  This tutorial can be done right after the RNA-seq SOM tutorial above.  If you would like to do the tutorial on the example data, use DNAdata for the SOM name instead of RNAdata and ../examples/example.matrix instead of trainingMatrix.

### Heirarical cluster on SOM Units

From base SOMatic folder

`./cd scripts`

`./getClusters.sh -SOMName RNAdata`

### Calculate Dimensionality

`./getDimensionality.sh -SOMName RNAdata -CutLevel 40`

### Metacluster

`./metaClusterSOM.sh -SOMName RNAdata -Rows 20 -Cols 30 -Metaclusters 5 -MetaclustersEnd 50 -Trials 10`
This step is multithreaded, so if you have the option to provide multiple cores to this program, it can use up to the number of Trials you specified.

### Generate Metacluster reports and build cluster heatmaps for website

### Requirements

R with the following packages installed are required:
reshape2
ggplot2
ggdendro
grid
RColorBrewer
plyr

Also, Rscript must be a runable application.

### Process

`./generateMetaclusterReports.sh -SOMName RNAdata -Rows 40 -Cols 60 -Matrix trainingMatrix -ShowSegments 0 -OutputPrefix ../RNAdata-`

Creates a number of pdf files with the following output names "RNAdata-#.pdf", where # is the metacluster.  It also sets up files in the website to draw metacluster clustering heatmaps.  You need to re-copy your website to the webserver to access them.

### Check Metaclusters on SOMatic viewer

It is important at this stage to make sure that your metaclustering was done properly.  

#### Things to inspect
Did you get a number of metaclusters on the edge of your search space?  Aka did you get 5 or 50 metaclusters in this case?
Solution: Run the Metacluster step above with a different range.

Did you get a lot of single unit metaclusters?
Solution: Your SOM is too small, and one unit is trying to cover a big cluster on its own.  Re-run the SOM at a larger size.

Did one of your metaclusters go all of the way across the rows or columns of your SOM?
Solution: Your SOM is too large, and you are overclustering the differences between your observations.  Re-run the SOM at a smaller size.

## Tutorial: Metacluster downstream analysis

### Requirements
A SOM needs to have been trained with the metaclustering step done.  In this tutorial, we will assume it is RNAdata from above.

Also, R with the following packages installed are required:
reshape2
ggplot2
ggdendro
grid
RColorBrewer
plyr

Also, Rscript must be a runable application.

### Trait analysis
A trait descriptor file needs to be made with the following tab-delimited format:
(tab)	Trait#1	Trait#2
sample1	1	0
sample2	0	1

With this file made, the following can be run:
`./SOMMeta.sh -SOMName RNAData -TraitFile traits -Output ../RNAdata-Traits.pdf

This creates a PDF graph with significantly enriched or de-enriched metaclusters for the traits you specified.  These metaclusters could be analyzed further.

### Convert Metacluster output into a geneID format to calculate GO terms

The contents of the metaclusters are stored in files in the data folder of your website with the format: Genes_(Metacluster #).  For RNA, using cut in the proper way, you can remove everything from each row except for the geneID or gene name.  This file can be uploaded to PantherDB or David for GO analysis.For DNA, these files can be transformed into Bed files to input into GREAT to find GO terms for nearby genes or be further transformed to fasta files for motif analysis.  

## Tutorial: Linking RNA and DNA SOMs (and downstream GO/motif analysis)

### Requirements

Two SOMs need to have been trained: 1 from RNA data (RNAdata) and 1 from DNA data (DNAdata).  Both need to have been metaclustered as well.  A GTF file for your organism needs to be downloaded as well.

### Perform link

./Link.sh -SOMName1 DNAdata -Row1 40 -Col1 60 -SOMName2 RNAData -Rows2 40 -Col2 60 -OutputFolder SOMLinkage -GTFFile mm10.gtf 
