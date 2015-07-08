Getting Started:

This short tutorial will explain how to use SOMatic to build a Self-Organizing
Map of your data and set up the visualization platform to explore your data
set.

Install:
To build SOMatic, be sure that the gcc installed on your system has a version
>2.8.2.  You can check your version of gcc by running gcc --version.  Once you
are sure move into the bin directory and run the following:
make

Building a SOM and website:
To easily build a SOM and it's cooresponding website, see the File Format
section below for a description of the required training matrix and sample list 
file formats and, then, move into the scripts folder and run the following script:

Usage: buildsite.sh [required options] 
Required Options:
-SOMName <SOM name>
-Matrix <Training Matrix File Location>
-Rows <Number of rows you'd like in your SOM>
-Cols <Number of Columns you'd like in your SOM>
-SampleList <File with list of samples>
-Timesteps <Number of timesteps for your SOM>
-Trials <Number of trials you'd like to run.  The best SOM will be chosen>

Populating your website with gene names if your segments are in genomic
coordinates (chrN:X-Y):

If your training file is partitioned over genomic coordinates and you'd like
to overlay gene information onto your SOM, see the File Format section below
for a description of the required GTF file, move into the scripts folder and
run the following script.  This script requires that the SOM viewer built
above is in the SOMatic top level folder (where is was orgininally placed).

Usage: getgenes.sh [required options] [options]
Required Options:
-SOMName: SOM name
-Rows: Number of rows you'd like in your SOM 
-Cols: Number of columns you'd like in your SOM 
-GTFFile: Gene annotations file.  See below for file format.
Options: [choices] <default>
-Method: GREAT algorithm of choice. [TwoClosest] <TwoClosest>
-AddToChrom: If your gtf file uses a different format for it's chromosomes,
			 this option allows you to add text to all the chromosomes in the 
			 gtf file. <>

Populating your website with GO Terms (If segments are in genomic
coordinates):

If your training file is partitioned over genomic coordinates and you'd like
to overlay GO term information onto your SOM, see the File Format section below
for a description of the required Gene2GO and GeneInfo files, move into the 
scripts folder and run the following script.  This script requires that the 
SOM viewer built above is in the SOMatic top level folder (where is was 
orgininally placed).

Usage: getGOGenomic.sh [required options] [options]
Required Options:
-SOMName: SOM name
-Rows: Number of rows you'd like in your SOM 
-Cols: Number of columns you'd like in your SOM 
-Gene2GO: Gene2GO file.  See below for file format 
-GeneInfo: gene_info file.  See below for file format
-GOFile: OBO file from geneontology.org. http://geneontology.org/ontology/go.obo
Options: [choices] <default>
Sanity: If set to true, only GO terms with 5 genes in the unit will be reported. 
[true, false] <true>

Populating your website with GO Terms (If segments are genes):

If your training file is partitioned by genes and you'd like to overlay GO term 
information onto your SOM, see the File Format section below for a description 
of the required Gene2GO and GeneInfo files, move into the scripts folder and 
run the following script.  This script requires that the SOM viewer built above 
is in the SOMatic top level folder (where is was orgininally placed).

Usage: getGOGene.sh [required options] [options]
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

File Formats:
Training Matrix:
This file contains all of the data that will be trained on to generate a SOM.

Format(seperated by tabs):
SegmentName	Sample1RPKM	Sample2RPKM	...	SampleNRPKM

Examples:
Genes:
TSPAN6  0.046060667	0.033884333 0.00147032  0  0   0   0.0966804   0
TNMD    0	0 0.014867333 0.00292713  0.0112718   0.00291146  0.001928637 0
DPM1    17.77493333 37.41216667 46.9975	32.53276667 40.19216667 43.59976667 0	0

Genome Coordinates:
chr4:1003422-1249035	0.046060667	0.033884333 0.00147032  0  0   0   0.0966804   0
chr4:1249036-1444424	0	0 0.014867333 0.00292713  0.0112718   0.00291146  0.001928637 0
chr4:1444425-1672431	17.77493333 37.41216667 46.9975 32.53276667 40.19216667 43.59976667 0

Warning:
Segments with RPKMs of 0 across the whole row should be removed before
running.  Failure to do this will cause the SOM to be not useful.

There is an example.matrix file included in the SOMatic top level folder.

Sample List:
Format(1 sample per line):
Sample1Name
Sample2Name
...
SampleNName

Example:
HL60
3hMac
6hMac
12hMac
24hMac

Warning:
These names will be the names refered to on your website.  Don't make the
names too long or too obscure, or navigation of your data may be difficult.

GTF File:
Format(Standard GTF File format:
http://www.ensembl.org/info/website/upload/gff.html

Warning:
The names in the gene_name field will be used for the gene overlay.  Be sure
that the gene names in your GeneInfo file are the same as these.

Gene2GO File:
NCBI gene2go File:
ftp://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz

Format:
tax_id	GeneID	GO_ID	Evidence	Qualifier	GO_term	PubMed	Category 
(tab is used as a separator, pound sign - start of a comment)

Used Fields:
If you're organism is not supported, you can create your own gene2go file.
The following fields are used by SOMatic:
tax_id:	Species ID.  Must match the ID in the GeneInfo file.
GeneID: Gene ID.  Must match the IDs in the GeneInfo file.
GO_ID: GO ID.  These must match the IDs in the GO file.
GO_term: GO term description.  These must match the IDs in the GO File.

Warning:
Only the child nodes should be included.  The getGO program will build the
rest of the heirarchy from the GO file.

GeneInfo File:
Location of NCBI gene_info files:
ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/

Format:
tax_id	GeneID	Symbol	LocusTag	Synonyms	dbXrefs	chromosome	map_location	description	type_of_gene	Symbol_from_nomenclature_authority	Full_name_from_nomenclature_authority	Nomenclature_status	Other_designations	Modification_date 
(tab is used as a separator, pound sign - start of a comment)

Used Fields:
If you're organism is not supported, you can create your own gene_info file.
The following fields are used by SOMatic:
tax_id: Species ID.  Must match the ID in the Gene2GO file
GeneID: Gene ID.  Must match the IDs in the Gene2GO file
Symbol: Gene Name.  Must match the gene names in the GTF file.

GO File:
geneontology.org file:
http://purl.obolibrary.org/obo/go.obo
