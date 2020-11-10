#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <vector>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>

using namespace std;

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

class TSSsite {
public:
	string chrom;
	int pos;
	string strand;
	bool protein;
	TSSsite() {
		chrom = "";
		pos = 0;
		strand = "";
	}
};

class genomicRegion {
public:
	string gene;
	string chrom;
	int start;
	int stop;
};

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
        std::string item;
        while (std::getline(ss, item, delim)) {
            elems.push_back(item);
        }
        return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
        split(s, delim, elems);
        return elems;
}


map<string, TSSsite>* parseGtfFile(string gtfFileName, string geneidtype, bool Xeno) {
	ifstream gtffile(gtfFileName.c_str());
    map<string, TSSsite>* TSSsites=new map<string, TSSsite>();
	string line;
	int counter=0;
	string temptext = "";
	while(getline(gtffile,line)) {
		if(!Xeno) {
			if(line[0]=='#') continue;
			vector<string> splitz = split(line, '\t');
			if(splitz[2].compare("exon")==0) {
				string geneName;
				vector<string> attributes = split(splitz[8], ';');
				bool protein_coding=false;
				bool protein_coding_gene=false;
				for(int i = 0; i < attributes.size(); i++) {
					vector<string> pairItems = split(attributes[i],' ');
//					cout<<geneidtype<<endl;
					if(geneidtype.compare("gene_id")==0) {
//						cout<<pairItems[0]<<'\t'<<pairItems[1]<<endl;
						if(pairItems[0].compare("gene_id")==0) {
							temptext = pairItems[1].substr(1,pairItems[1].size()-2);
							vector<string> temptextlist = split(temptext,'.');
							temptext = temptextlist[0];
						}
						if(pairItems[1].compare("gene_name")==0) {
							//geneName = temptext+"-"+pairItems[2].substr(1,pairItems[2].size()-2);
							geneName = temptext;
						}
					} else if(geneidtype.compare("gene_name")==0) {
					
						if(pairItems[1].compare("gene_name")==0) {
							vector<string> splitz2 = split(pairItems[2].substr(1,pairItems[2].size()-2),'.');
							geneName = pairItems[2].substr(1,pairItems[2].size()-2);
//						cout<<geneName<<endl;
//						int temp;
//						cin>>temp;
						}
					}
					if(pairItems[1].compare("gene_type")==0) {
						string splitz2 = pairItems[2];
						if(splitz2.compare("\"protein_coding\"")==0) protein_coding_gene=true;
					}
					if(pairItems[1].compare("transcript_type")==0) {
						string splitz2 = pairItems[2];
					//cout<<splitz2<<endl;
						if(splitz2.compare("\"protein_coding\"")==0) protein_coding=true;
					//if(protein_coding) cout<<"coding"<<endl;
					//int temp;
					//cin>>temp;
					}
				}
//							cout<<geneName<<endl;
//					int temp;
//					cin>>temp;
				TSSsite record = TSSsites->operator[](geneName);
				if(!protein_coding&&protein_coding_gene) continue;
				record.strand = splitz[6];
				record.chrom = splitz[0];
				int pos;
				if(record.strand.compare("+")==0) {
					istringstream(splitz[3])>>pos;
					if(pos < record.pos || record.pos == 0) record.pos = pos;
				} else {
					istringstream(splitz[4])>>pos;
					if(pos > record.pos || record.pos == 0) record.pos = pos;
				}
				TSSsites->operator[](geneName) = record;
				counter++;
			}
		} else {
			if(line[0]=='#') {
				continue;
			}
			vector<string> splitz = split(line,'\t');
			if(splitz[2].compare("gene")==0) {
				TSSsite temp;
				string strand = splitz[6];
				int start;
				int stop;
				istringstream(splitz[3])>>start;
				istringstream(splitz[4])>>stop;
				temp.chrom = splitz[0];
				//cout<<temp.chrom<<endl;
				//cout<<start<<endl;
				//cout<<stop<<endl;
				if(strand.compare("+")==0) {
					temp.pos = start;
				} else {
					temp.pos=stop;
				}
				vector<string> splitz1 = split(splitz[8],';');
				vector<string> splitz2 = split(splitz1[0],'=');
				vector<string> splitz3 = split(splitz1[1],'=');
				string genename =splitz2[1]+"-"+splitz3[1];
				cout<<genename<<endl;
				TSSsites->operator[](genename) = temp;
				
			}
			/*if(counter>0) {
				vector<string> splitz = split(line,'\t');
				TSSsite temp;
				string genename = splitz[0];
				if(genename.compare("unnamed")==0) {
					genename = splitz[0];
				}
				vector<string> chromsplitz = split(splitz[2],':');
				vector<string> startsplitz = split(chromsplitz[1],'-');
				temp.chrom = chromsplitz[0];
				istringstream(startsplitz[0])>>temp.pos;
				temp.strand = '+';
				TSSsites->operator[](genename) = temp;
			}
			counter++;	*/
		}
	}
	cout<<"GTF Genes: "<<TSSsites->size()<<endl;
	return TSSsites;
}

vector <vector<int> > hexSurround(vector<double> input, int radius, int numRows, int numCols) {
    vector<vector<int> > result;
    vector<int> xs;
    vector<int> ys;
    vector<int> zs;
    int x = input[1] - (input[0] - (abs((int)input[0])%2))/2;
    int z = input[0];
    int y = -x-z;
    for(int i = -1*radius; i <= radius; i++) {
        for(int j = max(-1*radius,-i-radius); j <= min(radius,-i+radius);j++) {
            int dz = -i-j;
            xs.push_back(x+i);
            ys.push_back(y+j);
            zs.push_back(z+dz);
            //cout<<i<<'\t'<<j<<'\t'<<dz<<endl;
        }
    }
    for(int i = 0; i < xs.size(); i++) {
        vector<int> temp;
        temp.push_back(zs[i]);
        if(temp[0]<0) temp[0]+=numRows;
        if(temp[0]>=numRows) temp[0]-=numRows;
        temp.push_back(xs[i]+(zs[i]-(abs(zs[i])%2))/2);
        if(temp[1]<0) temp[1]+=numCols;
        if(temp[1]>=numCols) temp[1]-=numCols;
        result.push_back(temp);
    }
    return result;
}

string createKey(int row1, int col1, int row2, int col2) {
	return SSTR(row1)+"_"+SSTR(col1)+"_"+SSTR(row2)+"_"+SSTR(col2);
}

vector<genomicRegion> GetRegRegions(map<string, TSSsite>* TSSsites, vector<string> genes, string CompareType, int SearchRange,bool underscore) {
	vector<genomicRegion> regions;
	//cout<<genes.size()<<" Passed in."<<endl;
	for(int i = 0; i < genes.size(); i++) {
		vector<string> splitz;
		if(!underscore) {
			splitz = split(genes[i],'-');
		} else {
			splitz = split(genes[i],'_');
		}
		map<string, TSSsite>::iterator it = TSSsites->find(splitz[0]);
		TSSsite TSS;
		if(it!=TSSsites->end()) {
			TSS = it->second;
		} else {
			cout<<"Could not find: "<<splitz[0]<<".  Skipping..."<<endl;
			continue;
		}
		//if(genes[i].compare("ENSMUSG00000075370")==0||genes[i].compare("ENSMUSG00000059305")==0) cout<<genes[i]<<'\t'<<TSS.pos<<endl;
		int start;
		int stop;
		if(CompareType.compare("TwoClosest")==0) {
			int prevdist = SearchRange;
			if(TSS.pos < SearchRange) start = 0;
			else start = TSS.pos - SearchRange;
			int nextdist = SearchRange;
			if(TSS.pos < SearchRange) stop = 0;
            else stop = TSS.pos + SearchRange;
			//cout<<start<<'\t'<<stop;
			for(map<string, TSSsite>::iterator looper = TSSsites->begin(); looper != TSSsites->end(); looper++) {
				TSSsite tester = looper->second;
				if(TSS.chrom.compare(tester.chrom)==0 && ((TSS.pos<tester.pos&&tester.pos-TSS.pos<nextdist)||(TSS.pos>tester.pos&&TSS.pos-tester.pos<prevdist))) {
					if(TSS.pos<tester.pos) {
						stop = tester.pos;
						nextdist = tester.pos-TSS.pos;
					} else {
						start = tester.pos;
						prevdist = TSS.pos - tester.pos;
					}
				}
			}
			genomicRegion temp;
			temp.gene = genes[i];
			temp.chrom = TSS.chrom;
			temp.start = start;
			temp.stop = stop;
			regions.push_back(temp);
		} else if (CompareType.compare("OneClosest")==0) {
			float distance1=SearchRange;
			float distance2=SearchRange;
			for(map<string, TSSsite>::iterator looper = TSSsites->begin(); looper != TSSsites->end(); looper++) {
				TSSsite tester = looper->second;
				if(TSS.chrom.compare(tester.chrom)!=0) continue;
				if(TSS.pos < tester.pos && tester.pos-TSS.pos < distance1) distance1=tester.pos-TSS.pos;
				if(TSS.pos > tester.pos && TSS.pos-tester.pos < distance2) distance2=TSS.pos-tester.pos;
			}
			if(distance2!=SearchRange)
				start = TSS.pos-distance2/2;
			else
				start = TSS.pos-distance2;
			if(distance1!=SearchRange)
				stop = TSS.pos+distance1/2;
			else
				stop = TSS.pos+distance1;
			genomicRegion temp;
			temp.gene = genes[i];
			temp.chrom = TSS.chrom;
			temp.start = start;
			temp.stop = stop;
			regions.push_back(temp);
//			cout<<temp.gene<<endl;
//			cout<<temp.chrom<<endl;
//			cout<<temp.start<<endl;
//			cout<<temp.stop<<endl;
		}
	}
//	int temp;
//	cin>>temp;
	return regions; 
}

int hexdist(int row1, int col1, int row2, int col2, int rows, int cols) {
    row1=row1-row2;
    col1=col1-col2;
    row2=0;
    col2=0;
    int x1 = col1 - (row1 - (abs(row1)%2))/2;
    int z1 = row1;
    int y1 = -x1-z1;
    int x2 = col2 - (row2 - (abs(row2)%2))/2;
    int z2 = row2;
    int y2 = -x2-z2;
    return (abs(x1-x2) + abs(y1-y2) + abs(z1-z2))/2;
}


int main(int argc, char* argv[]) {
	if(argc < 2) {
        cout << "Usage: ./LinkUnits -UnitPrefix1 <Prefix of Unit Files from First SOM> -Row1 <Rows from First SOM> -Col1 <Cols from First SOM> -UnitPrefix2 <Prefix of Unit Files from Second SOM> -Rows2 <Rows from Second SOM> -Col2 <Cols from Second SOM> -Output <Output File Location> [options]"<<endl;
        cout << "Options: <default>[options]" <<endl;
        cout << "-Type: What type of data is being compared.  First SOM comes first part of type. e.g. ATACxRNA means that the first SOM is ATAC and the second SOM is RNA.<ATACxRNA>[ATACxRNA]"<<endl;
        cout << "ATACxRNA options:\n\t-Algorithm: Which GREAT algorithm is used. <OneClosest>[OneClosest,TwoClosest]"<<endl;
        cout << "\t-GTFFile: GTF File for your organism."<<endl;
        cout << "\t-GeneIDType: Type of gene ID in RNA SOM unit files. <gene_id>[gene_id,gene_name]"<<endl;
	cout << "\t-SearchRange: The distance from a gene to use for GREAT <1000000>"<<endl;
        return 0;
    }
    string type = "ATACxRNA";
	string CompareType="OneClosest";
	string geneidtype="gene_id";
	string prefix1;
	int row1;
	int row2;
	string prefix2;
	int col1;
	int col2;
	string outfileName;
	string gtfFileName;
	bool Xeno = false;
	int SearchRange=1000000;
	bool addChr = false;
	bool underscore=false;
	for(int i = 0; i < argc; i++) {
        	string temp = argv[i];
	        if(temp.compare("-UnitPrefix1")==0)
			prefix1= argv[i+1];	
        	if(temp.compare("-Row1")==0)
	        	istringstream(argv[i+1])>>row1;
        	if(temp.compare("-Row2")==0)
            		istringstream(argv[i+1])>>row2;
	        if(temp.compare("-UnitPrefix2")==0)
			prefix2= argv[i+1];	
        	if(temp.compare("-Col1")==0)
        		istringstream(argv[i+1])>>col1;
        	if(temp.compare("-Col2")==0)
        		istringstream(argv[i+1])>>col2;
        	if(temp.compare("-Output")==0)
			outfileName = argv[i+1];	
        	if(temp.compare("-Algorithm")==0)
			CompareType = argv[i+1];	
        	if(temp.compare("-Type")==0)
			type = argv[i+1];	
        	if(temp.compare("-GTFFile")==0)
			gtfFileName = argv[i+1];	
        	if(temp.compare("-GeneIDType")==0)
			geneidtype = argv[i+1];	
		if(temp.compare("-Xeno")==0)
			Xeno=true;
		if(temp.compare("-SearchRange")==0) 
			istringstream(argv[i+1])>>SearchRange;
		if(temp.compare("-AddChr")==0)
			addChr = true;
		if(temp.compare("-Underscore")==0)
			underscore = true;
	}

	cout<<prefix1<<" rows: "<<row1<<" cols: "<<col1<<endl;
	cout<<prefix2<<" rows: "<<row2<<" cols: "<<col2<<endl;
	ofstream outfile(outfileName.c_str());
	if(type.compare("RNAxRNA")==0) {
		cout<<"RNAxRNA"<<endl;
		vector<vector<vector<vector<string> > > > RNAOverlaps;
		vector<vector<vector<string> > > genes1;
		for(int RNARow1 = 0; RNARow1 < row1; RNARow1++) {
			vector<vector<string> > geneRow;
			for(int RNACol1 = 0; RNACol1 < col1; RNACol1++) {
				ifstream RNAUnit1((prefix1+"_"+SSTR(RNARow1)+"_"+SSTR(RNACol1)+".unit").c_str());
				string line;
				vector<string> geneCol;
				while(getline(RNAUnit1,line)) {
					vector<string> splitz = split(line,'\t');
					geneCol.push_back(splitz[0]);
				}
				RNAUnit1.close();
				geneRow.push_back(geneCol);
			}
			genes1.push_back(geneRow);
		}
		vector<vector<vector<string> > > genes2;
                for(int RNARow2 = 0; RNARow2 < row2; RNARow2++) {
                        vector<vector<string> > geneRow;
                        for(int RNACol2 = 0; RNACol2 < col2; RNACol2++) {
                                ifstream RNAUnit2((prefix2+"_"+SSTR(RNARow2)+"_"+SSTR(RNACol2)+".unit").c_str());
                                string line;
                                vector<string> geneCol;
                                while(getline(RNAUnit2,line)) {
                                        vector<string> splitz = split(line,'\t');
                                        geneCol.push_back(splitz[0]);
                                }
                                RNAUnit2.close();
                                geneRow.push_back(geneCol);
                        }
                        genes2.push_back(geneRow);
                }

		for(int RNARow1 = 0; RNARow1 < row1; RNARow1++) {
			for(int RNACol1 = 0; RNACol1 < col1; RNACol1++) {
				outfile<<RNARow1<<'\t'<<RNACol1;
				for(int RNARow2 = 0; RNARow2 < row2; RNARow2++) {
					for(int RNACol2 = 0; RNACol2 < col2; RNACol2++) {
						vector<string> matches;
						for(int i = 0; i < genes1[RNARow1][RNACol1].size(); i++) {
							for(int j = 0; j < genes2[RNARow2][RNACol2].size(); j++) {
								if(genes1[RNARow1][RNACol1][i].compare(genes2[RNARow2][RNACol2][j])==0) {
									matches.push_back(genes1[RNARow1][RNACol1][i]);
								}
							}
						}	
						if(matches.size() >0) {
							outfile<<'\t'<<RNARow2<<'\t'<<RNACol2<<'\t'<<matches.size();
							for(int i = 0; i < matches.size(); i++) {
								outfile<<'\t'<<matches[i]<<'\t'<<matches[i];
							}
						}
					}
				}
				outfile<<endl;
			}
		}
		outfile.close();
	}
	if(type.compare("ATACxRNA")==0) {
		cout<<"ATACxRNA"<<endl;
		cout<<gtfFileName<<endl;
		map<string, TSSsite>* TSSsites = parseGtfFile(gtfFileName, geneidtype,Xeno);
		vector<vector<int> > AtacSizes;
		vector<vector<int> > RNASizes;
		vector<vector<genomicRegion> > regions;
		vector<vector<vector<vector<string> > > > TotalAtacOverlaps;
		vector<vector<vector<vector<string> > > > TotalRnaOverlaps;
		vector<vector<genomicRegion> > TotalAtacRegions;
		for(int AtacRow = 0; AtacRow < row1; AtacRow++) {
			//cout<<"Atac Row: "<<AtacRow<<endl;
			vector<int> AtacSizeRow;
			for(int AtacCol = 0; AtacCol < col1; AtacCol++) {
				//cout<<"Atac Col: "<<AtacCol<<endl;
				vector<genomicRegion> AtacRegions;
				cout<<"Opening "<<prefix1+"_"+SSTR(AtacRow)+"_"+SSTR(AtacCol)+".unit"<<endl;
				
				ifstream AtacUnit((prefix1+"_"+SSTR(AtacRow)+"_"+SSTR(AtacCol)+".unit").c_str());
				string line;
				while(getline(AtacUnit,line)) {
//					cout<<line<<endl;
					genomicRegion temp;
					vector<string> frist = split(line,'\t');
					vector<string> splitz = split(frist[0],':');
					if(splitz.size() == 1) continue;
					temp.chrom = splitz[0];
					vector<string> splitz2 = split(splitz[1],'-');
					if(splitz2.size() > 1) {
						istringstream(splitz2[0])>>temp.start;
						istringstream(splitz2[1])>>temp.stop;
					} else {
					//cout<<temp.start<<'\t'<<temp.stop<<'\t'<<line<<endl;
						istringstream(splitz[1])>>temp.start;
						istringstream(splitz[2])>>temp.stop;
					}
					AtacRegions.push_back(temp);
				}
				vector<vector<vector<string> > > AtacOverlaps;
				vector<vector<vector<string> > > RNAOverlaps;
				//cout<<"Regions: "<<AtacRegions.size()<<endl;
				AtacSizeRow.push_back(AtacRegions.size());
				for(int RnaRow = 0; RnaRow < row2; RnaRow++) {
					//cout<<"RNA row: "<<RnaRow<<endl;
					vector< vector<string> > AtacOverRow;
					vector< vector<string> > RNAOverRow;
					vector<int> RNASizeRow;
					for(int RnaCol = 0; RnaCol < col2; RnaCol++) {
						//cout<<"RNA col: "<<RnaCol<<endl;
				//		cout<<"Opening "<<prefix2+"_"+SSTR(RnaRow)+"_"+SSTR(RnaCol)+".unit"<<endl;
						ifstream RnaUnit((prefix2+"_"+SSTR(RnaRow)+"_"+SSTR(RnaCol)+".unit").c_str());
						//cout<<"Opened: "<<(prefix2+"_"+SSTR(RnaRow)+"_"+SSTR(RnaCol)+".unit")<<endl;
						string line;
						vector<string> genes;
						vector<string> genes2;
						while(getline(RnaUnit, line)) {
							vector<string> splitz=split(line,'\t');
						//	cout<<splitz[0]<<endl;
							vector<string> splitz2 = split(splitz[0],'-');
							vector<string> splitz3 = split(splitz2[0],'.');
						//	cout<<splitz2[0]<<endl;
							genes.push_back(splitz3[0]);
							
						}
						RNASizeRow.push_back(genes.size());
						//cout<<"Genes: "<<genes.size()<<endl;
						if(regions.size() < row2 * col2) {
							//cout<<"Getting Reg Regions"<<endl;
							regions.push_back(GetRegRegions(TSSsites, genes, CompareType,SearchRange,underscore));
							//cout<<"Gotten"<<endl;
						}
						vector<string> overlaps;
						vector<string> genomicOverlaps;
						//cout<<"Calculating overlaps"<<endl;
						//cout<<AtacRegions.size()<<'\t'<<regions[RnaRow*col2+RnaCol].size()<<endl;
						for(int Atacs = 0; Atacs < AtacRegions.size(); Atacs++) {
							for(int RNAs = 0; RNAs < regions[RnaRow*col2+RnaCol].size(); RNAs++) {
					//			cout<< AtacRegions[Atacs].start<<'\t'<<regions[RnaRow*col2+RnaCol][RNAs].start<<endl;
								int temp = 0;
								if(addChr)
									if(regions[RnaRow*col2+RnaCol][RNAs].chrom[0]!='c')
										regions[RnaRow*col2+RnaCol][RNAs].chrom="chr"+regions[RnaRow*col2+RnaCol][RNAs].chrom;
								//cout<<AtacRegions[Atacs].chrom<<'\t'<<regions[RnaRow*col2+RnaCol][RNAs].chrom<<endl;
								//cin>>temp;
								if(AtacRegions[Atacs].chrom.compare(regions[RnaRow*col2+RnaCol][RNAs].chrom)==0) {
									if((AtacRegions[Atacs].start>=regions[RnaRow*col2+RnaCol][RNAs].start && AtacRegions[Atacs].start<=regions[RnaRow*col2+RnaCol][RNAs].stop)||(AtacRegions[Atacs].start>=regions[RnaRow*col2+RnaCol][RNAs].stop && AtacRegions[Atacs].stop <= regions[RnaRow*col2+RnaCol][RNAs].stop)||(AtacRegions[Atacs].start<=regions[RnaRow*col2+RnaCol][RNAs].start && AtacRegions[Atacs].stop >= regions[RnaRow*col2+RnaCol][RNAs].stop)) {
										int midDist = abs((AtacRegions[Atacs].stop-AtacRegions[Atacs].start)-(regions[RnaRow*col2+RnaCol][RNAs].stop-regions[RnaRow*col2+RnaCol][RNAs].start));
										overlaps.push_back(regions[RnaRow*col2+RnaCol][RNAs].gene);
										genomicOverlaps.push_back(AtacRegions[Atacs].chrom+":"+SSTR(AtacRegions[Atacs].start)+"-"+SSTR(AtacRegions[Atacs].stop));
									}
								}
							}
						}
						//cout<<overlaps.size()<<endl;
						//outfile<<AtacRow<<"\t"<<AtacCol<<"\t"<<RnaRow<<"\t"<<RnaCol<<"\t";
						if(AtacRegions.size()!=0)
							AtacOverRow.push_back(genomicOverlaps);
						else {
							vector<string> temp;
							AtacOverRow.push_back(temp);
						}
						if(regions[RnaRow*col2+RnaCol].size()!=0)
							RNAOverRow.push_back(overlaps);
						else {
								vector<string> temp;
								RNAOverRow.push_back(temp);
							}
							
						//cout<<"Closing RNA unit"<<endl;
						RnaUnit.close();
						//cout<<"Closed"<<endl;
					}
					AtacOverlaps.push_back(AtacOverRow);
					RNAOverlaps.push_back(RNAOverRow);
					if(RNASizes.size() <= row2) {
						RNASizes.push_back(RNASizeRow);
					}
					//cout<<"Down with RNA Row"<<endl;
				}
			
				//int intemp;
				//cin>>intemp;
				TotalAtacOverlaps.push_back(AtacOverlaps);
				TotalRnaOverlaps.push_back(RNAOverlaps);
				TotalAtacRegions.push_back(AtacRegions);
			}
		}
		for(int AtacRow = 0; AtacRow < row1; AtacRow++) {
            for(int AtacCol = 0; AtacCol < col1; AtacCol++) {
				outfile<<AtacRow<<'\t'<<AtacCol;
				int printamount=0;
				int totalamount=0;
                for(int AtacOverRows = 0; AtacOverRows<TotalAtacOverlaps[AtacRow*col1+AtacCol].size(); AtacOverRows++)
                    for(int AtacOverCols = 0; AtacOverCols<TotalAtacOverlaps[AtacRow*col1+AtacCol][AtacOverRows].size(); AtacOverCols++) {
			if (TotalAtacOverlaps[AtacRow*col1+AtacCol][AtacOverRows][AtacOverCols].size()>0) {
							vector<double> input;
							input.push_back(AtacRow);
							input.push_back(AtacCol);
							vector <vector<int> > surroundinghexs = hexSurround(input, 1, row1, col1);
							int total = 0;
							for(int i = 0; i < surroundinghexs.size(); i++) {
								total+=TotalAtacOverlaps[surroundinghexs[i][0]*col1+surroundinghexs[i][1]][AtacOverRows][AtacOverCols].size();
							}
							bool success=false;
							if(total>1) {
								success=true;
							} else {
								for(int AtacOverRows2 = 0; AtacOverRows2<TotalAtacOverlaps[AtacRow*col1+AtacCol].size(); AtacOverRows2++) {
				                    for(int AtacOverCols2 = 0; AtacOverCols2<TotalAtacOverlaps[AtacRow*col1+AtacCol][AtacOverRows].size(); AtacOverCols2++) {
										if (TotalAtacOverlaps[AtacRow*col1+AtacCol][AtacOverRows][AtacOverCols].size()>0 && !(AtacOverRows2==AtacOverRows&&AtacOverCols2==AtacOverCols)) {
											int dist = hexdist(AtacOverRows2, AtacOverCols2, AtacOverRows, AtacOverCols, row2, col2);
											if(dist==1) {
												success=true;
												break;	
											}
										}
									}
									if(success) break;
								}
							}
							if(success) {
								outfile<<'\t'<<AtacOverRows<<'\t'<<AtacOverCols<<'\t'<<TotalAtacOverlaps[AtacRow*col1+AtacCol][AtacOverRows][AtacOverCols].size();
                                printamount+= TotalAtacOverlaps[AtacRow*col1+AtacCol][AtacOverRows][AtacOverCols].size();
                                for(int genes = 0; genes < TotalAtacOverlaps[AtacRow*col1+AtacCol][AtacOverRows][AtacOverCols].size(); genes++) {
                                    outfile<<'\t'<<TotalRnaOverlaps[AtacRow*col1+AtacCol][AtacOverRows][AtacOverCols][genes]<<'\t'<<TotalAtacOverlaps[AtacRow*col1+AtacCol][AtacOverRows][AtacOverCols][genes];
                                }
							}
							totalamount+= TotalAtacOverlaps[AtacRow*col1+AtacCol][AtacOverRows][AtacOverCols].size();
						}
                    }
                //cout<<"AtacUnit closing"<<endl;
			outfile<<'\t'<<printamount<<'\t'<<totalamount<<endl;
            }
		}        
		outfile.close();
	}
} 
