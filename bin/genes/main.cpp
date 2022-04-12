/* getgenes: converts genomic regions into genes using GREAT algorithms
    Copyright (C) 2015  Camden Sinclair Jansen

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<map>
#include<string>
#include<math.h>
#include<algorithm>


using namespace std;

class TSSsite {
public:
    string chrom;
    int pos;
    string strand;
    TSSsite() {
        chrom = "";
        pos = 0;
        strand = "";
    }
};

//#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
//        ( std::ostringstream() << std::dec << x ) ).str()
//#define SSTR( x ) dynamic_cast< std::ostringstream>(std::ostringstream()<< std::dec << x).str()
string SSTR(int a)
{
    ostringstream temp;
    temp << a;
    return temp.str();
}

string SSTRF(int a)
{
    ostringstream temp;
    temp << a;
    return temp.str();
}


class genomicRegion {
public:
    string gene;
    string chrom;
    int start;
    int stop;
};


template <typename T>
  string NumberToString ( T Number )
  {
     ostringstream ss;
     ss << Number;
     return ss.str();
  }

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

map<string, TSSsite>* parseGtfFile(string gtfFileName, string addtochrom) {
    ifstream gtffile(gtfFileName.c_str());
    map<string, TSSsite>* TSSsites=new map<string, TSSsite>();
    string line;
	
    while(getline(gtffile,line)) {
        if(line[0]=='#') continue;
        vector<string> splitz = split(line, '\t');
        if(splitz[2].compare("exon")==0) {
            string geneName;
            vector<string> attributes = split(splitz[8], ';');
            for(int i = 0; i < attributes.size(); i++) {
                vector<string> pairItems = split(attributes[i],' ');
                /*if(pairItems[1].compare("gene_name")==0) {
                    geneName = pairItems[2].substr(1,pairItems[2].size()-2);
                    break;
                }*/
				//cout<<pairItems[0]<<'\t'<<pairItems[1]<<endl;
				//cout<<pairItems[1].substr(1,pairItems[1].size()-2)<<endl;
				//int temp;
                    //cin>>temp;
				if(pairItems[0].compare("gene_id")==0) {
                    geneName = pairItems[1].substr(1,pairItems[1].size()-2);
					//cout<<geneName<<endl;
					//int temp;
					//cin>>temp;
                    break;
                }
            }

            TSSsite record = TSSsites->operator[](geneName);
            record.strand = splitz[6];
            record.chrom = addtochrom+splitz[0];
			//cout<<record.chrom<<'\t'<<splitz[0];
			//int temp;
			//cin>>temp;
            int pos;
            if(record.strand.compare("+")==0) {
                istringstream(splitz[3])>>pos;
                if(pos < record.pos || record.pos == 0) record.pos = pos;
            } else {
                istringstream(splitz[4])>>pos;
                if(pos > record.pos || record.pos == 0) record.pos = pos;
            }
            TSSsites->operator[](geneName) = record;
        }
    }
    return TSSsites;
}

vector<genomicRegion> GetRegRegions(map<string, TSSsite>* TSSsites, vector<string> genes, string CompareType) {
    vector<genomicRegion> regions;
    for(int i = 0; i < genes.size(); i++) {
		if(i%10000==0) cout<<i<<" Done"<<endl;
        map<string, TSSsite>::iterator it = TSSsites->find(genes[i]);
        TSSsite TSS;
        if(it!=TSSsites->end()) {
            TSS = it->second;
        } else {
            cout<<"Could not find: "<<genes[i]<<".  Skipping..."<<endl;
            continue;
        }
        int start;
        int stop;
        if(CompareType.compare("TwoClosest")==0) {
            int prevdist = 50000;
            if(TSS.pos < 50000) start = 0;
            else start = TSS.pos - 50000;
            int nextdist = 50000;
            if(TSS.pos < 50000) stop = 0;
            else stop = TSS.pos + 50000;

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
            float distance1=50000;
            float distance2=50000;
            for(map<string, TSSsite>::iterator looper = TSSsites->begin(); looper != TSSsites->end(); looper++) {
                TSSsite tester = looper->second;
                if(TSS.chrom.compare(tester.chrom)!=0) continue;
                if(TSS.pos < tester.pos && tester.pos-TSS.pos < distance1) distance1=tester.pos-TSS.pos;
                if(TSS.pos > tester.pos && TSS.pos-tester.pos < distance2) distance2=TSS.pos-tester.pos;
            }
            if(distance2!=50000)
                start = TSS.pos-distance2/2;
            else
                start = TSS.pos-distance2;
            if(distance1!=50000)
                stop = TSS.pos+distance1/2;
            else
                stop = TSS.pos+distance1;
            genomicRegion temp;
            temp.gene = genes[i];
            temp.chrom = TSS.chrom;
            temp.start = start;
            temp.stop = stop;
            regions.push_back(temp);
        }
    }
    return regions;
}


int main(int argc, char* argv[]) {
	if(argc < 2) {
        cout << "Usage: ./getgenes [options] -GTFFile <Gene Anotation (GTF) File Location> -InputPrefix <Input File Location> -OutputPrefix <Output File Location> -Rows <Number of rows in your SOM> -Cols <Number of columns in your SOM>" <<endl;
        cout << "Options: <default>" <<endl;
        cout << "-Method [OneClosest, TwoClosest] <TwoClosest>"<<endl;
        return 0;
    }
	string inputprefix;
	int row1;
    int col1;
	string CompareType = "TwoClosest";
    string gtfFileName;
	string outputprefix;
	string addtochrom="";
	cout<<"Reading Inputs"<<endl;
	for(int i = 0; i < argc; i++) {
        string temp = argv[i];
        if(temp.compare("-Rows")==0)
			istringstream(argv[i+1])>>row1;
		if(temp.compare("-Cols")==0)
            istringstream(argv[i+1])>>col1;
		if(temp.compare("-GTFFile")==0)
			gtfFileName = argv[i+1];
		if(temp.compare("-Method")==0)
            CompareType = argv[i+1];
		if(temp.compare("-OutputPrefix")==0)
            outputprefix = argv[i+1];
		if(temp.compare("-InputPrefix")==0)
            inputprefix = argv[i+1];
		if(temp.compare("-AddToChrom")==0)
			if(i+1 < argc)
	            addtochrom = argv[i+1];
	}

	
	cout<<"Parsing GTF file"<<endl;
	map<string, TSSsite>* TSSsites = parseGtfFile(gtfFileName, addtochrom);
	vector<string> genes;
    for(map<string,TSSsite>::iterator it = TSSsites->begin(); it != TSSsites->end(); ++it) {
		genes.push_back(it->first);
    }
	cout<<genes.size()<<" Genes Loaded."<<endl;
	cout<<"Getting Regulatory Regions"<<endl;
    vector<genomicRegion> TSSregions = GetRegRegions(TSSsites, genes, CompareType);

    for(int row = 0; row < row1; row++) {
        for(int col = 0; col < col1; col++) {
			cout<<"Row: "<<row<<" Col: "<<col<<endl;
            vector<genomicRegion> regions;
            ifstream unit((inputprefix+"_"+SSTR(row)+"_"+SSTR(col)+".unit").c_str());
            string line;
            while(getline(unit,line)) {
				if(line[0]=='#') continue;
				vector<string> checksplitz = split(line,'\t');
				string first = checksplitz[0];
                genomicRegion temp;
                vector<string> splitz = split(first,':');
                if(splitz.size() == 1) continue;
                temp.chrom = splitz[0];
                vector<string> splitz2 = split(splitz[1],'-');
                istringstream(splitz2[0])>>temp.start;
                istringstream(splitz2[1])>>temp.stop;
				regions.push_back(temp);
            }
			vector<string> foundgenes;
			if(regions.size()>0 && TSSregions.size()>0)
				cout<<"Compare the format of these: "<<regions[0].chrom<<'\t'<<TSSregions[0].chrom<<endl;
			for(int i = 0; i < regions.size(); i++) {
				vector<int> found;
				for(int j = 0; j < TSSregions.size(); j++) {
					//int temp;
					//cin>>temp;
					if(regions[i].chrom.compare(TSSregions[j].chrom)==0 && ((regions[i].start<=TSSregions[j].start && regions[i].stop >= TSSregions[j].start)||(regions[i].start <= TSSregions[j].stop && regions[i].stop >= TSSregions[j].stop)||(regions[i].start>=TSSregions[j].start && regions[i].stop <= TSSregions[j].stop)||(regions[i].start<=TSSregions[j].start&&regions[i].stop>=TSSregions[j].stop))) {
						found.push_back(j);
					}
				}
				for(int j = 0; j < found.size(); j++) {
					string tempgene = genes[found[j]];
					bool found2=false;
					for(int k = 0; k < foundgenes.size(); k++) {
						if(foundgenes[k].compare(tempgene) == 0) {
							found2=true;
							break;
						}
					}
					if(!found2) foundgenes.push_back(tempgene);
				}
			}
			cout<<"Found Genes: "<<foundgenes.size()<<endl;

			ofstream outfile((outputprefix+"_"+SSTR(row)+"_"+SSTR(col)+".unit").c_str());
			for(int i = 0; i < foundgenes.size(); i ++) {
				outfile<<foundgenes[i]<<endl;
			}
			outfile.close();
		}
	}
}
