#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <vector>
#include <map>

using namespace std;

class BedNode {
public:
	string chr;
	int start;
	int stop;
};

//Move to util header
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

bool bedSort(BedNode i, BedNode j) {
	return i.start < j.start;
}


int main(int argc, char *argv[]) {
    if(argc < 2) {
        cout << "Usage: ./RSemToTrainingMatrix [options] -RsemFileList <file list location> -Output <output file location>" <<endl;
        cout << "Options: <default>" <<endl;

		return 0;
    }
    int mergeRegion=0;
    int minFeature=25;
	int ignoreRandom = 0;
	int padRegion = 0;
	string rawDataFileName;
	string peakDataFileName;
    string trainingFileName;
    for(int i = 0; i < argc; i++) {
        string temp = argv[i];
        if(temp.compare("-Output")==0)
            trainingFileName=argv[i+1];
        if(temp.compare("-RsemFileList")==0)
            peakDataFileName=argv[i+1];
    }
	cout<<"Reading RSem Data"<<endl;
	ifstream RsemDataFiles(peakDataFileName.c_str());
	string line1;
	vector<string> geneList;
	vector<vector<double> > RsemData;

	while(getline(RsemDataFiles,line1)) {
		ifstream RsemDataFile((line1).c_str());
		string line;
		vector<double> rsem;
		while(getline(RsemDataFile, line)) {
			vector<string> splitz = split(line, '\t');
			if(splitz.size() < 6) continue;
			if(splitz[0].compare("gene_id")==0) continue;
			int found = -1;
			for(int i = 0; i < geneList.size(); i++) {
				if(splitz[0].compare(geneList[i])==0) {
					found = i;
					break;
				}
			}
			if(found == -1) {
				geneList.push_back(splitz[0]);
				for(int i = 0; i < RsemData.size(); i++) {
					RsemData[i].push_back(0);
				}
			}
			double temp;
			istringstream(splitz[5])>>temp;
			rsem.push_back(temp);
		}
		while(rsem.size() < geneList.size()) {
			rsem.push_back(0);
		}
		RsemDataFile.close();
		RsemData.push_back(rsem);
	}
	RsemDataFiles.close();
	
	cout<<"Creating Training Matrix"<<endl;
	
	ofstream outfile(trainingFileName.c_str());
	for(int i = 0; i < geneList.size(); i++) {
		bool found = false;
		for(int j = 0; j < RsemData.size(); j++) {
			if(RsemData[j][i]!=0) {
				found = true;
				break;
			}
		}
		if(!found) continue;
		string gene = geneList[i];
		outfile<<gene;
		for(int j = 0; j < RsemData.size(); j++) {
			double RPKM = RsemData[j][i];
			outfile<<'\t'<<RPKM;
		}
		outfile<<endl;
	}
	outfile.close();
}
