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
        cout << "Usage: ./regionCounts [options] -RawDataFile <raw file list location> -Output <output file location>" <<endl;
        cout << "Options: <default>" <<endl;
        cout << "-MergeRegion: Regions this close together will be merged <0>"<<endl;
        cout << "-MinFeature: Size of smallest partition. <25>"<<endl;
        cout << "-IgnoreRandom: Ignores random chromosomes.  <Off> [Off, On]"<<endl;
        cout << "-PadRegion: Pads regions to this minimum size. <0>"<<endl;

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
        if(temp.compare("-MergeRegion")==0)
            istringstream(argv[i+1])>>mergeRegion;
        if(temp.compare("-MinFeature")==0)
            istringstream(argv[i+1])>>minFeature;
        if(temp.compare("-Output")==0)
            trainingFileName=argv[i+1];
        if(temp.compare("-RawDataFile")==0)
            peakDataFileName=argv[i+1];
        if(temp.compare("-IgnoreRandom")==0) {
			string option = argv[i+1];
			if(option.compare("On") == 0) {
				ignoreRandom=1;
			}
		}   
        if(temp.compare("-PadRegion")==0)
            istringstream(argv[i+1])>>padRegion;
    }
	
	cout<<"Reading Partition Data"<<endl;
	
	ifstream partDataFiles(peakDataFileName.c_str());
	string line1;
	vector<BedNode> allPartData;
	vector<string> chrList;
		//vector<BedNode> peakData;
		while(getline(partDataFiles,line1)) {
			if(line1[0] == '#') continue;
			vector<string> splitz = split(line1, '\t');
			if(splitz.size() < 3) continue;
			BedNode temp;
			temp.chr = splitz[0];
			bool found = false;
			for(int i = 0; i < chrList.size(); i++) {
				if(chrList[i].compare(temp.chr)==0) {
					found = true;
					break;
				}
			}
			if(!found) {
				chrList.push_back(temp.chr);
			}
			istringstream(splitz[1])>>temp.start;
			istringstream(splitz[2])>>temp.stop;
			allPartData.push_back(temp);
		}
		peakDataFile.close();
		//allPeakData.push_back(peakData);

	cout<<"Reading Raw Data"<<endl;
	vector<string> fileNames;
	vector<vector<BedNode> > rawData;
	ifstream rawDataFiles(rawDataFileName.c_str());
	while(getline(rawDataFiles,line1)) {
		ifstream rawDataFile(line1.c_str());
		fileNames.push_back(line1);
		string line;
		vector<BedNode> raw;
		while(getline(rawDataFile, line)) {
			if(line[0]=='@') continue;
			vector<string> splitz = split(line, '\t');
			BedNode temp;
			temp.chr = splitz[2];
			if(temp.chr.compare("*")==0) continue;
			istringstream(splitz[3])>>temp.start;
			istringstream(splitz[8])>>temp.stop;
			temp.stop = temp.start+temp.stop;
			raw.push_back(temp);
		}
		rawDataFile.close();
		rawData.push_back(raw);
	}
	rawDataFiles.close();
	
	cout<<"Counting Regions"<<endl;
	vector<vector<int> > regionCounts;
	for(int i = 0; i < rawData.size(); i++) {
		cout<<fileNames[i]<<endl;
		vector<int> regionCount;
		for(int j = 0; j < allPartData.size(); j++) {
			regionCount.push_back(0);
		}
		for(int j = 0; j < rawData[i].size(); j++) {
			BedNode temp = rawData[i][j];
			for(int k = 0; k < allPartData.size(); k++) {
				BedNode part = allPartData[k];
				if(temp.chr.compare(part.chr)==0 && ((temp.start > part.start && temp.start < part.stop) || (temp.stop > part.start && temp.stop < part.stop) || (temp.start < part.start && temp.stop > part.stop))) {
					regionCount[k]++;
				}
			}
		}
		regionCounts.push_back(regionCount);
	}
	
	cout<<"Outputing Training Matrix"<<endl;
	ofstream outfile(trainingFileName.c_str());
	for(int i = 0; i < allPartData.size(); i++) {
		BedNode part = allPartData[i];
		outfile<<part.chr<<':'<<part.start<<'-'<<part.stop;
		for(int j = 0; j < regionCounts.size(); j++) {
			double RPKM = ((double)regionCounts[j][i])/((double)rawData[j].size())/(1000000.0);
			outfile<<'\t'<<RPKM;
		}
		outfile<<endl;
	}
	outfile.close();
}
