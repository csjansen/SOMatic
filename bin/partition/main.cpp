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
        cout << "Usage: ./partition [options] -PeakDataFile <peak file list location> -Output <output file location>" <<endl;
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
        if(temp.compare("-PeakDataFile")==0)
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
	
	cout<<"Reading Peak Data"<<endl;
	
	ifstream peakDataFiles(peakDataFileName.c_str());
	string line;
	vector<BedNode> allPeakData;
	vector<string> chrList;
	while(getline(peakDataFiles,line)) {
		ifstream peakDataFile(line.c_str());
		string line1;
		//vector<BedNode> peakData;
		while(getline(peakDataFile,line1)) {
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
			allPeakData.push_back(temp);
		}
		peakDataFile.close();
		//allPeakData.push_back(peakData);
	}
	peakDataFiles.close();

	cout<<"Sorting peaks"<<endl;
	sort(allPeakData.begin(), allPeakData.end(), bedSort);

	cout<<"Building Partitions"<<endl;
	vector<vector<BedNode> > partitions;
	for(int i = 0; i < chrList.size(); i++) {
		string chr = chrList[i];
		int position = 0;
		int last = 0;
		vector<BedNode> partition;
		for(int j = 0; j < allPeakData.size(); j++) {
			BedNode temp;
			BedNode peak = allPeakData[j];
			if(peak.chr.compare(chr)==0) {
				if(peak.start - position > minFeature) {
					BedNode temp;
					temp.chr = chr;
					temp.start = position;
					temp.stop = peak.start-1;
					position = peak.start;
					partition.push_back(temp);
				}
			}
		}
		partitions.push_back(partition);
	}
	
	cout<<"Outputing Partitions"<<endl;
	ofstream outfile(trainingFileName.c_str());
	for(int i = 0; i < partitions.size(); i++) {
		for(int j = 0; j < partitions[i].size(); j++) {
			outfile<<partitions[i][j].chr<<':'<<partitions[i][j].start<<'-'<<partitions[i][j].stop<<endl;
		}
	}
	outfile.close();
}
