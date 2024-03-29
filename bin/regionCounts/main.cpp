// Counts sorted sam files into training matrices
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <vector>
#include <map>
#include <thread>

using namespace std;

//#define SSTR( x ) static_cast< std::ostringstream & >( \
//        ( std::ostringstream() << std::dec << x ) ).str()
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

int binary_search(vector<BedNode>* chr, int first, int last, BedNode bed) {
	int index;
	//cout<<first<<endl;
	if(first > last) index=-1;
	else {
		int mid = (first+last)/2;
		BedNode temp = chr->at(mid);
		if((temp.start > bed.start && temp.start < bed.stop) || (temp.stop > bed.start && temp.stop < bed.stop) || (temp.start < bed.start && temp.stop > bed.stop) || (temp.start > bed.start && temp.stop < bed.stop)) {
			index = mid;
		} else {
			if(temp.stop > bed.start)
				index = binary_search(chr, first, mid-1, bed);
			else
				index = binary_search(chr, mid+1, last, bed);
		}
	}
	//cout<<"returning "<<index<<endl;
	return index;
}

void regionCount(map<string, vector<BedNode> >* allPartData, string filename, vector<vector<vector<int> > >* counts, int i, vector<string>* chrList, vector<int>* totalRegions) {
	while(true) {
	try {
	vector<string> splitzfile = split(filename,'/');
	cout<<"Started "<<filename<<endl;
	ifstream checkProcessed(("temp_"+splitzfile[splitzfile.size()-1]).c_str());
	cout<<"checking if temp_"<<splitzfile[splitzfile.size()-1]<<" exists"<<endl;
    	string line;
    	vector<vector<int> > regionCount;
	if(checkProcessed.good()) {
		cout<<"importing previous data"<<endl;	
		while(getline(checkProcessed,line)) {
			vector<int> temp;
			vector<string> splitz = split(line,'\t');
			for(int j = 0; j < splitz.size(); j++) {
				int store;
				istringstream(splitz[j])>>store;
				totalRegions->at(i)+=store;
				temp.push_back(store);
			}
			regionCount.push_back(temp);
		}
		counts->at(i) = regionCount;
		cout<<"Finished "<<filename<<endl;
	        return;
	} 
	ifstream rawDataFile(filename.c_str());
	cout<<filename<<endl;
	if(rawDataFile.fail()) {
		cout<<"Failed?"<<endl;
	}
    vector<BedNode> raw;
	vector<int> chrpos;
	totalRegions->at(i)=0;
	cout<<"Opening File"<<endl;
    while(getline(rawDataFile, line)) {
        if(line[0]=='@') continue;
        vector<string> splitz = split(line, '\t');
	if(splitz.size()>=10) {
        	BedNode temp;
	        temp.chr = splitz[2];
		//temp.chr="chr"+temp.chr;
		//cout<<temp.chr<<endl;
        	if(temp.chr.compare("*")==0) continue;
			bool found = false;
			for(int j = 0; j < chrList->size(); j++) {
				if(chrList->at(j).compare(temp.chr)==0) {
					found = true;
					chrpos.push_back(j);
					break;
				}
			}
			if(!found) continue;
	        istringstream(splitz[3])>>temp.start;
		temp.stop=splitz[9].size();
        	//istringstream(splitz[8])>>temp.stop;
	        temp.stop = temp.start+temp.stop;
		totalRegions->at(i)=totalRegions->at(i)+1;
		//cout<<temp.chr<<'\t'<<SSTR(temp.start)<<'\t'<<SSTR(temp.stop)<<'\t'<<endl;
        	raw.push_back(temp);
	}
    }
    rawDataFile.close();
	cout<<"Building region data structure"<<endl;
    for(int j = 0; j < chrList->size(); j++) {
		vector<int> temp;
		for(int k = 0; k < allPartData->at(chrList->at(j)).size(); k++) {
			temp.push_back(0);
		}
		regionCount.push_back(temp);
    }
	cout<<"Counting: "<<filename<<" "<<raw.size()<<endl;
    for(int j = 0; j < raw.size(); j++) {
	        BedNode temp = raw[j];
		//cout<<"binary"<<'\t'<<filename<<'\t'<<j<<endl;
		//cout<<"looking for "<<temp.chr<<'\t'<<temp.start<<'\t'<<temp.stop<<endl;
		int index = binary_search(&(allPartData->at(temp.chr)), 0, allPartData->at(temp.chr).size()-1, temp);
		//cout<<"returned "<<index<<"/"<<regionCount[chrpos[j]].size()<<endl;
		//int temper;
		//cin>>temper;
		if(index != -1) {
		//	cout<<"It was in: "<<allPartData->at(temp.chr)[index].chr<<'\t'<<allPartData->at(temp.chr)[index].start<<'\t'<<allPartData->at(temp.chr)[index].stop<<endl;
			regionCount[chrpos[j]][index]++;
		//int temper;
		//cin>>temper;
		}
		/*for(int k = 0; k < allPartData->at(temp.chr).size(); k++) {
            BedNode part = allPartData->at(temp.chr).at(k);
            if(temp.chr.compare(part.chr)==0 && ((temp.start > part.start && temp.start < part.stop) || (temp.stop > part.start && temp.stop < part.stop) || (temp.start < part.start && temp.stop > part.stop))) {
                regionCount[k]++;
            }
        }*/
    }
    counts->at(i) = regionCount;
	ofstream outfile(("temp_"+splitzfile[splitzfile.size()-1]).c_str());
	for(int j = 0; j < regionCount.size(); j++) {
		for(int k = 0; k < regionCount[j].size(); k++) {
			if(k==0) {
				outfile<<regionCount[j][k];
			} else {
				outfile<<'\t'<<regionCount[j][k];
			}
		}
		outfile<<endl;
	}
	cout<<"Finished "<<filename<<endl;
	return;
	} catch(...) {
		cout<<"Failure with "<<filename<<" - Restarting"<<endl;
	}
	}
}

int main(int argc, char *argv[]) {
    if(argc < 2) {
        cout << "Usage: ./regionCounts [options] -RawDataFile <raw sam file list location> -Partitions <partition file> -Output <output file location>" <<endl;
        cout << "Options: <default>" <<endl;
	cout << "-Log2: Log scale RPKM"<<endl;
		return 0;
    }
    int mergeRegion=0;
    int minFeature=25;
	int ignoreRandom = 0;
	int padRegion = 0;
	string rawDataFileName;
	string peakDataFileName;
    string trainingFileName;
	bool logScale = false;
    for(int i = 0; i < argc; i++) {
        string temp = argv[i];
        if(temp.compare("-Output")==0)
            trainingFileName=argv[i+1];
        if(temp.compare("-RawDataFile")==0)
            rawDataFileName=argv[i+1];
	if(temp.compare("-Partitions")==0)
	    peakDataFileName=argv[i+1];
	if(temp.compare("-Log2")==0)
            logScale=true;
    }
	
	cout<<"Reading Partition Data"<<endl;
	cout<<peakDataFileName<<endl;	
	ifstream partDataFiles(peakDataFileName.c_str());
	string line1;
	map<string, vector<BedNode> > allPartData;
	vector<string> chrList;
	int parts = 0;
		//vector<BedNode> peakData;
		while(getline(partDataFiles,line1)) {
			if(line1[0] == '#') continue;
			vector<string> splitz = split(line1, ':');
			if(splitz.size() < 2) continue;
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
			vector<string> splitz2 = split(splitz[1], '-');
			istringstream(splitz2[0])>>temp.start;
			istringstream(splitz2[1])>>temp.stop;
			allPartData[temp.chr].push_back(temp);
			parts++;
		}
		partDataFiles.close();
		//allPeakData.push_back(peakData);
	cout<<parts<<" partitions loaded."<<endl;
	cout<<"Reading Raw Data"<<endl;
	cout<<rawDataFileName<<endl;
	vector<string> fileNames;
	vector<vector<BedNode> > rawData;
	ifstream rawDataFiles(rawDataFileName.c_str());
	int count = 0;
	vector<thread> threads;
	vector<vector<vector<int> > > regionCounts;
	vector<int> totalRegions;
	while(getline(rawDataFiles,line1)) {
		regionCounts.push_back(vector<vector<int> >());
		totalRegions.push_back(0);
		threads.push_back(thread(regionCount, &allPartData, line1, &regionCounts, count, &chrList,&totalRegions));
		count++;
		if(count%10==0) {
			for (auto& th : threads) th.join();
			threads.clear();
		}
	}
	rawDataFiles.close();
	for (auto& th : threads) th.join();
	for(int i = 0; i<count; i++) {
		cout<<totalRegions[i]<<endl;
	}
	cout<<"Outputing Training Matrix"<<endl;
	ofstream outfile(trainingFileName.c_str());
	//cout<<trainingFileName<<'\t'<<regionCounts[0].size()<<endl;
	for(int i = 0; i < regionCounts[0].size(); i++) {
		//cout<<chrList[i]<<endl;
		for(int k = 0; k < allPartData[chrList[i]].size(); k++) {
			BedNode part = allPartData[chrList[i]][k];
		//	cout<<part.chr<<':'<<part.start<<'-'<<part.stop;
			outfile<<part.chr<<':'<<part.start<<'-'<<part.stop;
			for(int j = 0; j < regionCounts.size(); j++) {
				//for(int h = 0; h < regionCounts[j][i].size(); h++) {
				if(part.stop==part.start || totalRegions[j]==0) {
					outfile <<'\t'<<0;
				} else {
					double RPKM = ((double)regionCounts[j][i][k])/((double)(part.stop-part.start)/(1000.0)*(totalRegions[j]/1000000.0));
		//			cout<<'\t'<<RPKM;
					if(logScale) RPKM = log2(RPKM+1);
					outfile<<'\t'<<RPKM;
				}
				//}
			}
		//	cout<<endl;
			outfile<<endl;
		}
	}
	outfile.close();
}
