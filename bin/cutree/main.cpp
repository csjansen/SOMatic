#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

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


int main(int argc, char *argv[]) {
    if(argc < 2) {
        cout << "Usage: ./cutree [options] -ClusterFile <Cluster file location> -CutLevel <\% of max hieght to make cut at> -SampleList <Sample list used in SOM creation> -Output <Output File Location>" <<endl;
        cout << "Options: <default>" <<endl;
        return 0;
    }
    string clusterFileLocation;
	double cutLevel=70;
	string outputFileLocation;
	string sampleListFileLocation;
	for(int i = 0; i < argc; i++) {
        string temp = argv[i];
        if(temp.compare("-ClusterFile")==0)
			clusterFileLocation=argv[i+1];
		if(temp.compare("-Output")==0)
            outputFileLocation=argv[i+1];
		if(temp.compare("-CutLevel")==0)
            istringstream(argv[i+1])>>cutLevel;
		if(temp.compare("-SampleList")==0)
            sampleListFileLocation=argv[i+1];
	}
	int maxDistance = 0;
	ifstream clusterFile(clusterFileLocation.c_str());
	string line;
	vector<int> posOrder;
	vector<double> distanceOrder;
	cout<<"Open Cluster File"<<endl;
	while(getline(clusterFile,line)) {
		vector<string> splitz = split(line,'\t');
		int tempPos;
		istringstream(splitz[0])>>tempPos;
		double tempDistance;
		istringstream(splitz[1])>>tempDistance;
		posOrder.push_back(tempPos);
		distanceOrder.push_back(tempDistance);
		if(maxDistance < tempDistance) maxDistance=tempDistance;
	}
	cout<<"done"<<endl;
	double limit = maxDistance*cutLevel/100.0;
	cout<<"Limit: "<<limit<<endl;
	clusterFile.close();
	vector<string> sampleList;
	ifstream sampleListFile(sampleListFileLocation.c_str());
	cout<<"Open Sample List"<<endl;
	while(getline(sampleListFile,line)) {
		sampleList.push_back(line);
	}
	sampleListFile.close();
	cout<<"done"<<endl;
	vector<vector<int> > clusters;
	vector<int> temp;
	while(posOrder.size()>0) {
		cout<<posOrder.size()<<endl;
		// Find smallest value remaining that isn't 0
		int mini = 0;
		double minDist = 1000000;	
		for(int i = 0; i < distanceOrder.size(); i++) {
			if(distanceOrder[i] < minDist && distanceOrder[i]!=0) {
				minDist = distanceOrder[i];
				mini = i;
			}
		}
		temp.push_back(posOrder[mini]);
		posOrder.erase(posOrder.begin()+mini);
		distanceOrder.erase(distanceOrder.begin()+mini);
		//begin adding closest values till one is over the limit
		bool cont = true;
		while(cont && posOrder.size()>0) {
			cout<<posOrder.size()<<endl;
			double compDist = 100000000;
			int compi=0;
			if(mini>0) {
				compDist=distanceOrder[mini-1];
				compi=mini-1;
			}
			if(mini<distanceOrder.size()) {
				if(compDist > distanceOrder[mini]) {
					compDist = distanceOrder[mini];
					compi=mini;
				}
			}
			cout<<compDist<<'\t'<<compi<<endl;
			if(compDist <= limit) {
				temp.push_back(compi);
				posOrder.erase(posOrder.begin()+compi);
				distanceOrder.erase(distanceOrder.begin()+compi);
			} else {
				clusters.push_back(temp);
				temp.clear();
				cont=false;
			}	
			if(mini > distanceOrder.size()) {
				mini--;
			}
		}	
	}
	ofstream outputFile(outputFileLocation.c_str());
	outputFile<<"#\t"<<clusters.size()<<endl;
	for(int i = 0; i < clusters.size(); i++) {
		for(int j = 0; j < clusters[i].size(); j++) {
			outputFile<<sampleList[clusters[i][j]]<<endl;
		}
		outputFile<<endl;
	}
	outputFile.close();
}	
