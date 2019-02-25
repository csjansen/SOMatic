#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <vector>
#include <queue>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string>
#include <locale>
#include <ctime>
using namespace std;

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()
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

int main(int argc, char* argv[]) {
	if(argc < 2) {
        cout << "Usage: ./LinkMeta [options] -FusionFile <Output from Fusion Program> -ClusterFile1 <MetaCluster file from SOM1> -ClusterFile2 <MetaCluster file from SOM2> -Output <Output File Location>" <<endl;
        return 0;
    }
	string infilename;
	string ATACClusterName;
	string RNAClusterName;
	string outfilename;
    for(int i = 0; i < argc; i++) {
        string temp = argv[i];
        if(temp.compare("-FusionFile")==0)
            infilename=argv[i+1];
        if(temp.compare("-ClusterFile1")==0)
            ATACClusterName=argv[i+1];
        if(temp.compare("-ClusterFile2")==0)
            RNAClusterName=argv[i+1];
        if(temp.compare("-Output")==0)
            outfilename=argv[i+1];
	}

	cout<<infilename<<endl;
    ifstream infile(infilename.c_str());
	vector<vector<double> > allpoints;
    map<string,vector<string> > allgenes;
    map<string,vector<string> > allregions;
    cout<<"Opening file"<<endl;
    map<string,int> pointcounts;
	string line;
    while(getline(infile, line)) {
        //cout<<line<<endl;
        vector<string> splitz = split(line,'\t');
        int AtacRow;
        istringstream(splitz[0])>>AtacRow;
        int AtacCol;
        istringstream(splitz[1])>>AtacCol;
        //cout<<AtacRow<<'\t'<<AtacCol<<endl;
        int i = 2;
        while(i < splitz.size()-2) {
            //cout<<i<<'\t'<<splitz.size()<<endl;
            int RNARow;
            int RNACol;
            int AtacCount;
            istringstream(splitz[i])>>RNARow;
            istringstream(splitz[i+1])>>RNACol;
            istringstream(splitz[i+2])>>AtacCount;
            //cout<<RNARow<<'\t'<<RNACol<<endl;
            //for(int j = 0; j < AtacCount; j++) {
              //  allgenes[SSTR(AtacRow)+"."+SSTR(AtacCol)+"."+SSTR(RNARow)+"."+SSTR(RNACol)].push_back(splitz[i+j*2+3]);
                //allregions[SSTR(AtacRow)+"."+SSTR(AtacCol)+"."+SSTR(RNARow)+"."+SSTR(RNACol)].push_back(splitz[i+j*2+1+3]);
            //}
     //       for(int j = 0; j < AtacCount; j++) {
                vector<double> temp;
                temp.push_back(AtacRow);
                temp.push_back(AtacCol);
                temp.push_back(RNARow);
                temp.push_back(RNACol);
                temp.push_back(-1);
                temp.push_back(-1);
                temp.push_back(-1);
                temp.push_back(-1);
                pointcounts[SSTR(AtacRow)+"."+SSTR(AtacCol)+"."+SSTR(RNARow)+"."+SSTR(RNACol)]++;
                allpoints.push_back(temp);
       //     }
            i+=2*AtacCount+3;
        }
    }
    infile.close();
	ifstream ATACClusterFile(ATACClusterName.c_str());
	vector<vector<int> > ATACPoints;
	while(getline(ATACClusterFile, line)) {
		if(line[0]=='#') continue;
		vector<string> splitz = split(line,'\t');
		int row;
		int col;
		int cluster;
		vector<int> temp;
		istringstream(splitz[0])>>row;
		istringstream(splitz[1])>>col;
		istringstream(splitz[2])>>cluster;
		temp.push_back(row);
		temp.push_back(col);
		temp.push_back(cluster);
		ATACPoints.push_back(temp);
	}
	ATACClusterFile.close();
	ifstream RNAClusterFile(RNAClusterName.c_str());
    vector<vector<int> > RNAPoints;
    while(getline(RNAClusterFile, line)) {
        if(line[0]=='#') continue;
        vector<string> splitz = split(line,'\t');
        int row;
        int col;
        int cluster;
        vector<int> temp;
        istringstream(splitz[0])>>row;
        istringstream(splitz[1])>>col;
        istringstream(splitz[2])>>cluster;
        temp.push_back(row);
        temp.push_back(col);
        temp.push_back(cluster);
        RNAPoints.push_back(temp);
    }
    RNAClusterFile.close();
	ofstream outfile(outfilename.c_str());
	for(int i = 0; i < allpoints.size(); i++) {
		for(int j = 0; j < ATACPoints.size(); j++) {
			if(allpoints[i][0]==ATACPoints[j][0] && allpoints[i][1]==ATACPoints[j][1]) {
				allpoints[i][4]=ATACPoints[j][2];
				break;
			}
		}
		for(int j = 0; j < RNAPoints.size(); j++) {
            if(allpoints[i][2]==RNAPoints[j][0] && allpoints[i][3]==RNAPoints[j][1]) {
                allpoints[i][5]=RNAPoints[j][2];
                break;
            }
        }
		outfile<<allpoints[i][0]<<'\t'<<allpoints[i][1]<<'\t'<<allpoints[i][2]<<'\t'<<allpoints[i][3]<<'\t'<<allpoints[i][4]<<'\t'<<allpoints[i][5]<<endl;
	}
	outfile.close();
}
	
