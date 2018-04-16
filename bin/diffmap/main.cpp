#include<iostream>
#include<fstream>
#include<vector>
#include<sstream>
#include<math.h>
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

int main(int argc, char* argv[]) {
	cout<<argv[1]<<endl;
	ifstream firstmaps(argv[1]);
	cout<<argv[2]<<endl;
	ifstream secondmaps(argv[2]);
	ofstream output(argv[3]);
	string line;
	vector<string> firstMapFileNames;
	vector<string> secondMapFileNames;
	while(getline(firstmaps,line)) {
		firstMapFileNames.push_back(line+".map");
	}	
	while(getline(secondmaps,line)) {
        secondMapFileNames.push_back(line+".map");
    }
	firstmaps.close();
	secondmaps.close();
	vector<vector<double> > AverageMap1;
	vector<vector<double> > AverageMap2;
	for(int i = 0; i < firstMapFileNames.size(); i++) {
		ifstream mapFile(firstMapFileNames[i].c_str());
		int row = 0;
		while(getline(mapFile,line)) {
			vector<string> splitz = split(line,'\t');
			vector<double> temp;
			for(int j = 0; j < splitz.size(); j++) {
				double num;
				istringstream(splitz[j])>>num;
	//			num = exp(num)-1;
				if(i==0) {
					temp.push_back(num);
				} else {
					AverageMap1[row][j]+=num;
				}
			}
			if(i==0)
				AverageMap1.push_back(temp);
			row++;
		}
		mapFile.close();
	}
	for(int i = 0; i < secondMapFileNames.size(); i++) {
        ifstream mapFile(secondMapFileNames[i].c_str());
        int row = 0;
        while(getline(mapFile,line)) {
			//cout<<line<<endl;
            vector<string> splitz = split(line,'\t');
            vector<double> temp;
            for(int j = 0; j < splitz.size(); j++) {
                double num;
                istringstream(splitz[j])>>num;
	//			num = exp(num)-1;
                if(i==0) {
                    temp.push_back(num);
                } else {
                    AverageMap2[row][j]+=num;
                }
            }
            if(i==0)
                AverageMap2.push_back(temp);
            row++;
        }
        mapFile.close();
    }
	for(int i = 0; i < AverageMap1.size(); i++) {
		for(int j = 0; j < AverageMap1[i].size(); j++) {
			AverageMap1[i][j]/=(double)firstMapFileNames.size();
			AverageMap2[i][j]/=(double)secondMapFileNames.size();	
		}
	}
	
	for(int i = 0; i < AverageMap1.size(); i++) {
        for(int j = 0; j < AverageMap1[i].size(); j++) {
			//double outputNum = (AverageMap1[i][j]-AverageMap2[i][j]);
			double outputNum;
			outputNum = log((AverageMap1[i][j]+.01)/(AverageMap2[i][j]+.01))/log(2);
			if(j==0) {
				output<<outputNum;
			} else {
				output<<'\t'<<outputNum;
			}
		}
		output<<endl;
	}
	output.close();
}
