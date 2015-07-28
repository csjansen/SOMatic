/* mapsom: Creates slices of the self-organizing map from the som and score files 
    Copyright (C) 2015 Camden Sinclair Jansen

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<map>
#include<string>
#include<math.h>
#include<algorithm>

using namespace std;

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
        cout << "Usage: ./mapsom -SampleList <Sample List File Location> -SOMFile <SOM File Location> -Prefix <Output Map File Prefix>" <<endl;
        return 0;
    }

	string somFileName;
	string labelfile;
	string prefix;

	for(int i = 0; i < argc; i++) {
        string temp = argv[i];
        if(temp.compare("-SampleList")==0)
			labelfile = argv[i+1];
		if(temp.compare("-SOMFile")==0)
            somFileName = argv[i+1];
		if(temp.compare("-Prefix")==0)
            prefix = argv[i+1];
	}


	cout<<"Inputing SOM file"<<endl;
	ifstream somFile(somFileName.c_str());
    string line;
    vector<vector<vector<double> > > inputMap;
    int lastcol;
    vector<vector<double> > temp;
    while(getline(somFile, line)) {
        if(line[0] == '#') {
            continue;
        }
        vector<string> splitz = split(line,'\t');
        int row;
        int col;
        istringstream(splitz[0])>>row;
        istringstream(splitz[1])>>col;
        if(lastcol > col) {
            inputMap.push_back(temp);
            temp.erase(temp.begin(),temp.end());
            lastcol = col;
        } else {
            lastcol = col;
        }
        vector<double> temp2;
        for(int i = 2; i < splitz.size(); i++) {
            double entry;
            istringstream(splitz[i])>>entry;
            temp2.push_back(entry);
        }
        temp.push_back(temp2);
    }
    inputMap.push_back(temp);
	somFile.close();

	cout<<"Inputing Labels"<<endl;
	vector<string> labels;
	ifstream Labels(labelfile.c_str());
	while(getline(Labels, line)) {
		labels.push_back(line);
	}
	Labels.close();
	vector<vector<double> > summ;
	for(int row = 0; row < inputMap.size(); row++) {
		vector<double> temp;
        for(int col = 0; col < inputMap[row].size(); col++) {
			temp.push_back(0);
		}
		summ.push_back(temp);
	}
	cout<<"Outputing maps"<<endl;
	for(int i = 0; i < labels.size(); i++) {
		ofstream outmap((prefix+labels[i]+".map").c_str());
		cout<<prefix+labels[i]+".map"<<endl;
		for(int row = 0; row < inputMap.size(); row++) {
			for(int col = 0; col < inputMap[row].size()-1; col++) {
				outmap<<inputMap[row][col][i]<<'\t';
				summ[row][col]+=inputMap[row][col][i];
			}
			outmap<<inputMap[row][inputMap[row].size()-1][i]<<endl;
			summ[row][summ[row].size()-1]+=inputMap[row][inputMap[row].size()-1][i];
		}
		outmap.close();
	}
	
	cout<<"Outputing Summary"<<endl;
	ofstream summfile((prefix+"summery.map").c_str());
	for(int row = 0; row < inputMap.size(); row++) {
        for(int col = 0; col < inputMap[row].size()-1; col++) {
			summfile<<summ[row][col]<<'\t';
		}
		summfile<<summ[row][summ[row].size()-1]<<endl;
	}
	summfile.close();
}
