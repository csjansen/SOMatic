#include<iostream>
#include<fstream>
#include<vector>
#include<sstream>
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
	ifstream firstmap(argv[1]);
	ifstream secondmap(argv[2]);
	ofstream output(argv[3]);
	string line;
	vector<vector<double> > map1;
	while(getline(firstmap,line)) {
		vector<double> temp;
		vector<string> splitz = split(line,'\t');
		for(int i = 0; i < splitz.size(); i++) {
			double push;
			istringstream(splitz[i])>>push;
			temp.push_back(push);
		}
		map1.push_back(temp);
	}
	firstmap.close();
	vector<vector<double> > map2;
    while(getline(secondmap,line)) {
        vector<double> temp;
        vector<string> splitz = split(line,'\t');
        for(int i = 0; i < splitz.size(); i++) {
            double push;
            istringstream(splitz[i])>>push;
            temp.push_back(push);
        }
		map2.push_back(temp);
    }
	secondmap.close();
	for(int i = 0; i < map1.size(); i++) {
		for(int j = 0; j < map1[i].size()-1; j++) {
			output<<(map1[i][j]-map2[i][j])<<'\t';
		}
		output<<(map1[i][map1.size()-1]-map2[i][map2.size()-1])<<endl;
	}
	output.close();
	return 0;
}
