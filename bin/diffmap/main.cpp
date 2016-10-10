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
	cout<<argv[1]<<endl;
	ifstream firstmap1(argv[1]);
	cout<<argv[2]<<endl;
	ifstream firstmap2(argv[2]);
	//cout<<argv[3]<<endl;
	//ifstream firstmap3(argv[3]);
	cout<<argv[3]<<endl;
	ifstream secondmap1(argv[3]);
	cout<<argv[4]<<endl;
	ifstream secondmap2(argv[4]);
	//cout<<argv[6]<<endl;
	//ifstream secondmap3(argv[6]);
	ofstream output(argv[5]);
	string line;
	vector<vector<double> > map11;
	while(getline(firstmap1,line)) {
		vector<double> temp;
		vector<string> splitz = split(line,'\t');
		for(int i = 0; i < splitz.size(); i++) {
			double push;
			istringstream(splitz[i])>>push;
			temp.push_back(push);
		}
		map11.push_back(temp);
	}
	firstmap1.close();
	vector<vector<double> > map12;
    while(getline(firstmap2,line)) {
        vector<double> temp;
        vector<string> splitz = split(line,'\t');
        for(int i = 0; i < splitz.size(); i++) {
            double push;
            istringstream(splitz[i])>>push;
            temp.push_back(push);
        }
        map12.push_back(temp);
    }
    firstmap2.close();
	/*vector<vector<double> > map13;
    while(getline(firstmap3,line)) {
        vector<double> temp;
        vector<string> splitz = split(line,'\t');
        for(int i = 0; i < splitz.size(); i++) {
            double push;
            istringstream(splitz[i])>>push;
            temp.push_back(push);
        }
        map13.push_back(temp);
    }
    firstmap3.close();
*/
	vector<vector<double> > map21;
    while(getline(secondmap1,line)) {
        vector<double> temp;
        vector<string> splitz = split(line,'\t');
        for(int i = 0; i < splitz.size(); i++) {
            double push;
            istringstream(splitz[i])>>push;
            temp.push_back(push);
        }
		map21.push_back(temp);
    }
	secondmap1.close();

	vector<vector<double> > map22;
    while(getline(secondmap2,line)) {
        vector<double> temp;
        vector<string> splitz = split(line,'\t');
        for(int i = 0; i < splitz.size(); i++) {
            double push;
            istringstream(splitz[i])>>push;
            temp.push_back(push);
        }
        map22.push_back(temp);
    }
    secondmap2.close();

/*	vector<vector<double> > map23;
    while(getline(secondmap3,line)) {
        vector<double> temp;
        vector<string> splitz = split(line,'\t');
        for(int i = 0; i < splitz.size(); i++) {
            double push;
            istringstream(splitz[i])>>push;
            temp.push_back(push);
        }
        map23.push_back(temp);
    }
    secondmap3.close();
*/
	cout<<"Doing the math."<<endl;
	cout<<map11[0][0]<<endl;
	cout<<map12[0][0]<<endl;
	cout<<map21[0][0]<<endl;
	cout<<map22[0][0]<<endl;
	for(int i = 0; i < map11.size(); i++) {
		for(int j = 0; j < map11[i].size(); j++) {
			if(j!=0)
				output<<'\t';
			output<<(map11[i][j]+map12[i][j])/2.0-(map21[i][j]+map22[i][j])/2.0;
		}
		output<<endl;
	}
	output.close();
	return 0;
}
