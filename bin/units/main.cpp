/* getunits: Generates unit files from a training matrix and a self organizing map
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

template <typename T>
  string NumberToString ( T Number )
  {
     ostringstream ss;
     ss << Number;
     return ss.str();
  }

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
        cout << "Usage: ./getunits -Rows <Rows in your SOM> -Cols <Columns in your SOM> -ScoreFile <Score File Location> -Prefix <Output Prefix for Unit Files>" <<endl;
        return 0;
    }

	string scorefileName;
	int numRows;
	int numCols;
	string prefix;
	for(int i = 0; i < argc; i++) {
        string temp = argv[i];
        if(temp.compare("-Rows")==0)
			istringstream(argv[i+1])>>numRows;
		if(temp.compare("-Cols")==0)
            istringstream(argv[i+1])>>numCols;
		if(temp.compare("-ScoreFile")==0)
            scorefileName=argv[i+1];
		if(temp.compare("-Prefix")==0)
            prefix=argv[i+1];
	}

	ifstream scoreFile(scorefileName.c_str()); 
	string line;
	vector<vector<vector<string> > > Scores;
	vector<vector<string> > Scorerow;
	vector<vector<vector<double> > > Dists;
	vector<vector<double> > Distrow;
	vector<string> temp;
	vector<double> tempdist;
	int rowcoord=0;
	int colcoord=0;
	cout<<"Inputing Score file"<<endl;
	while(getline(scoreFile, line)) {
		vector<string> splitz = split(line, '\t');
		if(splitz[0].compare("unit")==0) {
			vector<string> coords = split(splitz[1],',');
			int row;
			int col;
			istringstream(coords[0])>>row;
			istringstream(coords[1])>>col;
			if(row==0 && col ==0) continue;
			Scorerow.push_back(temp);
			Distrow.push_back(tempdist);
			temp.erase(temp.begin(),temp.end());
			tempdist.erase(tempdist.begin(),tempdist.end());
			if(row!=rowcoord) {
				Scores.push_back(Scorerow);
				Dists.push_back(Distrow);
				Scorerow.erase(Scorerow.begin(), Scorerow.end());
				Distrow.erase(Distrow.begin(), Distrow.end());
			}
			rowcoord = row;
		} else {
			temp.push_back(splitz[1]);
			string lastvalstr = splitz[splitz.size()-1];
			double lastval;
			istringstream(lastvalstr)>>lastval;
			tempdist.push_back(lastval);
		}
	}
	scoreFile.close();
	Scorerow.push_back(temp);
	Distrow.push_back(tempdist);
	Scores.push_back(Scorerow);
	Dists.push_back(Distrow);
	
	for(int i = 0; i < numRows; i++) {
		for(int j = 0; j < numCols; j++) {
			cout<<"Creating "<<prefix+'_'+NumberToString(i)+'_'+NumberToString(j)+".unit"<<endl;
			ofstream outfile((prefix+'_'+NumberToString(i)+'_'+NumberToString(j)+".unit").c_str());
			cout<<Scores[i][j].size()<<endl;
			for(int k = 0; k < Scores[i][j].size(); k++) {
				outfile<<Scores[i][j][k]<<'\t'<<Dists[i][j][k]<<endl;
			}
			outfile.close();
		}
	}
}
