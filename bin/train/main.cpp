/* trainsom: trains a self-organizing map
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
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <vector>
#include <map>

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

vector <vector<int> > merge(vector <vector<int> > lhs, vector<vector<int> > rhs) {
	vector <vector<int> > result;
	result.insert(result.begin(),lhs.begin(), lhs.end());
	for(int i = 0; i < rhs.size(); i++) {
		bool found = false;
		for(int j = 0; j < lhs.size(); j++) {
			if(rhs[i] == lhs[j]){
				found = true;
				break;
			}
		}
		if(!found) result.push_back(rhs[i]);
	}
	return result;
}


vector <vector<int> > hexSurround(vector<int> input, int radius, int numRows, int numCols) {
	vector<vector<int> > result;
	int factorCol = numRows/numCols+1;
	int factorRow = numCols/numRows+1;
	for(int i = -1*radius; i <= radius; i++) {
		int row = (input[0]+i+factorRow*numRows)%numRows;
		vector<int> distance;
		for(int j = -1*radius; j<= radius; j++) {
			int col = (input[1]+j+factorCol*numCols)%numCols;
			if((i==-radius&&j==0) || (i==radius&&j==0)) {
				continue;
			}
			vector<int> temp;
			temp.push_back(row);
			temp.push_back(col);
			result.push_back(temp);
		}
	}
	double temp;
	return result;
	/*vector <vector <int> > result;
	input.push_back(radius);
	result.push_back(input);
	if(radius == 0) return result;
	int arow = input[0];
	int acol = input[1];
	
	int arowm = (arow - 1+numRows) % numRows;
	int arowp = (arow + 1) % numRows;
    int acolm = (acol - 1+numCols) % numCols;
    int acolp = (acol + 1) % numCols;
	vector<int> temp1;
	temp1.push_back(arowm);
	temp1.push_back(acol);
	vector<int> temp2;
	temp2.push_back(arowp);
	temp2.push_back(acol);
	vector<int> temp3;
	temp3.push_back(arow);
	temp3.push_back(acolm);
	vector<int> temp4;
	temp4.push_back(arow);
	temp4.push_back(acolp);
	vector<int> temp5;
	vector<int> temp6;
        if(arow % 2 == 0) {
		temp5.push_back(arowm);
		temp5.push_back(acolm);
		temp6.push_back(arowp);
		temp6.push_back(acolm);
	} else {
		temp5.push_back(arowm);
		temp5.push_back(acolp);
		temp6.push_back(arowp);
		temp6.push_back(acolp);
	}
	result = merge(result,hexSurround(temp1, radius-1, numRows, numCols));
	result = merge(result,hexSurround(temp2, radius-1, numRows, numCols));
	result = merge(result,hexSurround(temp3, radius-1, numRows, numCols));
	result = merge(result,hexSurround(temp4, radius-1, numRows, numCols));
	result = merge(result,hexSurround(temp5, radius-1, numRows, numCols));
	result = merge(result,hexSurround(temp6, radius-1, numRows, numCols));
	return result;*/
}

vector<int> propagate(vector<double>* trainingVector, int numRows, int numCols, vector<vector<vector<double> > >* trainingMap) {
	vector<int> winUnit;
    winUnit.push_back(0);
    winUnit.push_back(0);
    double smallest = -999;
    for(int row = 0; row < numRows; row++) {
      	for(int col = 0; col < numCols; col++) {
       		double magSquared=0;
       		for(int num = 0; num < trainingVector->size(); num++) {
          		magSquared += pow((*trainingVector)[num] - (*trainingMap)[row][col][num],2);
           	}
           	if(smallest == -999 || magSquared < smallest) {
           		(winUnit)[0] = row;
               	(winUnit)[1] = col;
               	smallest = magSquared;
           	}
       	}
    }
	return winUnit;
}

int main(int argc, char *argv[]) {
	if(argc < 2) {
		cout << "Usage: ./trainsom [options] -TrainingMatrix <Training Matrix File Location> -SOMFile <Output File Location>" <<endl;
		cout << "Options: <default>" <<endl;
		cout << "-Rows: Number of rows for your SOM. <20>"<<endl;
		cout << "-Cols: Number of columns for your SOM. <30>"<<endl;
		cout << "-Trials: Number of trials run for your SOM. The package will choose the best SOM after all trials.  <3>"<<endl;
		cout << "-Timesteps: Number of timesteps for each trial. <3000000>"<<endl;
		cout << "-Topology: Topology for your SOM. (Only toroid currently supported) <toroid>"<<endl;
		cout << "-Seed: Set seed for random initialization."<<endl; 
		return 0;
	}
	double bestScore = 1000000;
	int numRows=20;
	int numCols=30;
	string trainingFileName;
	string somFile="out.som";
	int trials=3;
	int timesteps=3000000;
	string topology="toroid";
	int seed=-1;
	for(int i = 0; i < argc; i++) {
		string temp = argv[i];
		if(temp.compare("-Rows")==0) 
			istringstream(argv[i+1])>>numRows;
		if(temp.compare("-Cols")==0)
			istringstream(argv[i+1])>>numCols;
		if(temp.compare("-TrainingMatrix")==0)
			trainingFileName=argv[i+1];
		if(temp.compare("-SOMFile")==0)
			somFile=argv[i+1];
		if(temp.compare("-Trials")==0)
			istringstream(argv[i+1])>>trials;
		if(temp.compare("-Timesteps")==0)
			istringstream(argv[i+1])>>timesteps;
		if(temp.compare("-Topology")==0)
			topology=argv[i+1];
		if(temp.compare("-Seed")==0)
			istringstream(argv[i+1])>>seed;
	}
	int dimension = -1;
	ifstream trainingFile(trainingFileName.c_str());
	string line;
	map<string, vector<double> > dataMap;
	vector<string> dataKeys;
	cout<<"Reading File: "<<trainingFileName<<endl;
	int count = 0;
	while(getline(trainingFile,line)) {
		if(line.size() < 2) continue;
		if(line[0]=='#') continue;
		vector<string> fields = split(line, '\t');
		dimension = fields.size()-1;
		vector<double> temp;
		for(int i = 1; i < fields.size(); i++) {
			double input;
			istringstream(fields[i])>>input;
			temp.push_back(input);
		}
		dataMap[fields[0]]=temp;
		dataKeys.push_back(fields[0]);
		count++;
		if(count%100000==0) 
		cout<<"line: "<<count<<endl;
	}	
	cout<<"Done!"<<endl;
	int radius;
        if(numCols > numRows)   radius = numCols/2;
        else            radius = numRows/2;

	float learningRate = 0.2;
	for(int i = 0; i < trials; i++) {
		cout<<"Trial: "<<i+1<<endl;
		if(seed==-1)
			srand(time(0));
		else
			srand(seed);
		random_shuffle(dataKeys.begin(), dataKeys.end());
		// make new map
		vector<vector <vector <double> > > trainingMap;
		cout<<"Building Training map"<<endl;
		cout<<dataMap.size()<<endl;
        for(int row = 0; row < numRows; row++) {
            vector<vector <double> > temp2;
			//cout<<row<<endl;
            for(int col = 0; col < numCols; col++) {
				//cout<<col<<endl;
                vector<double> temp;
				//cout<<dataKeys.size()<<endl;
				string randomID = dataKeys[row*numCols+col];
               	for(int k = 0; k < dataMap[randomID].size(); k++) {
					//cout<<k<<endl;
               		temp.push_back(dataMap[randomID][k]);
               	}
               	temp2.push_back(temp);
            }
            trainingMap.push_back(temp2);
		}
		cout<<"done"<<endl;
		// train
		double multiplier = -1 * log(radius) / float(timesteps);
		cout<<"Beginning Training"<<endl;
		random_shuffle(dataKeys.begin(), dataKeys.end());
		for(int j = 0; j < timesteps; j++) {
			string trainingID = dataKeys[j % dataKeys.size()];
			vector<double> trainingVector = dataMap[trainingID];
			double decay = exp(j * multiplier);
			double timeLearningRate = learningRate * decay;
			int timeRadius = int(radius * decay);
	
			// Finding Winning Unit
			//cout<<"Finding Winning Unit"<<endl;
			vector<int> winUnit = propagate(&(trainingVector),numRows,numCols, &trainingMap);
			//cout<<"Found... Getting Neighbors: "<<winUnit[0]<<" "<<winUnit[1]<<endl;
			vector<vector<int> > winUnitNeighbors = hexSurround(winUnit,timeRadius,numRows, numCols);
			//cout<<"Got them.  Updating winning unit."<<endl;
			for(int k = 0; k < trainingVector.size(); k++) {
				trainingMap[(winUnit)[0]][(winUnit)[1]][k] += timeLearningRate * (trainingVector[k] - trainingMap[(winUnit)[0]][(winUnit)[1]][k]);
			}
			//cout<<"Changing neighbors"<<endl;
//			for(int dist = 1; dist<timeRadius+1;dist++) {
//				float theWeight = timeLearningRate * exp(-0.5*pow(dist,2)/pow(radius+1,2));
				for(int i = 0; i < winUnitNeighbors.size(); i++) {
					if(winUnitNeighbors[i][0]!=winUnit[0]||winUnitNeighbors[i][1]!=winUnit[1]) {
						vector<int> unit = winUnitNeighbors[i];
						int xdist = (unit[0] - winUnit[0]);
						int ydist = (unit[1] - winUnit[1]);
						int ddist = abs(xdist-ydist);
						int dist = 0;
						if(abs(xdist) > abs(ydist) && abs(xdist) > abs(ddist)) dist = abs(xdist);
						if(abs(ydist) > abs(xdist) && abs(ydist) > abs(ddist)) dist = abs(ydist);
						if(abs(ddist) > abs(xdist) && abs(ddist) > abs(ydist)) dist = abs(ddist);
						float theWeight = timeLearningRate * exp(-0.5*pow(dist,2)/pow(radius+1,2));
						for(int k = 0; k < trainingVector.size(); k++) {
							trainingMap[unit[0]][unit[1]][k] += theWeight*(trainingVector[k]-trainingMap[unit[0]][unit[1]][k]);
						}
					}
				}
//			}
			//cout<<"done"<<endl;
			//Score
			// New winner map
			vector<vector<vector<string > > > winnerMap;
            for(int row = 0; row < numRows; row ++) {
				vector<vector<string > > temp2;
                   	for(int col = 0; col < numCols; col++) {
					vector<string> temp;
					temp2.push_back(temp);
				}
				winnerMap.push_back(temp2);
			}
			//Get winners	
			vector<int> winunit = propagate(&(dataMap[trainingID]),numRows,numCols, &trainingMap);
			winnerMap[(winunit)[0]][(winunit)[1]].push_back(trainingID);
			double totalScore = 0;
			//Count up score	
			for(int row = 0; row < numRows; row++) {
				for(int col = 0; col < numCols; col++) {
					
					for(int k = 0; k < winnerMap[row][col].size(); k++) {
						double tempScore = 0;
						for(int u = 0; u < dataMap[winnerMap[row][col][k]].size(); u++) {
							tempScore+=pow(dataMap[winnerMap[row][col][k]][u]-trainingMap[row][col][u],2);
						}
						totalScore += sqrt(tempScore);
					}
				}
			}
			totalScore/=(double)dataKeys.size();
			if(totalScore < bestScore && totalScore!=0) {
				ofstream outfile(somFile.c_str());
				outfile<<"# "<<numRows<<" rows\t"<<numCols<<" cols\t"<<trainingVector.size()<<" dimensions\ttoroid topology"<<endl;
				for(int row = 0; row < numRows; row++) {
					for(int col = 0; col <numCols; col++) {
						outfile<<row<<'\t'<<col<<'\t';
						for(int k = 0; k < trainingVector.size(); k++) {
							outfile.precision(5);
							outfile<<fixed<<trainingMap[row][col][k];
							if(k != trainingVector.size()-1)
								outfile<<'\t';
						}
						outfile<<endl;
					}
				}
				outfile.close();
				bestScore = totalScore;
			}		
			if(j % 10000 == 0) {
				cout << "time " << j << " radius " << timeRadius<< " learn " << timeLearningRate<<endl;
			}
		}		
	}
	cout<<"Best Score = "<<bestScore<<endl;
}
