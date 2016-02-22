/* scoresom: Generates the score file for a self-organizing map
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

vector<double> propagate(vector<double>* trainingVector, int numRows, int numCols, vector<vector<vector<double> > >* trainingMap, bool sparse,bool print) {
    vector<double> winUnit;
    winUnit.push_back(0);
    winUnit.push_back(0);
    winUnit.push_back(0);
    double smallest = -999;
	if(!sparse) {
	    for(int row = 0; row < numRows; row++) {
			    for(int col = 0; col < numCols; col++) {
				double magSquared=0;
		        for(int num = 0; num < trainingVector->size(); num++) {
			        magSquared += pow((*trainingVector)[num] - (*trainingMap)[row][col][num],2);
				}
				if(smallest == -999 || magSquared < smallest) {
					(winUnit)[0] = row;
					(winUnit)[1] = col;
					(winUnit)[2] = magSquared;
					smallest = magSquared;
				}
			}
		}
	} else {
		for(int row = 0; row < numRows; row++) {
            for(int col = 0; col < numCols; col++) {
                /*double magSquared = 0;
                for(int num = 0; num < trainingVector->size(); num++) {
                    if((*trainingVector)[num]!=0) {
                        magSquared += (*trainingVector)[num]*((*trainingVector)[num]-2*(*trainingMap)[row][col][num]);
                    }
                }
                for(int num = 0; num < trainingVector->size(); num++) {
                    magSquared += pow((*trainingMap)[row][col][num],2);
                }
				
				if(print) cout<<row<<'\t'<<col<<'\t'<<magSquared<<endl;*/
				double similarity=0;
                double mag1 = 0;
                double mag2 = 0;
                for(int num = 0; num < trainingVector->size(); num++) {
                    similarity += (*trainingVector)[num]*(*trainingMap)[row][col][num];
                    mag1 += pow((*trainingVector)[num],2);
                    mag2 += pow((*trainingMap)[row][col][num],2);
                }
				 if(mag1==0 || mag2==0) {
                    similarity = -1;
                } else {
                    similarity /= (sqrt(mag1) * sqrt(mag2));
                }

                double magSquared = 1-similarity;

                if(smallest == -999 || magSquared < smallest) {
                    (winUnit)[0] = row;
                    (winUnit)[1] = col;
					(winUnit)[2] = magSquared;
                    smallest = magSquared;
                }
            }
        }	
	}
    return winUnit;
}


int main(int argc, char *argv[]) {
	if(argc < 2) {
        cout << "Usage: ./scoresom -TrainingMatrix <Training Matrix File Location> -SOMFile <SOM File Location> -ScoreFile <Output Score File Location" <<endl;
        return 0;
    }
	string somFileName;
    string dataFileName;
    string scoreFileName;
	bool sparse = false;
	bool sub1 = false;
	for(int i = 0; i < argc; i++) {
        string temp = argv[i];
        if(temp.compare("-TrainingMatrix")==0)
			dataFileName=argv[i+1];
        if(temp.compare("-SOMFile")==0)
			somFileName=argv[i+1];
		if(temp.compare("-ScoreFile")==0)
			scoreFileName=argv[i+1];
		if(temp.compare("-Sub1")==0) 
			sub1=true;
		if(temp.compare("-Sparse")==0) 
			sparse=true;
	}
	int numRows = 0;
	int numCols = 0;

	cout<<"Opening SOM file"<<endl;
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
		numRows = row+1;
		istringstream(splitz[1])>>col;
		if(lastcol > col) {
			numCols = lastcol+1;
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
	cout<<"Opening data file"<<endl;
	int dimension = -1;
    map<string, vector<double> > dataMap;
    vector<string> dataKeys;
	ifstream trainingFile(dataFileName.c_str());
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
			if(sub1) {
                if(input==1) input=0;
                else input--;
            }
            temp.push_back(input);
        }
        dataMap[fields[0]]=temp;
        dataKeys.push_back(fields[0]);
        count++;
        if(count%100000==0)
        cout<<"line: "<<count<<endl;
    }
    cout<<"Done!"<<endl;
	cout<<"Scoring!"<<endl;
	// New winner map
    vector<vector<vector<string > > > winnerMap;
    vector<vector<vector<double > > > winnerDist;
    for(int row = 0; row < numRows; row ++) {
    	vector<vector<string > > temp2;
		vector<vector<double > > temp3;
        for(int col = 0; col < numCols; col++) {
        	vector<string> temp;
			vector<double> temp4;
            temp2.push_back(temp);
			temp3.push_back(temp4);
        }
        winnerMap.push_back(temp2);
		winnerDist.push_back(temp3);
    }
	ofstream scoreFile(scoreFileName.c_str());
    //Get winners
    cout<<"Getting winners"<<endl;
	for(int i = 0; i < dataKeys.size(); i++) {
		bool print = false;
		/*if(dataKeys[i].compare("ENSMUSG00000075370")==0) {
			print = true;
			cout<<"75370"<<endl;
		}
		if(dataKeys[i].compare("ENSMUSG00000059305")==0) {
			print = true;
			cout<<"59305"<<endl;
		}
		if(dataKeys[i].compare("chr16:16865355:16865480")==0) {
            print = true;
            cout<<"75370"<<endl;
        }
        if(dataKeys[i].compare("chr16:16869268:16869375")==0) {
            print = true;
            cout<<"59305"<<endl;
        }*/	
    	vector<double> winunit = propagate(&(dataMap[dataKeys[i]]),numRows,numCols, &inputMap,sparse,print);
    	winnerMap[(int)(winunit)[0]][(int)(winunit)[1]].push_back(dataKeys[i]);
		if(sparse) 
			winnerDist[(int)(winunit)[0]][(int)(winunit)[1]].push_back((winunit)[2]);
		else
			winnerDist[(int)(winunit)[0]][(int)(winunit)[1]].push_back((winunit)[2]/(double)dataMap[dataKeys[i]].size());
	}
    double totalScore = 0;
    //Count up score
    cout<<"Counting up score"<<endl;
    for(int row = 0; row < numRows; row++) {
    	for(int col = 0; col < numCols; col++) {
			scoreFile<<"unit\t"<<row<<','<<col;
			for(int i = 0; i < inputMap[row][col].size(); i++) {
				scoreFile<<'\t'<<inputMap[row][col][i];
			}
			scoreFile<<endl;
        	for(int k = 0; k < winnerMap[row][col].size(); k++) {
            	double tempScore = 0;
				if(!sparse) {
	                for(int u = 0; u < dataMap[winnerMap[row][col][k]].size(); u++) {
		            	tempScore+=pow(dataMap[winnerMap[row][col][k]][u]-inputMap[row][col][u],2);
			        }
					tempScore=sqrt(tempScore)/(double)dataMap[winnerMap[row][col][k]].size();
				} else {					 
					/*for(int u = 0; u < dataMap[winnerMap[row][col][k]].size(); u++) {
			            if(dataMap[winnerMap[row][col][k]][u]!=0) {
	                        tempScore += dataMap[winnerMap[row][col][k]][u]*(dataMap[winnerMap[row][col][k]][u]-2*inputMap[row][col][u]);
						}
					}
					for(int u = 0; u < inputMap[row][col].size(); u++) {
						tempScore += pow(inputMap[row][col][u],2);
					}*/
					double similarity=0;
	                double mag1 = 0;
		            double mag2 = 0;
			        for(int num = 0; num < dataMap[winnerMap[row][col][k]].size(); num++) {
				        similarity += (dataMap[winnerMap[row][col][k]])[num]*(inputMap)[row][col][num];
					    mag1 += pow((dataMap[winnerMap[row][col][k]])[num],2);
						mag2 += pow((inputMap)[row][col][num],2);
					}
					 if(mag1==0 || mag2==0) {
                    similarity = -1;
                } else {
                    similarity /= (sqrt(mag1) * sqrt(mag2));
                }

					double magSquared = 1-similarity;
					//cout<<tempScore;
					tempScore=magSquared;
				}
				scoreFile.precision(5);
				scoreFile<<tempScore<<'\t'<<winnerMap[row][col][k];
				for(int u = 0; u < dataMap[winnerMap[row][col][k]].size(); u++) {
					scoreFile<<'\t'<<dataMap[winnerMap[row][col][k]][u];
				}
				scoreFile<<'\t'<<winnerDist[row][col][k];
				scoreFile<<endl;
                totalScore += tempScore;
            }
			
        }
    }
	scoreFile.close();
    totalScore/=(double)dataKeys.size();
	cout<<"Total score: "<<totalScore<<endl;
		
}
