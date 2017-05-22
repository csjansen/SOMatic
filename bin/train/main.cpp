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
#include <chrono>

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


vector <vector<int> > hexSurround(vector<double> input, int radius, int numRows, int numCols) {
	vector<vector<int> > result;
	vector<int> xs;
	vector<int> ys;
	vector<int> zs;
	int x = input[1] - (input[0] - (abs((int)input[0])%2))/2;
	int z = input[0];
	int y = -x-z;
	for(int i = -1*radius; i <= radius; i++) {
		for(int j = max(-1*radius,-i-radius); j <= min(radius,-i+radius);j++) {
			int dz = -i-j;
			xs.push_back(x+i);
			ys.push_back(y+j);
			zs.push_back(z+dz);
			//cout<<i<<'\t'<<j<<'\t'<<dz<<endl;
		}
	}
	for(int i = 0; i < xs.size(); i++) {
		vector<int> temp;
		temp.push_back(zs[i]);
		if(temp[0]<0) temp[0]+=numRows;
		while(temp[0]>=numRows) temp[0]-=numRows;
		temp.push_back(xs[i]+(zs[i]-(abs(zs[i])%2))/2);
		if(temp[1]<0) temp[1]+=numCols;
		while(temp[1]>=numCols) temp[1]-=numCols;
		result.push_back(temp);
	}
	return result;
	/*int factorCol = numRows/numCols+1;
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
	return result;*/
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

vector<double> propagate(vector<double>* trainingVector, int numRows, int numCols, vector<vector<vector<double> > >* trainingMap, bool sparse) {
	vector<double> winUnit;
    winUnit.push_back(0);
    winUnit.push_back(0);
	winUnit.push_back(0);
	winUnit.push_back(0);
	winUnit.push_back(0);
	double secondsmallest = -999;
    double smallest = -999;
	if(!sparse) {
	    for(int row = 0; row < numRows; row++) {
		  	for(int col = 0; col < numCols; col++) {
				double magSquared=0;
				for(int num = 0; num < trainingVector->size(); num++) {
					magSquared += pow((*trainingVector)[num] - (*trainingMap)[row][col][num],2);
	           	}
				if(smallest == -999 || magSquared < smallest) {
                    (winUnit)[2] = (winUnit)[0];
                    (winUnit)[3] = (winUnit)[1];
                    secondsmallest = smallest;
                    (winUnit)[0] = row;
                    (winUnit)[1] = col;
                    smallest = magSquared;
                    (winUnit)[4] = smallest;
                } else if(secondsmallest == -999 | magSquared < secondsmallest) {
                    (winUnit)[2] = row;
                    (winUnit)[3] = col;
                    secondsmallest = magSquared;
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
				}*/
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
					(winUnit)[2] = (winUnit)[0];
					(winUnit)[3] = (winUnit)[1];
					secondsmallest = smallest;
				    (winUnit)[0] = row;
			        (winUnit)[1] = col;
					smallest = magSquared;
					(winUnit)[4] = smallest;
		        } else if(secondsmallest == -999 | magSquared < secondsmallest) {
					(winUnit)[2] = row;
					(winUnit)[3] = col;
					secondsmallest = magSquared;
				}
			}
		}
	}	
	return winUnit;
}

int hexdist(int row1, int col1, int row2, int col2, int rows, int cols) {
	//row1=row1-row2;
	//col1=col1-col2;
	//row2=0;
	//col2=0;
	//while(row1>=rows) row1-=rows;
	//while(col1>=cols) col1-=cols;
	//while(row2>=rows) row2-=rows;
	//while(col2>=cols) col2-=cols;
	if(abs(row1-row2)>=rows/2.0) {
		if(row1<row2) {
			row1+=rows;
		} else {
			row2+=rows;
		}
	}
	if(abs(col1-col2)>=cols/2.0) {
        if(col1<col2) {
            col1+=cols;
        } else {
            col2+=cols;
        }
    }
	//Convert odd-r to cube
	int x1 = col1 - (row1 - (row1&1))/2;
	int z1 = row1;
	int y1 = -x1-z1;
	int x2 = col2 - (row2 - (row2&1))/2;
    int z2 = row2;
    int y2 = -x2-z2;
	int dist1 = (abs(x1-x2) + abs(y1-y2) + abs(z1-z2))/2;
	//x1 = col1 - (row1 - (row1&1))/2;
	//	39      16      0       18      39
//	cout<<row1<<'\t'<<col1<<'\t'<<row2<<'\t'<<col2<<'\t'<<'\t'<<rows<<'\t'<<cols<<'\t'<<dist1<<endl;
//	int temp;
//	cin>>temp;
	return (abs(x1-x2) + abs(y1-y2) + abs(z1-z2))/2;
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
		cout << "-LearningRate: Set Learning Rate <.2>"<<endl; 
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
	bool sparse=false;
	bool sub1=false;
	float learningRate = 0.2;
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
		if(temp.compare("-Sub1")==0) 
			sub1=true;
		if(temp.compare("-Sparse")==0) 
			sparse=true;
		if(temp.compare("-LearningRate")==0)
			istringstream(argv[i+1])>>learningRate;
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
			if(sub1) {
				if(input==1) input=0;
				else input--;
			}
			temp.push_back(input);
		}
		dataMap[fields[0]]=temp;
		dataKeys.push_back(fields[0]);
		count++;
		if(count%10000==0) 
		cout<<"line: "<<count<<endl;
	}	
	cout<<"Done! "<<dataKeys.size()<<endl;
	int radius;
        if(numCols > numRows)   radius = numCols/2;
        else            radius = numRows/2;

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
		cout<<"datakeys "<<dataKeys.size()<<endl;
		random_shuffle(dataKeys.begin(), dataKeys.end());
		cout<<"Timesteps "<<timesteps<<endl;
		for(int j = 0; j < timesteps; j++) {
			cout<<j<<endl;
			std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
			string trainingID = dataKeys[j % dataKeys.size()];
			vector<double> trainingVector = dataMap[trainingID];
			double decay = exp(j * multiplier);
			double timeLearningRate = learningRate * decay;
			int timeRadius = int(radius * decay);
	std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
            std::cout << "1 - Time difference (sec) = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) /1000000.0 <<std::endl;
            begin = std::chrono::steady_clock::now();
			// Finding Winning Unit
			//cout<<"Finding Winning Unit"<<endl;
			vector<double> winUnit = propagate(&(trainingVector),numRows,numCols, &trainingMap, sparse);
			end= std::chrono::steady_clock::now();
            std::cout << "1 - Time difference (sec) = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) /1000000.0 <<std::endl;
            begin = std::chrono::steady_clock::now();
			//cout<<"Found... Getting Neighbors: "<<winUnit[0]<<" "<<winUnit[1]<<endl;
			//vector<vector<int> > winUnitNeighbors = hexSurround(winUnit,timeRadius,numRows, numCols);
			vector<vector<int> > winUnitNeighbors = hexSurround(winUnit,2,numRows, numCols);
			end= std::chrono::steady_clock::now();
            std::cout << "1 - Time difference (sec) = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) /1000000.0 <<std::endl;
            begin = std::chrono::steady_clock::now();
			//cout<<"Got them.  Updating winning unit."<<endl;
			for(int k = 0; k < trainingVector.size(); k++) {
				trainingMap[(int)(winUnit)[0]][(int)(winUnit)[1]][k] += timeLearningRate * (trainingVector[k] - trainingMap[(int)(winUnit)[0]][(int)(winUnit)[1]][k]);
			}
			end= std::chrono::steady_clock::now();
            std::cout << "1 - Time difference (sec) = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) /1000000.0 <<std::endl;
            begin = std::chrono::steady_clock::now();
			//cout<<"Changing neighbors"<<endl;
//			for(int dist = 1; dist<timeRadius+1;dist++) {
//				float theWeight = timeLearningRate * exp(-0.5*pow(dist,2)/pow(radius+1,2));
				for(int i = 0; i < winUnitNeighbors.size(); i++) {
					if(winUnitNeighbors[i][0]!=winUnit[0]||winUnitNeighbors[i][1]!=winUnit[1]) {
						vector<int> unit = winUnitNeighbors[i];
						/*int xdist = (unit[0] - winUnit[0]);
						int ydist = (unit[1] - winUnit[1]);
						int ddist = abs(xdist-ydist);
						int dist = 0;
						if(abs(xdist) > abs(ydist) && abs(xdist) > abs(ddist)) dist = abs(xdist);
						if(abs(ydist) > abs(xdist) && abs(ydist) > abs(ddist)) dist = abs(ydist);
						if(abs(ddist) > abs(xdist) && abs(ddist) > abs(ydist)) dist = abs(ddist);*/
						int dist = hexdist(unit[0],unit[1],winUnit[0],winUnit[1],numRows,numCols);
						//cout<<unit[0]<<'\t'<<unit[1]<<'\t'<<winUnit[0]<<'\t'<<winUnit[1]<<'\t'<<dist;
						//int temp;
						//cin>>temp;
						float theWeight = timeLearningRate * exp(-0.5*pow(dist,2)/pow(radius+1,2));
						for(int k = 0; k < trainingVector.size(); k++) {
							trainingMap[unit[0]][unit[1]][k] += theWeight*(trainingVector[k]-trainingMap[unit[0]][unit[1]][k]);
						}
					}
				}
			end= std::chrono::steady_clock::now();
            std::cout << "1 - Time difference (sec) = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) /1000000.0 <<std::endl;
//			}
			//cout<<"done"<<endl;
			//Score
			// New winner map
			/*vector<vector<vector<string > > > winnerMap;
            for(int row = 0; row < numRows; row ++) {
				vector<vector<string > > temp2;
                   	for(int col = 0; col < numCols; col++) {
					vector<string> temp;
					temp2.push_back(temp);
				}
				winnerMap.push_back(temp2);
			}
			
			//Get winners	
			winnerMap[(winunit)[0]][(winunit)[1]].push_back(trainingID);
			double totalScore = 0;
			//Count up score	
			for(int row = 0; row < numRows; row++) {
				for(int col = 0; col < numCols; col++) {
					
					for(int k = 0; k < winnerMap[row][col].size(); k++) {
						/*double tempScore = 0;
						for(int u = 0; u < dataMap[winnerMap[row][col][k]].size(); u++) {
							tempScore+=pow(dataMap[winnerMap[row][col][k]][u]-trainingMap[row][col][u],2);
						}*/
			int temp;
		cin>>temp;
			if((j+1) % 10000 == 0) {
			double totalScore = 0;
			double ta = 0;
            for(int i = 0; i < dataKeys.size(); i++) {
                vector<double> winUnit = propagate(&(dataMap[dataKeys[i]]),numRows,numCols, &trainingMap, sparse);
				//if(winUnit[4]<0)
				//	cout<<winUnit[4]<<endl;
				totalScore += winUnit[4];
				//cout<<totalScore<<endl;
                int dist = hexdist(winUnit[2],winUnit[3],winUnit[0],winUnit[1],numRows,numCols);
                if(dist<=1) ta++;
            }
            ta /= (double)dataKeys.size();

					/*double similarity=0;
					double mag1 = 0;
				    double mag2 = 0;
						for(int num = 0; num < dataMap[winnerMap[row][col][k]].size(); num++) {
							similarity += dataMap[winnerMap[row][col][k]][num]*trainingMap[row][col][num];
							mag1 += pow((dataMap[winnerMap[row][col][k]])[num],2);
							mag2 += pow((trainingMap)[row][col][num],2);
						}
						similarity /= (sqrt(mag1) * sqrt(mag2));
						double magSquared = 1-similarity;

						//totalScore += sqrt(tempScore);
						totalScore += magSquared;
					}
				}
			}*/
			totalScore/=(double)dataKeys.size();
			cout<<dataKeys.size()<<endl;
			if(!sparse) totalScore /= trainingVector.size();
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
				cout << "time " << j+1 << " radius " << timeRadius<< " learn " << timeLearningRate<<endl;
				cout << "Total Score: "<<totalScore<<endl;;	
			// Find Embeded feature
			double ea = 0;
			for(int fe = 0; fe < trainingVector.size(); fe++) {
				double meanMap=0;
				double varMap=0;
				double meanData=0;
				double varData=0;
				for(int row = 0; row < numRows; row++) {
					for(int col = 0; col < numCols; col++) {
						meanMap+=trainingMap[row][col][fe];
					}
				}
				meanMap /= (double)(numRows*numCols);
				for(int row = 0; row < numRows; row++) {
                    for(int col = 0; col < numCols; col++) {
                        varMap+=pow(trainingMap[row][col][fe]-meanMap,2);
                    }
                }
				varMap /= (double)(numRows*numCols);
				//string trainingID = dataKeys[j % dataKeys.size()];
	            //vector<double> trainingVector = dataMap[trainingID];

				for(int datait = 0; datait < dataKeys.size(); datait++) {
					vector<double> vec = dataMap[dataKeys[datait]];
					meanData+=vec[fe];
				}
				meanData /= (double)(dataKeys.size());
				for(int datait = 0; datait < dataKeys.size(); datait++) {
                    vector<double> vec = dataMap[dataKeys[datait]];
                    varData+=pow(vec[fe]-meanData,2);
                }
				varData /= (double)(dataKeys.size());
				if((meanMap-meanData)>(1.96*sqrt(varMap/(double)(numRows*numCols)+varData/(double)(dataKeys.size())))) {
					ea++;
				}
			}
			ea /= (double)trainingVector.size();
			cout<<"ea: "<<ea<<endl;
            cout<<"ta: "<<ta<<endl;
		}
		}
	}
	cout<<"Best Score = "<<bestScore<<endl;
}
