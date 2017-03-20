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


int** hexSurround(double input[], int radius, int numRows, int numCols, int* count) {
	int** result=0;
	(*count) = 0;
	for(int i = -1*radius; i <= radius; i++) {
        for(int j = max(-1*radius,-i-radius); j <= min(radius,-i+radius);j++) {
			(*count)++;
		}
	}
	result=new int*[(*count)];
	for(int i = 0; i < (*count); i++) {
		result[i]=new int[2];
	}
	int xs[(*count)];
	int ys[(*count)];
	int zs[(*count)];
	int x = input[1] - (input[0] - (abs((int)input[0])%2))/2;
	int z = input[0];
	int y = -x-z;
	int counter= 0;
	for(int i = -1*radius; i <= radius; i++) {
		for(int j = max(-1*radius,-i-radius); j <= min(radius,-i+radius);j++) {
			int dz = -i-j;
			xs[counter]=x+i;
			ys[counter]=y+j;
			zs[counter]=z+dz;
			counter++;
			//cout<<i<<'\t'<<j<<'\t'<<dz<<endl;
		}
	}
	for(int i = 0; i < (*count); i++) {
		result[i][0]=zs[i];
		if(result[i][0]<0) result[i][0]+=numRows;
		while(result[i][0]>=numRows) result[i][0]-=numRows;
		result[i][1]=xs[i]+(zs[i]-(abs(zs[i])%2))/2;
		if(result[i][1]<0) result[i][1]+=numCols;
		while(result[i][1]>=numCols) result[i][1]-=numCols;
	}
	return result;
}

double* propagate(double* trainingVector, int numRows, int numCols, int colsTraining, double*** trainingMap, bool sparse) {
	double* winUnit=new double[5];
    winUnit[0]=0;
    winUnit[1]=0;
    winUnit[2]=0;
    winUnit[3]=0;
    winUnit[4]=0;
	double secondsmallest = -999;
    double smallest = -999;
	if(!sparse) {
	    for(int row = 0; row < numRows; row++) {
		  	for(int col = 0; col < numCols; col++) {
				double magSquared=0;
				for(int num = 0; num < colsTraining; num++) {
					magSquared += pow(trainingVector[num] - trainingMap[row][col][num],2);
	           	}
				if(smallest == -999 || magSquared < smallest) {
                    winUnit[2] = winUnit[0];
                    winUnit[3] = winUnit[1];
                    secondsmallest = smallest;
                    winUnit[0] = row;
                    winUnit[1] = col;
                    smallest = magSquared;
                    winUnit[4] = smallest;
                } else if(secondsmallest == -999 | magSquared < secondsmallest) {
                    winUnit[2] = row;
                    winUnit[3] = col;
                    secondsmallest = magSquared;
                }

			}
		}
	} else {
		for(int row = 0; row < numRows; row++) {
            for(int col = 0; col < numCols; col++) {
				double similarity=0;
				double mag1 = 0;
				double mag2 = 0;
				for(int num = 0; num < colsTraining; num++) {
					similarity += trainingVector[num]*trainingMap[row][col][num];
					mag1 += pow(trainingVector[num],2);
					mag2 += pow(trainingMap[row][col][num],2);
				}
				if(mag1==0 || mag2==0) {
					similarity = -1;
				} else {
					similarity /= (sqrt(mag1) * sqrt(mag2));
				}
				double magSquared = 1-similarity;
				if(smallest == -999 || magSquared < smallest) {
					winUnit[2] = winUnit[0];
					winUnit[3] = winUnit[1];
					secondsmallest = smallest;
				    winUnit[0] = row;
			        winUnit[1] = col;
					smallest = magSquared;
					winUnit[4] = smallest;
		        } else if(secondsmallest == -999 | magSquared < secondsmallest) {
					winUnit[2] = row;
					winUnit[3] = col;
					secondsmallest = magSquared;
				}
			}
		}
	}	
	return winUnit;
}

int hexdist(int row1, int col1, int row2, int col2, int rows, int cols) {
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
		if(temp.compare("-Sparse")==0) 
			sparse=true;
		if(temp.compare("-LearningRate")==0)
			istringstream(argv[i+1])>>learningRate;
	}
	int dimension = -1;
	string line;
	int linesTraining=1;
	int colsTraining=0;
	cout<<"Getting Training Matrix Stats"<<endl;
	ifstream trainingFile(trainingFileName.c_str());
	getline(trainingFile,line);
	for(int i = 0; i < line.size(); i++) {
		if(line[i]=='\t') colsTraining++;
	}
	cout<<"Training Matrix data cols: "<<colsTraining<<endl;
	while(getline(trainingFile,line)) {
		if(line.size() < 2) continue;
        if(line[0]=='#') continue;
		linesTraining++;
	}
	trainingFile.close();
	cout<<"Training Matrix lines: "<<linesTraining<<endl;
	cout<<"Building data structures"<<endl;
	double** dataMap=0;
	dataMap = new double*[linesTraining];
	for(int i = 0; i < linesTraining; i++) {
		dataMap[i]=new double[colsTraining];
	}
	
	string dataKeys[linesTraining];
	int dataOrder[linesTraining];
	cout<<"Reading File: "<<trainingFileName<<endl;
	ifstream trainingFile2(trainingFileName.c_str());
	int lineCount = 0;
	while(getline(trainingFile2,line)) {
		if(line.size() < 2) continue;
		if(line[0]=='#') continue;
		vector<string> fields = split(line, '\t');
		for(int i = 1; i < fields.size(); i++) {
			istringstream(fields[i])>>dataMap[lineCount][i-1];
		}
		dataKeys[lineCount]=fields[0];
		dataOrder[lineCount]=lineCount;
		lineCount++;
		if(lineCount%10000==0) 
		cout<<"line: "<<lineCount<<endl;
	}	
	cout<<"Done!"<<endl;
	int radius;
    if(numCols > numRows)	radius = numCols/2;
    else					radius = numRows/2;

	for(int i = 0; i < trials; i++) {
		cout<<"Trial: "<<i+1<<endl;
		if(seed==-1)
			srand(time(0));
		else
			srand(seed);
		random_shuffle(dataOrder,dataOrder+linesTraining);
		// make new map
		double*** trainingMap=0;
		trainingMap=new double**[numRows];
		for(int j = 0; j < numRows; j++) 
			trainingMap[j]=new double*[numCols];
		for(int j = 0; j <numRows; j++) {
			for(int k = 0; k <numCols; k++) {
				trainingMap[j][k]=new double[colsTraining];
			}
		}
		cout<<"Building Training map"<<endl;
        for(int row = 0; row < numRows; row++) {
            for(int col = 0; col < numCols; col++) {
               	for(int k = 0; k < colsTraining; k++) {
               		trainingMap[row][col][k]=dataMap[dataOrder[row*numCols+col]][k];
               	}
            }
		}
		cout<<"done"<<endl;
		// train
		double multiplier = -1 * log(radius) / float(timesteps);
		cout<<"Beginning Training"<<endl;
		random_shuffle(dataOrder, dataOrder+linesTraining);
			 std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
		for(int j = 0; j < timesteps; j++) {
			int trainingNum = dataOrder[j%linesTraining];
			string trainingID = dataKeys[trainingNum];
			double trainingVector[colsTraining];
			for(int k = 0; k < colsTraining; k++) {
				trainingVector[k]=dataMap[trainingNum][k];
			}
			double decay = exp(j * multiplier);
			double timeLearningRate = learningRate * decay;
			int timeRadius = int(radius * decay);
//	std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
  //          std::cout << "1 - Time difference (sec) = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) /1000000.0 <<std::endl;
    //        begin = std::chrono::steady_clock::now();
			// Finding Winning Unit
			double* winUnit = propagate(trainingVector,numRows,numCols,colsTraining, trainingMap, sparse);
//end= std::chrono::steady_clock::now();
  //          std::cout << "2 - Time difference (sec) = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) /1000000.0 <<std::endl;
    //        begin = std::chrono::steady_clock::now();
			//cout<<"Found... Getting Neighbors: "<<winUnit[0]<<" "<<winUnit[1]<<endl;
			int neighborCount=0;
			int** winUnitNeighbors = hexSurround(winUnit,timeRadius,numRows, numCols,&neighborCount);
//end= std::chrono::steady_clock::now();
      //     std::cout << "3 - Time difference (sec) = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) /1000000.0 <<std::endl;
        //    begin = std::chrono::steady_clock::now();			
//cout<<"Got them.  Updating winning unit."<<endl;
			for(int k = 0; k < colsTraining; k++) {
				trainingMap[(int)(winUnit)[0]][(int)(winUnit)[1]][k] += timeLearningRate * (trainingVector[k] - trainingMap[(int)(winUnit)[0]][(int)(winUnit)[1]][k]);
			}
//end= std::chrono::steady_clock::now();
  //          std::cout << "4 - Time difference (sec) = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) /1000000.0 <<std::endl;
    //        begin = std::chrono::steady_clock::now();
			//cout<<"Changing neighbors"<<endl;
			for(int i = 0; i < neighborCount; i++) {
				if(winUnitNeighbors[i][0]!=winUnit[0]||winUnitNeighbors[i][1]!=winUnit[1]) {
					int unit[2];
					unit[0] = winUnitNeighbors[i][0];
					unit[1] = winUnitNeighbors[i][1];
					int dist = hexdist(unit[0],unit[1],winUnit[0],winUnit[1],numRows,numCols);
					float theWeight = timeLearningRate * exp(-0.5*pow(dist,2)/pow(radius+1,2));
					for(int k = 0; k < colsTraining; k++) {
						trainingMap[unit[0]][unit[1]][k] += theWeight*(trainingVector[k]-trainingMap[unit[0]][unit[1]][k]);
					}
				}
			}
			for(int k = 0; k < neighborCount; k++) {
				delete [] winUnitNeighbors[k];
			}
			delete [] winUnitNeighbors;
			delete [] winUnit;
  //          std::cout << "5 - Time difference (sec) = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) /1000000.0 <<std::endl;
    //        begin = std::chrono::steady_clock::now();
	//	int temp;
//	cin>>temp;
			if((j+1) % 100000 == 0) {
				std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
				std::cout << "Average timestep: (sec) = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) /(1000000.0*100000) <<std::endl;
            begin = std::chrono::steady_clock::now();
    //
				double totalScore = 0;
				double ta = 0;
		        for(int i = 0; i < linesTraining; i++) {
					double winVect[colsTraining];
					for(int k = 0; k < colsTraining; k++) {
					    winVect[k]=dataMap[i][k];
		            }

		            double* winUnit = propagate(winVect,numRows,numCols,colsTraining, trainingMap, sparse);
					totalScore += winUnit[4];
				    int dist = hexdist(winUnit[2],winUnit[3],winUnit[0],winUnit[1],numRows,numCols);
					if(dist<=1) ta++;
	            }
		        ta /= (double)linesTraining;
				totalScore/=(double)linesTraining;
				if(!sparse) totalScore /= colsTraining;
				if(totalScore < bestScore && totalScore!=0) {
					ofstream outfile(somFile.c_str());
					outfile<<"# "<<numRows<<" rows\t"<<numCols<<" cols\t"<<colsTraining<<" dimensions\ttoroid topology"<<endl;
					for(int row = 0; row < numRows; row++) {
						for(int col = 0; col <numCols; col++) {
							outfile<<row<<'\t'<<col<<'\t';
							for(int k = 0; k < colsTraining; k++) {
								outfile.precision(5);
								outfile<<fixed<<trainingMap[row][col][k];
								if(k != colsTraining-1)
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
/*				double ea = 0;
				for(int fe = 0; fe < colsTraining; fe++) {
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

					for(int datait = 0; datait < linesTraining; datait++) {
						meanData+=dataMap[datait][fe];
					}
					meanData /= (double)(linesTraining);
					for(int datait = 0; datait < linesTraining; datait++) {
				        varData+=pow(dataMap[datait][fe]-meanData,2);
					}
					varData /= (double)(linesTraining);
					if((meanMap-meanData)>(1.96*sqrt(varMap/(double)(numRows*numCols)+varData/(double)(linesTraining)))) {
						ea++;
					}
				}
				ea /= (double)colsTraining;
				cout<<"ea: "<<ea<<endl;
				cout<<"ta: "<<ta<<endl;*/
			}
		}
		for(int j = 0; j <numRows; j++) {
            for(int k = 0; k <numCols; k++) {
				delete [] trainingMap[j][k];
			}
		}
		for(int j = 0; j < numRows; j++) {
			delete [] trainingMap[j];
		}
		delete [] trainingMap;
	}
	cout<<"Best Score = "<<bestScore<<endl;
}
