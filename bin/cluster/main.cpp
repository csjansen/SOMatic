/* cluster: Does hierarchical clustering on SOMs.
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
#include <string>
#include <cmath>
using namespace std;

class clustervertex {
public: 
	vector<vector<int> > contents;
	vector<int> name;
	vector<double> level;
};

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

int main(int argc, char *argv[]) {
	if(argc < 2) {
        cout << "Usage: ./cluster -SOMFile <SOM File Location> -Clusters1 <Output file with clusters on SOM units> -Clusters2 <Output file with clusters on profiles> -NoNormalize -Col <Number of Cols in SOM" <<endl;
        return 0;
    }
	string somFileName;
    string outFileName1;
    string outFileName2;
	int cols=60;
	bool Normalize = true;
	string DistanceMetric="Euclid";
	for(int i = 0; i < argc; i++) {
        string temp = argv[i];
        if(temp.compare("-Clusters1")==0)
			outFileName1=argv[i+1];
        if(temp.compare("-SOMFile")==0)
			somFileName=argv[i+1];
		if(temp.compare("-Clusters2")==0)
			outFileName2=argv[i+1];
		if(temp.compare("-NoNormalize")==0)
			Normalize = false;
		if(temp.compare("-Col")==0)
			istringstream(argv[i+1])>>cols;
		if(temp.compare("-DistanceMetric")==0)
			DistanceMetric=argv[i+1];
	}
	int numRows = 0;
	int numCols = cols;

	cout<<"Opening SOM file"<<endl;
	ifstream somFile(somFileName.c_str());
	string line;
	vector<vector<vector<double> > > inputMap;
	int lastcol;
	vector<vector<double> > temp;
	int index = 0;
	while(getline(somFile, line)) {
		if(line[0] == '#') {
			continue;
		}
		vector<string> splitz = split(line,'\t');
		if(index >= cols) {
			inputMap.push_back(temp);
			temp.erase(temp.begin(),temp.end());
			index=0;
			numRows++;
		}
		vector<double> temp2;
		for(int i = 0; i < splitz.size(); i++) {
			double entry;
			istringstream(splitz[i])>>entry;
			temp2.push_back(entry);
		}
		temp.push_back(temp2);
		index++;
	}
	inputMap.push_back(temp);
	numRows++;
	if(Normalize) {
		cout<<"Normalizing"<<endl;
		for(int i = 0; i < inputMap.size(); i++) {
			for(int j = 0; j < inputMap[i].size(); j++) {
				double total;
				for(int k = 0; k < inputMap[i][j].size(); k++) {
					total+=inputMap[i][j][k];
				}
				total /= inputMap[i][j].size();
				for(int k = 0; k < inputMap[i][j].size(); k++) {
					inputMap[i][j][k]/=total;
				}
			}
		}
	}
	//Clustering on Units
	vector<clustervertex> Clusters;
	cout<<"Clustering on Units"<<endl;
	for(int i=0; i < inputMap.size(); i++) {
		for(int j=0; j < inputMap[i].size(); j++) {
			clustervertex temp;
			temp.name.push_back(i*numCols+j);
			vector<int> tempcoord;
			tempcoord.push_back(i);
			tempcoord.push_back(j);
			temp.contents.push_back(tempcoord);
			Clusters.push_back(temp);
			//temp.level.push_back(0);
		}
	}
	vector<vector<double> > distances;
	/*for(int i = 0; i < Clusters.size(); i++) {
		vector<double> temp;
		
        for(int j = i+1; j<Clusters.size(); j++) {
			double dist = 0;
			temp.push_back(dist);
		}
		distances.push_back(temp);
	}*/
	for(int i = 0; i < Clusters.size(); i++) {
		vector<double> temp;
		for(int j = 0; j<Clusters.size(); j++) {
			if(DistanceMetric == "Euclid") {
				double dist = 0;
				for(int q = 0; q<inputMap[Clusters[i].contents[0][0]][Clusters[i].contents[0][1]].size(); q++) {
                			dist += pow(inputMap[Clusters[i].contents[0][0]][Clusters[i].contents[0][1]][q]-inputMap[Clusters[j].contents[0][0]][Clusters[j].contents[0][1]][q],2);
				}
				dist /= inputMap[Clusters[i].contents[0][0]][Clusters[i].contents[0][1]].size();
				dist=sqrt(dist);
				temp.push_back(dist);
			} else if(DistanceMetric == "Pearson") {
                                double Ex = 0;
                                double Ey = 0;
				for(int q = 0; q<inputMap[Clusters[i].contents[0][0]][Clusters[i].contents[0][1]].size(); q++) {
                                
                                        Ex += inputMap[Clusters[i].contents[0][0]][Clusters[i].contents[0][1]][q];
                                        Ey += inputMap[Clusters[j].contents[0][0]][Clusters[j].contents[0][1]][q];
                                }
                                Ex /= inputMap[Clusters[i].contents[0][0]][Clusters[i].contents[0][1]].size();
                                Ey /= inputMap[Clusters[i].contents[0][0]][Clusters[i].contents[0][1]].size();
                                double stdDevX=0;
                                double stdDevY=0;
                                double cov = 0;
				for(int q = 0; q<inputMap[Clusters[i].contents[0][0]][Clusters[i].contents[0][1]].size(); q++) {
                                        stdDevX += pow(inputMap[Clusters[i].contents[0][0]][Clusters[i].contents[0][1]][q]-Ex,2);
                                        stdDevY += pow(inputMap[Clusters[j].contents[0][0]][Clusters[j].contents[0][1]][q]-Ey,2);
                                        cov += (inputMap[Clusters[i].contents[0][0]][Clusters[i].contents[0][1]][q]-Ex)*(inputMap[Clusters[j].contents[0][0]][Clusters[j].contents[0][1]][q]-Ey);
                                }
                                stdDevX = sqrt(stdDevX/inputMap[Clusters[i].contents[0][0]][Clusters[i].contents[0][1]].size());
                                stdDevY = sqrt(stdDevY/inputMap[Clusters[i].contents[0][0]][Clusters[i].contents[0][1]].size());
                                cov = cov/inputMap[Clusters[i].contents[0][0]][Clusters[i].contents[0][1]].size();

                                double magSquared = 1-(cov/(stdDevX*stdDevY));
				temp.push_back(magSquared);
			}
		}
		distances.push_back(temp);
	}

	while(Clusters.size() > 1) {
		cout<<Clusters.size()<<" Clusters"<<endl;
		int mini = -1;
		int minj = -1;
		double mindist = -1;
		for(int i = 0; i < Clusters.size(); i++) {
			for(int j = i+1; j<Clusters.size(); j++) {
				double averagedist = 0;
				int numdist = 0;
				int mink = -1;
				int minp = -1;	
				for(int k = 0; k<Clusters[i].name.size(); k++) {
					for(int p = 0; p<Clusters[j].name.size(); p++) {
						double dist = 0;
						if(Clusters[j].name[p]==Clusters[i].name[k]) {
							int temp;
							cin>>temp;
						}
						//if(Clusters[i].name[k] < Clusters[j].name[p]) 
						if(Clusters[i].name[k]>=distances.size())
							cout<<i<<'\t'<<j<<'\t'<<k<<'\t'<<p<<'\t'<<Clusters[i].name[k]<<endl;
							dist = distances[Clusters[i].name[k]][Clusters[j].name[p]];
						//else 
						//	dist = distances[Clusters[j].name[p]][Clusters[i].name[k]];
						averagedist+=dist;
						numdist++;
						/*if(dist < mintempdist) {
							mintempdist = dist;
							mink = k;
							minp = p;
						}*/
					}
				}
				if(numdist>0)
					averagedist/=numdist;
				if(averagedist < mindist || mindist==-1) {
					//cout<<i<<'\t'<<j<<'\t'<<mink<<'\t'<<minp<<endl;
					/*for(int q = 0; q<inputMap[Clusters[i].contents[0][0]][Clusters[i].contents[0][1]].size(); q++) {
						cout<<inputMap[Clusters[i].contents[0][0]][Clusters[i].contents[0][1]][q]<<'\t';
					}
					cout<<endl;
					for(int q = 0; q<inputMap[Clusters[j].contents[0][0]][Clusters[j].contents[0][1]].size(); q++) {
                        cout<<inputMap[Clusters[j].contents[0][0]][Clusters[j].contents[0][1]][q]<<'\t';
                    }
*/
					mindist = averagedist;
					mini=i;
					minj=j;
				}				
			}
		}
		if(mini==-1||minj==-1) {
			cout<<"-1s?"<<endl;
			cout<<Clusters.size()<<endl;
			cout<<mindist<<endl;
		}
		cout<<"merging "<<mini<<'\t'<<minj<<endl;	
		clustervertex temp;
		for(int i = 0; i < Clusters[mini].contents.size(); i++) {
			temp.name.push_back(Clusters[mini].name[i]);
			//cout<<Clusters[mini].name[i]<<endl;
			temp.contents.push_back(Clusters[mini].contents[i]);
		}
		for(int i = 0; i < Clusters[mini].level.size(); i++) {
			temp.level.push_back(Clusters[mini].level[i]);
		}
		temp.level.push_back(mindist);
		//cout<<endl;
		for(int i = 0; i < Clusters[minj].contents.size(); i++) {
			temp.name.push_back(Clusters[minj].name[i]);
			//cout<<Clusters[mini].name[i]<<endl;
			temp.contents.push_back(Clusters[minj].contents[i]);
		}
		for(int i = 0; i < Clusters[minj].level.size(); i++) {
			temp.level.push_back(Clusters[minj].level[i]);
		}
		Clusters.erase(Clusters.begin()+mini);
		Clusters.erase(Clusters.begin()+minj-1);
		Clusters.push_back(temp);
	}
	ofstream outFile1(outFileName1.c_str());
	cout<<"Writing "<<outFileName1<<endl;
	for(int i = 0; i < Clusters[0].name.size(); i++) {
		if(i != Clusters[0].name.size()-1)
			outFile1<<Clusters[0].name[i]<<'\t'<<Clusters[0].level[i]<<endl;
		else
			outFile1<<Clusters[0].name[i]<<'\t'<<0<<endl;
			
	}
	outFile1.close();
	Clusters.clear();
    cout<<"Clustering on Profiles"<<endl;
	for(int k = 0; k < inputMap[0][0].size(); k++) {
		clustervertex temp;
        temp.name.push_back(k);
        vector<int> tempcoord;
		tempcoord.push_back(k);
		temp.contents.push_back(tempcoord);
		Clusters.push_back(temp);
	}
    distances.clear();

for(int i = 0; i < Clusters.size(); i++) {
        vector<double> temp;
        for(int j = 0; j<Clusters.size(); j++) {
            double dist = 0;
			for(int p = 0; p<inputMap.size(); p++) {
				for(int q = 0; q<inputMap[p].size(); q++) {
					dist += pow(inputMap[p][q][Clusters[i].contents[0][0]]-inputMap[p][q][Clusters[j].contents[0][0]],2);
				}
			}
			dist=sqrt(dist);
            temp.push_back(dist);
        }
        distances.push_back(temp);
    }



    /*for(int i = 0; i < Clusters.size(); i++) {
        vector<double> temp;
        for(int j = i+1; j<Clusters.size(); j++) {
            double dist = 0;
            temp.push_back(dist);
        }
        distances.push_back(temp);
    }
    for(int i = 0; i < Clusters.size(); i++) {
        for(int j = i+1; j<Clusters.size(); j++) {
            double dist = 0;
            for(int p = 0; p<inputMap.size(); p++) {
				for(int q = 0; q < inputMap[p].size(); q++) {
                	dist += abs(inputMap[p][q][Clusters[i].contents[0][0]]-inputMap[p][q][Clusters[j].contents[0][0]]);
				}
            }
            distances[i][j-i-1]=dist;
        }
    }*/
	while(Clusters.size() > 1) {
        cout<<Clusters.size()<<" Clusters"<<endl;
        int mini = -1;
        int minj = -1;
        double mindist = -1;
        for(int i = 0; i < Clusters.size(); i++) {
            for(int j = i+1; j<Clusters.size(); j++) {
                double averagedist = 0;
                int numdist = 0;
                int mink = -1;
                int minp = -1;
                for(int k = 0; k<Clusters[i].name.size(); k++) {
                    for(int p = 0; p<Clusters[j].name.size(); p++) {
                        double dist = 0;
                        if(Clusters[j].name[p]==Clusters[i].name[k]) {
                            cout<<j<<'\t'<<p<<'\t'<<i<<'\t'<<k<<'\t'<<Clusters[j].name[p]<<endl;
                            int temp;
                            cin>>temp;
                        }
                        //if(Clusters[i].name[k] < Clusters[j].name[p])
                            dist = distances[Clusters[i].name[k]][Clusters[j].name[p]];
                        //else
                        //    dist = distances[Clusters[j].name[p]][Clusters[i].name[k]];
                        averagedist+=dist;
                        numdist++;
                        /*if(dist < mintempdist) {
                            mintempdist = dist;
                            mink = k;
                            minp = p;
                        }*/
                    }
                }
				if(numdist>0)
                averagedist/=numdist;
                if(averagedist < mindist || mindist == -1) {
                    //cout<<mintempdist<<'\t'<<i<<'\t'<<j<<'\t'<<mink<<'\t'<<minp<<endl;
                    /*for(int q = 0; q<inputMap[Clusters[i].contents[0][0]][Clusters[i].contents[0][1]].size(); q++) {
                        cout<<inputMap[Clusters[i].contents[0][0]][Clusters[i].contents[0][1]][q]<<'\t';
                    }
                    cout<<endl;
                    for(int q = 0; q<inputMap[Clusters[j].contents[0][0]][Clusters[j].contents[0][1]].size(); q++) {
                        cout<<inputMap[Clusters[j].contents[0][0]][Clusters[j].contents[0][1]][q]<<'\t';
                    }
*/
                    mindist = averagedist;
                    mini=i;
                    minj=j;
                }
            }
        }
	cout<<"merging "<<mini<<'\t'<<minj<<endl;
        clustervertex temp;
        for(int i = 0; i < Clusters[mini].contents.size(); i++) {
            temp.name.push_back(Clusters[mini].name[i]);
            //cout<<Clusters[mini].name[i]<<endl;
            temp.contents.push_back(Clusters[mini].contents[i]);
        }
        for(int i = 0; i < Clusters[mini].level.size(); i++) {
            temp.level.push_back(Clusters[mini].level[i]);
        }
        temp.level.push_back(mindist);
        //cout<<endl;
        for(int i = 0; i < Clusters[minj].contents.size(); i++) {
            temp.name.push_back(Clusters[minj].name[i]);
            //cout<<Clusters[mini].name[i]<<endl;
            temp.contents.push_back(Clusters[minj].contents[i]);
        }
        for(int i = 0; i < Clusters[minj].level.size(); i++) {
            temp.level.push_back(Clusters[minj].level[i]);
        }
        Clusters.erase(Clusters.begin()+mini);
        Clusters.erase(Clusters.begin()+minj-1);
        Clusters.push_back(temp);
    }
    ofstream outFile2(outFileName2.c_str());
	cout<<"Writing "<<outFileName2<<endl;
    for(int i = 0; i < Clusters[0].name.size(); i++) {
        if(i != Clusters[0].name.size()-1)
            outFile2<<Clusters[0].name[i]<<'\t'<<Clusters[0].level[i]<<endl;
        else
            outFile2<<Clusters[0].name[i]<<'\t'<<0<<endl;

    }
    outFile2.close();

}
