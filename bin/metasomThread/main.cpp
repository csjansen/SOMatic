#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <vector>
#include <array>
#include <queue>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string>
#include <locale>
#include <ctime>
#include <thread>
using namespace std;

//#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
//        ( std::ostringstream() << std::dec << x ) ).str()
string SSTR(int a)
{
    ostringstream temp;
    temp << a;
    return temp.str();
}

string SSTRF(int a)
{
    ostringstream temp;
    temp << a;
    return temp.str();
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

void ClusterAndScore(vector<vector<float> >* allpoints, vector<vector<float> >* Som1, int Dimensionality, int clusternum, string DistanceMetric, vector<int>* clusterIndex, float* score, int row1, int col1, int trial) {
	vector<vector<float> > allpoints2=*allpoints;
	// Make random shuffle
	vector<int> index;
	for(int i = 0; i < allpoints->size(); i++) {
		//cout<<i<<endl;
		index.push_back(i);
	}
	random_shuffle(index.begin(),index.end());
	vector<vector<float> > MidPoints;
	vector<float> Rads;
        //cout<<"Assignment: "<<trial<<endl;
    // Initial Mid Points
    int uniqs = 0;
    for(int k = 0; k < clusternum+uniqs; k++) {
    	vector<float> temp;
        for(int m = 0; m < Som1->at(0).size(); m++) {
        	temp.push_back(Som1->at((int)(allpoints->at(index[k])[0])*col1+(int)(allpoints->at(index[k])[1]))[m]);
        }
        bool found = false;
        for(int m = 0; m < MidPoints.size(); m++) {
        	bool inMidpoint = true;
            for(int n = 0; n < temp.size(); n++) {
            	if(temp[n]!=MidPoints[m][n]) {
                	inMidpoint=false;
                    break;
                }
            }
            if(inMidpoint) {
            	found = true;
                break;
            }
        }
        if(!found) {
        	MidPoints.push_back(temp);
            Rads.push_back(0);
        } else {
        	uniqs++;
		}
    }
	bool finished = false;
	bool finished2 = false;
    int iter = 0;
    bool cont = false;
	vector<vector<vector<float> > > points;
    vector<vector<vector<float> > > points2;
    vector<vector<float> > points3;
	int pointlist[row1][col1];
    int Nc[allpoints->size()];
	int K = clusternum;
    int N = allpoints->size();
    int P = Som1->at(0).size();
	vector<vector<float> > DistMatrix;
	while(!finished2) {
    while(!finished) {
        int count=0;
        //Assignment
        for(int k = 0; k < allpoints->size(); k++) {
        	float lowestdist1=999999999;
            float lowestdist=99999999;
            float lowestdist2=999999999;
            int lowestk1=-1;
            int lowestk2=-1;
            for(int m = 0; m < clusternum; m++) {
		float dist = 0;
		if(DistanceMetric == "Euclid") {
	                for(int num = 0; num < Som1->at(0).size(); num++) {
                        	dist += pow(MidPoints[m][num]-Som1->at((int)(allpoints->at(k)[0])*col1+(int)(allpoints->at(k)[1]))[num],2);
                	}
			dist = dist/Som1->at(0).size();
                	dist = sqrt(dist);	
		} else if(DistanceMetric == "Pearson") {
			double Ex = 0;
                        double Ey = 0;
			for(int num = 0; num < Som1->at(0).size(); num++) {
                        	Ex += MidPoints[m][num];
                                Ey += Som1->at((int)(allpoints->at(k)[0])*col1+(int)(allpoints->at(k)[1]))[num];
                        }
                        Ex /= Som1->at(0).size();
                        Ey /= Som1->at(0).size();
                        double stdDevX=0;
                        double stdDevY=0;
                        double cov = 0;
			for(int num = 0; num < Som1->at(0).size(); num++) {
	                        stdDevX += pow(MidPoints[m][num]-Ex,2);
                                stdDevY += pow(Som1->at((int)(allpoints->at(k)[0])*col1+(int)(allpoints->at(k)[1]))[num]-Ey,2);
                                cov += (MidPoints[m][num]-Ex)*(Som1->at((int)(allpoints->at(k)[0])*col1+(int)(allpoints->at(k)[1]))[num]-Ey);
                        }
                        stdDevX = sqrt(stdDevX/Som1->at(0).size());
                        stdDevY = sqrt(stdDevY/Som1->at(0).size());
                        cov = cov/Som1->at(0).size();

                        dist = 1-(cov/(stdDevX*stdDevY));

		}
                if(lowestdist1>dist) {
                        lowestdist1 = dist;
                    lowestk1 = m;
                    lowestdist=dist;
                }

		
            	/*float similarity=0;
                float mag1 = 0;
                float mag2 = 0;
                float dist = 0;
                for(int num = 0; num < Som1->at(0).size(); num++) {
                	dist += pow(MidPoints[m][num]-Som1->at((int)(allpoints->at(k)[0])*col1+(int)(allpoints->at(k)[1]))[num],2);
                    similarity += MidPoints[m][num]*Som1->at((int)(allpoints->at(k)[0])*col1+(int)(allpoints->at(k)[1]))[num];
                    mag1 += pow(MidPoints[m][num],2);
                    mag2 += pow(Som1->at((int)(allpoints->at(k)[0])*col1+(int)(allpoints->at(k)[1]))[num],2);
                }
                dist = sqrt(dist);

                if(mag1==0 || mag2==0) {
                	similarity = -1;
                } else {
                    similarity /= (sqrt(mag1) * sqrt(mag2));
                }
                float dist1;
                if(sparse)
                	dist1 = 1-similarity;
                else
                    dist1 = dist;
                if(lowestdist1>dist1) {
                	lowestdist1 = dist1;
                    lowestk1 = m;
                    lowestdist=dist1;
                }*/
            }
            if(allpoints2.at(k)[2] != lowestk1) {
            	count++;
                allpoints2.at(k)[2] = lowestk1;
                allpoints2.at(k)[3] = lowestdist1;
                finished = false;
            }
            if(lowestdist > Rads[lowestk1]) Rads[lowestk1]=lowestdist;
       	}
		for(int k = 0; k < clusternum; k++) {
			int count = 0;
			for(int m = 0; m < allpoints->size(); m++) {
                if((int)(allpoints2.at(m)[2])==k) {
                    count++;
                }
            }
			if(count==0) {
				int index = rand() % (allpoints->size());
				for(int n = 0; n < Som1->at(0).size(); n++) {
					MidPoints[k][n]=Som1->at((int)(allpoints->at(index)[0])*col1+(int)(allpoints->at(index)[1]))[n];
				}
				allpoints2.at(index)[2]=k;
				finished=false;
			}
		}
    //      cout<<count<<endl;
       	if(finished) break;
		for(int k = 0; k < Rads.size(); k++) {
        	Rads[k]=0;
        }
		//cout<<"Update: "<<trial<<endl;
        cont=true;
		for(int k = 0; k < clusternum; k++) {
		//	cout<<k<<endl;
        	vector<float> totals1;
            int count = 0;
            for(int m = 0; m < Som1->at(0).size(); m++) totals1.push_back(0);
            for(int m = 0; m < allpoints->size(); m++) {
            	if((int)(allpoints2.at(m)[2])==k) {
                	for(int n = 0; n < Som1->at(0).size(); n++) totals1[n]+=Som1->at((int)(allpoints->at(m)[0])*col1+(int)(allpoints->at(m)[1]))[n];
                   	count++;
                }
            }
		//	cout<<count<<endl;
            for(int m = 0; m < totals1.size(); m++) {
				if(count==0)
					MidPoints[k][m]=0;
				else
	                MidPoints[k][m] = totals1[m]/(float)count;
                //      cout<<MidPoints[k][m]<<endl;
            }
        }
		//cout<<"Done"<<endl;
    	finished = true;
    }
//    cout<<"Cluster: "<<trial<<endl;
        //rebuild clusters to be adjacient
        //Calculating Distance Matrix
//	cout<<"Dist: "<<trial<<endl;
	DistMatrix.clear();
    //float DistMatrix[N][K];
//	cout<<"Dist: "<<trial<<endl;
    for(int n = 0; n < N; n++) {
		vector<float> DistMatrixRow;
    	for(int k = 0; k < K; k++) {
		//	DistMatrixRow.push_back(0);
/*           	if(!sparse) {
            	for(int p = 0; p < P; p++) {
                	DistMatrixRow[k]+=pow(Som1->at((int)(allpoints->at(n)[0])*col1+(int)(allpoints->at(n)[1]))[p]-MidPoints[k][p],2);
            	}
				DistMatrixRow[k]=sqrt(DistMatrixRow[k]);
            } else {
                float similarity=0;
                float mag1 = 0;
                float mag2 = 0;
                for(int p = 0; p < P; p++) {
                	similarity += MidPoints[k][p]*Som1->at((int)(allpoints->at(n)[0])*col1+(int)(allpoints->at(n)[1]))[p];
                    mag1 += pow(MidPoints[k][p],2);
                    mag2 += pow(Som1->at((int)(allpoints->at(n)[0])*col1+(int)(allpoints->at(n)[1]))[p],2);
                }

                if(mag1==0 || mag2==0) {
                	similarity = -1;
                } else {
                	similarity /= (sqrt(mag1) * sqrt(mag2));
                }
                float dist1 = 1-similarity;
                DistMatrixRow.push_back(dist1);

            }*/
		float dist = 0;
                if(DistanceMetric == "Euclid") {
                        for(int p = 0; p < P; p++) {
                                dist += pow(Som1->at((int)(allpoints->at(n)[0])*col1+(int)(allpoints->at(n)[1]))[p]-MidPoints[k][p],2);
                        }
                        dist = dist/Som1->at(0).size();
                        dist = sqrt(dist);
			
                } else if(DistanceMetric == "Pearson") {
                        double Ex = 0;
                        double Ey = 0;
                        for(int p = 0; p < P; p++) {
                                Ex += Som1->at((int)(allpoints->at(n)[0])*col1+(int)(allpoints->at(n)[1]))[p];
                                Ey += MidPoints[k][p];
                        }
                        Ex /= P;
                        Ey /= P;
                        double stdDevX=0;
                        double stdDevY=0;
                        double cov = 0;
                        for(int p = 0; p < P; p++) {
                                stdDevX += pow(Som1->at((int)(allpoints->at(n)[0])*col1+(int)(allpoints->at(n)[1]))[p]-Ex,2);
                                stdDevY += pow(MidPoints[k][p]-Ey,2);
                                cov += (Som1->at((int)(allpoints->at(n)[0])*col1+(int)(allpoints->at(n)[1]))[p]-Ex)*(MidPoints[k][p]-Ey);
                        }
                        stdDevX = sqrt(stdDevX/P);
                        stdDevY = sqrt(stdDevY/P);
                        cov = cov/P;

                        dist = 1-(cov/(stdDevX*stdDevY));

                }
		DistMatrixRow.push_back(dist);

        }
	DistMatrix.push_back(DistMatrixRow);
    }
    for(int k = 0; k < row1; k++) {
    	for(int n = 0; n < col1; n++) {
        	pointlist[k][n]=-1;
        }
    }
	//cout<<"Distance calculated. "<<trial<<endl;
    for(int k = 0; k < K; k++) {
    	vector<vector<float> > temp;
        vector<vector<float> > temp2;
        vector<float> temp3;
        points.push_back(temp);
        points2.push_back(temp2);
        points3.push_back(temp3);
        Nc[k]=0;
    }
    for(int k = 0; k < K; k++) {
    	float lowest = 99999999;
        float lowestindex = -1;
        for(int n = 0; n < N; n++) {
        	if(DistMatrix[(int)allpoints->at(n)[0]*col1+(int)allpoints->at(n)[1]][k]<lowest || lowestindex==-1) {
            	lowest = DistMatrix[(int)allpoints->at(n)[0]*col1+(int)allpoints->at(n)[1]][k];
                lowestindex=n;
            }
      	}
        points[k].push_back(Som1->at((int)(allpoints->at(lowestindex)[0])*col1+(int)(allpoints->at(lowestindex)[1])));
        points2[k].push_back(allpoints2.at(lowestindex));
        points3[k]=allpoints2.at(lowestindex);
        pointlist[(int)allpoints->at(lowestindex)[0]][(int)allpoints->at(lowestindex)[1]]=k;
    }
	finished2=true;
	for(int k = 0; k < clusternum; k++) {
        int count = 0;
        for(int m = 0; m < row1; m++) {
			for(int n = 0; n < col1; n++) {
				if(pointlist[m][n]==k) {
                    count++;
                }
            }
		}
        if(count==0) {
			int index = rand() % (allpoints->size());
            for(int n = 0; n < Som1->at(0).size(); n++) {
				MidPoints[k][n]=Som1->at((int)(allpoints->at(index)[0])*col1+(int)(allpoints->at(index)[1]))[n];
            }
            allpoints2.at(index)[2]=k;
            finished2=false;
        }
	}
	}
	//build adjecency list
	//cout<<"Build Adjecency List: "<<trial<<endl;
    bool done = false;
    bool recalc=true;
    vector<vector<vector<float> > > adj;
    for(float r = .01; !done; r+=.01) {
    	if(recalc) {
        	recalc=false;
            adj.clear();
            for(int k = 0; k < K; k++) {
            	vector<vector<float> > temp;
                adj.push_back(temp);
            }
            for(int k = 0; k < K; k++) {
            	for(int n = 0; n < points[k].size(); n++) {
                	vector<vector<float> > possadj;
                    if((int)points2[k][n][0]%2==0) {
                    	vector<float> point1;
                        point1.push_back(points2[k][n][0]-1);
                        point1.push_back(points2[k][n][1]-1);
                        possadj.push_back(point1);
                        vector<float> point2;
                        point2.push_back(points2[k][n][0]-1);
                        point2.push_back(points2[k][n][1]);
                        possadj.push_back(point2);
                        vector<float> point3;
                        point3.push_back(points2[k][n][0]);
                        point3.push_back(points2[k][n][1]-1);
                        possadj.push_back(point3);
                        vector<float> point4;
                        point4.push_back(points2[k][n][0]);
                        point4.push_back(points2[k][n][1]+1);
                        possadj.push_back(point4);
                        vector<float> point5;
                        point5.push_back(points2[k][n][0]+1);
                        point5.push_back(points2[k][n][1]-1);
                        possadj.push_back(point5);
                        vector<float> point6;
                        point6.push_back(points2[k][n][0]+1);
                        point6.push_back(points2[k][n][1]);
                        possadj.push_back(point6);
                    } else {
						vector<float> point1;
                        point1.push_back(points2[k][n][0]-1);
                        point1.push_back(points2[k][n][1]);
                        possadj.push_back(point1);
                        vector<float> point2;
                        point2.push_back(points2[k][n][0]-1);
                        point2.push_back(points2[k][n][1]+1);
                        possadj.push_back(point2);
                        vector<float> point3;
                        point3.push_back(points2[k][n][0]);
                        point3.push_back(points2[k][n][1]-1);
                        possadj.push_back(point3);
                        vector<float> point4;
                        point4.push_back(points2[k][n][0]);
                        point4.push_back(points2[k][n][1]+1);
                        possadj.push_back(point4);
                        vector<float> point5;
                        point5.push_back(points2[k][n][0]+1);
                        point5.push_back(points2[k][n][1]);
                        possadj.push_back(point5);
                        vector<float> point6;
                        point6.push_back(points2[k][n][0]+1);
                        point6.push_back(points2[k][n][1]+1);
                        possadj.push_back(point6);
                    }
                    for(int m = 0; m < possadj.size(); m++) {
                    	if(possadj[m][0]>=row1) possadj[m][0]-=row1;
                        if(possadj[m][0]<0) possadj[m][0]+=row1;
                        if(possadj[m][1]>=col1) possadj[m][1]-=col1;
                        if(possadj[m][1]<0) possadj[m][1]+=col1;
                    }
                    for(int m = 0; m < possadj.size(); m++) {
                    	bool found = false;
                        if(pointlist[(int)possadj[m][0]][(int)possadj[m][1]]!=-1) {
                        	found = true;
                        }
                        for(int p = 0; p < adj[k].size() && !found; p++) {
                        	if(adj[k][p][0]==possadj[m][0]&&adj[k][p][1]==possadj[m][1]) {
                            	found = true;
                                break;
                            }
                        }
                        if(!found) {
                        	adj[k].push_back(possadj[m]);
                        }
                    }
                }
            }
        }
		for(int k = 0; k < K; k++) {
        	for(int n = 0; n < adj[k].size(); n++) {
            	if(DistMatrix[(int)adj[k][n][0]*col1+(int)adj[k][n][1]][k]<=r) {
                	if(pointlist[(int)adj[k][n][0]][(int)adj[k][n][1]]==-1) {
                    	points[k].push_back(Som1->at((int)(adj[k][n][0])*col1+(int)(adj[k][n][1])));
                        points2[k].push_back(adj[k][n]);
                        pointlist[(int)adj[k][n][0]][(int)adj[k][n][1]]=k;
                        recalc=true;
                        r=.01;
                    }
                }
            }
        }
        done = true;
        for(int k = 0; k < row1; k++) {
        	for(int n = 0; n < col1; n++) {
            	if(pointlist[k][n]==-1) {
                	done = false;
                }
            }
        }
    }
    for(int n = 0; n < N; n++) {
    	clusterIndex->push_back(pointlist[(int)allpoints->at(n)[0]][(int)allpoints->at(n)[1]]);
    }
	for(int i = 0; i < points2.size(); i++) {
                for(int j = 0; j < points2[i].size(); j++) {
                        points2[i][j].clear();
                }
                points2[i].clear();
        }
        points2.clear();
        for(int i = 0; i < points3.size(); i++) {
                points3[i].clear();
        }
        points3.clear();
        for(int i = 0; i < DistMatrix.size(); i++) {
                DistMatrix[i].clear();
        }
        DistMatrix.clear();
	for(int i = 0; i < MidPoints.size();i++) {
		MidPoints[i].clear();
	}
	MidPoints.clear();
	Rads.clear();
	//cout<<"Calculate AIC: "<<trial<<endl;
	// Calculate AIC
    // Nc contains number of objects in each cluster
    for(int k = 0; k < K; k++) {
    	Nc[k]=points[k].size();
    }
    //Calculate Midpoint of points
    //Vc is a P x K matrix that contains variances by cluster
	float** Vc = new float*[P];
	for(int i = 0; i < P; i++)
		Vc[i] = new float[K];
    	if(!Vc) {
		cout<<"Vc failed to allocate."<<endl;
    	}
    	for(int k = 0; k < K; k++) {
		float* Mids = new float[P];
		if(!Mids)
			cout<<"Mids failed to allocate: "<<k<<endl;
        	for(int p = 0; p < P; p++) {
        		Mids[p]=0;
			Vc[p][k] = 0;
        	}
        	for(int p = 0; p < P; p++) {
        		for(int n = 0; n < points[k].size(); n++) {
            			Mids[p]+=points[k][n][p];
            		}
            		if(points[k].size()>0)
            			Mids[p]/=points[k].size();
			for(int n = 0; n < points[k].size(); n++) {
                		Vc[p][k]+=pow(points[k][n][p]-Mids[p],2);
            		}
        	}
		delete [] Mids;
    	}
    	// Mid is the mid point of the whole SOM
    	float* Mid = new float[P];
    	if(!Mid)
		cout<<"Mid failed to allocate."<<endl;

    	for(int p = 0; p < P; p++) {
    		Mid[p]=0;
    	}
    	for(int n = 0; n < N; n++) {
    		for(int p = 0; p < P; p++) {
        		Mid[p]+=Som1->at(allpoints->at(n)[0]*col1+allpoints->at(n)[1])[p];
        	}
    	}
	if(N>0)
    		for(int p = 0; p < P; p++) {
        		Mid[p]=Mid[p]/N;
        	}

    	// V is a P x 1 matrix that contains variances for whole sample
    	float* V = new float[P];
	if(!V)
		cout<<"V failed to allocate."<<endl;
    	for(int p = 0; p < P; p++) {
        	float Var = 0;
        	for(int n = 0; n < N; n++) {
        		Var += pow(Som1->at(allpoints->at(n)[0]*col1+allpoints->at(n)[1])[p]-Mid[p],2);
        	}
		if(N>1)
	        	V[p]=Var/(N-1);
		else
			V[p]=Var;
    }
    //Compute log-like LL, 1 x K
    float* LL = new float[K];
    for(int k = 0; k < K; k++) {
    	float csum = 0;
        for(int p = 0; p < P; p++) {
        	csum+=log(Vc[p][k]+V[p])/2.0;
        }
        LL[k]=-1*Nc[k] * csum;
    }
    //Compute AIC and BIC
    float rsum = 0;
    for(int k = 0; k < K; k++) {
    	rsum+=LL[k];
    }
    float AIC;
    if(Dimensionality==-1)
    	AIC = -2*rsum+4*K*P;
    else {
    	AIC=-2*rsum+4*K*Dimensionality;
    }
	*score = AIC;
}
/*class distorder {
	int pos;
	float dist;
	bool operator<(distorder d1, distorder d2) {
		return d1.dist<d2.dist;
	}
}*/

int main(int argc, char* argv[]) {
	if(argc < 2) {
        cout << "Usage: ./metasom [options] -SOMFile <SOM File Location> -Outfile <Output File Location>" <<endl;
        cout << "Options: <default>" <<endl;
        cout << "-Rows: Number of rows for your SOM. <20>"<<endl;
        cout << "-Cols: Number of columns for your SOM. <30>"<<endl;
		cout << "-Metaclusters: Number of metaclusters to start with. <20>"<<endl;
		cout << "-MetaclustersEnd: Number of metaclusters to end with. <20>"<<endl;
        cout << "-Trials: Number of trials for your metaclustering. The best clustering will be chosen.  <20>"<<endl;
		cout << "-Dimensionality: Effective number of dimensions in your sample set. Defaults to the total number of samples."<<endl;
        return 0;
    }

    srand ( unsigned ( std::time(0) ) );
    int row1=20;
    int col1=30;
	int kmeans1 = 20;
	int kmeans2 = 20;
	int numberOfTrials=20;
	string ATACSom;
	string BestClusterName;
	string DistanceMetric="Euclid";
	string geneFilePrefix="";
	int dimensionality=-1;
	for(int i = 0; i < argc; i++) {
        string temp = argv[i];
        if(temp.compare("-Rows")==0)
            istringstream(argv[i+1])>>row1;
		if(temp.compare("-Cols")==0)
            istringstream(argv[i+1])>>col1;
		if(temp.compare("-Metaclusters")==0)
            istringstream(argv[i+1])>>kmeans1;
		if(temp.compare("-MetaclustersEnd")==0)
            istringstream(argv[i+1])>>kmeans2;
		if(temp.compare("-Trials")==0)
            istringstream(argv[i+1])>>numberOfTrials;
		if(temp.compare("-SOMFile")==0)
            ATACSom = argv[i+1];
		if(temp.compare("-Outfile")==0)
            BestClusterName = argv[i+1];
		if(temp.compare("-genePrefix")==0) 
			geneFilePrefix=argv[i+1];
		if(temp.compare("-Dimensionality")==0) 
			istringstream(argv[i+1])>>dimensionality;
		if(temp.compare("-DistanceMetric")==0)
			DistanceMetric=argv[i+1];
	}
	cout<<"Opening Gene Files"<<endl;
	vector<vector<int> > geneCounts;
	for(int i = 0; i < row1; i++) {
		vector<int> temp;
		for(int j = 0; j < col1; j++) {
			int count = 0;
			cout<<geneFilePrefix+"_"+SSTR(i)+"_"+SSTR(j)+".unit"<<endl;
			ifstream geneFile((geneFilePrefix+"_"+SSTR(i)+"_"+SSTR(j)+".unit").c_str());
			string line;
			while(getline(geneFile,line)) {
				count++;
			}
			temp.push_back(count);
			cout<<count<<endl;
		}
		geneCounts.push_back(temp);
	}
	cout<<"Opening SOM file"<<endl;	
    ifstream som1file(ATACSom.c_str());
    vector<vector<float> > Som1;
    string line;
	vector<vector<float> > allpoints;
	int tempcol = 0;
        int temprow = 0;
        while(getline(som1file,line)) {
                if(line[0]=='#') continue;
                cout<<tempcol<<'\t'<<temprow<<endl;
                vector<float> Som1Row;
                vector<string> splitz = split(line, '\t');
                vector<float> temp;
                temp.push_back(temprow);
                temp.push_back(tempcol);
                temp.push_back(0);
                temp.push_back(0);
                temp.push_back(geneCounts[temprow][tempcol]+1);
                allpoints.push_back(temp);
                for(int i = 0; i < splitz.size(); i++) {
                        float temp;
                        istringstream(splitz[i])>>temp;
                        Som1Row.push_back(temp);
                }
                Som1.push_back(Som1Row);
                tempcol++;
                if(tempcol == col1) {
                        tempcol=0;
                        temprow++;
                }
    //    if(Som1Row.size()==0) cout<<line<<endl;
    	}
    //
	
	vector<vector<int> > Indexes;
	vector<int> bestClusters;
	vector<float> bestScores;
	for(int i = kmeans1; i <= kmeans2; i++) {
		vector<thread> threads;
		cout<<"Cluster Num: "<<i<<endl;
		vector<vector<int> > IndexRow;
		vector<float> ScoreRow;
		// Run numberOfTrials
		for(int j = 0; j < numberOfTrials; j++) {
			ScoreRow.push_back(0.0);
			vector<int> temp;
			IndexRow.push_back(temp);
		}
		int count = 0;
		cout<<"Pushing threads"<<endl;
		for(int j = 0; j < numberOfTrials; j++) {
			threads.push_back(thread(ClusterAndScore,&allpoints, &Som1, dimensionality, i, DistanceMetric, &(IndexRow[j]), &(ScoreRow[j]), row1, col1,j));
		}
		for (auto& th : threads) th.join();
		float minScore = -999;
                int bestTrial = -1;
                for(int j = 0; j < ScoreRow.size(); j++) {
                        if(ScoreRow[j]<minScore || bestTrial == -1) {
                                minScore = ScoreRow[j];
                                bestTrial = j;
                        }
                }
		cout<<minScore<<endl;
                bestScores.push_back(minScore);
		Indexes.push_back(IndexRow[bestTrial]);
	}

	/*for(int i = 0; i < Scores.size(); i++) {
		float minScore = -999;
		int bestTrial = -1;
		for(int j = 0; j < Scores[i].size(); j++) {
			if(Scores[i][j]<minScore || bestTrial == -1) {
				minScore = Scores[i][j];
				bestTrial = j;
			}
		}
		bestClusters.push_back(bestTrial);
		bestScores.push_back(minScore);
	}*/
	float minScore = -999;
	int bestClusterNum = -1;
	for(int i = 0; i < bestScores.size(); i++) {
		if(bestScores[i] < minScore || bestClusterNum==-1) {
			minScore = bestScores[i];
			bestClusterNum = i;
		}
	}
	
	ofstream outfile2(BestClusterName.c_str());
	outfile2<<"# Cluster Number: "<<bestClusterNum+kmeans1<<endl;
	vector<int> BestIndex = Indexes[bestClusterNum];
	for(int j = 0; j < BestIndex.size(); j++) {
		outfile2<<(int)j/col1<<'\t'<<j%col1<<'\t'<<BestIndex[j]<<endl;
	}		
	outfile2.close();
}	
