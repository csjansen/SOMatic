#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <vector>
#include <queue>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string>
#include <locale>
#include <ctime>
using namespace std;

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()
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

/*class distorder {
	int pos;
	double dist;
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
		cout << "-Sparse: Show this option to use Cos Distance"<<endl;
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
	bool Sparse = false;
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
		if(temp.compare("-Sparse")==0)
			Sparse = true;
		if(temp.compare("-genePrefix")==0) 
			geneFilePrefix=argv[i+1];
		if(temp.compare("-Dimensionality")==0) 
			istringstream(argv[i+1])>>dimensionality;
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
		}
		geneCounts.push_back(temp);
	}
	cout<<"Opening SOM file"<<endl;	
    ifstream som1file(ATACSom.c_str());
    vector<vector<double> > Som1;
    string line;
	vector<vector<double> > allpoints;
	while(getline(som1file,line)) {
        if(line[0]=='#') continue;
        vector<double> Som1Row;
        vector<string> splitz = split(line, '\t');
		vector<double> temp;
		int temprow;
		int tempcol;
		istringstream(splitz[0])>>temprow;
		istringstream(splitz[1])>>tempcol;
		temp.push_back(temprow);
		temp.push_back(tempcol);
		temp.push_back(0);
		temp.push_back(0);
		temp.push_back(geneCounts[temprow][tempcol]+1);
		allpoints.push_back(temp);
        for(int i = 2; i < splitz.size(); i++) {
            double temp;
            istringstream(splitz[i])>>temp;
            Som1Row.push_back(temp);
        }
        Som1.push_back(Som1Row);
    //    if(Som1Row.size()==0) cout<<line<<endl;
    }
	
	vector<double> finalAverageRads;
	vector<double> RSSs;
	vector<double> RSSsb;
	vector<int> RSSindexs;
	vector<vector<double> > AllPointsStorage;
	
	vector<double> averageRads;
	int BestKmeans = 0;
	double BestRad=9999999999;
	int BestRadIndex=-1;
	for(int i = kmeans1; i <= kmeans2; i++) {
		cout<<"Cluster Num: "<<i<<endl;
	// Run numberOfTrials
	for(int j = 0; j < numberOfTrials; j++) {
		double RSS=0;
		cout<<"Trial: "<<j<<endl;
		vector<double> Rads;
    	random_shuffle(allpoints.begin(),allpoints.end());
		cout<<"Midpoints created"<<endl;
		vector<vector<double> > MidPoints;
		// Initial Mid Points
		int uniqs = 0;	
		for(int k = 0; k < i+uniqs; k++) {
			vector<double> temp;
			int sizer = allpoints.size();
			for(int m = 0; m < Som1[0].size(); m++) {
				temp.push_back(Som1[(int)(allpoints[k%sizer][0])*col1+(int)(allpoints[k%sizer][1])][m]);
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
			}
			else
				uniqs++;
		}
			
		cout<<"Clustering"<<endl;
		bool finished = false;
		int iter = 0;
		bool cont = false;
		while(!finished) {
			finished = true;
		    int count=0;
        	//Assignment
		    for(int k = 0; k < allpoints.size(); k++) {
		        double lowestdist1=999999999;
				double lowestdist=99999999;
        		double lowestdist2=999999999;
		        int lowestk1=-1;
        		int lowestk2=-1;
		        for(int m = 0; m <i; m++) {
					double similarity=0;
		            double mag1 = 0;
        		    double mag2 = 0;
					double dist = 0;
                	for(int num = 0; num < Som1[0].size(); num++) {
						dist += pow(MidPoints[m][num]-Som1[(int)(allpoints[k][0])*col1+(int)(allpoints[k][1])][num],2);
                    	similarity += MidPoints[m][num]*Som1[(int)(allpoints[k][0])*col1+(int)(allpoints[k][1])][num];
                    	mag1 += pow(MidPoints[m][num],2);
		                mag2 += pow(Som1[(int)(allpoints[k][0])*col1+(int)(allpoints[k][1])][num],2);
                	}
					dist = sqrt(dist);
						
					if(mag1==0 || mag2==0) {
                    	similarity = -1;
                	} else {
                    	similarity /= (sqrt(mag1) * sqrt(mag2));
                	}
                	double dist1; 
					if(Sparse)
						dist1 = 1-similarity;
					else
						dist1 = dist;
                	if(lowestdist1>dist1) {
                    	lowestdist1 = dist1;
                    	lowestk1 = m;
						lowestdist=dist1;
                	}
            	}
				if(allpoints[k][2] != lowestk1) {
					count++;
        		    allpoints[k][2] = lowestk1;
                	allpoints[k][3] = lowestdist1;
		            finished = false;
        		}
				if(lowestdist > Rads[lowestk1]) Rads[lowestk1]=lowestdist;
			}
	//		cout<<count<<endl;
			if(finished) break;
			for(int k = 0; k < Rads.size(); k++) {
				Rads[k]=0;
			}
			cont=true;
        	//Update
        	for(int k = 0; k < i; k++) {
            	vector<double> totals1;
			    int count = 0;
            	for(int m = 0; m < Som1[0].size(); m++) totals1.push_back(0);
			    for(int m = 0; m < allpoints.size(); m++) {
            		if((int)(allpoints[m][2])==k) {
			        	for(int n = 0; n < Som1[0].size(); n++) totals1[n]+=Som1[(int)(allpoints[m][0])*col1+(int)(allpoints[m][1])][n];
            			count++;
                	}
				}
			    for(int m = 0; m < totals1.size(); m++) {
					MidPoints[k][m] = totals1[m]/(double)count;
				//		cout<<MidPoints[k][m]<<endl;
				}
        	}
		}
		int K = i;
        int N = allpoints.size();
        int P = Som1[0].size();
		cout<<K<<endl;
		cout<<N<<endl;
		cout<<P<<endl;
		//rebuild clusters to be adjacient
		//Calculating Distance Matrix
		cout<<"Calculating Distance Matrix..."<<endl;
		double DistMatrix[N][K];
		cout<<"Allocated"<<endl;
		for(int n = 0; n < N; n++) {
			for(int k = 0; k < K; k++) {
				DistMatrix[(int)allpoints[n][0]*col1+(int)allpoints[n][1]][k]=0;
				if(!Sparse) {
					for(int p = 0; p < P; p++) {
						DistMatrix[(int)allpoints[n][0]*col1+(int)allpoints[n][1]][k]+=pow(Som1[(int)(allpoints[n][0])*col1+(int)(allpoints[n][1])][p]-MidPoints[k][p],2);
					}
					DistMatrix[(int)allpoints[n][0]*col1+(int)allpoints[n][1]][k]=sqrt(DistMatrix[(int)allpoints[n][0]*col1+(int)allpoints[n][1]][k]);
				} else {
					double similarity=0;
                    double mag1 = 0;
                    double mag2 = 0;
                    for(int p = 0; p < P; p++) {
                        similarity += MidPoints[k][p]*Som1[(int)(allpoints[n][0])*col1+(int)(allpoints[n][1])][p];
                        mag1 += pow(MidPoints[k][p],2);
                        mag2 += pow(Som1[(int)(allpoints[n][0])*col1+(int)(allpoints[n][1])][p],2);
                    }

                    if(mag1==0 || mag2==0) {
                        similarity = -1;
                    } else {
                        similarity /= (sqrt(mag1) * sqrt(mag2));
                    }
                    double dist1 = 1-similarity;
                    DistMatrix[(int)allpoints[n][0]*col1+(int)allpoints[n][1]][k]=dist1;

				}
			}
		}
		cout<<"Assigning Closest Points"<<endl;
		//assign closest points
		int pointlist[row1][col1];
		for(int k = 0; k < row1; k++) {
        	for(int n = 0; n < col1; n++) {
            	pointlist[k][n]=-1;
            }
        }
		vector<vector<vector<double> > > points;
		vector<vector<vector<double> > > points2;
		vector<vector<double> > points3;
		int Nc[allpoints.size()];
        for(int k = 0; k < K; k++) {
        	vector<vector<double> > temp;
            vector<vector<double> > temp2;
            vector<double> temp3;
            points.push_back(temp);
            points2.push_back(temp2);
            points3.push_back(temp3);
            Nc[k]=0;
        }
		for(int k = 0; k < K; k++) {
			double lowest = 99999999;
			//double seclowest = 99999999;
			double lowestindex = -1;
			//double seclowestind = -1;
			for(int n = 0; n < N; n++) {
				//cout<<N<<'\t'<<n<<'\t'<<k<<endl;
				//cout<<(int)allpoints[n][0]*col1+(int)allpoints[n][1];
				if(DistMatrix[(int)allpoints[n][0]*col1+(int)allpoints[n][1]][k]<lowest || lowestindex==-1) {
			//		seclowest=lowest;
					lowest = DistMatrix[(int)allpoints[n][0]*col1+(int)allpoints[n][1]][k];
			//		seclowestind=lowestindex;
					lowestindex=n;
			//	} else if(DistMatrix[(int)allpoints[n][0]*col1+(int)allpoints[n][1]][k]<seclowest) {
			//		seclowest = DistMatrix[n][k];
			//		seclowestind=n;
				}
			}
			//cout<<"after"<<endl;
			//cout<<lowestindex<<endl;
			//cout<<allpoints.size()<<endl;
			//cout<<(int)(allpoints[lowestindex][0])*col1+(int)(allpoints[lowestindex][1])<<endl;
			points[k].push_back(Som1[(int)(allpoints[lowestindex][0])*col1+(int)(allpoints[lowestindex][1])]);
			points2[k].push_back(allpoints[lowestindex]);
			points3[k]=allpoints[lowestindex];
			pointlist[(int)allpoints[lowestindex][0]][(int)allpoints[lowestindex][1]]=k;
		}
		//build adjecency list
		cout<<"Building Adj List"<<endl;
		bool done = false;
		bool recalc=true;
		vector<vector<vector<double> > > adj;
		for(double r = .01; !done; r+=.01) {
			if(recalc) {
				recalc=false;
				adj.clear();
				for(int k = 0; k < K; k++) {
					vector<vector<double> > temp;
					adj.push_back(temp);
				}
				for(int k = 0; k < K; k++) {
					for(int n = 0; n < points[k].size(); n++) {
						vector<vector<double> > possadj;
						if((int)points2[k][n][0]%2==0) {
							vector<double> point1;
							point1.push_back(points2[k][n][0]-1);
							point1.push_back(points2[k][n][1]-1);
							possadj.push_back(point1);
							vector<double> point2;
						    point2.push_back(points2[k][n][0]-1);
							point2.push_back(points2[k][n][1]);
			                possadj.push_back(point2);
							vector<double> point3;
							point3.push_back(points2[k][n][0]);
			                point3.push_back(points2[k][n][1]-1);
						    possadj.push_back(point3);
							vector<double> point4;
			                point4.push_back(points2[k][n][0]);
						    point4.push_back(points2[k][n][1]+1);
							possadj.push_back(point4);
							vector<double> point5;
						    point5.push_back(points2[k][n][0]+1);
							point5.push_back(points2[k][n][1]-1);
			                possadj.push_back(point5);
							vector<double> point6;
							point6.push_back(points2[k][n][0]+1);
			                point6.push_back(points2[k][n][1]);
						    possadj.push_back(point6);
						} else {
							vector<double> point1;
			                point1.push_back(points2[k][n][0]-1);
						    point1.push_back(points2[k][n][1]);
							possadj.push_back(point1);
			                vector<double> point2;
						    point2.push_back(points2[k][n][0]-1);
							point2.push_back(points2[k][n][1]+1);
			                possadj.push_back(point2);
						    vector<double> point3;
							point3.push_back(points2[k][n][0]);
			                point3.push_back(points2[k][n][1]-1);
						    possadj.push_back(point3);
							vector<double> point4;
			                point4.push_back(points2[k][n][0]);
						    point4.push_back(points2[k][n][1]+1);
							possadj.push_back(point4);
			                vector<double> point5;
						    point5.push_back(points2[k][n][0]+1);
							point5.push_back(points2[k][n][1]);
			                possadj.push_back(point5);
						    vector<double> point6;
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
							points[k].push_back(Som1[(int)(adj[k][n][0])*col1+(int)(adj[k][n][1])]);
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
			allpoints[n][2]=pointlist[(int)allpoints[n][0]][(int)allpoints[n][1]];
		}
		cout<<"Calculating AIC"<<endl;
		//Calculate AIC
	    //Nc contains number of objects in each cluster
    	for(int k = 0; k < K; k++) {
        	Nc[k]=points[k].size();
	    }
		//for(int k = 0; k < K; k++) {
		//	cout<<k<<'\t'<<Nc[k]<<endl;
		//}
		//Calculate Midpoint of points
		vector<vector<double> > Mids;
		for(int k = 0; k < K; k++) {
			vector<double> temp;
			for(int p = 0; p < P; p++) {
				temp.push_back(0);
			}
			for(int p = 0; p < P; p++) {
				for(int n = 0; n < points[k].size(); n++) {
					temp[p]+=points[k][n][p];
				}
				if(points[k].size()>0)
					temp[p]/=points[k].size();
			}
			Mids.push_back(temp);
		}
		//cout<<"MidPoint"<<endl;
		//for(int p = 0; p < P; p++) {
		//	cout<<Mids[0][p]<<endl;
		//}
    	//Vc is a P x K matrix that contains variances by cluster
	    double Vc[P][K];
    	for(int p = 0; p < P; p++) {
        	for(int k = 0; k < K; k++) {
                double Var=0;
               	for(int n = 0; n < points[k].size(); n++) {
                   	Var+=pow(points[k][n][p]-Mids[k][p],2);
                }
				if(points[k].size() >0)
					Var = Var/(points[k].size());
				Vc[p][k]=Var;
			}
		}
		//cout<<"Variances in first cluster"<<endl;
		//for(int p = 0; p < P; p++) {
        //    cout<<Vc[p][0]<<endl;
        //}
		// Mid is the mid point of the whole SOM
		double Mid[P];
		for(int p = 0; p < P; p++) {
			Mid[p]=0;
		}
		for(int n = 0; n < N; n++) {
			for(int p = 0; p < P; p++) {
				Mid[p]+=Som1[allpoints[n][0]*col1+allpoints[n][1]][p];
			}
		}
		//cout<<"Midpoint of SOM"<<endl;
		if(N>0)
			for(int p = 0; p < P; p++) {
				Mid[p]=Mid[p]/N;
		//	cout<<Mid[p]<<endl;
			}
		//V is a P x 1 matrix that contains variances for whole sample
		double V[P];
		for(int p = 0; p < P; p++) {
			double Var = 0;
			for(int n = 0; n < N; n++) {
				Var += pow(Som1[allpoints[n][0]*col1+allpoints[n][1]][p]-Mid[p],2);
			}
			V[p]=Var/(N-1);
		}
		//cout<<"Variances in whole SOM"<<endl;
        //for(int p = 0; p < P; p++) {
        //    cout<<V[p]<<endl;
        //}
		//Compute log-like LL, 1 x K 
		//cout<<"LL"<<endl;
		double LL[K];
		for(int k = 0; k < K; k++) {
			double csum = 0;
			for(int p = 0; p < P; p++) {
				csum+=log(Vc[p][k]+V[p])/2.0;
			}
			LL[k]=-1*Nc[k] * csum;
		//	cout<<LL[k]<<endl;
		}
		//Compute AIC and BIC
		double rsum = 0;
		for(int k = 0; k < K; k++) {
			rsum+=LL[k];
		}
		//cout<<rsum<<endl;
		double AIC;
		if(dimensionality==-1) 
			AIC = -2*rsum+4*K*P;
		else { 
			AIC=-2*rsum+4*K*dimensionality;
			cout<<"Dimensionality: "<<dimensionality<<endl;
		}
		double BIC = -2*rsum+2*K*P*log(N);
		cout<<AIC<<endl;
		if(BestRad>AIC) {
			BestKmeans=i;
			BestRad=AIC;
			BestRadIndex=j;
			AllPointsStorage=allpoints;
		}
	}
	}
	cout<<"Cluster Number: "<<BestKmeans<<endl;
	ofstream outfile2(BestClusterName.c_str());
	outfile2<<"# Cluster Number: "<<BestKmeans<<endl;
	vector<vector<double> > uniq;
	for(int j = 0; j < AllPointsStorage.size(); j++) {
		bool found = false;
		for(int i = 0; i < uniq.size(); i++) {
			if(uniq[i][0]==AllPointsStorage[j][0]&&uniq[i][1]==AllPointsStorage[j][1]) {
				found = true;
				break;
			}
		}
		if(!found) {
			outfile2<<AllPointsStorage[j][0]<<'\t'<<AllPointsStorage[j][1]<<'\t'<<AllPointsStorage[j][2]<<endl;
			uniq.push_back(AllPointsStorage[j]);
		}
	}		
	outfile2.close();
}	
