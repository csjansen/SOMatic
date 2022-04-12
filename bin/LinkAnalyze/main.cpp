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

class treeNode {
public:
    string GO_ID;
    string GO_Desc;
    double pval;
    vector<treeNode*> parents;
    vector<treeNode*> children;
};


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
}


void treeadd(vector<treeNode*>* Tree, string Id, double pval, map<string, vector<string> > GoHeir) {
	// If the node already exists, update pval and return
	bool found = false;
	for(int i = 0; i < Tree->size(); i++) {
		if(Tree->at(i)->GO_ID.compare(Id)==0) {
			Tree->at(i)->pval = pval;
			found = true;
			break;
		}
	}
	if(found) return;

	//Otherwise, add as child node and build into tree
	treeNode* newNode = new treeNode;
	newNode->GO_ID = Id;
	newNode->pval = pval;
	
	queue<treeNode*> addqueue;
	addqueue.push(newNode);

	while(addqueue.size()>0) {
		treeNode* present = addqueue.front();
		addqueue.pop();
		Tree->push_back(present);
		vector<string> parents = GoHeir[present->GO_ID];
		for(int i = 0; i < parents.size(); i++) {
			bool found = false;
			for(int j = 0; j < Tree->size(); j++) {
				if(Tree->at(j)->GO_ID.compare(parents[i])==0) {
					Tree->at(j)->children.push_back(present);
					present->parents.push_back(Tree->at(j));
					found = true;
					break;
				}
			}
			if(!found) {
				treeNode* newNode = new treeNode;
				present->parents.push_back(newNode);
				newNode->children.push_back(present);
				newNode->GO_ID = parents[i];
				newNode->pval = 1;
				addqueue.push(newNode);
			}
		}
	}
	return;
}

void filltree(vector<treeNode*>* Tree, vector<string> Ids, vector<double> pvals, map<string, vector<string> > GoHeir) {
	cout<<Ids.size()<<endl;
	for(int i = 0; i < Ids.size(); i++) {
		treeadd(Tree, Ids[i], pvals[i], GoHeir);
	}
	cout<<Tree->size()<<endl;
}


#define MAXIT 10000
#define EPS 3.0e-7
#define FPMIN 1.0e-30


double betacf(double a, double b, double x)
// Used by betai: Evaluates continued fraction for incomplete beta function
// by modified Lentz s method (§5.2).
{
    int m,m2;
    double aa,c,d,del,h,qab,qam,qap;

    qab=a+b;
    // These q s will be used in factors that occur in the coefficients (6.4.6).
    qap=a+1.0;
    qam=a-1.0;
    c=1.0; // First step of Lentz s method.
    d=1.0-qab*x/qap;
    if (fabs(d) < FPMIN)
        d=FPMIN;
    d=1.0/d;
    h=d;
    for (m=1; m<=MAXIT; m++) {
        m2=2*m;
        aa=m*(b-m)*x/((qam+m2)*(a+m2));
        d=1.0+aa*d; // One step (the even one) of the recurrence.
        if (fabs(d) < FPMIN)
            d=FPMIN;
        c=1.0+aa/c;
        if (fabs(c) < FPMIN)
            c=FPMIN;
        d=1.0/d;
        h *= d*c;
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
        d=1.0+aa*d; // Next step of the recurrence (the odd one).
        if (fabs(d) < FPMIN)
            d=FPMIN;
        c=1.0+aa/c;
        if (fabs(c) < FPMIN)
            c=FPMIN;
        d=1.0/d;
        del=d*c;
        h *= del;
        if (fabs(del-1.0) < EPS)
            break;
    }
    //if (m > MAXIT)
        //errAbort("a or b too big, or MAXIT too small in betacf");
    return h;
}

double betai(double a, double b, double x)
// Returns the incomplete beta function Ix(a, b).
{
    double bt;

    //if (x < 0.0 || x > 1.0)
        //errAbort("Bad x in routine betai");
    if(x==0.0)
		return 1;
    if (x == 0.0 || x == 1.0)
        bt=0.0;
    else // Factors in front of the continued fraction.
        bt=exp(lgamma(a+b)-lgamma(a)-lgamma(b)+a*log(x)+b*log(1.0-x));
    if (x < (a+1.0)/(a+b+2.0)) // Use continued fraction directly.
        return bt*betacf(a,b,x)/a;
    else // Use continued fraction after making the symmetry transformation.
        return 1.0-bt*betacf(b,a,1.0-x)/b;
}


double getBinomPval(int n, int k, double p)
{
    if (k == 0) return 1;
    else return betai(k, n-k+1, p);
}



int hexdist(int row1, int col1, int row2, int col2, int rows, int cols) {
    row1=row1-row2;
    col1=col1-col2;
    row2=0;
    col2=0;
    int x1 = col1 - (row1 - (abs(row1)%2))/2;
    int z1 = row1;
    int y1 = -x1-z1;
    int x2 = col2 - (row2 - (abs(row2)%2))/2;
    int z2 = row2;
    int y2 = -x2-z2;
    return (abs(x1-x2) + abs(y1-y2) + abs(z1-z2))/2;
}


double lngamm(double z)
// Reference: "Lanczos, C. 'A precision approximation 
// of the gamma function', J. SIAM Numer. Anal., B, 1, 86-96, 1964."
// Translation of  Alan Miller's FORTRAN-implementation
// See http://lib.stat.cmu.edu/apstat/245
{
  double x = 0;
  x += 0.1659470187408462e-06/(z+7);
  x += 0.9934937113930748e-05/(z+6);
  x -= 0.1385710331296526    /(z+5);
  x += 12.50734324009056     /(z+4);
  x -= 176.6150291498386     /(z+3);
  x += 771.3234287757674     /(z+2);
  x -= 1259.139216722289     /(z+1);
  x += 676.5203681218835     /(z);
  x += 0.9999999999995183;
  return(log(x)-5.58106146679532777-z+(z-0.5)*log(z+6.5));
}

double fakt(int n) {
if(n<=1) return 0;
else return lngamm((double)n+1);
}

double lnbico(int n, int k) {
	return (fakt(n)-fakt(k)-fakt(n-k));
}

int sn11;
int sn10;
int sn01;
int sn;
double sprob;

double hyper0(int n11i, int n10i, int n01i, int ni) 
{
    if(n10i == 0 && n01i == 0 && ni==0){
		if(!(n11i % 10 == 0))
		{
			if(n11i==sn11+1)  
			{
				sprob *= ((sn10-sn11)/(double)(n11i))*((sn01-sn11)/(double)(n11i+sn-sn10-sn01));
				sn11 = n11i;
				return sprob;
			}
			if(n11i==sn11-1)
			{
				sprob *= ((sn11)/(double)(sn10-n11i))*((sn11+sn-sn10-sn01)/(double)(sn01-n11i));
				sn11 = n11i;
				return sprob;
			}
		}
		sn11 = n11i;	
	} else {
		sn11=n11i;
		sn10=n10i;
		sn01=n01i;
		sn=ni;
	}
	sprob = exp(lnbico(sn10,sn11)+lnbico(sn-sn10,sn01-sn11)-lnbico(sn,sn01));//hyper_323(sn11,sn10,sn01,sn);
	return sprob;
}

double hyper(int n11) {
	return hyper0(n11,0,0,0);
}

double exact(int n11, int n10, int n01, int n)
{
	double sleft;
	double sright;
	double sless;
	double slarg;
    double p;
	double prob;
	double max=n10;
	int i;
	int j;
	if(n01<max) max=n01;
	double min = n10+n01-n;
	if(min<0) min=0;
	if(min==max)
	{
		sless = 1;
	    sright= 1;
		sleft = 1;
	    slarg = 1;
		return 1;
	}
	prob=hyper0(n11, n10, n01, n);
	sleft=0;
	p=hyper(min);
	//cout<<n11<<'\t'<<n10<<'\t'<<n01<<'\t'<<n<<endl;
	//cout<<p<<'\t'<<min<<endl;
	//int temp;
	//cin>>temp;
	for(i=min+1; p<=0.99999999*prob; i++)
	{
		sleft += p;
		p=hyper(i);
	//	cout<<i<<'\t'<<p<<endl;
	}
	i--;
	if(p<=1.00000001*prob) sleft += p;
	else i--;
	sright=0;
	p=hyper(max);
	for(j=max-1; p<=0.99999999*prob; j--)
	{
		sright += p;
		p=hyper(j);
	}
	j++;
	if(p<=1.00000001*prob) sright += p;
	else j++;
	if(abs(i-n11)<abs(j-n11)) 
	{
	    sless = sleft;
		slarg = 1 - sleft + prob;
	} else {
	    sless = 1 - sright + prob;
		slarg = sright;
	}
	return slarg;
}


double Fishers(int GoCount, int TotalGoCount, int OverallGo, int TotalOverallGo) {
    int a = GoCount;
    int c = TotalGoCount;
    int b = OverallGo;
    int d = TotalOverallGo;
     //return(Math.exp(lnbico(n1_,n11)+lnbico(n-n1_,n_1-n11)-lnbico(n,n_1)));

    double pval=exact(a, a+b, a+c, a+b+c+d);
    //cout<<a<<'\t'<<b<<'\t'<<c<<'\t'<<d<<'\t'<<pval<<endl;
    //int temp;
    //cin>>temp;
    //cout<<pval<<endl;
	//cout<<endl;
    return pval;
}

struct bioms {
public:
	string GO;
	string GODesc;
	double biom;
	int GoCounts;
	int TotalGoCounts;
	double annot;
};
	bool operator<(bioms const &a, bioms const &b) {
		return a.biom<b.biom;
	}

double MWUtest(double n1, double n2, double R1, double R2) {
	double U1 = n1*n2+(n1*(n1+1))/2-R1;
	double U2 = n2*n2+(n1*(n2+1))/2-R2;
	double U = U1;
	if(U<U2) U=U2;
	return (U-(n1*n2)/2)/(double)sqrt(n1*n2*(n1+n2+1)/12);
}

int main(int argc, char* argv[]) {
	 if(argc < 2) {
        cout << "Usage: ./LinkAnalyze [options] -FusionFile <Output from Fusion> -FusionClusterFile <Output from FusionCluster> -ClusterNum1 <Cluster Number from SOM 1> -ClusterNum2 <Cluster Number from SOM 2> -OutputPrefix <Prefix for Output Files> -Row1 -Col1 -Row2 -Col2" <<endl;
        return 0;
    }
	string clusterfile;
	int kmeans1;
	int kmeans2;
	string outputprefix;
	string ConnectionsFileName;
	bool Fuse=false;
    for(int i = 0; i < argc; i++) {
        string temp = argv[i];
        if(temp.compare("-FusionClusterFile")==0)
            clusterfile=argv[i+1];
		if(temp.compare("-FusionFile")==0)
			ConnectionsFileName=argv[i+1];
        if(temp.compare("-ClusterNum1")==0)
            istringstream(argv[i+1])>>kmeans1;
        if(temp.compare("-ClusterNum2")==0)
            istringstream(argv[i+1])>>kmeans2;
        if(temp.compare("-OutputPrefix")==0)
            outputprefix=argv[i+1];
	if(temp.compare("-Fuse")==0)
	    Fuse=true;
	}
	ifstream ClusterFile(clusterfile.c_str());
	ifstream ConnectionsFile(ConnectionsFileName.c_str());
	vector<vector<vector<string> > > ConnectionGenes;
    string line;
	vector<vector<double> > allpoints;
    map<string,vector<string> > allgenes;
    map<string,vector<string> > allregions;
	map<string,vector<string> > regionToGene;
    cout<<"Opening file"<<endl;
    map<string,int> pointcounts;
    while(getline(ConnectionsFile, line)) {
        vector<string> splitz = split(line,'\t');
        int AtacRow;
        istringstream(splitz[0])>>AtacRow;
        int AtacCol;
        istringstream(splitz[1])>>AtacCol;
        int i = 2;
        while(i < splitz.size()-2) {
            int RNARow;
            int RNACol;
            int AtacCount;
            istringstream(splitz[i])>>RNARow;
            istringstream(splitz[i+1])>>RNACol;
            istringstream(splitz[i+2])>>AtacCount;
            for(int j = 0; j < AtacCount; j++) {
                allgenes[SSTR(AtacRow)+"\t"+SSTR(AtacCol)+"\t"+SSTR(RNARow)+"\t"+SSTR(RNACol)].push_back(splitz[i+j*2+3]);
                allregions[SSTR(AtacRow)+"\t"+SSTR(AtacCol)+"\t"+SSTR(RNARow)+"\t"+SSTR(RNACol)].push_back(splitz[i+j*2+1+3]);
				if(regionToGene.find(splitz[i+j*2+1+3])!=regionToGene.end())
					regionToGene[splitz[i+j*2+1+3]].push_back(splitz[i+j*2+3]);	
				else {
					vector<string> temp;
					temp.push_back(splitz[i+j*2+3]);
					regionToGene[splitz[i+j*2+1+3]]=temp;
				}
            }
            for(int j = 0; j < AtacCount; j++) {
                vector<double> temp;
                temp.push_back(AtacRow);
                temp.push_back(AtacCol);
                temp.push_back(RNARow);
                temp.push_back(RNACol);
                temp.push_back(-1);
                temp.push_back(-1);
                temp.push_back(-1);
                temp.push_back(-1);
                pointcounts[SSTR(AtacRow)+"\t"+SSTR(AtacCol)+"\t"+SSTR(RNARow)+"\t"+SSTR(RNACol)]++;
                allpoints.push_back(temp);
            }
            i+=2*AtacCount+3;
        }
    }
	
	cout<<"Gathering Units in each cluster"<<endl;
	vector<vector<vector<string> > > UnitsInEachCluster;
	vector<vector<vector<string> > > GenesInEachCluster;
	vector<vector<string> > TotalGenesInEachCluster;
	vector<vector<vector<string> > > RegionsInEachCluster;
	vector<string> FullListOfGenes;
	
	for(int i = 0; i < kmeans1; i++) {
		vector<vector<string> > tempUnits;
		vector<vector<string> > tempGenes;
		vector<vector<string> > tempRegion;
		for(int j = 0; j < kmeans2; j++) {
			//cout<<j<<endl;
			vector<string> temp;
			vector<string> temp2;
			vector<string> temp3;
			vector<string> temp6;
			vector<string> temp7;
			vector<string> temp8;
			vector<int> temp4;
			vector<int> temp10;
			int temp9=0;
			int temp5=0;
			tempUnits.push_back(temp);
			tempGenes.push_back(temp2);
			tempRegion.push_back(temp6);
			TotalGenesInEachCluster.push_back(temp7);
		}
		UnitsInEachCluster.push_back(tempUnits);
		GenesInEachCluster.push_back(tempGenes);
		RegionsInEachCluster.push_back(tempRegion);
	}
	while(getline(ClusterFile,line)) {
		vector<string> splitz = split(line,'\t');
		cout<<line<<endl;
		if(splitz.size() < 5) {
			cout<<splitz.size()<<endl;
			cout<<line<<endl;
			continue;
		}
		int ATACCluster;
		int RNACluster;
		istringstream(splitz[4])>>ATACCluster;
		istringstream(splitz[5])>>RNACluster;
		string unit = splitz[0]+'\t'+splitz[1]+'\t'+splitz[2]+'\t'+splitz[3];
		//cout<<unit<<endl;	
		//cout<<ATACCluster<<endl;
		//cout<<RNACluster<<endl;
		//cout<<UnitsInEachCluster.size()<<endl;
		//cout<<UnitsInEachCluster[ATACCluster].size()<<endl;
		UnitsInEachCluster[ATACCluster][RNACluster].push_back(unit);
	}
	cout<<"Genes"<<endl;
	ofstream GeneTotalFile((outputprefix+"/GeneTotal.txt").c_str());
	ofstream RegionTotalFile((outputprefix+"/RegionTotal.txt").c_str());
	ofstream ComboTotalFile((outputprefix+"/SOMFusionTrack.bed").c_str());
	ComboTotalFile<<"track name=SOMFusion description=\"SOM Fusion Regions/Metaclusters\" visibility=2"<<endl;
	for(int i = 0; i < kmeans1; i++) {
		for(int j = 0; j < kmeans2; j++) {
			for(int k = 0; k < UnitsInEachCluster[i][j].size(); k++) {
				vector<string> splitz= split(UnitsInEachCluster[i][j][k],'\t');
				map<string,vector<string> >::iterator loc= allgenes.find(UnitsInEachCluster[i][j][k]);
				map<string,vector<string> >::iterator loc2= allregions.find(UnitsInEachCluster[i][j][k]);
				if(loc!=allgenes.end()) {
					vector<string> temp=vector<string>(loc->second);
					vector<string> temp2=vector<string>(loc2->second);
					for(int m = 0; m < temp2.size(); m++) {
						bool found = false;
						for(int u = 0; u < RegionsInEachCluster[i][j].size(); u++) {
                            if(RegionsInEachCluster[i][j][u].compare(temp2[m])==0) {
                                found = true;
                                break;
                            }
                        }
						if(!found)
							RegionsInEachCluster[i][j].push_back(temp2[m]);
					}

					for(int m = 0; m < temp.size(); m++) {
						bool found = false;
						for(int u = 0; u < GenesInEachCluster[i][j].size(); u++) {
							if(GenesInEachCluster[i][j][u].compare(temp[m])==0) {
								found = true;
								break;
							}
						}
						if(!found) {
							bool found2=false;
							GenesInEachCluster[i][j].push_back(temp[m]);
							for(int u = 0; u < FullListOfGenes.size(); u++) {
								if(temp[m].compare(FullListOfGenes[u]) == 0) {
									found=true;
									break;
								}
							}
							if(!found2) {
								FullListOfGenes.push_back(temp[m]);
							}
						}
						//if(i==9&&j==12) cout<<temp[m]<<endl;
						found = false;
						for(int u = 0; u < TotalGenesInEachCluster[j].size(); u++) {
							if(temp[m].compare(TotalGenesInEachCluster[j][u])==0) {
								found = true;
								break;
							}
						}
						if(!found) {
							TotalGenesInEachCluster[j].push_back(temp[m]);
						}
					}
				}
			}
			cout<<i<<'\t'<<j<<'\t'<<GenesInEachCluster[i][j].size()<<endl;
			ofstream outfile((outputprefix+"/Genes_"+SSTR(i)+"_"+SSTR(j)).c_str());
			ofstream outfile2((outputprefix+"/Regions_"+SSTR(i)+"_"+SSTR(j)).c_str());
			ofstream outfile3((outputprefix+"/Combo_"+SSTR(i)+"_"+SSTR(j)).c_str());
			for(int k = 0; k < GenesInEachCluster[i][j].size(); k++) {
				outfile<<GenesInEachCluster[i][j][k]<<endl;
			}
			if(j==0)
				GeneTotalFile<<GenesInEachCluster[i][j].size();
			else
				GeneTotalFile<<'\t'<<GenesInEachCluster[i][j].size();
			for(int k = 0; k < RegionsInEachCluster[i][j].size(); k++) {
				vector<string> splitz = split(RegionsInEachCluster[i][j][k],':');
				if(splitz.size() < 2) continue;
				vector<string> splitz2 = split(splitz[1],'-');
				if(Fuse)
					for(int m = k+1; m < RegionsInEachCluster[i][j].size(); m++) {
						vector<string> splitzcomp = split(RegionsInEachCluster[i][j][m],':');
						vector<string> splitzcomp2 = split(splitzcomp[1],'-');
						int start;
						int end;
						int compstart;
						int compend;
						istringstream(splitz2[0])>>start;
						istringstream(splitz2[1])>>end;
						istringstream(splitzcomp2[0])>>compstart;
						istringstream(splitzcomp2[1])>>compend;
						if(abs(start-compend)<=5 && splitz[0].compare(splitzcomp[0])==0) {
							splitz2[0]=SSTR(compstart);
							RegionsInEachCluster[i][j].erase(RegionsInEachCluster[i][j].begin()+m);
							m=k;
							continue;
						}
						if(abs(end-compstart)<=5 && splitz[0].compare(splitzcomp[0])==0) {
							splitz2[1]=SSTR(compend);
                                                        RegionsInEachCluster[i][j].erase(RegionsInEachCluster[i][j].begin()+m);
                                                        m=k;
                                                        continue;
						}
					}
				outfile2<<splitz[0]<<'\t'<<splitz2[0]<<'\t'<<splitz2[1]<<endl;
				for(int n = 0; n < regionToGene[RegionsInEachCluster[i][j][k]].size(); n++)
				{
					int found = -1;
					for(int m = 0; m < GenesInEachCluster[i][j].size(); m++) {
						if(GenesInEachCluster[i][j][m].compare(regionToGene[RegionsInEachCluster[i][j][k]][n])==0)
							found = m;
					}
					if(found >-1)
					{
						outfile3<<splitz[0]<<'\t'<<splitz2[0]<<'\t'<<splitz2[1]<<'\t'<<regionToGene[RegionsInEachCluster[i][j][k]][n]<<endl;
						ComboTotalFile<<splitz[0]<<'\t'<<splitz2[0]<<'\t'<<splitz2[1]<<'\t'<<SSTR(i)+"_"+SSTR(j)<<endl;
					}
				}
			}
			if(j==0)
                RegionTotalFile<<RegionsInEachCluster[i][j].size();
            else
                RegionTotalFile<<'\t'<<RegionsInEachCluster[i][j].size();

			outfile3.close();
			outfile2.close();
			outfile.close();
		}
		GeneTotalFile<<endl;
		RegionTotalFile<<endl;
	}
	GeneTotalFile.close();
	RegionTotalFile.close();
	ComboTotalFile.close();
	cout<<"Making Confusion Matrix."<<endl;
	vector<vector<double> > ConfusMatrix;
	for(int i = 0; i < kmeans1; i++) {
		int max=0;
		vector<double> row;
		for(int j = 0; j < kmeans2; j++) {
			row.push_back(GenesInEachCluster[i][j].size());
		}
		ConfusMatrix.push_back(row);
	}
	vector<vector<double> > ConfusMatrix2;
    for(int i = 0; i < kmeans1; i++) {
        int max=0;
        vector<double> row;
        for(int j = 0; j < kmeans2; j++) {
            row.push_back(RegionsInEachCluster[i][j].size());
        }
        ConfusMatrix2.push_back(row);
    }

	vector<int> AtacOrder;
	vector<int> RNAOrder;
	for(int i = 0; i < max(kmeans1,kmeans2); i++) {
		double max = -1;
		int maxj = -1;
		int maxk = -1;
		for(int j = 0; j < kmeans1; j++) {
			bool found = false;
			for(int k = 0; k < AtacOrder.size(); k++) {
				if(AtacOrder[k]==j) {
					found = true;
					break;
				}
			}
			if(!found) {
				for(int k = 0; k < kmeans2; k++) {
					bool found2 = false;
					for(int m = 0; m < RNAOrder.size(); m++) {
						if(RNAOrder[m]==k) {
							found2=true;
							break;
						}
					}
					if(!found2) {
						if(ConfusMatrix[j][k] > max) {
							max = ConfusMatrix[j][k];
							maxj=j;
							maxk=k;
						}
					}
				}
			}
		}
		if(maxj!=-1) {
			AtacOrder.push_back(maxj);
			RNAOrder.push_back(maxk);
		}
	}
	for(int i = 0; i < kmeans1; i++) {
        bool found2=false;
        for(int m = 0; m < AtacOrder.size(); m++) {
            if(AtacOrder[m]==i) {
                found2=true;
                break;
            }
        }
        if(!found2) AtacOrder.push_back(i);
    }

	for(int i = 0; i < kmeans2; i++) {
		bool found2=false;
		for(int m = 0; m < RNAOrder.size(); m++) {
			if(RNAOrder[m]==i) {
				found2=true;
                break;
            }
        }
		if(!found2) RNAOrder.push_back(i);
	}
	
	vector<int> AtacOrder2;
    vector<int> RNAOrder2;
    for(int i = 0; i < max(kmeans1,kmeans2); i++) {
        double max = -1;
        int maxj = -1;
        int maxk = -1;
        for(int j = 0; j < kmeans1; j++) {
            bool found = false;
            for(int k = 0; k < AtacOrder2.size(); k++) {
                if(AtacOrder2[k]==j) {
                    found = true;
                    break;
                }
            }
            if(!found) {
                for(int k = 0; k < kmeans2; k++) {
                    bool found2 = false;
                    for(int m = 0; m < RNAOrder2.size(); m++) {
                        if(RNAOrder2[m]==k) {
                            found2=true;
                            break;
                        }
                    }
                    if(!found2) {
                        if(ConfusMatrix2[j][k] > max) {
                            max = ConfusMatrix2[j][k];
                            maxj=j;
                            maxk=k;
                        }
                    }
                }
            }
        }
        if(maxj!=-1) {
            AtacOrder2.push_back(maxj);
            RNAOrder2.push_back(maxk);
        }
    }
	for(int i = 0; i < kmeans1; i++) {
        bool found2=false;
        for(int m = 0; m < AtacOrder2.size(); m++) {
            if(AtacOrder2[m]==i) {
                found2=true;
                break;
            }
        }
        if(!found2) AtacOrder2.push_back(i);
    }

    for(int i = 0; i < kmeans2; i++) {
        bool found2=false;
        for(int m = 0; m < RNAOrder2.size(); m++) {
            if(RNAOrder2[m]==i) {
                found2=true;
                break;
            }
        }
        if(!found2) RNAOrder2.push_back(i);
    }


	cout<<"Outputing Confusion."<<endl;
	ofstream GenesConfusion((outputprefix+"/GenesConfus.txt").c_str());	
	for(int i = 0; i < kmeans2; i++) {
		GenesConfusion<<'\t'<<RNAOrder[i];
	}
	GenesConfusion<<endl;
	for(int i = 0; i < kmeans1; i++) {
		cout<<i<<'\t'<<AtacOrder.size()<<'\t'<<RNAOrder.size()<<endl;
		GenesConfusion<<AtacOrder[i];
		for(int j = 0; j < kmeans2; j++) {
			cout<<j<<'\t'<<AtacOrder[i]<<'\t'<<RNAOrder[j]<<'\t'<<ConfusMatrix.size()<<'\t'<<ConfusMatrix[AtacOrder[i]].size()<<endl;
			GenesConfusion<<'\t'<<ConfusMatrix[AtacOrder[i]][RNAOrder[j]];
		}
		GenesConfusion<<endl;
	}
	GenesConfusion.close();

	cout<<"Outputing Confusion."<<endl;
    ofstream GenesConfusion2((outputprefix+"/RegionsConfus.txt").c_str());
    for(int i = 0; i < kmeans2; i++) {
        GenesConfusion2<<'\t'<<RNAOrder2[i];
    }
    GenesConfusion2<<endl;
    for(int i = 0; i < kmeans1; i++) {
        cout<<i<<'\t'<<AtacOrder2.size()<<'\t'<<RNAOrder2.size()<<endl;
        GenesConfusion2<<AtacOrder2[i];
        for(int j = 0; j < kmeans2; j++) {
            cout<<j<<'\t'<<AtacOrder2[i]<<'\t'<<RNAOrder2[j]<<'\t'<<ConfusMatrix2.size()<<'\t'<<ConfusMatrix2[AtacOrder[i]].size()<<endl;
            GenesConfusion2<<'\t'<<ConfusMatrix2[AtacOrder2[i]][RNAOrder2[j]];
        }
        GenesConfusion2<<endl;
    }
    GenesConfusion2.close();
}	
