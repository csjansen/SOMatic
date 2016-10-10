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
        cout << "Usage: ./metasomgene [options] -UnitPrefix <Unit files location> -MetaclusterFile <Metacluster File Location> -OutPrefix <Output File Location>" <<endl;
        cout << "Options: <default>" <<endl;
        cout << "-Rows: Number of rows for your SOM. <20>"<<endl;
        cout << "-Cols: Number of columns for your SOM. <30>"<<endl;
        cout << "-Metaclusters: Number of metaclusters. <20>"<<endl;
        return 0;
    }

	int kmeans1=0;
	int row1 = 20;
	int col1 = 30;
	string MetaClusterFileName;
	string outputprefix;
	string UnitPrefix;
	 for(int i = 0; i < argc; i++) {
        string temp = argv[i];
        if(temp.compare("-Rows")==0)
            istringstream(argv[i+1])>>row1;
        if(temp.compare("-Cols")==0)
            istringstream(argv[i+1])>>col1;
        if(temp.compare("-Metaclusters")==0)
            istringstream(argv[i+1])>>kmeans1;
		if(temp.compare("-MetaclusterFile")==0)
            MetaClusterFileName = argv[i+1];
		if(temp.compare("-OutPrefix")==0)
            outputprefix = argv[i+1];
		if(temp.compare("-UnitPrefix")==0)
            UnitPrefix = argv[i+1];
	} 
	
	ifstream MetaClusterFile(MetaClusterFileName.c_str());
	vector<vector<vector<int> > > MetaClusters;
	for(int i = 0; i < kmeans1; i++) {
		vector<vector<int> > temp;
		MetaClusters.push_back(temp);
	}
    cout<<"Opening file"<<endl;
	string line;
	while(getline(MetaClusterFile,line)) {
		if(line[0]=='#') continue;
		vector<string> splitz = split(line,'\t');
		vector<int> temp;
		int rows;
		int cols;
		int cluster;
		istringstream(splitz[0])>>rows;
		istringstream(splitz[1])>>cols;
		istringstream(splitz[2])>>cluster;
		temp.push_back(rows);
		temp.push_back(cols);
		MetaClusters[cluster].push_back(temp);
	}
	
	cout<<"Getting Genes"<<endl;
	vector<vector<string> > GenesInEachCluster;
	vector<vector<string> > RowsInEachCluster;
	vector<vector<string> > ColsInEachCluster;
	vector<string> FullListOfGenes;
	
	for(int i = 0; i < kmeans1; i++) {
		vector<string> tempGenes;
		vector<string> tempRows;
		vector<string> tempCols;
		GenesInEachCluster.push_back(tempGenes);
		RowsInEachCluster.push_back(tempRows);
		ColsInEachCluster.push_back(tempCols);
	}
	for(int i = 0; i < kmeans1; i++) {
		for(int k = 0; k < MetaClusters[i].size(); k++) {
			stringstream row;
			row<<MetaClusters[i][k][0];
			stringstream col;
			col<<MetaClusters[i][k][1];
			ifstream UnitFile((UnitPrefix+"_"+row.str()+"_"+col.str()+".unit").c_str());
			string line;
			while(getline(UnitFile,line)) {
				vector<string> splitz2=split(line,'\t');
				GenesInEachCluster[i].push_back(splitz2[0]);
				RowsInEachCluster[i].push_back(row.str());
				ColsInEachCluster[i].push_back(col.str());
				int found = -1;
				for(int j = 0; j < FullListOfGenes.size(); j++) {
					if(FullListOfGenes[j].compare(splitz2[0])==0) {
						found = j;
						break;
					}
				}
				if(found == -1) {
					FullListOfGenes.push_back(splitz2[0]);
				}
			}
			UnitFile.close();
		}
		cout<<i<<'\t'<<GenesInEachCluster[i].size()<<endl;
		ofstream outfile2((outputprefix+"/Genes_"+SSTR(i)).c_str());
	    for(int k = 0; k < GenesInEachCluster[i].size(); k++) {
		    outfile2<<GenesInEachCluster[i][k]<<'\t'<<RowsInEachCluster[i][k]<<'\t'<<ColsInEachCluster[i][k]<<endl;
		}
		outfile2.close();
	}
}
