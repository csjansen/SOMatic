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
        cout << "Usage: ./FusionBreakup [options] -FusionFile <Output from Fusion> -FusionClusterFile <Output from FusionCluster> -ClusterNum1 <Cluster Number from SOM 1> -ClusterNum2 <Cluster Number from SOM 2> -OutputPrefix <Prefix for Output Files> -Row1 -Col1 -Row2 -Col2" <<endl;
        return 0;
    }
	string clusterfile;
	int kmeans1;
	int kmeans2;
	string outputprefix;
	string ConnectionsFileName;
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
	}
	ifstream ClusterFile(clusterfile.c_str());
	//string GOTermsFileName = argv[4];//"../Kmeans_GO/gene2go";
    //string MusTermsFileName = argv[5];//"../Kmeans_GO/Homo_sapiens.gene_info";
    /*ifstream GOTermsFile(GOTermsFileName.c_str());
    ifstream MusTermsFile(MusTermsFileName.c_str());
    string tax_id = argv[6];// = "9606";
	string UnitPrefix = argv[7];
	string GOFileName=argv[8];
	string ConnectionsFileName=argv[9];
	int row1;
	int row2;
	int col1;
	int col2;*/
	//istringstream(argv[10])>>row1;
	/*istringstream(argv[11])>>col1;
	istringstream(argv[12])>>row2;
	istringstream(argv[13])>>col2;
	string outputprefix=argv[14];*/
	ifstream ConnectionsFile(ConnectionsFileName.c_str());
	vector<vector<vector<string> > > ConnectionGenes;
    string line;
	vector<vector<double> > allpoints;
    map<string,vector<string> > allgenes;
    map<string,vector<string> > allregions;
	map<string,string> regionToGene;
    cout<<"Opening file"<<endl;
    map<string,int> pointcounts;
    while(getline(ConnectionsFile, line)) {
      //  cout<<line<<endl;
        vector<string> splitz = split(line,'\t');
        int AtacRow;
        istringstream(splitz[0])>>AtacRow;
        int AtacCol;
        istringstream(splitz[1])>>AtacCol;
    //    cout<<AtacRow<<'\t'<<AtacCol<<endl;
        int i = 2;
        while(i < splitz.size()-2) {
  //          cout<<i<<'\t'<<splitz.size()<<endl;
            int RNARow;
            int RNACol;
            int AtacCount;
            istringstream(splitz[i])>>RNARow;
            istringstream(splitz[i+1])>>RNACol;
            istringstream(splitz[i+2])>>AtacCount;
//            cout<<RNARow<<'\t'<<RNACol<<endl;
            for(int j = 0; j < AtacCount; j++) {
				//cout<<j<<endl;
				//cout<<SSTR(AtacRow)+"\t"+SSTR(AtacCol)+"\t"+SSTR(RNARow)+"\t"+SSTR(RNACol)<<endl;
                allgenes[SSTR(AtacRow)+"\t"+SSTR(AtacCol)+"\t"+SSTR(RNARow)+"\t"+SSTR(RNACol)].push_back(splitz[i+j*2+3]);
				//cout<<i+j*2+1+3<<endl;
                allregions[SSTR(AtacRow)+"\t"+SSTR(AtacCol)+"\t"+SSTR(RNARow)+"\t"+SSTR(RNACol)].push_back(splitz[i+j*2+1+3]);
				//cout<<"3"<<endl;
				regionToGene[splitz[i+j*2+1+3]]=splitz[i+j*2+3];
            }
			//cout<<"second for"<<endl;
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
	
	/*ifstream GOFile(GOFileName.c_str());
cout<<"Getting Go terms"<<endl;
    map<string, vector<string> > Goterms;
    vector<string> AllGoTerms;
    map<string,int> AllGoTermsCounts;

    int lines = 0;
	map<string, vector<string> > GoHeir;
    cout<<"Building GO Heirarchy"<<endl;
    int mode = 0;
    string GO_ID;
    int GOTermNumber=0;
    vector<string> is_a;
    map<string, string> GOIDconversion;
    while(getline(GOFile,line)) {
        vector<string> splitz = split(line, ' ');
        if(mode == 0) {
            if(splitz.size() == 0) continue;
            if(splitz[0].compare("id:")==0) {
                GO_ID = splitz[1];
                mode = 1;
            }
        } else if(mode == 1) {
            if(splitz.size() == 0) continue;
            if(splitz[0].compare("name:")==0) {
                vector<string> splitz2 = split(line, ':');
                splitz2[1].erase(0,1);
                GOIDconversion[GO_ID]=splitz2[1];
            //  cout<<splitz2[1]<<endl;
                mode=2;
            }
        } else if(mode == 2) {
            if(splitz.size()==0) {
                GoHeir[GO_ID]=is_a;
                GOTermNumber++;
                mode = 0;
                is_a.clear();
            } else if(splitz[0].compare("is_obsolete:")==0) {
                mode = 0;
            } else if(splitz[0].compare("is_a:")==0) {
                string ID = splitz[1];
                is_a.push_back(ID);
            //  cout<<ID<<endl;
            }
        }
    }
    cout<<"Go Heirarchy Number: "<<GOTermNumber<<endl;
    map<string, string> geneIds;
    map<string, string> geneNames;
	vector<string> AllGenes;
	int AllGenesSize=0;
    cout<<"Getting gene ids"<<endl;
    int totalGenes=0;
	locale loc;
    //string tax_id="";
    //ifstream MusTermsFile(MusTermsFileName.c_str());
    while(getline(MusTermsFile,line)) {
        if(line[0]=='#') continue;
        vector<string> splitz = split(line,'\t');
        if(tax_id.compare("")==0) tax_id = splitz[0];
        string temp = splitz[2];
        string upper="";
        locale loc;
        for(string::size_type i = 0; i < temp.length(); i++)
            upper+=toupper(temp[i],loc);
        geneIds[upper]=splitz[1];
        geneNames[splitz[1]]=upper;
		//cout<<upper<<endl;
		//int temp2;
		//cin>>temp2;
		AllGenes.push_back(upper);
		AllGenesSize++;
        totalGenes++;
    }
    cout<<"Genes Loaded:" << totalGenes<<endl;
	//map<string, vector<string> > Goterms;
    //vector<string> AllGoTerms;
    //map<string,int> AllGoTermsCounts;
    cout<<"Getting Go terms"<<endl;

	while(getline(GOTermsFile,line)) {
    	vector<string> splitz = split(line, '\t');

            if(tax_id.compare(splitz[0])!=0&&tax_id.compare("")!=0) continue;
            string GoID = splitz[2];
            bool found = 0;
            for(int i = 0; i < Goterms[splitz[1]].size(); i++) {
                if(Goterms[splitz[1]][i].compare(GoID)==0) {
                    found = 1;
                    break;
                }
            }
            if(found) continue;
            lines++;
            Goterms[splitz[1]].push_back(GoID);
			//cout<<splitz[1]<<'\t'<<GoID<<'\t'<<GOIDconversion[GoID]<<endl;
            AllGoTerms.push_back(GoID);
            AllGoTermsCounts[GoID]++;
            vector<string> queue=GoHeir[GoID];
            while(queue.size() != 0) {
                GoID = queue[queue.size()-1];
                queue.pop_back();
                bool found = 0;
                for(int i = 0; i < Goterms[splitz[1]].size(); i++) {
                    if(Goterms[splitz[1]][i].compare(GoID)==0) {
                        found = 1;
                        break;
                    }
                }
                if(found) continue;

                Goterms[splitz[1]].push_back(GoID);
				//cout<<splitz[1]<<'\t'<<GoID<<'\t'<<GOIDconversion[GoID]<<endl;
                vector<string> moreGO = GoHeir[GoID];
                for(int i = 0; i < moreGO.size(); i++) {
                    queue.push_back(moreGO[i]);
                }
            }
			//int temp;
			//cin>>temp;
    }
	cout<<"Computing genes in each GO term"<<endl;
	map<string,vector<string> > GO2Genes;
	for(int i = 0; i < AllGenes.size(); i++) {
		map<string,string>::iterator it1 = geneIds.find(AllGenes[i]);
        if(it1 != geneIds.end()) {
			map<string, vector<string> >::iterator it2;
            it2 = Goterms.find((string)it1->second);
			if(it2 != Goterms.end()) {
				vector<string> temp = (vector<string>)it2->second;
				for(int u = 0 ; u < temp.size(); u++) {
					map<string, vector<string> >::iterator it3;
					it3 = GO2Genes.find(temp[u]);
					if(it3 == GO2Genes.end()) {
						vector<string> newlist;
						newlist.push_back(AllGenes[i]);
						GO2Genes[temp[u]]=newlist;
					} else {
						GO2Genes[temp[u]].push_back(AllGenes[i]);
					}
				}
			}
		} else {
			AllGenesSize--;
		}
	}*/
	cout<<"Gathering Units in each cluster"<<endl;
	vector<vector<vector<string> > > UnitsInEachCluster;
	vector<vector<vector<string> > > GenesInEachCluster;
	vector<vector<string> > TotalGenesInEachCluster;
	//vector<vector<vector<string> > > GOInEachCluster;
	//vector<vector<string> > TotalGOInEachCluster;
	vector<vector<vector<string> > > RegionsInEachCluster;
	//vector<vector<vector<int> > > GOCountsInEachCluster;
	//vector<vector<int> > TotalGOCountsInEachCluster;
	//vector<vector<int> > GOCountsTotalInEachCluster;
	//vector<int> TotalGOCountsTotalInEachCluster;
	vector<string> FullListOfGenes;
	//vector<string> FullListOfGO;
	//vector<int> FullListOfGOCounts;
	
	for(int i = 0; i < kmeans1; i++) {
		vector<vector<string> > tempUnits;
		vector<vector<string> > tempGenes;
	//	vector<vector<string> > tempGO;
		vector<vector<string> > tempRegion;
	//	vector<vector<int> > tempGOCounts;
	//	vector<int> tempGOCountsTotal;
		for(int j = 0; j < kmeans2; j++) {
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
	//		tempGO.push_back(temp3);
	//		tempGOCounts.push_back(temp4);
	//		tempGOCountsTotal.push_back(temp5);
			tempRegion.push_back(temp6);
			TotalGenesInEachCluster.push_back(temp7);
	//		TotalGOInEachCluster.push_back(temp8);
	//		TotalGOCountsTotalInEachCluster.push_back(temp9);
	//		TotalGOCountsInEachCluster.push_back(temp10);
		}
		UnitsInEachCluster.push_back(tempUnits);
		GenesInEachCluster.push_back(tempGenes);
	//	GOInEachCluster.push_back(tempGO);
	//	GOCountsInEachCluster.push_back(tempGOCounts);
	//	GOCountsTotalInEachCluster.push_back(tempGOCountsTotal);
		RegionsInEachCluster.push_back(tempRegion);
	}
	while(getline(ClusterFile,line)) {
		vector<string> splitz = split(line,'\t');
		int ATACCluster;
		int RNACluster;
		istringstream(splitz[4])>>ATACCluster;
		istringstream(splitz[5])>>RNACluster;
		string unit = splitz[0]+'\t'+splitz[1]+'\t'+splitz[2]+'\t'+splitz[3];
		UnitsInEachCluster[ATACCluster][RNACluster].push_back(unit);	
		//cout<<unit<<endl;
		//int temp2;
		//cin>>temp2;
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
				/*cout<<(UnitPrefix+"_"+splitz[2]+"_"+splitz[3]+".unit")<<endl;
				int temp;
				cin>>temp;*/
				/*ifstream UnitFile((UnitPrefix+"_"+splitz[2]+"_"+splitz[3]+".unit").c_str());
				while(getline(UnitFile,line)) {
					vector<string> splitz2=split(line,'\t');
					GenesInEachCluster[i][j].push_back(splitz2[0]);
				}
				UnitFile.close();*/
				map<string,vector<string> >::iterator loc= allgenes.find(UnitsInEachCluster[i][j][k]);
				map<string,vector<string> >::iterator loc2= allregions.find(UnitsInEachCluster[i][j][k]);
				if(loc!=allgenes.end()) {
					vector<string> temp=vector<string>(loc->second);
					vector<string> temp2=vector<string>(loc2->second);
					for(int m = 0; m < temp.size(); m++) {
						RegionsInEachCluster[i][j].push_back(temp2[m]);
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
				vector<string> splitz2 = split(splitz[1],'-');
				outfile2<<splitz[0]<<'\t'<<splitz2[0]<<'\t'<<splitz2[1]<<endl;
				outfile3<<splitz[0]<<'\t'<<splitz2[0]<<'\t'<<splitz2[1]<<'\t'<<regionToGene[RegionsInEachCluster[i][j][k]]<<endl;
				ComboTotalFile<<splitz[0]<<'\t'<<splitz2[0]<<'\t'<<splitz2[1]<<'\t'<<SSTR(i)+"_"+SSTR(j)<<endl;
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

	/*cout<<"Finding GO terms."<<endl;
	ofstream GOTotalFile((outputprefix+"/GOTotal.txt").c_str());
	for(int i = 0; i < AllGenes.size(); i++) {
		string gene=AllGenes[i];
		string upper="";
        for(string::size_type n = 0; n < gene.length(); n++)
			upper+=toupper(gene[n],loc);
        map<string,string>::iterator it1 = geneIds.find(upper);
		if(it1 != geneIds.end()) {
        	map<string, vector<string> >::iterator it2;
            it2 = Goterms.find((string)it1->second);
            if(it2 != Goterms.end()) {
            	vector<string> temp = (vector<string>)it2->second;
                //cout<<GenesInEachCluster[i][j][k]<<endl;
                for(int u = 0; u < temp.size(); u++) {
                	if(GO2Genes[temp[u]].size()<2000 && GO2Genes[temp[u]].size()>10) {
                    //cout<<temp[u]<<'\t'<<GOIDconversion[temp[u]]<<endl;
                    	int foundp = -1;
                        for(int p = 0; p < FullListOfGO.size(); p++) {
                        	if(FullListOfGO[p].compare(temp[u])==0) {
                            	foundp=p;
                                break;
                            }
                        }
                        if(foundp>-1) {
                        	FullListOfGOCounts[foundp]++;
                        } else {
                        	FullListOfGO.push_back(temp[u]);
                           	FullListOfGOCounts.push_back(1);
                        }
                    }
                }
                        //int temp2;
                        //cin>>temp2;
            }
        }
	}
	for(int i = 0; i < kmeans1; i++) {
		for(int j = 0; j < kmeans2;j++) {
			if(i==0) {
				for(int k = 0; k < TotalGenesInEachCluster[j].size(); k++) {
					string gene = TotalGenesInEachCluster[j][k];
					string upper="";
	                for(string::size_type n = 0; n < gene.length(); n++)
    	                upper+=toupper(gene[n],loc);
        	        map<string,string>::iterator it1 = geneIds.find(upper);
                //map<string,string>::iterator first = (geneIds.begin());
                    //for(;first!=geneIds.end();first++)
                    //  cout<<upper<<'\t'<<(string)first->first<<endl;
                    //int temp;
                    //cin>>temp;
            	    if(it1 != geneIds.end()) {
                	    map<string, vector<string> >::iterator it2;
                    	it2 = Goterms.find((string)it1->second);
	                    if(it2 != Goterms.end()) {
    	                    vector<string> temp = (vector<string>)it2->second;
            	            for(int u = 0; u < temp.size(); u++) {
                	            if(GO2Genes[temp[u]].size()<2000 && GO2Genes[temp[u]].size()>10) {
                    	        //cout<<temp[u]<<'\t'<<GOIDconversion[temp[u]]<<endl;
                        	        int foundp = -1;
                            	    for(int p = 0; p < GOInEachCluster[i][j].size(); p++) {
                                	    if(TotalGOInEachCluster[j][p].compare(temp[u])==0) {
                                    	    foundp=p;
                                        	break;
	                                    }
    	                            }
        	                        if(foundp>-1) {
            	                        TotalGOCountsInEachCluster[j][foundp]++;
                	                    TotalGOCountsTotalInEachCluster[j]++;
                    	            } else {
                        	            TotalGOInEachCluster[j].push_back(temp[u]);
                            	        TotalGOCountsInEachCluster[j].push_back(1);
                                	    TotalGOCountsTotalInEachCluster[j]++;
	                                }
    	                        }
        	                }
                        //int temp2;
                        //cin>>temp2;
            	        }
                	}
				}
			}
			for(int k = 0; k < GenesInEachCluster[i][j].size();k++) {
				string gene=GenesInEachCluster[i][j][k];
				string upper="";
                for(string::size_type n = 0; n < gene.length(); n++)
                    upper+=toupper(gene[n],loc);
				map<string,string>::iterator it1 = geneIds.find(upper);
				//map<string,string>::iterator first = (geneIds.begin());
					//for(;first!=geneIds.end();first++)
					//	cout<<upper<<'\t'<<(string)first->first<<endl;
					//int temp;
					//cin>>temp;
				if(it1 != geneIds.end()) {
                    map<string, vector<string> >::iterator it2;
                    it2 = Goterms.find((string)it1->second);
                    if(it2 != Goterms.end()) {
                        vector<string> temp = (vector<string>)it2->second;
						//cout<<GenesInEachCluster[i][j][k]<<endl;
                        for(int u = 0; u < temp.size(); u++) {
							if(GO2Genes[temp[u]].size()<2000 && GO2Genes[temp[u]].size()>10) {
							//cout<<temp[u]<<'\t'<<GOIDconversion[temp[u]]<<endl;
								int foundp = -1;
								for(int p = 0; p < GOInEachCluster[i][j].size(); p++) {
									if(GOInEachCluster[i][j][p].compare(temp[u])==0) {
										foundp=p;
										break;
									}
								}
								if(foundp>-1) {
									GOCountsInEachCluster[i][j][foundp]++;
									GOCountsTotalInEachCluster[i][j]++;
								} else {
									GOInEachCluster[i][j].push_back(temp[u]);
									GOCountsInEachCluster[i][j].push_back(1);
									GOCountsTotalInEachCluster[i][j]++;
								}
							}
						}
						//int temp2;
						//cin>>temp2;
					}
				}
			}
			cout<<i<<'\t'<<j<<'\t'<<GOInEachCluster[i][j].size()<<endl;	
			//cout<<j<<'\t'<<GOInEachCluster[i][j].size()<<endl;
		}
	}*/
/*	vector<vector<vector<bioms> > > Allbiom;
	vector<vector<bioms> > TotalBioms;
	for(int i = 0; i < kmeans1; i++) {
		vector<vector<bioms> > biorow;
		if(i==0) {
			for(int j = 0; j < kmeans2; j++) {
	            vector<bioms> biostore;
    	        for(int k = 0; k < TotalGOInEachCluster[j].size(); k++) {
        	        string GOTerm = TotalGOInEachCluster[j][k];
            	    int index=-1;
                	for(int u = 0; u < FullListOfGO.size(); u++) {
                    	if(FullListOfGO[u].compare(GOTerm)==0) {
                        	index=u;
	                        break;
    	                }
        	        }
            	    double Totalnum=FullListOfGOCounts[index];

                //double biom = getBinomPval(GOCountsTotalInEachCluster[i][j], GOCountsInEachCluster[i][j][k], annotationFactor);
                	double biom = Fishers(TotalGOCountsInEachCluster[j][k],TotalGenesInEachCluster[j].size(),Totalnum,AllGenes.size());
					cout<<FullListOfGO[index]<<'\t'<<biom<<endl;;
     	           bioms temp;
        	        temp.biom = biom;
            	    temp.GO = TotalGOInEachCluster[j][k];
                	temp.GODesc = GOIDconversion[TotalGOInEachCluster[j][k]];
	                temp.GoCounts=TotalGOCountsInEachCluster[j][k];
    	            temp.TotalGoCounts=TotalGenesInEachCluster[j].size();
        	        temp.annot=Totalnum;

            	    biostore.push_back(temp);
                //if(biom < .05) {///GOInEachCluster[i][j].size()) {
                //}
	            }
    	        TotalBioms.push_back(biostore);
        	}

		}
		for(int j = 0; j < kmeans2; j++) {
			vector<bioms> biostore;
			for(int k = 0; k < GOInEachCluster[i][j].size(); k++) {
				string GOTerm = GOInEachCluster[i][j][k];
				int index=-1;
				for(int u = 0; u < FullListOfGO.size(); u++) {
					if(FullListOfGO[u].compare(GOTerm)==0) {
						index=u;
						break;
					}
				}
				double Totalnum=FullListOfGOCounts[index];

                //double biom = getBinomPval(GOCountsTotalInEachCluster[i][j], GOCountsInEachCluster[i][j][k], annotationFactor);
                double biom = Fishers(GOCountsInEachCluster[i][j][k],GenesInEachCluster[i][j].size(),Totalnum,AllGenes.size());
				bioms temp;
				temp.biom = biom;
				temp.GO = GOInEachCluster[i][j][k];
				temp.GODesc = GOIDconversion[GOInEachCluster[i][j][k]];
				temp.GoCounts=GOCountsInEachCluster[i][j][k];
				temp.TotalGoCounts=GenesInEachCluster[i][j].size();
				temp.annot=Totalnum;

				biostore.push_back(temp);
				//if(biom < .05) {///GOInEachCluster[i][j].size()) {
				//}
			}
			biorow.push_back(biostore);
		}
		Allbiom.push_back(biorow);
	}
	for(int i = 0; i < kmeans1; i++) {
		if(i==0) {
			for(int j = 0; j < kmeans2; j++) {
         	   sort(TotalBioms[j].begin(), TotalBioms[j].end());
            	double alpha = .05;
	            int counter=0;
    	        ofstream outfile((outputprefix+"/GOTerms_"+SSTR(j)).c_str());
        	    for(int k = 0; k < TotalBioms[j].size(); k++) {
            	    if(TotalBioms[j][k].biom > alpha/(double)(TotalBioms[j].size()+1-k)) {
                	    break;
	                } else {
    	                outfile<<TotalBioms[j][k].GO<<'\t'<<TotalBioms[j][k].GODesc<<'\t'<<TotalBioms[j][k].biom<<'\t'<<TotalBioms[j][k].GoCounts<<'\t'<<TotalBioms[j][k].TotalGoCounts<<'\t'<<TotalBioms[j][k].annot<<'\t'<<AllGenes.size()<<endl;
        	            counter++;
            	    }
            	}	
	            if(j==0)
    	            GOTotalFile<<counter;
        	    else
            	    GOTotalFile<<'\t'<<counter;

	            outfile.close();
    	    }
		}
		for(int j = 0; j < kmeans2; j++) {
			sort(Allbiom[i][j].begin(), Allbiom[i][j].end());
			double alpha = .05;
			int counter=0;
			ofstream outfile((outputprefix+"/GOTerms_"+SSTR(i)+"_"+SSTR(j)).c_str());
			for(int k = 0; k < Allbiom[i][j].size(); k++) {
				if(Allbiom[i][j][k].biom > alpha/(double)(Allbiom[i][j].size()+1-k)) {
					break;
				} else {
					outfile<<Allbiom[i][j][k].GO<<'\t'<<Allbiom[i][j][k].GODesc<<'\t'<<Allbiom[i][j][k].biom<<'\t'<<Allbiom[i][j][k].GoCounts<<'\t'<<Allbiom[i][j][k].TotalGoCounts<<'\t'<<Allbiom[i][j][k].annot<<'\t'<<AllGenes.size()<<endl;
					counter++;
				}
			}
			if(j==0)
                GOTotalFile<<counter;
            else
                GOTotalFile<<'\t'<<counter;

			outfile.close();
		}
		GOTotalFile<<endl;
	}
	GOTotalFile.close();*/
}
	//srand ( unsigned ( std::time(0) ) );
    /*string infilename = argv[1];
    int row1;
    int col1;
    int row2;
    int col2;
    istringstream(argv[2])>>row1;
    istringstream(argv[3])>>col1;
    istringstream(argv[4])>>row2;
    istringstream(argv[5])>>col2;
	vector<int> maxs;
	maxs.push_back(row1);
	maxs.push_back(col1);
	maxs.push_back(row2);
	maxs.push_back(col2);
	string ATACSom = argv[8];
	string RNASom = argv[9];	
	
	ifstream som1file(ATACSom.c_str());
	vector<vector<double> > Som1;
    string line;
	double average1=0;
	vector<double> averages1;
	double stddev1;
	vector<double> stddevs1;
	vector<double> values1;
	int count = 0;
	int counts = 0;
	while(getline(som1file,line)) {
		if(line[0]=='#') continue;
		counts++;
		vector<double> Som1Row;
		vector<string> splitz = split(line, '\t');
		for(int i = 2; i < splitz.size(); i++) {
			if(averages1.size() == 0) {
				for(int j = 2; j < splitz.size();j++) {
					averages1.push_back(0);
				}
			}
			double temp;
			istringstream(splitz[i])>>temp;
			//if(temp>2) temp=2;
			average1+=temp;
			averages1[i-2]+=temp;
			values1.push_back(temp);
			count++;
			Som1Row.push_back(temp);
		}
		Som1.push_back(Som1Row);
		if(Som1Row.size()==0) cout<<line<<endl;
	}
	average1/=(double)count;
	for(int i = 0; i < averages1.size(); i++) {
		averages1[i]/=(double)counts;
	}
	for(int i = 0; i < values1.size(); i++) {
		stddev1+=pow(values1[i]-average1,2);
	}
	stddev1/=(double)values1.size();
	stddev1 = sqrt(stddev1);
	for(int i = 0; i < averages1.size(); i++) {
		stddevs1.push_back(0);
		for(int j=0; j < Som1.size(); j++) {
			stddevs1[i]+=pow(Som1[j][i]-averages1[i],2);
		}
		stddevs1[i]/=(double)counts;
		stddevs1[i]=sqrt(stddevs1[i]);
	}
	ifstream som2file(RNASom.c_str());
	vector<vector<double> > Som2;
    double average2=0;
    vector<double> averages2;
    double stddev2;
    vector<double> stddevs2;
    vector<double> values2;
    count = 0;
    counts = 0;
    while(getline(som2file,line)) {
        if(line[0]=='#') continue;
        counts++;
        vector<double> Som2Row;
        vector<string> splitz = split(line, '\t');
        for(int i = 2; i < splitz.size(); i++) {
            if(averages2.size() == 0) {
                for(int j = 2; j < splitz.size();j++) {
                    averages2.push_back(0);
                }
            }
            double temp;
            istringstream(splitz[i])>>temp;
            //if(temp>2) temp=2;
            average2+=temp;
            averages2[i-2]+=temp;
            values2.push_back(temp);
            count++;
            Som2Row.push_back(temp);
        }
        Som2.push_back(Som2Row);
        if(Som2Row.size()==0) cout<<line<<endl;
    }
    average2/=(double)count;
    for(int i = 0; i < averages2.size(); i++) {
        averages2[i]/=(double)counts;
    }
    for(int i = 0; i < values2.size(); i++) {
        stddev2+=pow(values2[i]-average2,2);
    }
    stddev2/=(double)values2.size();
    stddev2 = sqrt(stddev2);
    for(int i = 0; i < averages2.size(); i++) {
        stddevs2.push_back(0);
        for(int j=0; j < Som2.size(); j++) {
            stddevs2[i]+=pow(Som2[j][i]-averages2[i],2);
        }
        stddevs2[i]/=(double)counts;
        stddevs2[i]=sqrt(stddevs2[i]);
    }



//	cout<<averages1[0]<<'\t'<<stddevs1[0]<<endl;
//	cout<<averages2[0]<<'\t'<<stddevs2[0]<<endl;

	cout<<infilename<<endl;	
    ifstream infile(infilename.c_str());
    /*vector<vector<vector<vector<double> > > > points;
    for(int i = 0; i < row1; i++) {
        vector<vector<vector<double> > > temp3;
        for(int j = 0; j < col1; j++) {
            vector<vector<double> > temp2;
            for(int k = 0; k < row2; k++) {
                vector<double> temp;
                for(int o = 0; o < col2; o++) {
                    temp.push_back(0);
                }
                temp2.push_back(temp);
            }
            temp3.push_back(temp2);
        }
        points.push_back(temp3);
    }*/
	/*vector<vector<double> > allpoints;
	map<string,vector<string> > allgenes;
	map<string,vector<string> > allregions;
	cout<<"Opening file"<<endl;
	map<string,int> pointcounts;
    while(getline(infile, line)) {
		//cout<<line<<endl;
        vector<string> splitz = split(line,'\t');
        int AtacRow;
        istringstream(splitz[0])>>AtacRow;
        int AtacCol;
        istringstream(splitz[1])>>AtacCol;
		//cout<<AtacRow<<'\t'<<AtacCol<<endl;
		int i = 2;
		while(i < splitz.size()-2) {
			//cout<<i<<'\t'<<splitz.size()<<endl;
            int RNARow;
            int RNACol;
            int AtacCount;
            istringstream(splitz[i])>>RNARow;
            istringstream(splitz[i+1])>>RNACol;
            istringstream(splitz[i+2])>>AtacCount;
			//cout<<RNARow<<'\t'<<RNACol<<endl;
			for(int j = 0; j < AtacCount; j++) {
				allgenes[SSTR(AtacRow)+"."+SSTR(AtacCol)+"."+SSTR(RNARow)+"."+SSTR(RNACol)].push_back(splitz[i+j*2+3]);
				allregions[SSTR(AtacRow)+"."+SSTR(AtacCol)+"."+SSTR(RNARow)+"."+SSTR(RNACol)].push_back(splitz[i+j*2+1+3]);
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
				pointcounts[SSTR(AtacRow)+"."+SSTR(AtacCol)+"."+SSTR(RNARow)+"."+SSTR(RNACol)]++;
				allpoints.push_back(temp);
			}	
			i+=2*AtacCount+3;
        }
    }
    infile.close();
	cout<<"Shuffling"<<endl;
	std::srand ( unsigned ( std::time(0) ) );
  std::vector<int> myvector;	
	random_shuffle(allpoints.begin(),allpoints.end());
	//int kmeans=(int)round(sqrt(allpoints.size()/2));
	int kmeans=25;//round(sqrt(row1*col1/2.0));
	cout<<"Kmeans size 1: "<<kmeans<<endl;
	//cout<<Som1[0].size()<<'\t'<<Som2[0].size()<<endl;
	vector<vector<double> > means1;
	int counter1=0;
	for(int i = 0; i < kmeans; i++) {
		bool done = false;
		while(!done) {
			done=true;
			vector<double> temp1;
			for(int j = 0; j < Som1[0].size(); j++) temp1.push_back(Som1[allpoints[i+counter1][0]*col1+allpoints[i+counter1][1]][j]);
			for(int k = 0; k < means1.size();k++) {
				double dist = 0;
				for(int u = 0; u < temp1.size(); u++) {
					dist+=means1[k][u]-temp1[u];
				}
				if(dist == 0) {
					done = false;
					counter1++;
					break;
				}
			}
			if(done)
				means1.push_back(temp1);
		}
	}
	random_shuffle(allpoints.begin(),allpoints.end());
	
	vector<vector<double> > means2;
	int counter2=0;
    for(int i = 0; i < kmeans; i++) {
        bool done = false;
        while(!done) {
            done=true;
            vector<double> temp2;
            for(int j = 0; j < Som2[0].size(); j++) temp2.push_back(Som2[allpoints[i+counter2][2]*col2+allpoints[i+counter2][3]][j]);
            for(int k = 0; k < means2.size();k++) {
                double dist = 0;
                for(int u = 0; u < temp2.size(); u++) {
                    dist+=means2[k][u]-temp2[u];
                }
                if(dist == 0) {
                    done = false;
					counter2++;
                    break;
                }
            }
            if(done)
                means2.push_back(temp2);
        }
    }

	cout<<"Clustering"<<endl;
	bool finished = false;
	int iter = 0;
	bool cont = false;
	while(!finished) {
		finished = true;
		int count=0;
		//Assignment
		for(int i = 0; i < allpoints.size(); i++) {
			double lowestdist1=999999999;
			double lowestdist2=999999999;
			int lowestk1=-1;
			int lowestk2=-1;
			for(int j = 0; j < kmeans; j++) {
				//double dist = pow(hexdist(row1,col1,means[j][0],means[j][1],allpoints[i][0],allpoints[i][1]),2)/(double)(row1*col1)+pow(hexdist(row2,col2,means[j][2],means[j][3],allpoints[i][2],allpoints[i][3]),2)/(double)(row2*col2);
				//double dist1 = 0;
				//double maxdist1 = 0;
				//cout<<stddev1<<endl;
				//for(int k = 0; k < Som2[0].size(); k++) {
					//dist1 += pow(means1[j][k]-Som2[allpoints[i][2]*col2+allpoints[i][3]][k],2);
						//dist1 += pow((means[j][k]-average1)/stddev1-(Som1[allpoints[i][0]*col1+allpoints[i][1]][k]-average1)/stddev1,2);
						//double tempdist1 = pow((means1[j][k]-averages1[k])/stddevs1[k]-(Som1[allpoints[i][0]*col1+allpoints[i][1]][k]-averages1[k])/stddevs1[k],2)/(double)Som1[0].size();
						//dist1+=tempdist1;
						//if(maxdist1<pow(means1[j][k]-Som1[allpoints[i][0]*col1+allpoints[i][1]][k],2)) maxdist1=pow(means1[j][k]-Som1[allpoints[i][0]*col1+allpoints[i][1]][k],2);
				//}
//				dist1-=maxdist1;
				//cout<<dist<<endl;
				/*double dist2 = 0;
				double maxdist2=0;
				for(int k = 0; k < Som2[0].size(); k++) {
                    //dist2 += pow(means[j][k+Som1[0].size()]-Som2[allpoints[i][2]*col2+allpoints[i][3]][k],2);
						double tempdist2 = pow((means2[j][k]-averages2[k])/stddevs2[k]-(Som2[allpoints[i][2]*col2+allpoints[i][3]][k]-averages2[k])/stddevs2[k],2)/(double)Som2[0].size();
						dist2 += tempdist2;
						if(maxdist2<tempdist2) maxdist2=tempdist2;
                }
				dist2-=maxdist2;*/
				//cout<<dist1<<'\t'<<dist2<<endl;
				//double dist=pow(dist1,2)+pow(dist2,2); 
				//double dist = max(dist1/Som1[0].size(),dist2/Som2[0].size());
				//double dist = dist1;
				
				//cout<<means[j][0]<<'\t'<<means[j][1]<<'\t'<<allpoints[i][0]<<'\t'<<allpoints[i][1]<<'\t'<<hexdist(row1,col1,means[j][0],means[j][1],allpoints[i][0],allpoints[i][1])<<endl;
				//pow(means[j][0]-allpoints[i][0],2)+pow(means[j][1]-allpoints[i][1],2)+pow(means[j][2]-allpoints[i][2],2)+pow(means[j][3]-allpoints[i][3],2);
/*				double similarity=0;
                double mag1 = 0;
                double mag2 = 0;
				
                for(int num = 0; num < Som1[0].size(); num++) {	
                    similarity += means1[j][num]*Som1[allpoints[i][0]*col1+allpoints[i][1]][num];
                    mag1 += pow(means1[j][num],2);
                    mag2 += pow(Som1[allpoints[i][0]*col1+allpoints[i][1]][num],2);
				}
				if(mag1==0 || mag2==0) {
                    similarity = -1;
                } else {
                    similarity /= (sqrt(mag1) * sqrt(mag2));
                }
                double dist1 = 1-similarity;
				/*double dist1 = 0;
				for(int num = 0; num < Som1[0].size(); num++) {
					dist1+=pow(means1[j][num]-Som1[allpoints[i][0]*col1+allpoints[i][1]][num],2);
				}*/
/*				if(lowestdist1>dist1) {
                    lowestdist1 = dist1;
                    lowestk1 = j;
                }
				similarity=0;
				mag1=0;
				mag2=0;
				//cout<<similarity<<'\t'<<1<<endl;
				//cout<<Som1[0].size()<<'\t'<<Som2[0].size()<<endl;
				//cout<<means1[j].size()<<'\t'<<endl;
				/*double dist2 = 0;
                for(int num = 0; num < Som2[0].size(); num++) {
                    dist2+=pow(means2[j][num]-Som2[allpoints[i][2]*col2+allpoints[i][3]][num],2);
                }*/

  /*              for(int num = 0; num < Som2[0].size(); num++) {
                    similarity += means2[j][num]*Som2[allpoints[i][2]*col2+allpoints[i][3]][num];
                    mag1 += pow(means2[j][num],2);
                    mag2 += pow(Som2[allpoints[i][2]*col2+allpoints[i][3]][num],2);
                }
				//cout<<similarity<<'\t'<<2<<endl;
                if(mag1==0 || mag2==0) {
                    similarity = -1;
                } else {
                    similarity /= (sqrt(mag1) * sqrt(mag2));
                }
                double dist2 = 1-similarity;
				//double dist = (dist1+dist2)/2.0;
				if(lowestdist2>dist2) {
					lowestdist2 = dist2;
					lowestk2 = j;
				}
			}
			if(allpoints[i][4] != lowestk1) {
				count++;
				allpoints[i][4] = lowestk1;
				allpoints[i][5] = lowestdist1;
				finished = false;
			}
			if(allpoints[i][6] != lowestk2) {
                count++;
                allpoints[i][6] = lowestk2;
				allpoints[i][7] = lowestdist2;
                finished = false;
            }
	
		}
		cout<<count<<endl;
		/*if(finished) {
			for(int i = 0; i < kmeans; i++) {
				int total = 0;
				for(int j = 0; j < allpoints.size(); j++) {
					if(allpoints[j][4]==i) total++;
				}
				if(total < 40) {
					finished = false;
					kmeans--;
					//cout<<means1[1][0]<<endl;
					means1.erase(means1.begin()+i);
					//means2.erase(means2.begin()+i);
					//cout<<means1[1][0]<<endl;
					//cout<<"Cluster "<<i<<" is too small (Size "<<total<<"), deleting and re-running.  Cluster Number:"<<kmeans<<endl;
					break;
				}
			}
			if(!finished) continue;
		}*/
	/*	if(finished) break;
			
		cont=true;
		//Update
		for(int i = 0; i < kmeans; i++) {
			vector<double> totals1;
			int count = 0;
			for(int j = 0; j < Som1[0].size(); j++) totals1.push_back(0);
			for(int j = 0; j < allpoints.size(); j++) {
				if(allpoints[j][4]==i) {
					for(int k = 0; k < Som1[0].size(); k++) totals1[k]+=Som1[allpoints[j][0]*col1+allpoints[j][1]][k];
					count++;
				}
			}
			for(int j = 0; j < totals1.size(); j++) means1[i][j] = totals1[j]/(double)count;
		}
		for(int i = 0; i < kmeans; i++) {
            vector<double> totals2;
            int count = 0;
            for(int j = 0; j < Som2[0].size(); j++) totals2.push_back(0);
            for(int j = 0; j < allpoints.size(); j++) {
                if(allpoints[j][6]==i) {
                    for(int k = 0; k < Som2[0].size(); k++) totals2[k]+=Som2[allpoints[j][2]*col2+allpoints[j][3]][k];
                    count++;
                }
            }
            for(int j = 0; j < totals2.size(); j++) means2[i][j] = totals2[j]/(double)count;
        }
		iter++;
		if(iter%1000000==0) cout<<iter<<" iterations"<<endl;
	}

	vector<vector<int> > kmeanscounts;
	for(int i = 0; i < kmeans; i++) {
		vector<int> temp;
		for(int j = 0; j < kmeans; j++) {
			temp.push_back(0);
		}
		kmeanscounts.push_back(temp);
	}
	string outfilename = argv[6];
	ofstream outfile(outfilename.c_str());
	for(int i = 0; i < allpoints.size(); i++) {
		outfile<<allpoints[i][0]<<'\t'<<allpoints[i][1]<<'\t'<<allpoints[i][2]<<'\t'<<allpoints[i][3]<<'\t'<<allpoints[i][4]<<'\t'<<allpoints[i][6];
		vector<string> regions=allregions[SSTR(allpoints[i][0])+"."+SSTR(allpoints[i][1])+"."+SSTR(allpoints[i][2])+"."+SSTR(allpoints[i][3])];
		vector<string> genes=allgenes[SSTR(allpoints[i][0])+"."+SSTR(allpoints[i][1])+"."+SSTR(allpoints[i][2])+"."+SSTR(allpoints[i][3])];
		kmeanscounts[allpoints[i][4]][allpoints[i][6]]+=genes.size();
		outfile<<'\t'<<pointcounts[SSTR(allpoints[i][0])+"."+SSTR(allpoints[i][1])+"."+SSTR(allpoints[i][2])+"."+SSTR(allpoints[i][3])];

		outfile<<endl;
	}
	outfile.close();
	string countsout = "counts.out";
	ofstream countsoutfile(countsout.c_str());
	for(int i = 0; i < kmeans; i++) {
		countsoutfile<<kmeanscounts[i][0];
		for(int j = 1; j < kmeans; j++) {
			countsoutfile<<'\t'<<kmeanscounts[i][j];
		}
		countsoutfile<<endl;
	}
	for(int i = 0; i < kmeans; i++) {
		for(int j = 0; j < kmeans; j++) {
			ofstream genesoutputfile(("Genes/Genes"+SSTR(i)+"."+SSTR(j)).c_str());
			for(int k = 0; k < allpoints.size(); k++) {
				if(allpoints[k][4]==i && allpoints[k][6]==j) {
					genesoutputfile<<allpoints[k][0]<<'\t'<<allpoints[k][1]<<'\t'<<allpoints[k][2]<<'\t'<<allpoints[k][3]<<'\t'<<allpoints[k][4]<<'\t'<<allpoints[k][5]<<'\t'<<allpoints[k][6]<<'\t'<<allpoints[k][7]<<endl;
					vector<string> regions=allregions[SSTR(allpoints[k][0])+"."+SSTR(allpoints[k][1])+"."+SSTR(allpoints[k][2])+"."+SSTR(allpoints[k][3])];
			        vector<string> genes=allgenes[SSTR(allpoints[k][0])+"."+SSTR(allpoints[k][1])+"."+SSTR(allpoints[k][2])+"."+SSTR(allpoints[k][3])];
					for(int u = 0; u < genes.size(); u++) {
						genesoutputfile<<regions[u]<<'\t'<<genes[u]<<endl;
					}
				}
			}
			genesoutputfile.close();
		}
    }
}
	/*string clusterpoints = "clusterpoints.txt";
	ofstream clusterfile(clusterpoints.c_str());
	for(int i = 0; i < means1.size(); i++) {
		clusterfile<<means1[i][0]<<'\t'<<means1[i][1]<<'\t'<<means[i][2]<<'\t'<<means[i][3]<<'\t'<<i<<endl;
	}
	clusterfile.close();*/
    /*map<string, vector<string> > GoHeir;
	
    cout<<"Building GO Heirarchy"<<endl;
	string GOFileName="go.obo";
    ifstream GOFile(GOFileName.c_str());
    int mode = 0;
    string GO_ID;
    int GOTermNumber=0;
    vector<string> is_a;
    map<string, string> GOIDconversion;
    while(getline(GOFile,line)) {
        vector<string> splitz = split(line, ' ');
        if(mode == 0) {
            if(splitz.size() == 0) continue;
            if(splitz[0].compare("id:")==0) {
                GO_ID = splitz[1];
                mode = 1;
            }
        } else if(mode == 1) {
            if(splitz.size() == 0) continue;
            if(splitz[0].compare("name:")==0) {
                vector<string> splitz2 = split(line, ':');
                splitz2[1].erase(0,1);
                GOIDconversion[GO_ID]=splitz2[1];
            //  cout<<splitz2[1]<<endl;
                mode=2;
            }
        } else if(mode == 2) {
            if(splitz.size()==0) {
                GoHeir[GO_ID]=is_a;
                GOTermNumber++;
                mode = 0;
                is_a.clear();
            } else if(splitz[0].compare("is_obsolete:")==0) {
                mode = 0;
            } else if(splitz[0].compare("is_a:")==0) {
                string ID = splitz[1];
                is_a.push_back(ID);
            //  cout<<ID<<endl;
            }
        }
    }

	
	string GOTermsFileName = "../Kmeans_GO/gene2go";
	string MusTermsFileName = "../Kmeans_GO/Homo_sapiens.gene_info";
	ifstream GOTermsFile(GOTermsFileName.c_str());
	ifstream MusTermsFile(MusTermsFileName.c_str());
	string tax_id = "9606";
	/*cout<<"Getting Go terms"<<endl;	
	while(getline(GOTermsFile,line)) {
		vector<string> splitz = split(line, '\t');
		int test;
		istringstream(splitz[0])>>test;
		if(test!=tax_id) continue;
		int geneid;
		istringstream(splitz[1])>>geneid;
		string GoTerm = splitz[5];
		Goterms[geneid].push_back(GoTerm);	
		AllGoTerms.push_back(GoTerm);
	}*/
	/*cout<<"Getting Go terms"<<endl;
	map<string, vector<string> > Goterms;
    vector<string> AllGoTerms;
    map<string,int> AllGoTermsCounts;

    int lines = 0;
    while(getline(GOTermsFile,line)) {
        /*if(NameType.compare("Xeno")==0) {
            vector<string> splitz = split(line, '\t');
            vector<string> splitz2 = split(splitz[0],'|');
            Goterms[splitz2[0]].push_back(splitz[3]);
            //cout<<splitz2[0]<<'\t'<<splitz[3]<<endl;
            AllGoTerms.push_back(splitz[3]);
        } else {*/
        /*    vector<string> splitz = split(line, '\t');
            //cout<<tax_id<<" "<<splitz[0]<<endl;

            if(tax_id.compare(splitz[0])!=0&&tax_id.compare("")!=0) continue;
            //int geneid;
            //istringstream(splitz[1])>>geneid;
            string GoID = splitz[2];
            bool found = 0;
            for(int i = 0; i < Goterms[splitz[1]].size(); i++) {
                if(Goterms[splitz[1]][i].compare(GoID)==0) {
                    found = 1;
                    break;
                }
            }
            if(found) continue;
            lines++;
            Goterms[splitz[1]].push_back(GoID);
            //cout<<geneNames[splitz[1]]<<endl;
            //cout<<GOIDconversion[GoID]<<endl;
            AllGoTerms.push_back(GoID);
            vector<string> queue=GoHeir[GoID];
            while(queue.size() != 0) {
                GoID = queue[queue.size()-1];
                queue.pop_back();
                bool found = 0;
                for(int i = 0; i < Goterms[splitz[1]].size(); i++) {
                    if(Goterms[splitz[1]][i].compare(GoID)==0) {
                        found = 1;
                        break;
                    }
                }
                if(found) continue;

                Goterms[splitz[1]].push_back(GoID);
                //cout<<GOIDconversion[GoID]<<endl;
                AllGoTerms.push_back(GoID);
                vector<string> moreGO = GoHeir[GoID];
                for(int i = 0; i < moreGO.size(); i++) {
                    queue.push_back(moreGO[i]);
                }
            }
            //int stop;
            //cin>>stop;
        //}
    }

	cout<<AllGoTerms.size()<<endl;
	map<string, string> geneIds;
	cout<<"Getting gene ids"<<endl;
	int totalGenes=0;
	while(getline(MusTermsFile,line)) {
		if(line[0]=='#') continue;
		vector<string> splitz = split(line,'\t');
		string geneid;
		geneid=splitz[1];//istringstream(splitz[1])>>geneid;
		string temp = splitz[2];
		string upper="";
		locale loc;
		for(string::size_type i = 0; i < temp.length(); i++)
			upper+=toupper(temp[i],loc);
		geneIds[upper]=geneid;
		totalGenes++;
		map<string, vector<string> >::iterator it2;
                //if(NameType.compare("Xeno")==0) {
                //  it2 = Goterms.find(genes[k]);
                //} else {
        it2 = Goterms.find(geneid);
                //}
        if(it2 != Goterms.end()) {
        	vector<string> temp = (vector<string>)it2->second;
                    //cout<<genes[k]<<endl;
            for(int k = 0; k < temp.size(); k++) {
                        //cout<<GOIDconversion[temp[k]]<<endl;
            	AllGoTermsCounts[temp[k]]++;
			}
                        //int stop;
                        //cin>>stop;
		}
	}
	cout<<totalGenes<<endl;
	cout<<"Doing enrichments"<<endl;
	string prefix = argv[7];
	string outprefix = argv[10];
	for(int i = 0; i < kmeans; i++) {
	for(int k = 0; k < kmeans; k++) {
		//int i = 2;
		vector<string> genes;
		vector<string> genomicCoords;
		vector<vector<int> > coords;
		for(int j = 0; j < allpoints.size(); j++) {
			if(allpoints[j][4]!=i || allpoints[j][5]!=k) continue;
			/*ifstream unitfile((prefix+SSTR(allpoints[j][2])+"_"+SSTR(allpoints[j][3])+".unit").c_str());
			string line;
			vector<string> uniques;
			while(getline(unitfile, line)) {
				bool found = false;
				for(int k = 0; k < genes.size(); k++) {
					if(genes[k].compare(line)==0) {
						found = true;
						break;
					}
				}
				if(!found) 
					genes.push_back(line);
			}*/
/*				vector<string> genestorun = allgenes[SSTR(allpoints[j][0])+"."+SSTR(allpoints[j][1])+"."+SSTR(allpoints[j][2])+"."+SSTR(allpoints[j][3])];
				vector<string> regionstorun = allregions[SSTR(allpoints[j][0])+"."+SSTR(allpoints[j][1])+"."+SSTR(allpoints[j][2])+"."+SSTR(allpoints[j][3])];

				for(int q = 0; q < genestorun.size(); q++) {
					bool found = false;
					for(int p = 0; p < genes.size(); p++) {
						if((genes[p]+genomicCoords[p]).compare(genestorun[q]+regionstorun[q])==0) {
							found = true;
							break;
						}
					}
                
					if(!found) {
						genes.push_back(genestorun[q]);
						genomicCoords.push_back(regionstorun[q]);
						coords.push_back(allpoints[j]);
					}
				}
		}
		cout<<"genes: "<<i<<'\t'<<k<<'\t'<<genes.size()<<endl;
		//for(int j = 0; j < genes.size(); j++) {
		//	cout<<genes[j]<<endl;
		//}
		ofstream genesoutputfile(("Genes/Genes"+SSTR(i)+"."+SSTR(k)).c_str());
		for(int j = 0; j < genes.size(); j++) {
			genesoutputfile<<genes[j]<<'\t'<<genomicCoords[j]<<'\t'<<coords[j][0]<<'\t'<<coords[j][1]<<'\t'<<coords[j][2]<<'\t'<<coords[j][3]<<'\t'<<coords[j][4]<<'\t'<<coords[j][5]<<endl;
		}	
		vector<string> enrGOTerms;
		vector<string> uniqueGOTerms;
		map<string, int> GOTermNumber;
		for(int j = 0; j < genes.size(); j++) {
			/*map<string,int>::iterator it1 = geneIds.find(genes[j]);
			if(it1 != geneIds.end()) {
				map<int, vector<string> >::iterator it2 = Goterms.find((int)it1->second);
				if(it2 != Goterms.end()) {
					vector<string> temp = (vector<string>)it2->second;
					for(int k = 0; k < temp.size(); k++) {
						enrGOTerms.push_back(temp[k]);
						bool found = false;
						for(int p = 0; p < uniqueGOTerms.size(); p++) {
							if(uniqueGOTerms[p].compare(temp[k])==0) {
								found = true;
								break;
							}
						}
						if(!found) uniqueGOTerms.push_back(temp[k]);
					}
				}
			}*/
			/*string upper="";
                locale loc;
                for(string::size_type u = 0; u < genes[j].length(); u++)
                    upper+=toupper(genes[j][u],loc);

			map<string,string>::iterator it1 = geneIds.find(upper);
            //cout<<upper<<endl;
            //cout<<genes[k]<<'\t'<<it1->second<<endl;
            //  int temp;
            //  cin>>temp;
            if(it1 != geneIds.end()) {
            	map<string, vector<string> >::iterator it2;
                //if(NameType.compare("Xeno")==0) {
                //  it2 = Goterms.find(genes[k]);
                //} else {
                it2 = Goterms.find((string)it1->second);
                //}
                if(it2 != Goterms.end()) {
                	vector<string> temp = (vector<string>)it2->second;
                    //cout<<genes[k]<<endl;
                    for(int u = 0; u < temp.size(); u++) {
                    	//cout<<GOIDconversion[temp[k]]<<endl;
                        enrGOTerms.push_back(temp[u]);
                        GOTermNumber[temp[u]]++;
                        bool found = false;
                        for(int p = 0; p < uniqueGOTerms.size(); p++) {
                        	if(uniqueGOTerms[p].compare(temp[u])==0) {
                            	found = true;
                                break;
                            }
                        }
                        if(!found) uniqueGOTerms.push_back(temp[u]);
                    }
                        //int stop;
                        //cin>>stop;
                }
            }
		}
		ofstream outfile2((outprefix+SSTR(i)+"."+SSTR(k)).c_str());
		/*vector<string> GOTermNames;
		vector<double> biomv;
		for(int j = 0; j < uniqueGOTerms.size(); j++) {
			int Clusternum = 0;
			int Totalnum = 0;
			for(int k = 0; k < enrGOTerms.size(); k++)
				if(uniqueGOTerms[j].compare(enrGOTerms[k])==0)
					Clusternum++;
			for(int k = 0; k < AllGoTerms.size(); k++) 
				if(uniqueGOTerms[j].compare(AllGoTerms[k])==0)
					Totalnum++;
			double annotationFactor = Totalnum/(double)totalGenes;
			double biom = getBinomPval(genes.size(), Clusternum, annotationFactor);
			GOTermNames.push_back(uniqueGOTerms[j]);
			biomv.push_back(biom);
		}*/
		/*vector<string> GOTermNames;
        vector<double> biomv;
        vector<int> Totalnums;
        double biomvnum = 0;
        for(int m = 0; m < uniqueGOTerms.size(); m++) {
        	//cout<<m<<'\t'<<enrGOTerms.size()<<'\t'<<AllGoTerms.size()<<endl;
            int Clusternum = 0;
            int Totalnum = 0;
            /*for(int k = 0; k < enrGOTerms.size(); k++)
                    if(uniqueGOTerms[m].compare(enrGOTerms[k])==0)
                        Clusternum++;*/
          /*  Clusternum=GOTermNumber[uniqueGOTerms[m]];
            //Totalnum++;
            /*for(int k = 0; k < AllGoTerms.size(); k++) {
                    //cout<<uniqueGOTerms[m]<<'\t'<<AllGoTerms[k]<<endl;
                    if(uniqueGOTerms[m].compare(AllGoTerms[k])==0)
                        Totalnum++;
                }*/
            /*Totalnum=AllGoTermsCounts[uniqueGOTerms[m]];
			string Sanity = "true";
            if(Sanity.compare("true")!=0 || Clusternum >= 5) {
            	double annotationFactor = Totalnum/(double)totalGenes;

                double biom = getBinomPval(genes.size(), Clusternum, annotationFactor);
                    //cout<<uniqueGOTerms[m]<<'\t'<<Clusternum<<'\t'<<'\t'<<genes.size()<<'\t'<<endl;
					//cout<<Totalnum<<'\t'<<totalGenes<<'\t'<<annotationFactor<<'\t'<<Clusternum/(double)genes.size()<<'\t'<<biom<<endl;
                    //int temp;
                    //cin>>temp;
                GOTermNames.push_back(uniqueGOTerms[m]);
                biomv.push_back(biom);
                biomvnum++;
                Totalnums.push_back(Clusternum);
            } else {
            	GOTermNames.push_back(uniqueGOTerms[m]);
                biomv.push_back(1);
                Totalnums.push_back(Clusternum);
            }
       }
		int biocount = 0;
		vector<string> Enriched;
		vector<double> pvals;
			for(int m = 0; m < biomv.size(); m++) {
                if(biomv[m] < .05/biomvnum && biomv[m] < .05 && biomv[m]!=0) {
                    biocount++;
                //cout<<GOTermNames[m]<<'\t'<<biomv[m]<<'\t'<<GOIDconversion[GOTermNames[m]]<<endl;
                    outfile2<<biomv[m]<<'\t'<<GOTermNames[m]<<'\t'<<GOIDconversion[GOTermNames[m]]<<'\t'<<Totalnums[m]<<'/'<<genes.size()<<endl;
					Enriched.push_back(GOTermNames[m]);
					pvals.push_back(biomv[m]);
					
                }
            }
		for(int j = 0; j < biomv.size(); j++) {
			if(biomv[j] < .05/biomv.size()) {
				outfile2<<biomv[j]<<'\t'<<GOTermNames[j]<<endl;
			}
		}
		outfile2.close();
		/*vector<treeNode*> Tree;
		filltree(&Tree, Enriched, pvals, GoHeir);
		ofstream outfile3((outprefix+SSTR(i)+"."+SSTR(k)+".dot").c_str());
		cout<<Tree.size()<<endl;
		outfile3<<"digraph GO_"<<i<<" {"<<endl;
		outfile3<<"size=\"6,6\";"<<endl;
		outfile3<<"node [shape=box];"<<endl;
		for(int j = 0; j < Tree.size(); j++) {
			if(Tree[j]->pval!=1)
				outfile3<<Tree[j]->GO_ID<<" [color=crimson, style=filled, label=\""<<GOIDconversion[Tree[j]->GO_ID]<<"\"];"<<endl;
			else
				outfile3<<Tree[j]->GO_ID<<" [color=white, style=filled, label=\""<<GOIDconversion[Tree[j]->GO_ID]<<"\"];"<<endl;
		}
		for(int j = 0; j < Tree.size(); j++) {
			for(int u = 0; u < Tree[u]->children.size(); u++) {
				outfile3<<Tree[j]->GO_ID<<"->"<<Tree[j]->children[u]->GO_ID<<endl;
			}
		}
		outfile3<<"}"<<endl;*/
//	}
//	}
//}	
