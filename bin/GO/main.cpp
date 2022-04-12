/* getGO: Converts gene names into GO term enrichments
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
    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<map>
#include<string>
#include<math.h>
#include<algorithm>


using namespace std;

//#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
//        ( std::ostringstream() << std::dec << x ) ).str()
//
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


void myReplace(std::string& str, const std::string& oldStr, const std::string& newStr){
  size_t pos = 0;
  while((pos = str.find(oldStr, pos)) != std::string::npos){
     str.replace(pos, oldStr.length(), newStr);
     pos += newStr.length();
  }
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

struct GO {
	string name;
	vector<vector<int> > coords;
	int total;
};

struct by_total { 
    bool operator()(GO const &a, GO const &b) { 
        return a.total > b.total;
    }
};

struct by_name {
    bool operator()(GO const &a, GO const &b) {
        return a.name.compare(b.name)<0;
    }
};

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
    //  cout<<i<<'\t'<<p<<endl;
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





int main(int argc, char* argv[]) {
	if(argc < 2) {
        cout << "Usage: ./getGO [options] -Rows <Number of rows in your SOM> -Cols <Number of cols in your SOM> -GenePrefix <Input File Location> -Gene2GO <Gene2GO File.  See main README.txt for format> -GeneInfo <GeneInfo File. see main README.txt for format> -GOFile <GO term file. See main README.txt for format> -OutputPrefix <Output File Location>" <<endl;
        cout << "Options: <default>" <<endl;
        cout << "-Sanity: If set to true, only GO terms with 5 genes in the unit will be reported. [true, false] <true>"<<endl;
		cout << "-OutputTotals: If set to a file location, will output total genes and GO terms in each unit. <>"<<endl;
        return 0;
    }

	string inputprefix;
    int row;
    int col;
	string GOTermsFileName="";
	string GOFileName="";
    string MusTermsFileName="";
	string outputprefix;
	string Sanity = "true";
	string outputTotalsName="";
	for(int i = 0; i < argc; i++) {
        string temp = argv[i];
        if(temp.compare("-Rows")==0)
			istringstream(argv[i+1])>>row;
		if(temp.compare("-Cols")==0)
            istringstream(argv[i+1])>>col;
		if(temp.compare("-GenePrefix")==0)
			inputprefix=argv[i+1];
		if(temp.compare("-Gene2GO")==0)
            GOTermsFileName=argv[i+1];
		if(temp.compare("-GeneInfo")==0)
            MusTermsFileName=argv[i+1];
		if(temp.compare("-GOFile")==0)
            GOFileName=argv[i+1];
		if(temp.compare("-OutputPrefix")==0)
            outputprefix=argv[i+1];
		if(temp.compare("-Sanity")==0)
            Sanity=argv[i+1];
		if(temp.compare("-OutputTotals")==0)
			outputTotalsName=argv[i+1];
	}
	ifstream GOTermsFile(GOTermsFileName.c_str());

	ifstream GOFile(GOFileName.c_str());
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
	string line;
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
    string tax_id="";
    ifstream MusTermsFile(MusTermsFileName.c_str());
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
    }
	vector<GO> GOs;
	vector<vector<int> > GOtotals;
    vector<vector<int> > genetotals;
	for(int i = 0; i < row; i++) {
        vector<int> Gototals;
        vector<int> genetotal;
        for(int j = 0; j < col; j++) {
            cout<<"Row: "<<i<<" Col: "<<j<<'\t'<<(inputprefix+"_"+SSTR(i)+"_"+SSTR(j)+".unit")<<endl;
            ifstream genefile((inputprefix+"_"+SSTR(i)+"_"+SSTR(j)+".unit").c_str());
            vector<string> genes;
            while(getline(genefile, line)) {
                bool found = false;
                vector<string> splitz = split(line,'\t');
                //cout<<splitz[0]<<endl;
                //int temp;
                //cin>>temp;
                for(int k = 0; k < genes.size(); k++) {
                    if(genes[k].compare(splitz[0])==0) {
                        found = true;
                        break;
                    }
                }
                if(!found)
                    genes.push_back(splitz[0]);
            }
            cout<<"Genes: "<<genes.size()<<endl;
			vector<string> GOlist;
			genetotal.push_back(genes.size());
			vector<int> GOcounts;
			int GoTotal=0;
			for(int k = 0; k < genes.size();k++) {
                string gene=genes[k];
                string upper="";
                for(string::size_type n = 0; n < gene.length(); n++)
                    upper+=toupper(gene[n],loc);
                map<string,string>::iterator it1 = geneIds.find(upper);
                //map<string,string>::iterator first = (geneIds.begin());
                  //  for(;first!=geneIds.end();first++) {
                    //  cout<<upper<<'\t'<<(string)first->first<<endl;
				//		cout<<upper<<endl;
				//		cout<<(string)first->first<<endl;
                	//}
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
                                for(int p = 0; p < GOlist.size(); p++) {
                                    if(GOlist[p].compare(temp[u])==0) {
                                        foundp=p;
                                        break;
                                    }
                                }
                                if(foundp>-1) {
                                    GOcounts[foundp]++;
                                    GoTotal++;
                                } else {
                                    GOlist.push_back(temp[u]);
                                    GOcounts.push_back(1);
                                    GoTotal++;
                                }
                            }
                        }
                    }
                }
            }
			vector<bioms>  Allbiom;
            for(int k = 0; k < GOlist.size(); k++) {
                double Totalnum=GO2Genes[GOlist[k]].size();

                double biom = Fishers(GOcounts[k],genes.size()-GOcounts[k],Totalnum,AllGenesSize-Totalnum);
                bioms temp;
                temp.biom = biom;
                temp.GO = GOlist[k];
                temp.GODesc = GOIDconversion[GOlist[k]];
                temp.GoCounts=GOcounts[k];
                temp.TotalGoCounts=genes.size()-GOcounts[k];
                temp.annot=Totalnum;

                Allbiom.push_back(temp);
            }
			sort(Allbiom.begin(), Allbiom.end());
            double alpha = .05;
			ofstream outfile2((outputprefix+"_"+SSTR(i)+"_"+SSTR(j)+".unit").c_str());
			Gototals.push_back(0);
            for(int k = 0; k < Allbiom.size(); k++) {
                if(Allbiom[k].biom > alpha/(double)(Allbiom.size()+1-k)) {
					 break;
                } else {
                    outfile2<<Allbiom[k].GO<<'\t'<<Allbiom[k].GODesc<<'\t'<<Allbiom[k].biom<<'\t'<<Allbiom[k].GoCounts<<'\t'<<Allbiom[k].TotalGoCounts<<'\t'<<Allbiom[k].annot<<'\t'<<AllGenesSize-Allbiom[k].annot<<endl;
					Gototals[Gototals.size()-1]++;
					int coordsfound = -1;
                    for(int m = 0; m < GOs.size(); m++) {
                        if(GOs[m].name.compare(Allbiom[k].GO)==0) {
                            coordsfound = m;
                            break;
                        }
                    }
                    if(coordsfound>=0) {
                        GOs[coordsfound].coords[i][j]++;
                        GOs[coordsfound].total++;
                    } else {
                        GO temp;
                        temp.name=Allbiom[k].GO;
                        vector<vector<int> > temp2;
                        for(int o = 0; o < row; o++) {
                            vector<int> temp3;
                            for(int p = 0; p < col; p++) {
                                temp3.push_back(0);
                            }
                            temp2.push_back(temp3);
                        }
                        temp2[i][j]++;
                        temp.coords=temp2;
                        temp.total=1;
                        GOs.push_back(temp);
                    }

                }
			}
			outfile2.close();
		}
		GOtotals.push_back(Gototals);
        genetotals.push_back(genetotal);

	}
	sort(GOs.begin(), GOs.end(), by_name());
    ofstream outfile4((outputprefix+"_List").c_str());
    for(int i = 0; i < GOs.size(); i++) {
        myReplace(GOs[i].name,"/"," ");
        ofstream outfile3((outputprefix+"_"+GOs[i].name+".map").c_str());
        for(int j = 0; j < row; j++) {
            for(int k = 0; k < col-1; k++) {
                outfile3<<GOs[i].coords[j][k]<<'\t';
            }
            outfile3<<GOs[i].coords[j][col-1]<<endl;
        }
        outfile4<<GOs[i].name<<'\t'<<GOs[i].total<<'\t'<<GOIDconversion[GOs[i].name]<<endl;
        outfile3.close();
    }
    outfile4.close();
    sort(GOs.begin(), GOs.end(), by_total());
    ofstream outfile5((outputprefix+"_List_Totals").c_str());
    for(int i = 0; i < GOs.size(); i++) {
        outfile5<<GOs[i].name<<'\t'<<GOs[i].total<<endl;
    }
    outfile5.close();
    if(outputTotalsName.compare("")!=0) {
        ofstream outfile6(outputTotalsName.c_str());
        for(int i = 0; i < row; i++) {
            for(int j = 0; j < col; j++) {
                outfile6<<i<<'\t'<<j<<'\t'<<genetotals[i][j]<<'\t'<<GOtotals[i][j]<<endl;
                cout<<i<<'\t'<<j<<'\t'<<genetotals[i][j]<<'\t'<<GOtotals[i][j]<<endl;
            }
        }
        outfile6.close();
    }
}
    /*ifstream GOTermsFile(GOTermsFileName.c_str());

	map<string, vector<string> > GoHeir;
	cout<<"Building GO Heirarchy"<<endl;
	ifstream GOFile(GOFileName.c_str());
	string line;
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
			//	cout<<splitz2[1]<<endl;
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
			//	cout<<ID<<endl;
			}
		}
	}
	cout<<"Go Heirarchy Number: "<<GOTermNumber<<endl;
	map<string, string> geneIds;
	map<string, string> geneNames;
    cout<<"Getting gene ids"<<endl;
    int totalGenes=0;
	string tax_id="";
	if(MusTermsFileName.compare("-GOFile")!=0) {
    ifstream MusTermsFile(MusTermsFileName.c_str());
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
        totalGenes++;
        /*} else if (NameType.compare("Ensembl")==0) {
            vector<string> splitz2 = split(splitz[5],'|');
            string temp;
            for(int i = 0; i < splitz2.size(); i++) {
                size_t found = splitz2[i].find("Ensembl:");
                if(found != string::npos) {
                    vector<string> splitz3=split(splitz2[i],':');
                    temp = splitz3[1];
                    break;
                }
            }
            string upper="";
            locale loc;
            for(string::size_type i = 0; i < temp.length(); i++)
                upper+=toupper(temp[i],loc);
            geneIds[upper]=splitz[1];
            /*cout<<temp<<endl;
            int temp2;
            cin>>temp2;*/
          /*  totalGenes++;
        } else if(NameType.compare("Xeno")==0) {
            string temp = splitz[0];
            string upper="";
            locale loc;
            for(string::size_type i = 0; i < temp.length(); i++)
                upper+=toupper(temp[i],loc);
            geneIds[upper]=splitz[2];
            //cout<<upper<<'\t'<<splitz[2]<<endl;
            totalGenes++;
        }*/

    /*}
    cout<<"Genes Loaded:" << totalGenes<<endl;
	}
	
	map<string, vector<string> > Goterms;
    vector<string> AllGoTerms;
	map<string,int> AllGoTermsCounts;
    cout<<"Getting Go terms"<<endl;
	int lines = 0;
    while(getline(GOTermsFile,line)) {
		/*if(NameType.compare("Xeno")==0) {
			vector<string> splitz = split(line, '\t');
			vector<string> splitz2 = split(splitz[0],'|');
			Goterms[splitz2[0]].push_back(splitz[3]);
			//cout<<splitz2[0]<<'\t'<<splitz[3]<<endl;
			AllGoTerms.push_back(splitz[3]);
		} else {*/
	/*		vector<string> splitz = split(line, '\t');
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
				//cout<<GOIDconversion[GoID]<<endl;
				AllGoTerms.push_back(GoID);
				AllGoTermsCounts[GoID]++;
				vector<string> moreGO = GoHeir[GoID];
				for(int i = 0; i < moreGO.size(); i++) {
					queue.push_back(moreGO[i]);
				}
			}
			//int stop;
			//cin>>stop;
		//}
    }
	cout<<"Lines loaded: "<<lines<<endl;
	vector<GO> GOs;
	vector<vector<int> > GOtotals;
	vector<vector<int> > genetotals;
	//vector<vector<vector<int> > > GOCoords;
	//vector<string> GONames;
	//vector<int> GOTotal;
	for(int i = 0; i < row; i++) {
		vector<int> Gototal;
		vector<int> genetotal;
		for(int j = 0; j < col; j++) {
			cout<<"Row: "<<i<<" Col: "<<j<<'\t'<<(inputprefix+"_"+SSTR(i)+"_"+SSTR(j)+".unit")<<endl;
			ifstream genefile((inputprefix+"_"+SSTR(i)+"_"+SSTR(j)+".unit").c_str());
			vector<string> genes;
            while(getline(genefile, line)) {
                bool found = false;
				vector<string> splitz = split(line,' ');
				//cout<<splitz[0]<<endl;
				//int temp;
				//cin>>temp;
                for(int k = 0; k < genes.size(); k++) {
                    if(genes[k].compare(splitz[0])==0) {
                        found = true;
                        break;
                    }
                }
                if(!found)
                    genes.push_back(splitz[0]);
            }
			cout<<"Genes: "<<genes.size()<<endl;
			vector<string> enrGOTerms;
	        vector<string> uniqueGOTerms;
			map<string, int> GOTermNumber;
    	    for(int k = 0; k < genes.size(); k++) {
				string upper="";
				locale loc;
	            for(string::size_type i = 0; i < genes[k].length(); i++)
					upper+=toupper(genes[k][i],loc);
				if(MusTermsFileName.compare("-GOFile")!=0) {
					map<string,string>::iterator it1 = geneIds.find(upper);
				//cout<<upper<<endl;
				//cout<<genes[k]<<'\t'<<it1->second<<endl;
				//	int temp;
				//	cin>>temp;
            	if(it1 != geneIds.end()) {
                	map<string, vector<string> >::iterator it2;
					//if(NameType.compare("Xeno")==0) {
					//	it2 = Goterms.find(genes[k]);
					//} else {
						it2 = Goterms.find((string)it1->second);
					//}
	                if(it2 != Goterms.end()) {
    	                vector<string> temp = (vector<string>)it2->second;
						//cout<<genes[k]<<endl;
        	            for(int k = 0; k < temp.size(); k++) {
							//cout<<GOIDconversion[temp[k]]<<endl;
            	            enrGOTerms.push_back(temp[k]);
							GOTermNumber[temp[k]]++;
                	        bool found = false;
                    	    for(int p = 0; p < uniqueGOTerms.size(); p++) {
                        	    if(uniqueGOTerms[p].compare(temp[k])==0) {
                            	    found = true;
                                	break;
	                            }
    	                    }
        	                if(!found) uniqueGOTerms.push_back(temp[k]);
            	        }
						//int stop;
						//cin>>stop;
                	}
	            }
				} else {
					map<string, vector<string> >::iterator it2;
					it2 = Goterms.find(genes[k]);
					if(it2 != Goterms.end()) {
						vector<string> temp = (vector<string>)it2->second;
						for(int k = 0; k < temp.size(); k++) {
							enrGOTerms.push_back(temp[k]);
                            GOTermNumber[temp[k]]++;
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
				}
				
    	    }
			cout<<"GO terms before enrichment: "<<uniqueGOTerms.size()<<endl;
	        ofstream outfile2((outputprefix+"_"+SSTR(i)+"_"+SSTR(j)+".unit").c_str());
			vector<string> GOTermNames;
			vector<vector<vector<bioms> > > Allbiom;
		    for(int i = 0; i < kmeans1; i++) {
        		vector<vector<bioms> > biorow;
        		for(int j = 0; j < kmeans2; j++) {
            		vector<bioms> biostore;
            		for(int k = 0; k < GOInEachCluster[i][j].size(); k++) {
                		double Totalnum=GO2Genes[GOInEachCluster[i][j][k]].size();

                		double biom = Fishers(GOCountsInEachCluster[i][j][k],GenesInEachCluster[i][j].size()-GOCountsInEachCluster[i][j][k],Totalnum,AllGenesSize-Totalnum);
                		bioms temp;
                		temp.biom = biom;
                		temp.GO = GOInEachCluster[i][j][k];
                		temp.GODesc = GOIDconversion[GOInEachCluster[i][j][k]];
                		temp.GoCounts=GOCountsInEachCluster[i][j][k];
                		temp.TotalGoCounts=GenesInEachCluster[i][j].size()-GOCountsInEachCluster[i][j][k];
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
        		for(int j = 0; j < kmeans2; j++) {
            		sort(Allbiom[i][j].begin(), Allbiom[i][j].end());
            		double alpha = .05;
            		ofstream outfile((outputprefix+"/GOTerms_"+SSTR(i)+"_"+SSTR(j)).c_str());
            		for(int k = 0; k < Allbiom[i][j].size(); k++) {
                		if(Allbiom[i][j][k].biom > alpha/(double)(Allbiom[i][j].size()+1-k)) {
                    		break;
                		} else {
                    		outfile<<Allbiom[i][j][k].GO<<'\t'<<Allbiom[i][j][k].GODesc<<'\t'<<Allbiom[i][j][k].biom<<'\t'<<Allbiom[i][j][k].GoCounts<<'\t'<<Allbiom[i][j][k].TotalGoCounts<<'\t'<<Allbiom[i][j][k].annot<<'\t'<<AllGenesSize-Allbiom[i][j][k].annot<<endl;

                		}
            		}
            		outfile.close();
        		}
    		}

	        /*vector<double> biomv;
			vector<int> Totalnums;
			double biomvnum = 0;
			for(int m = 0; m < uniqueGOTerms.size(); m++) {
				//cout<<m<<'\t'<<enrGOTerms.size()<<'\t'<<AllGoTerms.size()<<endl;
           		int Clusternum = 0;
            	int Totalnum = 0;
            	/*for(int k = 0; k < enrGOTerms.size(); k++)
                	if(uniqueGOTerms[m].compare(enrGOTerms[k])==0)
                    	Clusternum++;*/
			/*	Clusternum=GOTermNumber[uniqueGOTerms[m]];
				//Totalnum++;                     
            	/*for(int k = 0; k < AllGoTerms.size(); k++) {
					//cout<<uniqueGOTerms[m]<<'\t'<<AllGoTerms[k]<<endl;
                	if(uniqueGOTerms[m].compare(AllGoTerms[k])==0)
                    	Totalnum++;
				}*/
			/*	Totalnum=AllGoTermsCounts[uniqueGOTerms[m]];
				if(Sanity.compare("true")!=0 || Clusternum >= 5) {
	            	double annotationFactor = Totalnum/(double)AllGoTerms.size();
	            	
					double biom = getBinomPval(enrGOTerms.size(), Clusternum, annotationFactor);
					//cout<<uniqueGOTerms[m]<<'\t'<<enrGOTerms.size()<<'\t'<<Clusternum<<'\t'<<annotationFactor<<'\t'<<biom<<endl;
					//int temp;
					//cin>>temp;
		        	//cout<<genes.size()<<'\t'<<Clusternum<<'\t'<<annotationFactor<<'\t'<<biom<<endl;
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
        	for(int m = 0; m < biomv.size(); m++) {
            	if(biomv[m] < .05/biomvnum && biomv[m] < .05 && biomv[m]!=0) {
					biocount++;
				//cout<<GOTermNames[m]<<'\t'<<biomv[m]<<'\t'<<GOIDconversion[GOTermNames[m]]<<endl;
                	outfile2<<biomv[m]<<'\t'<<GOTermNames[m]<<'\t'<<GOIDconversion[GOTermNames[m]]<<'\t'<<Totalnums[m]<<'/'<<genes.size()<<endl;
					int coordsfound = -1;
	                for(int k = 0; k < GOs.size(); k++) {
		                if(GOs[k].name.compare(uniqueGOTerms[m])==0) {
			                coordsfound = k;
				            break;
					    }
					}
					if(coordsfound>=0) {
						GOs[coordsfound].coords[i][j]++;
	                    GOs[coordsfound].total++;
		            } else {
						GO temp;
			            temp.name=uniqueGOTerms[m];
				        vector<vector<int> > temp2;
					    for(int o = 0; o < row; o++) {
						    vector<int> temp3;
							for(int p = 0; p < col; p++) {
								temp3.push_back(0);
	                        }
		                    temp2.push_back(temp3);
			            }
				        temp2[i][j]++;
					    temp.coords=temp2;
						temp.total=1;
						GOs.push_back(temp);
					}
            	}
        	}*/
/*			Gototal.push_back(biocount);
			genetotal.push_back(genes.size());
			cout<<"Final amount: "<<biocount<<endl;
        	outfile2.close();
		}
		GOtotals.push_back(Gototal);
		genetotals.push_back(genetotal);
	}

	sort(GOs.begin(), GOs.end(), by_name());
	ofstream outfile4((outputprefix+"_List").c_str());
	for(int i = 0; i < GOs.size(); i++) {
		myReplace(GOs[i].name,"/"," ");
		ofstream outfile3((outputprefix+"_"+GOs[i].name+".map").c_str());
		for(int j = 0; j < row; j++) {
			for(int k = 0; k < col-1; k++) {
				outfile3<<GOs[i].coords[j][k]<<'\t';
			}
			outfile3<<GOs[i].coords[j][col-1]<<endl;
		}
		outfile4<<GOs[i].name<<'\t'<<GOs[i].total<<'\t'<<GOIDconversion[GOs[i].name]<<endl;
		outfile3.close();
	}
	outfile4.close();
	sort(GOs.begin(), GOs.end(), by_total());
	ofstream outfile5((outputprefix+"_List_Totals").c_str());
	for(int i = 0; i < GOs.size(); i++) {
        outfile5<<GOs[i].name<<'\t'<<GOs[i].total<<endl;
    }
	outfile5.close();
	if(outputTotalsName.compare("")!=0) {
		ofstream outfile6(outputTotalsName.c_str());
		for(int i = 0; i < row; i++) {
			for(int j = 0; j < col; j++) {
				outfile6<<i<<'\t'<<j<<'\t'<<genetotals[i][j]<<'\t'<<GOtotals[i][j]<<endl;
				cout<<i<<'\t'<<j<<'\t'<<genetotals[i][j]<<'\t'<<GOtotals[i][j]<<endl;
			}
		}
		outfile6.close();
	}
}*/
