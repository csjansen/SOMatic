/*convertEnsemblIDs: Used to convert Ensembl IDs to gene names
    Copyright (C) 2015  Camden Sinclair Jansen

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


int main(int argc, char* argv[]) {
 if(argc < 2) {
        cout << "Usage: ./convertEnsemblGeneIDs -Rows <Number of rows in your SOM> -Cols <Number of cols in your SOM> -InputPrefix <Input File Location> -GeneInfo <ENSEMBL geneID file> -OutputPrefix <Output File Location>" <<endl;
        return 0;
    }

	string inputprefix;
    int row;
    int col;
    string MusTermsFileName;
	string outputprefix;

	for(int i = 0; i < argc; i++) {
        string temp = argv[i];
        if(temp.compare("-Rows")==0)
			istringstream(argv[i+1])>>row;
		if(temp.compare("-Cols")==0)
            istringstream(argv[i+1])>>col;
		if(temp.compare("-InputPrefix")==0)
			inputprefix = argv[i+1];
		if(temp.compare("-GeneInfo")==0)
            MusTermsFileName = argv[i+1];
		if(temp.compare("-OutputPrefix")==0)
            outputprefix = argv[i+1];
	}
    ifstream MusTermsFile(MusTermsFileName.c_str());
    map<string, int> geneIds;
	map<string, string> genenames;
    cout<<"Getting gene ids"<<endl;
    int totalGenes=0;
	string line;
    while(getline(MusTermsFile,line)) {
        if(line[0]=='#') continue;
        vector<string> splitz = split(line,'\t');
        int geneid;
        /*istringstream(splitz[1])>>geneid;
			vector<string> splitz2 = split(splitz[5],'|');
			string temp;
			for(int i = 0; i < splitz2.size(); i++) {
				size_t found = splitz2[i].find("Ensembl:");
				if(found != string::npos) {
					vector<string> splitz3=split(splitz2[i],':');
					temp = splitz3[1];
					break;
				}
			}*/
			genenames[splitz[0]]=splitz[2];
			cout<<splitz[2]<<'\t'<<genenames[splitz[2]]<<endl;
			totalGenes++;
    }
	cout<<"Genes Loaded:" << totalGenes<<endl;
			//cout<<"Row: "<<i<<" Col: "<<j<<" "<<(inputprefix+"_"+SSTR(i)+"_"+SSTR(j)+".unit")<<endl;
			ifstream genefile((inputprefix).c_str());
			vector<string> genes;
            while(getline(genefile, line)) {
               /* bool found = false;
                for(int k = 0; k < genes.size(); k++) {
                    if(genes[k].compare(line)==0) {
                        found = true;
                        break;
                    }
                }
                if(!found)*/
					vector<string> splitz = split(line,'\t');
                    genes.push_back(splitz[0]);
					cout<<genenames[splitz[0]]<<endl;
            }
			cout<<genes.size()<<endl;

	        ofstream outfile2((outputprefix).c_str());
	        for(int k = 0 ; k < genes.size(); k++) {
				cout<<genes[k]<<'\t'<<genenames[genes[k]]<<endl;
				if(genenames[genes[k]].compare("")==0) {
					outfile2<<genes[k]<<endl;
				} else {
					outfile2<<genenames[genes[k]]<<endl;
				}
			}
}
