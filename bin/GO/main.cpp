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

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

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
	string GOTermsFileName;
	string GOFileName;
    string MusTermsFileName;
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
    ifstream MusTermsFile(MusTermsFileName.c_str());

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

    }
    cout<<"Genes Loaded:" << totalGenes<<endl;

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
			lines++;
			vector<string> splitz = split(line, '\t');
			//cout<<tax_id<<" "<<splitz[0]<<endl;
			if(tax_id.compare(splitz[0])!=0) continue;
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
    	    }
			cout<<"GO terms before enrichment: "<<uniqueGOTerms.size()<<endl;
	        ofstream outfile2((outputprefix+"_"+SSTR(i)+"_"+SSTR(j)+".unit").c_str());
			vector<string> GOTermNames;
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
				Clusternum=GOTermNumber[uniqueGOTerms[m]];
				//Totalnum++;                     
            	/*for(int k = 0; k < AllGoTerms.size(); k++) {
					//cout<<uniqueGOTerms[m]<<'\t'<<AllGoTerms[k]<<endl;
                	if(uniqueGOTerms[m].compare(AllGoTerms[k])==0)
                    	Totalnum++;
				}*/
				Totalnum=AllGoTermsCounts[uniqueGOTerms[m]];
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
        	}
			Gototal.push_back(biocount);
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
}
