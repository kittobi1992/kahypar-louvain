#include <iostream>
#include <fstream>
#include <bits/stdc++.h>

using namespace std;

typedef pair<int,int> ii;


class SuffixArray {

public:
  #define MAXN 5000010
  
  SuffixArray(vector<int>& T) : n(T.size()+1), T(T.size()+1), sa(T.size()+1) {
    T.push_back(-1);
    memset(ra,0,sizeof(ra)); memset(temp_ra,0,sizeof(temp_ra)); memset(sa_idx,0,sizeof(sa_idx));
    for(int i = 0; i < sa.size(); i++)	{ this->T[i] = T[i]; sa[i] = i; sa_idx[i] = i; }
    prefixDoublingSorting();
    //computeLCP();
  }
  
  int operator[](const int i) {
    return sa[i];
  }

  
  vector<int> getSuffixArray() {
    return sa;
  }
  
  void printSuffixArray(int t) {
    for(int i = 0; i < n; i++) {
      cout << "SA["<<i<<"]=" << sa[i] << " -> ";
      for(int j = sa[i]; j < sa[i] + t; j++)
	cout << T[j] << " ";
      cout << endl;
    }
  }
  
private:
  
  //Construct a Suffix Array for a string s (Based on the prefix doubling algorithm). Time complexity O(n*log^2(n))
  void prefixDoublingSorting() {
    vector<ii> ra_p(n);
    //Initiliaze ra with ascii values of chars in s
    for(int i = 0; i < n; i++)	ra[i] = T[i];
    for(int k = 1; k < n; k <<= 1) {
      for(int j = 0; j < n; j++) { 
	ra_p[j] = make_pair(ra[sa[j]],(sa[j] + k < n ? ra[sa[j] + k] : 0)); 
	sa_idx[sa[j]] = j;
      }
      //If we want to construct Suffix-Array in O(n*log(n)) we have to replace this naive sorting
      //stradegy with Counting-Sort
      sort(sa.begin(),sa.end(),[&](const int i, const int j) {
	return (ra_p[sa_idx[i]] < ra_p[sa_idx[j]]) || (ra_p[sa_idx[i]] == ra_p[sa_idx[j]] && sa_idx[i] < sa_idx[j]);
      });
      //Re-rank after sorting
      temp_ra[sa[0]] = 0; int r = 0;
      for(int j = 1; j < n; j++)
	temp_ra[sa[j]] = (ra[sa[j]] == ra[sa[j-1]] && ra[sa[j] + (sa[j] + k < n ? k : 0)] == ra[sa[j-1] + (sa[j-1] + k < n ? k : 0)]) ?
			 r : ++r;
      for(int j = 0; j < n; j++) ra[j] = temp_ra[j];
      if(ra[sa[n-1]] == n-1) break;
    }
  }

  
  int n;
  vector<int> T, sa;
  int ra[MAXN], temp_ra[MAXN], sa_idx[MAXN];
};



void sop(int i, int l, int j, int *lps, int *prev_occ, int bot){
	if( j == 0 and l ==0 and i==0 )
		return;
	assert(i>j);
	assert(i!=j);
	if( lps[i] == bot ){
		lps[i] = l;
		prev_occ[i] = j;
	}else{
		if( lps[i] < l ){
			if( prev_occ[i] > j )
				sop(prev_occ[i], lps[i], j, lps, prev_occ, bot);
			else
				sop(j, lps[i],prev_occ[i], lps, prev_occ, bot);
			lps[i] = l;
			prev_occ[i] = j;
		}else{
			if( prev_occ[i] > j )
				sop(prev_occ[i], l, j, lps, prev_occ, bot);
			else
				sop(j, l, prev_occ[i],lps, prev_occ, bot);			
		}
	}
}


int main(int argc, char* argv[])
{
    if (argc > 1) {
        string file  = string(argv[1]);
        string ofile = file + ".lz";
        if (argc > 2) {
            ofile = argv[2];
        }
      
	std::ifstream rfile(file);
	int n; rfile >> n;
	vector<int> x(n);
	for(int i = 0; i < n; i++) {
	  rfile >> x[i];
	}
	SuffixArray *sa_array = new SuffixArray(x);
	vector<int> sa_v = sa_array->getSuffixArray();
	int *sa = new int [sa_v.size()];
	for(int i = 0; i < sa_v.size(); i++)
	  sa[i] = sa_v[i];
	//sa->printSuffixArray(10);
        
        
	//compute LZ factorization
	// compute PHI
	bool printlz = false;
	n++;
	int *phi = new int [n];
	int *prev_occ = new int [n];
	for(int i=0; i<n; ++i) {
		prev_occ[i] = x[i];
	}
	int to_add[2] = {-1,n-1};
	for(int i=0; i<n; ++i){
	      phi[sa[i]] = sa[i+to_add[i==0]];
	}
	// sa holds now LPS
	for(int i=0; i<n; ++i)
	  sa[i] = -1;	   

	int maxFactorLength = 0;
	int l = 0;
	for(int i=0; i<n; ++i){
	  int j = phi[i];
	      while( x[i+l] == x[j+l] ) ++l;
	      if( i>j ){
		      sop(i,l,j,sa,prev_occ,-1);
	      }else{
		      sop(j,l,i,sa,prev_occ,-1);
	      }
	  if( l > 0 ) --l;	 
	}
	int numfactors = 0;
	int longFactors = 0;
	float averageFactorLength = 0; 

	sa[0] = 0;
	std::ofstream out_stream(ofile.c_str());
	out_stream << n << endl;
	for(int i=0; i<n; ){
	      ++numfactors;
	      if( sa[i] > maxFactorLength )
		      maxFactorLength = sa[i];
	      if( sa[i] > 8 )
		      ++longFactors;
	  if( sa[i] < 1 ){
		      out_stream << sa[i] << " " << -(x[i]+1) << " ";
		      averageFactorLength += 1;
		      ++i; 
	      }else{
		      out_stream << sa[i] << " " << prev_occ[i] << " ";
		      averageFactorLength += sa[i];
		      i += sa[i];
	      }
	}
	out_stream << "0 0" << endl;
	out_stream.close();
	printf("maxFactorLength = %d\n",maxFactorLength);
	printf("numfactors = %d\n",numfactors);
	printf("longFactors (len>8) = %d\n", longFactors);
	averageFactorLength /= numfactors;
	printf("averageFactorLength = %.2f\n", averageFactorLength);

	delete [] sa;
	delete [] phi;
	delete [] prev_occ;
	
    } else {
        cout << "Usage: " << argv[0] << " file [ofile]" << endl;
        cout << " Computes the SA from an array of 64-bit integers." << endl;
        cout << " Result is stored in `ofile`, or `file`.sa if `ofile`" << endl;
        cout << " is not specified." << endl;
    }
}

