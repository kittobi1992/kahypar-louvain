#include <iostream>
#include <fstream>
#include <bits/stdc++.h>

using namespace std;



int main(int argc, char* argv[])
{
    if (argc > 1) {
        string file  = string(argv[1]);
        string ofile = file + ".lz";
      
	std::ifstream rfile(file);
	int n; rfile >> n;
	vector<int> x(n);
	for(int i = 0; i < n; i++) {
	  rfile >> x[i];
	}
	rfile.close();
	
	std::ifstream rofile(ofile);
	int lz_n = 0; rofile >> lz_n;
	lz_n = 0;
	vector<int> lps(n), prev_occ(n);
	while(rofile >> lps[lz_n] >> prev_occ[lz_n]) {
	  lz_n++;
	}
	rofile.close();
	
	vector<int> x_lz(n);
	int pos = 0;
	for(int i = 0; i < lz_n; i++) {
	    if(prev_occ[i] < 0) {
	      x_lz[pos++] = -prev_occ[i] - 1;
	    } else {
	      for(int j = prev_occ[i]; j < prev_occ[i]+lps[i]; j++)
		x_lz[pos++] = x_lz[j];
	    }
	}
	
	for(int i = 0; i < n; i++) {
	  if(x[i] != x_lz[i]) {
	    cout << "Mismatch at position " << i << " with x["<<i<<"]="<<x[i] << " and x_lz["<<i<<"]="<<x_lz[i] << endl;
	  }
	}
	
	
    } 
}

