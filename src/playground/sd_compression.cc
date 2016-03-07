#include <sdsl/bit_vectors.hpp>
#include <sdsl/vectors.hpp>
#include <random>
#include <iostream>
#include <chrono>
#include <sdsl/util.hpp>

using namespace sdsl;
using namespace std;

template<typename T>
void print_stats(vector<T> hg, string desc, int n) {
    float int_v_mb = 0;
    for(int i = 0; i < n; i++) {
      int_v_mb += size_in_mega_bytes(hg[i]);
    }
    cout << desc << " size in MB: " << int_v_mb << " MB" << endl;
}

int main(int argc, char* argv[])
{
  
    string file  = string(argv[1]);
    string ofile = file + ".lz";
    if (argc > 2) {
      ofile = argv[2];
    }
     
    std::ifstream rfile(file);
    int n; rfile >> n;
    vector<int_vector<>> hg(n,int_vector<>());
    for(int i = 0; i < n; i++) {
      int num; rfile >> num;
      hg[i].resize(num);
      for(int j = 0; j < num; j++) {
	int elem; rfile >> elem;
	hg[i].set_int(j,elem);
      }
    }
    rfile.close();
    
    print_stats<int_vector<>>(hg,"Uncompressed Int-Vector", n);
    for(int i = 0; i < n; i++) {
      util::bit_compress(hg[i]);
    }
    print_stats<int_vector<>>(hg,"Compressed Int-Vector", n);
  

    /*std::mt19937_64 rng;
    std::uniform_int_distribution<uint64_t> distribution(0, bv.size()-1);
    auto dice = bind(distribution, rng);
    // populate vectors with some other bits
    for (uint64_t i=0; i < bv.size()/25; ++i) {
        uint64_t x = dice();
        bv[x] = !default_value;
    }
    for(int i = 0; i < 20; i++)
      cout << bv[i] << endl;
    cout << "size in MiB: " << size_in_mega_bytes(bv) << endl;

    sd_vector<> bv_sd(bv);
    {
        bit_vector().swap(bv);
    }
    cout << "size in MiB: " << size_in_mega_bytes(bv_sd) << endl;
    cout << "wl = " << (size_t) bv_sd.wl << endl;
    cout << "n = " << bv_sd.size() << endl;
    cout << "2*m = " << bv_sd.high.size()<<endl;
    cout <<"n/m=" << (2.0*bv_sd.size())/bv_sd.high.size()<<endl;*/


}