#include "stdio.h" 
#include <iostream>
#include <vector> 
#include <omp.h>

using namespace std;

int main(int argc,char * argv[])  
{  
    omp_set_num_threads(6); 
    vector<int> vt(100000000, 1000);
    long long tot = 0;
    //#pragma omp parallel for
    for(long long i = 0; i < 100000000; i ++) {
        tot += vt[i];
    }
    cout << tot << endl;
    return 0;
}  