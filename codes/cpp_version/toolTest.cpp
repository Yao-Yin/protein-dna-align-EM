#include <iostream>
#include <vector>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <numeric>
//#include "DataTool.h"

using namespace std;
//DataTool dt;

template <typename Iter>
typename std::iterator_traits<Iter>::value_type
log_sum_exp(Iter begin, Iter end)
{
  using VT = typename std::iterator_traits<Iter>::value_type;
  if (begin==end) return VT{};
  using std::exp;
  using std::log;
  auto max_elem = *std::max_element(begin, end);
  auto sum = std::accumulate(begin, end, VT{}, 
     [max_elem](VT a, VT b) { return a + exp(b - max_elem); });
  return max_elem + log(sum);
}

struct num {
    double fp; // float part 
    double ep; // exponential part
    bool zf; // zeroflag
    num(double x) {
        fp = x;
        ep = 0;
        autofix();
    }
    num() {
        fp = 0.0;
        ep = 0;
    }
    void autofix() {
        double exx = floorf(log10f(fp));
        ep += exx;
        fp /= (double)pow(10, exx);
    }
    num operator+ (const num & other) const {
        if(other.ep - ep >= 30) return other;
        if(ep - other.ep >= 30) return *this;
        num res = other;
        res.fp *= pow(10, res.ep - ep);
        res.fp += fp;
        res.ep = ep;
        res.autofix();
        return res;
    }
    num operator* (const num & other) const {
        num res;
        res.fp = fp*other.fp;
        res.ep = ep + other.ep;
        res.autofix();
        return res;
    }
    void operator+= (const num & other) {
        if(ep - other.ep >= 100) return;
        if(other.ep - ep >= 100) {
            ep = other.ep;
            fp = other.fp;
            return;
        }
        fp *= pow(10, other.ep - ep);
        fp += other.fp;
        ep = other.ep;
        autofix();
    }

    void operator*= (const num & other) {
        fp *= other.fp;
        ep += other.ep;
        autofix();
    }
};

int main() {
    int n;
    cin >> n;
    vector<double> testdata(n);
    vector<double> testlog(n);
    vector<num> testmult(n);
    for(int i = 0; i < n; i ++) testdata[i] = (double)rand() / INT32_MAX;
    for(int i = 0; i < n; i ++) {
        testlog[i] = log(testdata[i]);
        testmult[i] = num(testdata[i]);
    }
    time_t curr = clock();
    double rres = 0;
    
    rres += log_sum_exp(testlog.begin(),testlog.end());
    cout << clock() - curr << endl;
    num res = testmult[0];
    curr = clock();
    for(int i = 1; i < n; i ++) {
        res += testmult[i];
    }
    
    cout<<clock()-curr<<" "<< res.fp << " 10^" << res.ep  << endl;
    double t = 0;
    for(int i = 0; i < n; i ++) t += testdata[i];
    cout << t << endl;
    /*num t(0.001);
    num k(0.006);
    cout << (t+k).fp << (t+k).ep << " " << (t+k).fp * exp((t+k).ep) << endl; */
    return 0;
}