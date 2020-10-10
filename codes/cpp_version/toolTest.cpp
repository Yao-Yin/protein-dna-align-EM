#include <iostream>
#include <math.h>
#include <numeric>
#include <vector>
#include <algorithm>
using namespace std;

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

int main() {
    int t;
    cin >> t;
    double r = log(0);
    vector<double> test;
    while (t --) {
        double curr;
        cin >> curr;
        test.push_back(log(curr));
    }
    cout << log_sum_exp(test.begin(), test.end()) << endl;
    return 0;
}