#include "NumType.h"

/*
template <typename Iter>
typename std::iterator_traits<Iter>::value_type
log_sum_exp(Iter begin, Iter end)
{
  using VT = typename std::iterator_traits<Iter>::value_type;
  if (begin==end) return log(0.0);
  //using std::exp;
  //using std::log;
  auto max_elem = *std::max_element(begin, end);
  auto sum = std::accumulate(begin, end, VT{}, 
     [max_elem](VT a, VT b) { return a + exp(b - max_elem); });
  return max_elem + log(sum);
}

template <typename Iter>
typename std::iterator_traits<Iter>::value_type
_log_sum_exp(Iter begin, Iter end)
{
  using VT = typename std::iterator_traits<Iter>::value_type;
  if (begin==end) return log(0.0);
  //using std::exp;
  //using std::log;
  std::vector<VT> curr;
  for(auto iter = begin; iter != end; iter ++) {
      if(*iter != log(0.0)) curr.push_back(*iter);
  }
  if(curr.empty()) return log(0.0);
  return _log_sum_exp(curr.begin(), curr.end());
}

/*LogNumType LogSumExp(const LogNumType & a, const LogNumType & b) {
    std::vector<LogNumType> curr {a, b};
    //std::cout << a << " & " << b << std::endl;
    return log_sum_exp(curr.begin(), curr.end());
}

LogNumType LogSumExp(const LogNumType & a, const LogNumType & b, const LogNumType & c, const LogNumType & d, const LogNumType & e,const LogNumType & f, const LogNumType & g) {
    std::vector<LogNumType> curr {a, b, c, d, e, f, g};
    return log_sum_exp(curr.begin(), curr.end());
}

/*LogNumType Log1m(const LogNumType & x) {
    // to calculate log(1-exp(x))
    if (x < -log(2.0)) {
      return log1p(-exp(x));
    }
    else return log(-expm1(x));
}

LogNumType LogSumExp(const LogNumType & a, const LogNumType & b, const LogNumType & c) {
    std::vector<LogNumType> curr {a, b, c};
    return log_sum_exp(curr.begin(), curr.end());
}

LogNumType LogSumExp(const LogNumType & a, const LogNumType & b, const LogNumType & c, const LogNumType & d) {
    std::vector<LogNumType> curr {a, b, c, d};
    return log_sum_exp(curr.begin(), curr.end());
}

*/