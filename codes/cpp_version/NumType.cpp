#include "NumType.h"
#include <math.h>
#include <numeric>
#include <iostream>

template <typename Iter>
typename std::iterator_traits<Iter>::value_type
_log_sum_exp(Iter begin, Iter end)
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
log_sum_exp(Iter begin, Iter end)
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

LogNumType LogSumExp(const LogNumType & a, const LogNumType & b) {
    std::vector<LogNumType> curr {a, b};
    //std::cout << a << " & " << b << std::endl;
    return log_sum_exp(curr.begin(), curr.end());
}

LogNumType LogSumExp(const LogNumType & a, const LogNumType & b, const LogNumType & c, const LogNumType & d, const LogNumType & e,const LogNumType & f, const LogNumType & g) {
    std::vector<LogNumType> curr {a, b, c, d, e, f, g};
    return log_sum_exp(curr.begin(), curr.end());
}

LogNumType Log1m(const LogNumType & x) {
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
/*
LogNumType LogSumExp(const LogNumType & a, const LogNumType & b) {
    LogNumType max_elem = a > b ? a : b;
    return max_elem == log(0.0) ? max_elem : max_elem + exp(a-max_elem) + exp(b-max_elem);  
}

LogNumType LogSumExp(const LogNumType & a, const LogNumType & b, const LogNumType & c, const LogNumType & d, const LogNumType & e,const LogNumType & f, const LogNumType & g) {
    using VT = LogNumType;
    LogNumType curr[7] {a, b, c, d, e, f, g};
    LogNumType max_elem = *std::max_element(curr, curr + 7);
    if (max_elem == log(0.0)) return log(0.0);
    auto sum = std::accumulate(curr, curr + 7, VT{}, 
     [max_elem](VT a, VT b) { return a + exp(b - max_elem); });
    return max_elem + log(sum);
    return log(0);
}

LogNumType Log1m(const LogNumType & x) {
    // to calculate log(1-exp(x))
    /*if (x < -log(2.0)) {
      return log1p(-exp(x));
    }
    else return log(-expm1(x));
    return log(1.0-exp(x));
}

LogNumType LogSumExp(const LogNumType & a, const LogNumType & b, const LogNumType & c) {
    /*using VT = LogNumType;
    VT curr[3] {a, b, c};
    VT max_elem = *std::max_element(curr, curr + 3);
    if (max_elem == log(0.0)) return log(0.0);
    auto sum = std::accumulate(curr, curr + 3, VT{}, 
     [max_elem](VT a, VT b) { return a + exp(b - max_elem); });
    return max_elem + log(sum);
    return log(0);
}

LogNumType LogSumExp(const LogNumType & a, const LogNumType & b, const LogNumType & c, const LogNumType & d) {
   /*using VT = LogNumType;
    VT curr[4] {a, b, c, d};
    VT max_elem = *std::max_element(curr, curr + 4);
    if (max_elem == log(0.0)) return log(0.0);
    auto sum = std::accumulate(curr, curr + 4, VT{}, 
     [max_elem](VT a, VT b) { return a + exp(b - max_elem); });
    return max_elem + log(sum);
    return log(0);
}
*/