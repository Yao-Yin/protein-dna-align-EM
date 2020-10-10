#include "NumType.h"


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

inline LogNumType LogSumExp(const LogNumType & a, const LogNumType & b) {
    std::vector<LogNumType> curr {a, b};
    return log_sum_exp(curr.begin(), curr.end());
}

inline LogNumType LogSumExp(const LogNumType & a, const LogNumType & b, const LogNumType & c, const LogNumType & d, const LogNumType & e,const LogNumType & f, const LogNumType & g) {
    std::vector<LogNumType> curr {a, b, c, d, e, f, g};
    return log_sum_exp(curr.begin(), curr.end());
}

inline LogNumType Log1m(const LogNumType & x) {
    // to calculate log(1-exp(x))
    if (x < -log(2.0)) {
      return log1p(-exp(x));
    }
    else return log(-expm1(x));
}

inline LogNumType LogSumExp(const LogNumType & a, const LogNumType & b, const LogNumType & c) {
    std::vector<LogNumType> curr {a, b, c};
    return log_sum_exp(curr.begin(), curr.end());
}

inline LogNumType LogSumExp(const LogNumType & a, const LogNumType & b, const LogNumType & c, const LogNumType & d) {
    std::vector<LogNumType> curr {a, b, c, d};
    return log_sum_exp(curr.begin(), curr.end());
}