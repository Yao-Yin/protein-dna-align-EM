#ifndef NUMTYPE_H
#define NUMTYPE_H
#include <numeric>
#include <math.h>
#include <algorithm>
#include <vector>
#include <iostream>

typedef long double NumType; // To represent Probability numerically;
typedef long double LogNumType;
/*template <typename Iter>
typename std::iterator_traits<Iter>::value_type
_log_sum_exp(Iter begin, Iter end);
template <typename Iter>
typename std::iterator_traits<Iter>::value_type
log_sum_exp(Iter begin, Iter end);*/
LogNumType log_sum_exp(std::vector<LogNumType>::iterator begin, std::vector<LogNumType>::iterator end);
LogNumType LogSumExp(const LogNumType & a, const LogNumType & b);
LogNumType LogSumExp(const LogNumType & a, const LogNumType & b, const LogNumType & c);
LogNumType LogSumExp(const LogNumType & a, const LogNumType & b, const LogNumType & c, const LogNumType & d);
LogNumType LogSumExp(const LogNumType & a, const LogNumType & b, const LogNumType & c, const LogNumType & d, const LogNumType & e,const LogNumType & f, const LogNumType & g);
LogNumType Log1m(const LogNumType & a);

inline LogNumType log_sum_exp(std::vector<LogNumType>::iterator begin, std::vector<LogNumType>::iterator end)
{
    using VT = LogNumType;
    if (begin==end) return log(0.0);
    //using std::exp;
    //using std::log;
    auto max_elem = *std::max_element(begin, end);
    if(max_elem == log(0.0)) return log(0.0);
    auto sum = std::accumulate(begin, end, VT{}, 
     [max_elem](VT a, VT b) { return a + exp(b - max_elem); });
    return max_elem + log(sum);
}

inline LogNumType LogSumExp(const LogNumType & a, const LogNumType & b) {
    LogNumType max_elem = a > b ? a : b;
    return max_elem == log(0.0) ? max_elem : max_elem + log(exp(a-max_elem) + exp(b-max_elem));  
}

inline LogNumType LogSumExp(const LogNumType & a, const LogNumType & b, const LogNumType & c) {
    using VT = LogNumType;
    VT max_elem = std::max(a, std::max(b, c));
    if (max_elem == log(0.0)) return log(0.0);
    return max_elem + log(exp(a-max_elem)+exp(b-max_elem)+exp(c-max_elem));
}

inline LogNumType LogSumExp(const LogNumType & a, const LogNumType & b, const LogNumType & c, const LogNumType & d) {
    using VT = LogNumType;
    VT max_elem = std::max(a, std::max(b, std::max(c, d)));
    if (max_elem == log(0.0)) return log(0.0);
    return max_elem + log(exp(a-max_elem)+exp(b-max_elem)+exp(c-max_elem)+exp(d-max_elem));
}

inline LogNumType LogSumExp(const LogNumType & a, const LogNumType & b, const LogNumType & c, const LogNumType & d, const LogNumType & e,const LogNumType & f, const LogNumType & g) {
    using VT = LogNumType;
    VT max_elem = std::max(a, std::max(b, std::max(c, std::max(d, std::max(e, std::max(f, g))))));
    if (max_elem == log(0.0)) return log(0.0);
    return max_elem + log(exp(a-max_elem)+exp(b-max_elem)+exp(c-max_elem)+exp(d-max_elem)+exp(e-max_elem)+exp(f-max_elem)+exp(g-max_elem));
}

inline LogNumType Log1m(const LogNumType & x) {
    return log(1.0-exp(x));
    /*if (x < -log(2.0)) {
      return log1p(-exp(x));
    }
    else return log(-expm1(x));*/
}

#endif