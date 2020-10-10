#ifndef NUMTYPE_H
#define NUMTYPE_H
#include <numeric>
#include <math.h>
#include <algorithm>
#include <vector>

typedef double NumType; // To represent Probability numerically;
typedef double LogNumType;
template <typename Iter>
typename std::iterator_traits<Iter>::value_type
log_sum_exp(Iter begin, Iter end);
inline LogNumType LogSumExp(const LogNumType & a, const LogNumType & b);
inline LogNumType LogSumExp(const LogNumType & a, const LogNumType & b, const LogNumType & c);
inline LogNumType LogSumExp(const LogNumType & a, const LogNumType & b, const LogNumType & c, const LogNumType & d);
inline LogNumType LogSumExp(const LogNumType & a, const LogNumType & b, const LogNumType & c, const LogNumType & d, const LogNumType & e,const LogNumType & f, const LogNumType & g);
inline LogNumType Log1m(const LogNumType & a);
/*struct NumType {
NumType(double num);
NumType();
double fp;
int ep;
bool zf; // zeroflag
void rescaling();
void operator=(const NumType & other);
NumType operator*(const NumType & other);
NumType operator+(const NumType & other);
};"*/
#endif