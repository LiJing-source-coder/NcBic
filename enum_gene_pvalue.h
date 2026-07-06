#ifndef ENUM_GENE_PVALUE
#define ENUM_GENE_PVALUE
#include<boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
using namespace boost::multiprecision;
using namespace boost::math;
cpp_dec_float_100 compute_p(const int *conditionL,int m,double tolerance,int &cons);
cpp_dec_float_100 pbinom_at_least(int k, int n, cpp_dec_float_100 p);
int compute_confidence(const int *conditionL,int m,double repvalue,int &cons);
#endif
