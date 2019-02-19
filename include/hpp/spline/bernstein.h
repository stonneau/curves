/**
* \file bezier_curve.h
* \brief class allowing to create a Bezier curve of dimension 1 <= n <= 3.
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*/


#ifndef _CLASS_BERNSTEIN
#define _CLASS_BERNSTEIN

#include "curve_abc.h"

#include "MathDefs.h"

#include <math.h>
#include <vector>
#include <stdexcept>

namespace spline
{
///
/// \brief Computes a binomal coefficient
///
inline unsigned int bin(const unsigned  int n, const unsigned  int k)
{
    if(k >  n)
        throw std::runtime_error("binomial coefficient higher than degree");
    if(k == 0)
        return 1;
    if(k > n/2)
        return bin(n,n-k);
    return n * bin(n-1,k-1) / k;
}


/// \class Bernstein
/// \brief Computes a Bernstein polynome
///
template <typename Numeric = double>
struct Bern{
Bern(const long unsigned int m, const long unsigned int i)
    :m_minus_i(m - i)
    ,i_(i)
    ,bin_m_i_(bin(m,i)) {}

~Bern(){}

Numeric operator()(const Numeric u) const
{
    assert(u >= 0. && u <= 1.);
    return bin_m_i_*(pow(u, i_)) *pow((1-u),m_minus_i);
}

Numeric m_minus_i;
Numeric i_;
Numeric bin_m_i_;
};


///
/// \brief Computes all Bernstein polynomes for a certain degree
///
template <typename Numeric>
std::vector<Bern<Numeric> > makeBernstein(const long unsigned int n)
{
    std::vector<Bern<Numeric> > res;
    for(long unsigned int i = 0; i<= n; ++i)
        res.push_back(Bern<Numeric>(n, i));
    return res;
}
} // namespace spline
#endif //_CLASS_BERNSTEIN

