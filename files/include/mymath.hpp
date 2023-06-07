#pragma once
#ifndef LO_MYMATH_HPP
#define LO_MYMATH_HPP

#include "constants.hpp"
#include "gauss_legendre/gauss_legendre.hpp"
#include <cmath>
#include <complex>
#include <vector>

namespace util
{

inline double delta(int m, int n) { return (m == n) ? 1.0 : 0.0; }

inline double sign_function(int m) { return (m % 2 == 0) ? 1.0 : -1.0; }

const std::complex<double> im_unit(0.0, 1.0); // unit imaginary number i.

std::vector<double> make_table(double xstart, double xend, int inter)
{
    std::vector<double> temp;
    for (int i = 1; i <= inter; i = i + 1)
    {
        double tempp = xstart + i * (xend - xstart) / inter;
        temp.push_back(tempp);
    }
    return temp;
}

} // namespace util

#endif // LO_MYMATH_HPP