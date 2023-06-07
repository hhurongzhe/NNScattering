#pragma once
#ifndef INTERACTION_aPWD_HPP
#define INTERACTION_aPWD_HPP

#include "lib_define.hpp"
#include "mymath.hpp"
#include "state.hpp"

namespace interaction_aPWD
{

using basic_math::quick_pow;
using std::string;

// Error message written and all processes stop
void error_message_print_abort(const string &error_message) { std::cout << error_message << std::endl; }

constexpr double ppi = 3.14159265358979323846;

// partial-wave projection by aPWD method under SYM-LSJ basis convention.
double V_kernal_SYM(const states::LSJ_State &s1, const states::LSJ_State &s2, const double &f1, const double &f2,
                    const double &f3, const double &f4, const double &f5, const double &f6, const double &x)
{
    const double pmag = s1.momentum;
    const double ppmag = s2.momentum;

    if (!(s1.j == s2.j && s1.s == s2.s && s1.j >= 0 && s2.j >= 0 && s1.l >= 0 && s2.l >= 0 && s1.s >= 0 && s2.s >= 0 &&
          std::abs(s1.l - s2.l) <= 2))
    {
        return 0;
    }
    else if (s2.l == 0 && s1.l == 0 && s1.s == 0 && s1.j == 0)
    {
        return 2 * ppi *
               (f1 - 3 * f2 - quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 - quick_pow(pmag, 2) * f5 -
                quick_pow(ppmag, 2) * f5 - 2 * pmag * ppmag * x * f5 - quick_pow(pmag, 2) * f6 -
                quick_pow(ppmag, 2) * f6 + 2 * pmag * ppmag * x * f6);
    }
    else if (s2.l == 1 && s1.l == 1 && s1.s == 1 && s1.j == 0)
    {
        return -2 * ppi *
               (-2 * pmag * ppmag * quick_pow(x, 2) * f3 +
                quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 3) * f4 + 2 * pmag * ppmag * (f3 + f5 - f6) +
                x * (-f1 - f2 - quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 + quick_pow(pmag, 2) * f5 +
                     quick_pow(ppmag, 2) * f5 + quick_pow(pmag, 2) * f6 + quick_pow(ppmag, 2) * f6));
    }
    else if (s2.l == 1 && s1.l == 1 && s1.s == 0 && s1.j == 1)
    {
        return 2 * ppi * x *
               (f1 - 3 * f2 - quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 - quick_pow(pmag, 2) * f5 -
                quick_pow(ppmag, 2) * f5 - 2 * pmag * ppmag * x * f5 - quick_pow(pmag, 2) * f6 -
                quick_pow(ppmag, 2) * f6 + 2 * pmag * ppmag * x * f6);
    }
    else if (s2.l == 1 && s1.l == 1 && s1.s == 1 && s1.j == 1)
    {
        return 2 * ppi *
               (pmag * ppmag * quick_pow(x, 2) * (f3 + f5 - f6) - pmag * ppmag * (f3 - f5 + f6) +
                x * (f1 + f2 + (quick_pow(pmag, 2) + quick_pow(ppmag, 2)) * (f5 + f6)));
    }
    else if (s2.l == 0 && s1.l == 0 && s1.s == 1 && s1.j == 1)
    {
        return (2 * ppi *
                (3 * f1 + 3 * f2 + quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                 quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 + quick_pow(pmag, 2) * f5 +
                 quick_pow(ppmag, 2) * f5 + 2 * pmag * ppmag * x * f5 + quick_pow(pmag, 2) * f6 +
                 quick_pow(ppmag, 2) * f6 - 2 * pmag * ppmag * x * f6)) /
               3.;
    }
    else if (s2.l == 0 && s1.l == 2 && s1.s == 1 && s1.j == 1)
    {
        return -(2 * sqrt(2) * ppi *
                 (4 * pmag * ppmag * x * (f5 - f6) + quick_pow(ppmag, 2) * (-1 + 3 * quick_pow(x, 2)) * (f5 + f6) +
                  quick_pow(pmag, 2) * (quick_pow(ppmag, 2) * (-1 + quick_pow(x, 2)) * f4 + 2 * (f5 + f6)))) /
               3.;
    }
    else if (s2.l == 2 && s1.l == 0 && s1.s == 1 && s1.j == 1)
    {
        return -(2 * sqrt(2) * ppi *
                 (4 * pmag * ppmag * x * (f5 - f6) + 2 * quick_pow(ppmag, 2) * (f5 + f6) +
                  quick_pow(pmag, 2) *
                      (quick_pow(ppmag, 2) * (-1 + quick_pow(x, 2)) * f4 + (-1 + 3 * quick_pow(x, 2)) * (f5 + f6)))) /
               3.;
    }
    else if (s2.l == 2 && s1.l == 2 && s1.s == 1 && s1.j == 1)
    {
        return (ppi *
                ((-3 + 9 * quick_pow(x, 2)) * f1 + (-3 + 9 * quick_pow(x, 2)) * f2 - 18 * pmag * ppmag * x * f3 +
                 18 * pmag * ppmag * quick_pow(x, 3) * f3 - 5 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                 14 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 -
                 9 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 4) * f4 + quick_pow(pmag, 2) * f5 +
                 quick_pow(ppmag, 2) * f5 - 4 * pmag * ppmag * x * f5 - 3 * quick_pow(pmag, 2) * quick_pow(x, 2) * f5 -
                 3 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f5 + quick_pow(pmag, 2) * f6 + quick_pow(ppmag, 2) * f6 +
                 4 * pmag * ppmag * x * f6 - 3 * quick_pow(pmag, 2) * quick_pow(x, 2) * f6 -
                 3 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f6)) /
               3.;
    }
    else if (s2.l == 2 && s1.l == 2 && s1.s == 0 && s1.j == 2)
    {
        return ppi * (-1 + 3 * quick_pow(x, 2)) *
               (f1 - 3 * f2 - quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 - quick_pow(pmag, 2) * f5 -
                quick_pow(ppmag, 2) * f5 - 2 * pmag * ppmag * x * f5 - quick_pow(pmag, 2) * f6 -
                quick_pow(ppmag, 2) * f6 + 2 * pmag * ppmag * x * f6);
    }
    else if (s2.l == 2 && s1.l == 2 && s1.s == 1 && s1.j == 2)
    {
        return ppi *
               ((-1 + 3 * quick_pow(x, 2)) * f1 + (-1 + 3 * quick_pow(x, 2)) * f2 - 2 * pmag * ppmag * x * f3 +
                2 * pmag * ppmag * quick_pow(x, 3) * f3 + quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                2 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 +
                quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 4) * f4 - quick_pow(pmag, 2) * f5 -
                quick_pow(ppmag, 2) * f5 + 3 * quick_pow(pmag, 2) * quick_pow(x, 2) * f5 +
                3 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f5 + 4 * pmag * ppmag * quick_pow(x, 3) * f5 -
                quick_pow(pmag, 2) * f6 - quick_pow(ppmag, 2) * f6 + 3 * quick_pow(pmag, 2) * quick_pow(x, 2) * f6 +
                3 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f6 - 4 * pmag * ppmag * quick_pow(x, 3) * f6);
    }
    else if (s2.l == 1 && s1.l == 1 && s1.s == 1 && s1.j == 2)
    {
        return (-2 * ppi *
                (2 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 3) * f4 +
                 pmag * ppmag * (-5 * f3 + f5 - f6) + pmag * ppmag * quick_pow(x, 2) * (5 * f3 - 3 * f5 + 3 * f6) -
                 x * (5 * f1 + 5 * f2 + 2 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 + quick_pow(pmag, 2) * f5 +
                      quick_pow(ppmag, 2) * f5 + quick_pow(pmag, 2) * f6 + quick_pow(ppmag, 2) * f6))) /
               5.;
    }
    else if (s2.l == 1 && s1.l == 3 && s1.s == 1 && s1.j == 2)
    {
        return -(2 * sqrt(6) * ppi *
                 (2 * pmag * ppmag * (-1 + 3 * quick_pow(x, 2)) * (f5 - f6) +
                  quick_pow(ppmag, 2) * x * (-3 + 5 * quick_pow(x, 2)) * (f5 + f6) +
                  quick_pow(pmag, 2) * x * (quick_pow(ppmag, 2) * (-1 + quick_pow(x, 2)) * f4 + 2 * (f5 + f6)))) /
               5.;
    }
    else if (s2.l == 3 && s1.l == 1 && s1.s == 1 && s1.j == 2)
    {
        return -(2 * sqrt(6) * ppi *
                 (2 * pmag * ppmag * (-1 + 3 * quick_pow(x, 2)) * (f5 - f6) + 2 * quick_pow(ppmag, 2) * x * (f5 + f6) +
                  quick_pow(pmag, 2) * x *
                      (quick_pow(ppmag, 2) * (-1 + quick_pow(x, 2)) * f4 + (-3 + 5 * quick_pow(x, 2)) * (f5 + f6)))) /
               5.;
    }
    else if (s2.l == 3 && s1.l == 3 && s1.s == 1 && s1.j == 2)
    {
        return -0.2 *
               (ppi *
                (-50 * pmag * ppmag * quick_pow(x, 4) * f3 +
                 25 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 5) * f4 -
                 2 * pmag * ppmag * (5 * f3 + f5 - f6) + 6 * pmag * ppmag * quick_pow(x, 2) * (10 * f3 + f5 - f6) +
                 x * (15 * f1 + 15 * f2 + 19 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                      3 * quick_pow(pmag, 2) * f5 - 3 * quick_pow(ppmag, 2) * f5 - 3 * quick_pow(pmag, 2) * f6 -
                      3 * quick_pow(ppmag, 2) * f6) +
                 quick_pow(x, 3) * (-25 * f1 - 25 * f2 - 44 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                                    5 * quick_pow(pmag, 2) * f5 + 5 * quick_pow(ppmag, 2) * f5 +
                                    5 * quick_pow(pmag, 2) * f6 + 5 * quick_pow(ppmag, 2) * f6)));
    }
    else if (s2.l == 3 && s1.l == 3 && s1.s == 0 && s1.j == 3)
    {
        return ppi * x * (-3 + 5 * quick_pow(x, 2)) *
               (f1 - 3 * f2 - quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 - quick_pow(pmag, 2) * f5 -
                quick_pow(ppmag, 2) * f5 - 2 * pmag * ppmag * x * f5 - quick_pow(pmag, 2) * f6 -
                quick_pow(ppmag, 2) * f6 + 2 * pmag * ppmag * x * f6);
    }
    else if (s2.l == 3 && s1.l == 3 && s1.s == 1 && s1.j == 3)
    {
        return (ppi * (5 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 5) * f4 +
                       5 * pmag * ppmag * quick_pow(x, 4) * (f3 + 3 * f5 - 3 * f6) -
                       6 * pmag * ppmag * quick_pow(x, 2) * (f3 + f5 - f6) + pmag * ppmag * (f3 - f5 + f6) +
                       10 * quick_pow(x, 3) *
                           (f1 + f2 - quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 + quick_pow(pmag, 2) * f5 +
                            quick_pow(ppmag, 2) * f5 + quick_pow(pmag, 2) * f6 + quick_pow(ppmag, 2) * f6) -
                       x * (6 * f1 + 6 * f2 - 5 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            6 * quick_pow(pmag, 2) * f5 + 6 * quick_pow(ppmag, 2) * f5 + 6 * quick_pow(pmag, 2) * f6 +
                            6 * quick_pow(ppmag, 2) * f6))) /
               2.;
    }
    else if (s2.l == 2 && s1.l == 2 && s1.s == 1 && s1.j == 3)
    {
        return -0.14285714285714285 *
               (ppi * ((7 - 21 * quick_pow(x, 2)) * f1 + (7 - 21 * quick_pow(x, 2)) * f2 - 28 * pmag * ppmag * x * f3 +
                       28 * pmag * ppmag * quick_pow(x, 3) * f3 + 5 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                       16 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 +
                       11 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 4) * f4 + quick_pow(pmag, 2) * f5 +
                       quick_pow(ppmag, 2) * f5 + 6 * pmag * ppmag * x * f5 -
                       3 * quick_pow(pmag, 2) * quick_pow(x, 2) * f5 - 3 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f5 -
                       10 * pmag * ppmag * quick_pow(x, 3) * f5 + quick_pow(pmag, 2) * f6 + quick_pow(ppmag, 2) * f6 -
                       6 * pmag * ppmag * x * f6 - 3 * quick_pow(pmag, 2) * quick_pow(x, 2) * f6 -
                       3 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f6 + 10 * pmag * ppmag * quick_pow(x, 3) * f6));
    }
    else if (s2.l == 2 && s1.l == 4 && s1.s == 1 && s1.j == 3)
    {
        return -(sqrt(3) * ppi *
                 (8 * pmag * ppmag * x * (-3 + 5 * quick_pow(x, 2)) * (f5 - f6) +
                  quick_pow(ppmag, 2) * (3 - 30 * quick_pow(x, 2) + 35 * quick_pow(x, 4)) * (f5 + f6) +
                  quick_pow(pmag, 2) * (quick_pow(ppmag, 2) * (1 - 6 * quick_pow(x, 2) + 5 * quick_pow(x, 4)) * f4 +
                                        4 * (-1 + 3 * quick_pow(x, 2)) * (f5 + f6)))) /
               7.;
    }
    else if (s2.l == 4 && s1.l == 2 && s1.s == 1 && s1.j == 3)
    {
        return -(sqrt(3) * ppi *
                 (8 * pmag * ppmag * x * (-3 + 5 * quick_pow(x, 2)) * (f5 - f6) +
                  4 * quick_pow(ppmag, 2) * (-1 + 3 * quick_pow(x, 2)) * (f5 + f6) +
                  quick_pow(pmag, 2) * (quick_pow(ppmag, 2) * (1 - 6 * quick_pow(x, 2) + 5 * quick_pow(x, 4)) * f4 +
                                        (3 - 30 * quick_pow(x, 2) + 35 * quick_pow(x, 4)) * (f5 + f6)))) /
               7.;
    }
    else if (s2.l == 4 && s1.l == 4 && s1.s == 1 && s1.j == 3)
    {
        return (ppi *
                (7 * (3 - 30 * quick_pow(x, 2) + 35 * quick_pow(x, 4)) * f1 +
                 7 * (3 - 30 * quick_pow(x, 2) + 35 * quick_pow(x, 4)) * f2 + 210 * pmag * ppmag * x * f3 -
                 700 * pmag * ppmag * quick_pow(x, 3) * f3 + 490 * pmag * ppmag * quick_pow(x, 5) * f3 +
                 27 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                 267 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 +
                 485 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 4) * f4 -
                 245 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 6) * f4 - 3 * quick_pow(pmag, 2) * f5 -
                 3 * quick_pow(ppmag, 2) * f5 + 24 * pmag * ppmag * x * f5 +
                 30 * quick_pow(pmag, 2) * quick_pow(x, 2) * f5 + 30 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f5 -
                 40 * pmag * ppmag * quick_pow(x, 3) * f5 - 35 * quick_pow(pmag, 2) * quick_pow(x, 4) * f5 -
                 35 * quick_pow(ppmag, 2) * quick_pow(x, 4) * f5 - 3 * quick_pow(pmag, 2) * f6 -
                 3 * quick_pow(ppmag, 2) * f6 - 24 * pmag * ppmag * x * f6 +
                 30 * quick_pow(pmag, 2) * quick_pow(x, 2) * f6 + 30 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f6 +
                 40 * pmag * ppmag * quick_pow(x, 3) * f6 - 35 * quick_pow(pmag, 2) * quick_pow(x, 4) * f6 -
                 35 * quick_pow(ppmag, 2) * quick_pow(x, 4) * f6)) /
               28.;
    }
    else if (s2.l == 4 && s1.l == 4 && s1.s == 0 && s1.j == 4)
    {
        return (ppi * (3 - 30 * quick_pow(x, 2) + 35 * quick_pow(x, 4)) *
                (f1 - 3 * f2 - quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                 quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 - quick_pow(pmag, 2) * f5 -
                 quick_pow(ppmag, 2) * f5 - 2 * pmag * ppmag * x * f5 - quick_pow(pmag, 2) * f6 -
                 quick_pow(ppmag, 2) * f6 + 2 * pmag * ppmag * x * f6)) /
               4.;
    }
    else if (s2.l == 4 && s1.l == 4 && s1.s == 1 && s1.j == 4)
    {
        return (ppi *
                ((3 - 30 * quick_pow(x, 2) + 35 * quick_pow(x, 4)) * f1 +
                 (3 - 30 * quick_pow(x, 2) + 35 * quick_pow(x, 4)) * f2 + 6 * pmag * ppmag * x * f3 -
                 20 * pmag * ppmag * quick_pow(x, 3) * f3 + 14 * pmag * ppmag * quick_pow(x, 5) * f3 -
                 3 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                 27 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 -
                 45 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 4) * f4 +
                 21 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 6) * f4 + 3 * quick_pow(pmag, 2) * f5 +
                 3 * quick_pow(ppmag, 2) * f5 - 30 * quick_pow(pmag, 2) * quick_pow(x, 2) * f5 -
                 30 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f5 - 40 * pmag * ppmag * quick_pow(x, 3) * f5 +
                 35 * quick_pow(pmag, 2) * quick_pow(x, 4) * f5 + 35 * quick_pow(ppmag, 2) * quick_pow(x, 4) * f5 +
                 56 * pmag * ppmag * quick_pow(x, 5) * f5 + 3 * quick_pow(pmag, 2) * f6 + 3 * quick_pow(ppmag, 2) * f6 -
                 30 * quick_pow(pmag, 2) * quick_pow(x, 2) * f6 - 30 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f6 +
                 40 * pmag * ppmag * quick_pow(x, 3) * f6 + 35 * quick_pow(pmag, 2) * quick_pow(x, 4) * f6 +
                 35 * quick_pow(ppmag, 2) * quick_pow(x, 4) * f6 - 56 * pmag * ppmag * quick_pow(x, 5) * f6)) /
               4.;
    }
    else if (s2.l == 3 && s1.l == 3 && s1.s == 1 && s1.j == 4)
    {
        return -0.05555555555555555 *
               (ppi * (55 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 5) * f4 +
                       3 * pmag * ppmag * (9 * f3 - f5 + f6) -
                       6 * pmag * ppmag * quick_pow(x, 2) * (27 * f3 - 5 * f5 + 5 * f6) +
                       5 * pmag * ppmag * quick_pow(x, 4) * (27 * f3 - 7 * f5 + 7 * f6) -
                       2 * quick_pow(x, 3) *
                           (45 * f1 + 45 * f2 + 47 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            5 * quick_pow(pmag, 2) * f5 + 5 * quick_pow(ppmag, 2) * f5 + 5 * quick_pow(pmag, 2) * f6 +
                            5 * quick_pow(ppmag, 2) * f6) +
                       x * (54 * f1 + 54 * f2 + 39 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            6 * quick_pow(pmag, 2) * f5 + 6 * quick_pow(ppmag, 2) * f5 + 6 * quick_pow(pmag, 2) * f6 +
                            6 * quick_pow(ppmag, 2) * f6)));
    }
    else if (s2.l == 3 && s1.l == 5 && s1.s == 1 && s1.j == 4)
    {
        return -(sqrt(5) * ppi *
                 (2 * pmag * ppmag * (3 - 30 * quick_pow(x, 2) + 35 * quick_pow(x, 4)) * (f5 - f6) +
                  quick_pow(ppmag, 2) * x * (15 - 70 * quick_pow(x, 2) + 63 * quick_pow(x, 4)) * (f5 + f6) +
                  quick_pow(pmag, 2) * x *
                      (quick_pow(ppmag, 2) * (3 - 10 * quick_pow(x, 2) + 7 * quick_pow(x, 4)) * f4 +
                       4 * (-3 + 5 * quick_pow(x, 2)) * (f5 + f6)))) /
               9.;
    }
    else if (s2.l == 5 && s1.l == 3 && s1.s == 1 && s1.j == 4)
    {
        return -(sqrt(5) * ppi *
                 (2 * pmag * ppmag * (3 - 30 * quick_pow(x, 2) + 35 * quick_pow(x, 4)) * (f5 - f6) +
                  4 * quick_pow(ppmag, 2) * x * (-3 + 5 * quick_pow(x, 2)) * (f5 + f6) +
                  quick_pow(pmag, 2) * x *
                      (quick_pow(ppmag, 2) * (3 - 10 * quick_pow(x, 2) + 7 * quick_pow(x, 4)) * f4 +
                       (15 - 70 * quick_pow(x, 2) + 63 * quick_pow(x, 4)) * (f5 + f6)))) /
               9.;
    }
    else if (s2.l == 5 && s1.l == 5 && s1.s == 1 && s1.j == 4)
    {
        return -0.027777777777777776 *
               (ppi *
                (-1134 * pmag * ppmag * quick_pow(x, 6) * f3 +
                 567 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 7) * f4 -
                 30 * pmag * ppmag * quick_pow(x, 2) * (27 * f3 + 2 * f5 - 2 * f6) +
                 6 * pmag * ppmag * (9 * f3 + f5 - f6) + 70 * pmag * ppmag * quick_pow(x, 4) * (27 * f3 + f5 - f6) +
                 quick_pow(x, 3) * (630 * f1 + 630 * f2 + 845 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                                    70 * quick_pow(pmag, 2) * f5 - 70 * quick_pow(ppmag, 2) * f5 -
                                    70 * quick_pow(pmag, 2) * f6 - 70 * quick_pow(ppmag, 2) * f6) -
                 7 * quick_pow(x, 5) *
                     (81 * f1 + 81 * f2 + 179 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                      9 * quick_pow(pmag, 2) * f5 - 9 * quick_pow(ppmag, 2) * f5 - 9 * quick_pow(pmag, 2) * f6 -
                      9 * quick_pow(ppmag, 2) * f6) -
                 3 * x *
                     (45 * f1 + 45 * f2 + 53 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                      5 * quick_pow(pmag, 2) * f5 - 5 * quick_pow(ppmag, 2) * f5 - 5 * quick_pow(pmag, 2) * f6 -
                      5 * quick_pow(ppmag, 2) * f6)));
    }
    else if (s2.l == 5 && s1.l == 5 && s1.s == 0 && s1.j == 5)
    {
        return (ppi * x * (15 - 70 * quick_pow(x, 2) + 63 * quick_pow(x, 4)) *
                (f1 - 3 * f2 - quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                 quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 - quick_pow(pmag, 2) * f5 -
                 quick_pow(ppmag, 2) * f5 - 2 * pmag * ppmag * x * f5 - quick_pow(pmag, 2) * f6 -
                 quick_pow(ppmag, 2) * f6 + 2 * pmag * ppmag * x * f6)) /
               4.;
    }
    else if (s2.l == 5 && s1.l == 5 && s1.s == 1 && s1.j == 5)
    {
        return (ppi * (42 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 7) * f4 +
                       21 * pmag * ppmag * quick_pow(x, 6) * (f3 + 5 * f5 - 5 * f6) -
                       35 * pmag * ppmag * quick_pow(x, 4) * (f3 + 3 * f5 - 3 * f6) +
                       15 * pmag * ppmag * quick_pow(x, 2) * (f3 + f5 - f6) - pmag * ppmag * (f3 - f5 + f6) -
                       70 * quick_pow(x, 3) *
                           (f1 + f2 - quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 + quick_pow(pmag, 2) * f5 +
                            quick_pow(ppmag, 2) * f5 + quick_pow(pmag, 2) * f6 + quick_pow(ppmag, 2) * f6) +
                       7 * quick_pow(x, 5) *
                           (9 * f1 + 9 * f2 - 14 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            9 * quick_pow(pmag, 2) * f5 + 9 * quick_pow(ppmag, 2) * f5 + 9 * quick_pow(pmag, 2) * f6 +
                            9 * quick_pow(ppmag, 2) * f6) +
                       x * (15 * f1 + 15 * f2 - 14 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            15 * quick_pow(pmag, 2) * f5 + 15 * quick_pow(ppmag, 2) * f5 +
                            15 * quick_pow(pmag, 2) * f6 + 15 * quick_pow(ppmag, 2) * f6))) /
               4.;
    }
    else if (s2.l == 4 && s1.l == 4 && s1.s == 1 && s1.j == 5)
    {
        return (ppi *
                (11 * (3 - 30 * quick_pow(x, 2) + 35 * quick_pow(x, 4)) * f1 +
                 11 * (3 - 30 * quick_pow(x, 2) + 35 * quick_pow(x, 4)) * f2 - 264 * pmag * ppmag * x * f3 +
                 880 * pmag * ppmag * quick_pow(x, 3) * f3 - 616 * pmag * ppmag * quick_pow(x, 5) * f3 +
                 27 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                 273 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 +
                 505 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 4) * f4 -
                 259 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 6) * f4 + 3 * quick_pow(pmag, 2) * f5 +
                 3 * quick_pow(ppmag, 2) * f5 + 30 * pmag * ppmag * x * f5 -
                 30 * quick_pow(pmag, 2) * quick_pow(x, 2) * f5 - 30 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f5 -
                 140 * pmag * ppmag * quick_pow(x, 3) * f5 + 35 * quick_pow(pmag, 2) * quick_pow(x, 4) * f5 +
                 35 * quick_pow(ppmag, 2) * quick_pow(x, 4) * f5 + 126 * pmag * ppmag * quick_pow(x, 5) * f5 +
                 3 * quick_pow(pmag, 2) * f6 + 3 * quick_pow(ppmag, 2) * f6 - 30 * pmag * ppmag * x * f6 -
                 30 * quick_pow(pmag, 2) * quick_pow(x, 2) * f6 - 30 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f6 +
                 140 * pmag * ppmag * quick_pow(x, 3) * f6 + 35 * quick_pow(pmag, 2) * quick_pow(x, 4) * f6 +
                 35 * quick_pow(ppmag, 2) * quick_pow(x, 4) * f6 - 126 * pmag * ppmag * quick_pow(x, 5) * f6)) /
               44.;
    }
    else if (s2.l == 4 && s1.l == 6 && s1.s == 1 && s1.j == 5)
    {
        return -(sqrt(7.5) * ppi *
                 (4 * pmag * ppmag * x * (15 - 70 * quick_pow(x, 2) + 63 * quick_pow(x, 4)) * (f5 - f6) +
                  quick_pow(ppmag, 2) * (-5 + 105 * quick_pow(x, 2) - 315 * quick_pow(x, 4) + 231 * quick_pow(x, 6)) *
                      (f5 + f6) +
                  quick_pow(pmag, 2) *
                      (quick_pow(ppmag, 2) * (-1 + 15 * quick_pow(x, 2) - 35 * quick_pow(x, 4) + 21 * quick_pow(x, 6)) *
                           f4 +
                       2 * (3 - 30 * quick_pow(x, 2) + 35 * quick_pow(x, 4)) * (f5 + f6)))) /
               22.;
    }
    else if (s2.l == 6 && s1.l == 4 && s1.s == 1 && s1.j == 5)
    {
        return -(sqrt(7.5) * ppi *
                 (4 * pmag * ppmag * x * (15 - 70 * quick_pow(x, 2) + 63 * quick_pow(x, 4)) * (f5 - f6) +
                  2 * quick_pow(ppmag, 2) * (3 - 30 * quick_pow(x, 2) + 35 * quick_pow(x, 4)) * (f5 + f6) +
                  quick_pow(pmag, 2) *
                      (quick_pow(ppmag, 2) * (-1 + 15 * quick_pow(x, 2) - 35 * quick_pow(x, 4) + 21 * quick_pow(x, 6)) *
                           f4 +
                       (-5 + 105 * quick_pow(x, 2) - 315 * quick_pow(x, 4) + 231 * quick_pow(x, 6)) * (f5 + f6)))) /
               22.;
    }
    else if (s2.l == 6 && s1.l == 6 && s1.s == 1 && s1.j == 5)
    {
        return (ppi *
                (11 * (-5 + 105 * quick_pow(x, 2) - 315 * quick_pow(x, 4) + 231 * quick_pow(x, 6)) * f1 +
                 11 * (-5 + 105 * quick_pow(x, 2) - 315 * quick_pow(x, 4) + 231 * quick_pow(x, 6)) * f2 -
                 770 * pmag * ppmag * x * f3 + 5390 * pmag * ppmag * quick_pow(x, 3) * f3 -
                 9702 * pmag * ppmag * quick_pow(x, 5) * f3 + 5082 * pmag * ppmag * quick_pow(x, 7) * f3 -
                 65 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                 1360 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 -
                 4970 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 4) * f4 +
                 6216 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 6) * f4 -
                 2541 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 8) * f4 + 5 * quick_pow(pmag, 2) * f5 +
                 5 * quick_pow(ppmag, 2) * f5 - 60 * pmag * ppmag * x * f5 -
                 105 * quick_pow(pmag, 2) * quick_pow(x, 2) * f5 - 105 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f5 +
                 280 * pmag * ppmag * quick_pow(x, 3) * f5 + 315 * quick_pow(pmag, 2) * quick_pow(x, 4) * f5 +
                 315 * quick_pow(ppmag, 2) * quick_pow(x, 4) * f5 - 252 * pmag * ppmag * quick_pow(x, 5) * f5 -
                 231 * quick_pow(pmag, 2) * quick_pow(x, 6) * f5 - 231 * quick_pow(ppmag, 2) * quick_pow(x, 6) * f5 +
                 5 * quick_pow(pmag, 2) * f6 + 5 * quick_pow(ppmag, 2) * f6 + 60 * pmag * ppmag * x * f6 -
                 105 * quick_pow(pmag, 2) * quick_pow(x, 2) * f6 - 105 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f6 -
                 280 * pmag * ppmag * quick_pow(x, 3) * f6 + 315 * quick_pow(pmag, 2) * quick_pow(x, 4) * f6 +
                 315 * quick_pow(ppmag, 2) * quick_pow(x, 4) * f6 + 252 * pmag * ppmag * quick_pow(x, 5) * f6 -
                 231 * quick_pow(pmag, 2) * quick_pow(x, 6) * f6 - 231 * quick_pow(ppmag, 2) * quick_pow(x, 6) * f6)) /
               88.;
    }
    else if (s2.l == 6 && s1.l == 6 && s1.s == 0 && s1.j == 6)
    {
        return (ppi * (-5 + 105 * quick_pow(x, 2) - 315 * quick_pow(x, 4) + 231 * quick_pow(x, 6)) *
                (f1 - 3 * f2 - quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                 quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 - quick_pow(pmag, 2) * f5 -
                 quick_pow(ppmag, 2) * f5 - 2 * pmag * ppmag * x * f5 - quick_pow(pmag, 2) * f6 -
                 quick_pow(ppmag, 2) * f6 + 2 * pmag * ppmag * x * f6)) /
               8.;
    }
    else if (s2.l == 6 && s1.l == 6 && s1.s == 1 && s1.j == 6)
    {
        return (ppi *
                ((-5 + 105 * quick_pow(x, 2) - 315 * quick_pow(x, 4) + 231 * quick_pow(x, 6)) * f1 +
                 (-5 + 105 * quick_pow(x, 2) - 315 * quick_pow(x, 4) + 231 * quick_pow(x, 6)) * f2 -
                 10 * pmag * ppmag * x * f3 + 70 * pmag * ppmag * quick_pow(x, 3) * f3 -
                 126 * pmag * ppmag * quick_pow(x, 5) * f3 + 66 * pmag * ppmag * quick_pow(x, 7) * f3 +
                 5 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                 100 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 +
                 350 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 4) * f4 -
                 420 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 6) * f4 +
                 165 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 8) * f4 - 5 * quick_pow(pmag, 2) * f5 -
                 5 * quick_pow(ppmag, 2) * f5 + 105 * quick_pow(pmag, 2) * quick_pow(x, 2) * f5 +
                 105 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f5 + 140 * pmag * ppmag * quick_pow(x, 3) * f5 -
                 315 * quick_pow(pmag, 2) * quick_pow(x, 4) * f5 - 315 * quick_pow(ppmag, 2) * quick_pow(x, 4) * f5 -
                 504 * pmag * ppmag * quick_pow(x, 5) * f5 + 231 * quick_pow(pmag, 2) * quick_pow(x, 6) * f5 +
                 231 * quick_pow(ppmag, 2) * quick_pow(x, 6) * f5 + 396 * pmag * ppmag * quick_pow(x, 7) * f5 -
                 5 * quick_pow(pmag, 2) * f6 - 5 * quick_pow(ppmag, 2) * f6 +
                 105 * quick_pow(pmag, 2) * quick_pow(x, 2) * f6 + 105 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f6 -
                 140 * pmag * ppmag * quick_pow(x, 3) * f6 - 315 * quick_pow(pmag, 2) * quick_pow(x, 4) * f6 -
                 315 * quick_pow(ppmag, 2) * quick_pow(x, 4) * f6 + 504 * pmag * ppmag * quick_pow(x, 5) * f6 +
                 231 * quick_pow(pmag, 2) * quick_pow(x, 6) * f6 + 231 * quick_pow(ppmag, 2) * quick_pow(x, 6) * f6 -
                 396 * pmag * ppmag * quick_pow(x, 7) * f6)) /
               8.;
    }
    else if (s2.l == 5 && s1.l == 5 && s1.s == 1 && s1.j == 6)
    {
        return -0.019230769230769232 *
               (ppi * (588 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 7) * f4 -
                       5 * pmag * ppmag * (13 * f3 - f5 + f6) +
                       15 * pmag * ppmag * quick_pow(x, 2) * (65 * f3 - 7 * f5 + 7 * f6) -
                       35 * pmag * ppmag * quick_pow(x, 4) * (65 * f3 - 9 * f5 + 9 * f6) -
                       5 * x *
                           (39 * f1 + 39 * f2 + 32 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            3 * quick_pow(pmag, 2) * f5 + 3 * quick_pow(ppmag, 2) * f5 + 3 * quick_pow(pmag, 2) * f6 +
                            3 * quick_pow(ppmag, 2) * f6) +
                       10 * quick_pow(x, 3) *
                           (91 * f1 + 91 * f2 + 86 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            7 * quick_pow(pmag, 2) * f5 + 7 * quick_pow(ppmag, 2) * f5 + 7 * quick_pow(pmag, 2) * f6 +
                            7 * quick_pow(ppmag, 2) * f6) -
                       7 * quick_pow(x, 5) *
                           (117 * f1 + 117 * f2 + 184 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            9 * quick_pow(pmag, 2) * f5 + 9 * quick_pow(ppmag, 2) * f5 + 9 * quick_pow(pmag, 2) * f6 +
                            9 * quick_pow(ppmag, 2) * f6) +
                       21 * pmag * ppmag * quick_pow(x, 6) * (65 * f3 + 11 * (-f5 + f6))));
    }
    else if (s2.l == 5 && s1.l == 7 && s1.s == 1 && s1.j == 6)
    {
        return -(sqrt(10.5) * ppi *
                 (2 * pmag * ppmag * (-5 + 105 * quick_pow(x, 2) - 315 * quick_pow(x, 4) + 231 * quick_pow(x, 6)) *
                      (f5 - f6) +
                  quick_pow(ppmag, 2) * x *
                      (-35 + 315 * quick_pow(x, 2) - 693 * quick_pow(x, 4) + 429 * quick_pow(x, 6)) * (f5 + f6) +
                  quick_pow(pmag, 2) * x *
                      (quick_pow(ppmag, 2) * (-5 + 35 * quick_pow(x, 2) - 63 * quick_pow(x, 4) + 33 * quick_pow(x, 6)) *
                           f4 +
                       2 * (15 - 70 * quick_pow(x, 2) + 63 * quick_pow(x, 4)) * (f5 + f6)))) /
               26.;
    }
    else if (s2.l == 7 && s1.l == 5 && s1.s == 1 && s1.j == 6)
    {
        return -(sqrt(10.5) * ppi *
                 (2 * pmag * ppmag * (-5 + 105 * quick_pow(x, 2) - 315 * quick_pow(x, 4) + 231 * quick_pow(x, 6)) *
                      (f5 - f6) +
                  2 * quick_pow(ppmag, 2) * x * (15 - 70 * quick_pow(x, 2) + 63 * quick_pow(x, 4)) * (f5 + f6) +
                  quick_pow(pmag, 2) * x *
                      (quick_pow(ppmag, 2) * (-5 + 35 * quick_pow(x, 2) - 63 * quick_pow(x, 4) + 33 * quick_pow(x, 6)) *
                           f4 +
                       (-35 + 315 * quick_pow(x, 2) - 693 * quick_pow(x, 4) + 429 * quick_pow(x, 6)) * (f5 + f6)))) /
               26.;
    }
    else if (s2.l == 7 && s1.l == 7 && s1.s == 1 && s1.j == 6)
    {
        return -0.009615384615384616 *
               (ppi *
                (-11154 * pmag * ppmag * quick_pow(x, 8) * f3 +
                 5577 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 9) * f4 +
                 70 * pmag * ppmag * quick_pow(x, 2) * (52 * f3 + 3 * f5 - 3 * f6) -
                 10 * pmag * ppmag * (13 * f3 + f5 - f6) - 630 * pmag * ppmag * quick_pow(x, 4) * (26 * f3 + f5 - f6) +
                 462 * pmag * ppmag * quick_pow(x, 6) * (52 * f3 + f5 - f6) -
                 33 * quick_pow(x, 7) *
                     (169 * f1 + 169 * f2 + 454 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                      13 * quick_pow(pmag, 2) * f5 - 13 * quick_pow(ppmag, 2) * f5 - 13 * quick_pow(pmag, 2) * f6 -
                      13 * quick_pow(ppmag, 2) * f6) -
                 35 * quick_pow(x, 3) *
                     (117 * f1 + 117 * f2 + 142 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                      9 * quick_pow(pmag, 2) * f5 - 9 * quick_pow(ppmag, 2) * f5 - 9 * quick_pow(pmag, 2) * f6 -
                      9 * quick_pow(ppmag, 2) * f6) +
                 5 * x *
                     (91 * f1 + 91 * f2 + 103 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                      7 * quick_pow(pmag, 2) * f5 - 7 * quick_pow(ppmag, 2) * f5 - 7 * quick_pow(pmag, 2) * f6 -
                      7 * quick_pow(ppmag, 2) * f6) +
                 693 * quick_pow(x, 5) *
                     (13 * f1 + 13 * f2 + 20 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 - quick_pow(pmag, 2) * f5 -
                      quick_pow(ppmag, 2) * f5 - quick_pow(pmag, 2) * f6 - quick_pow(ppmag, 2) * f6)));
    }
    else if (s2.l == 7 && s1.l == 7 && s1.s == 0 && s1.j == 7)
    {
        return (ppi * x * (-35 + 315 * quick_pow(x, 2) - 693 * quick_pow(x, 4) + 429 * quick_pow(x, 6)) *
                (f1 - 3 * f2 - quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                 quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 - quick_pow(pmag, 2) * f5 -
                 quick_pow(ppmag, 2) * f5 - 2 * pmag * ppmag * x * f5 - quick_pow(pmag, 2) * f6 -
                 quick_pow(ppmag, 2) * f6 + 2 * pmag * ppmag * x * f6)) /
               8.;
    }
    else if (s2.l == 7 && s1.l == 7 && s1.s == 1 && s1.j == 7)
    {
        return (ppi * (1287 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 9) * f4 +
                       429 * pmag * ppmag * quick_pow(x, 8) * (f3 + 7 * f5 - 7 * f6) -
                       924 * pmag * ppmag * quick_pow(x, 6) * (f3 + 5 * f5 - 5 * f6) +
                       630 * pmag * ppmag * quick_pow(x, 4) * (f3 + 3 * f5 - 3 * f6) -
                       140 * pmag * ppmag * quick_pow(x, 2) * (f3 + f5 - f6) + 5 * pmag * ppmag * (f3 - f5 + f6) +
                       1260 * quick_pow(x, 3) *
                           (f1 + f2 - quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 + quick_pow(pmag, 2) * f5 +
                            quick_pow(ppmag, 2) * f5 + quick_pow(pmag, 2) * f6 + quick_pow(ppmag, 2) * f6) +
                       132 * quick_pow(x, 7) *
                           (13 * f1 + 13 * f2 - 27 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            13 * quick_pow(pmag, 2) * f5 + 13 * quick_pow(ppmag, 2) * f5 +
                            13 * quick_pow(pmag, 2) * f6 + 13 * quick_pow(ppmag, 2) * f6) -
                       126 * quick_pow(x, 5) *
                           (22 * f1 + 22 * f2 - 27 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            22 * quick_pow(pmag, 2) * f5 + 22 * quick_pow(ppmag, 2) * f5 +
                            22 * quick_pow(pmag, 2) * f6 + 22 * quick_pow(ppmag, 2) * f6) -
                       5 * x *
                           (28 * f1 + 28 * f2 - 27 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            28 * quick_pow(pmag, 2) * f5 + 28 * quick_pow(ppmag, 2) * f5 +
                            28 * quick_pow(pmag, 2) * f6 + 28 * quick_pow(ppmag, 2) * f6))) /
               32.;
    }
    else if (s2.l == 6 && s1.l == 6 && s1.s == 1 && s1.j == 7)
    {
        return (ppi *
                (15 * (-5 + 105 * quick_pow(x, 2) - 315 * quick_pow(x, 4) + 231 * quick_pow(x, 6)) * f1 +
                 15 * (-5 + 105 * quick_pow(x, 2) - 315 * quick_pow(x, 4) + 231 * quick_pow(x, 6)) * f2 +
                 900 * pmag * ppmag * x * f3 - 6300 * pmag * ppmag * quick_pow(x, 3) * f3 +
                 11340 * pmag * ppmag * quick_pow(x, 5) * f3 - 5940 * pmag * ppmag * quick_pow(x, 7) * f3 -
                 65 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                 1370 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 -
                 5040 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 4) * f4 +
                 6342 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 6) * f4 -
                 2607 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 8) * f4 - 5 * quick_pow(pmag, 2) * f5 -
                 5 * quick_pow(ppmag, 2) * f5 - 70 * pmag * ppmag * x * f5 +
                 105 * quick_pow(pmag, 2) * quick_pow(x, 2) * f5 + 105 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f5 +
                 630 * pmag * ppmag * quick_pow(x, 3) * f5 - 315 * quick_pow(pmag, 2) * quick_pow(x, 4) * f5 -
                 315 * quick_pow(ppmag, 2) * quick_pow(x, 4) * f5 - 1386 * pmag * ppmag * quick_pow(x, 5) * f5 +
                 231 * quick_pow(pmag, 2) * quick_pow(x, 6) * f5 + 231 * quick_pow(ppmag, 2) * quick_pow(x, 6) * f5 +
                 858 * pmag * ppmag * quick_pow(x, 7) * f5 - 5 * quick_pow(pmag, 2) * f6 -
                 5 * quick_pow(ppmag, 2) * f6 + 70 * pmag * ppmag * x * f6 +
                 105 * quick_pow(pmag, 2) * quick_pow(x, 2) * f6 + 105 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f6 -
                 630 * pmag * ppmag * quick_pow(x, 3) * f6 - 315 * quick_pow(pmag, 2) * quick_pow(x, 4) * f6 -
                 315 * quick_pow(ppmag, 2) * quick_pow(x, 4) * f6 + 1386 * pmag * ppmag * quick_pow(x, 5) * f6 +
                 231 * quick_pow(pmag, 2) * quick_pow(x, 6) * f6 + 231 * quick_pow(ppmag, 2) * quick_pow(x, 6) * f6 -
                 858 * pmag * ppmag * quick_pow(x, 7) * f6)) /
               120.;
    }
    else if (s2.l == 6 && s1.l == 8 && s1.s == 1 && s1.j == 7)
    {
        return -(sqrt(3.5) * ppi *
                 (16 * pmag * ppmag * x *
                      (-35 + 315 * quick_pow(x, 2) - 693 * quick_pow(x, 4) + 429 * quick_pow(x, 6)) * (f5 - f6) +
                  quick_pow(ppmag, 2) *
                      (35 - 1260 * quick_pow(x, 2) + 6930 * quick_pow(x, 4) - 12012 * quick_pow(x, 6) +
                       6435 * quick_pow(x, 8)) *
                      (f5 + f6) +
                  quick_pow(pmag, 2) *
                      (quick_pow(ppmag, 2) *
                           (5 - 140 * quick_pow(x, 2) + 630 * quick_pow(x, 4) - 924 * quick_pow(x, 6) +
                            429 * quick_pow(x, 8)) *
                           f4 +
                       8 * (-5 + 105 * quick_pow(x, 2) - 315 * quick_pow(x, 4) + 231 * quick_pow(x, 6)) * (f5 + f6)))) /
               120.;
    }
    else if (s2.l == 8 && s1.l == 6 && s1.s == 1 && s1.j == 7)
    {
        return -(sqrt(3.5) * ppi *
                 (16 * pmag * ppmag * x *
                      (-35 + 315 * quick_pow(x, 2) - 693 * quick_pow(x, 4) + 429 * quick_pow(x, 6)) * (f5 - f6) +
                  8 * quick_pow(ppmag, 2) *
                      (-5 + 105 * quick_pow(x, 2) - 315 * quick_pow(x, 4) + 231 * quick_pow(x, 6)) * (f5 + f6) +
                  quick_pow(pmag, 2) * (quick_pow(ppmag, 2) *
                                            (5 - 140 * quick_pow(x, 2) + 630 * quick_pow(x, 4) - 924 * quick_pow(x, 6) +
                                             429 * quick_pow(x, 8)) *
                                            f4 +
                                        (35 - 1260 * quick_pow(x, 2) + 6930 * quick_pow(x, 4) -
                                         12012 * quick_pow(x, 6) + 6435 * quick_pow(x, 8)) *
                                            (f5 + f6)))) /
               120.;
    }
    else if (s2.l == 8 && s1.l == 8 && s1.s == 1 && s1.j == 7)
    {
        return (ppi *
                (15 *
                     (35 - 1260 * quick_pow(x, 2) + 6930 * quick_pow(x, 4) - 12012 * quick_pow(x, 6) +
                      6435 * quick_pow(x, 8)) *
                     f1 +
                 15 *
                     (35 - 1260 * quick_pow(x, 2) + 6930 * quick_pow(x, 4) - 12012 * quick_pow(x, 6) +
                      6435 * quick_pow(x, 8)) *
                     f2 +
                 9450 * pmag * ppmag * x * f3 - 113400 * pmag * ppmag * quick_pow(x, 3) * f3 +
                 374220 * pmag * ppmag * quick_pow(x, 5) * f3 - 463320 * pmag * ppmag * quick_pow(x, 7) * f3 +
                 193050 * pmag * ppmag * quick_pow(x, 9) * f3 + 595 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                 21385 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 +
                 131670 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 4) * f4 -
                 297066 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 6) * f4 +
                 282711 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 8) * f4 -
                 96525 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 10) * f4 -
                 35 * quick_pow(pmag, 2) * f5 - 35 * quick_pow(ppmag, 2) * f5 + 560 * pmag * ppmag * x * f5 +
                 1260 * quick_pow(pmag, 2) * quick_pow(x, 2) * f5 + 1260 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f5 -
                 5040 * pmag * ppmag * quick_pow(x, 3) * f5 - 6930 * quick_pow(pmag, 2) * quick_pow(x, 4) * f5 -
                 6930 * quick_pow(ppmag, 2) * quick_pow(x, 4) * f5 + 11088 * pmag * ppmag * quick_pow(x, 5) * f5 +
                 12012 * quick_pow(pmag, 2) * quick_pow(x, 6) * f5 +
                 12012 * quick_pow(ppmag, 2) * quick_pow(x, 6) * f5 - 6864 * pmag * ppmag * quick_pow(x, 7) * f5 -
                 6435 * quick_pow(pmag, 2) * quick_pow(x, 8) * f5 - 6435 * quick_pow(ppmag, 2) * quick_pow(x, 8) * f5 -
                 35 * quick_pow(pmag, 2) * f6 - 35 * quick_pow(ppmag, 2) * f6 - 560 * pmag * ppmag * x * f6 +
                 1260 * quick_pow(pmag, 2) * quick_pow(x, 2) * f6 + 1260 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f6 +
                 5040 * pmag * ppmag * quick_pow(x, 3) * f6 - 6930 * quick_pow(pmag, 2) * quick_pow(x, 4) * f6 -
                 6930 * quick_pow(ppmag, 2) * quick_pow(x, 4) * f6 - 11088 * pmag * ppmag * quick_pow(x, 5) * f6 +
                 12012 * quick_pow(pmag, 2) * quick_pow(x, 6) * f6 +
                 12012 * quick_pow(ppmag, 2) * quick_pow(x, 6) * f6 + 6864 * pmag * ppmag * quick_pow(x, 7) * f6 -
                 6435 * quick_pow(pmag, 2) * quick_pow(x, 8) * f6 -
                 6435 * quick_pow(ppmag, 2) * quick_pow(x, 8) * f6)) /
               960.;
    }
    else if (s2.l == 8 && s1.l == 8 && s1.s == 0 && s1.j == 8)
    {
        return (ppi *
                (35 - 1260 * quick_pow(x, 2) + 6930 * quick_pow(x, 4) - 12012 * quick_pow(x, 6) +
                 6435 * quick_pow(x, 8)) *
                (f1 - 3 * f2 - quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                 quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 - quick_pow(pmag, 2) * f5 -
                 quick_pow(ppmag, 2) * f5 - 2 * pmag * ppmag * x * f5 - quick_pow(pmag, 2) * f6 -
                 quick_pow(ppmag, 2) * f6 + 2 * pmag * ppmag * x * f6)) /
               64.;
    }
    else if (s2.l == 8 && s1.l == 8 && s1.s == 1 && s1.j == 8)
    {
        return (ppi * (70 * quick_pow(ppmag, 2) * quick_pow(-1 + quick_pow(x, 2), 2) *
                           (-1 + 33 * quick_pow(x, 2) - 143 * quick_pow(x, 4) + 143 * quick_pow(x, 6)) *
                           (quick_pow(pmag, 2) * f4 - f5 - f6) +
                       2 *
                           (35 - 1260 * quick_pow(x, 2) + 6930 * quick_pow(x, 4) - 12012 * quick_pow(x, 6) +
                            6435 * quick_pow(x, 8)) *
                           (f1 + f2 + quick_pow(pmag + ppmag * x, 2) * f5 + quick_pow(pmag - ppmag * x, 2) * f6) +
                       4 * ppmag * x * (1 - quick_pow(x, 2)) *
                           (-35 + 385 * quick_pow(x, 2) - 1001 * quick_pow(x, 4) + 715 * quick_pow(x, 6)) *
                           (-(pmag * (f3 - f5 + f6)) + ppmag * x * (f5 + f6)))) /
               128.;
    }
    else if (s2.l == 7 && s1.l == 7 && s1.s == 1 && s1.j == 8)
    {
        return -0.001838235294117647 *
               (ppi * (22737 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 9) * f4 +
                       35 * pmag * ppmag * (17 * f3 - f5 + f6) -
                       140 * pmag * ppmag * quick_pow(x, 2) * (119 * f3 - 9 * f5 + 9 * f6) +
                       35 * x *
                           (68 * f1 + 68 * f2 + 59 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            4 * quick_pow(pmag, 2) * f5 + 4 * quick_pow(ppmag, 2) * f5 + 4 * quick_pow(pmag, 2) * f6 +
                            4 * quick_pow(ppmag, 2) * f6) -
                       140 * quick_pow(x, 3) *
                           (153 * f1 + 153 * f2 + 143 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            9 * quick_pow(pmag, 2) * f5 + 9 * quick_pow(ppmag, 2) * f5 + 9 * quick_pow(pmag, 2) * f6 +
                            9 * quick_pow(ppmag, 2) * f6) -
                       132 * quick_pow(x, 7) *
                           (221 * f1 + 221 * f2 + 461 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            13 * quick_pow(pmag, 2) * f5 + 13 * quick_pow(ppmag, 2) * f5 +
                            13 * quick_pow(pmag, 2) * f6 + 13 * quick_pow(ppmag, 2) * f6) +
                       126 * quick_pow(x, 5) *
                           (374 * f1 + 374 * f2 + 445 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            22 * quick_pow(pmag, 2) * f5 + 22 * quick_pow(ppmag, 2) * f5 +
                            22 * quick_pow(pmag, 2) * f6 + 22 * quick_pow(ppmag, 2) * f6) +
                       630 * pmag * ppmag * quick_pow(x, 4) * (119 * f3 + 11 * (-f5 + f6)) -
                       924 * pmag * ppmag * quick_pow(x, 6) * (119 * f3 + 13 * (-f5 + f6)) +
                       429 * pmag * ppmag * quick_pow(x, 8) * (119 * f3 + 15 * (-f5 + f6))));
    }
    else if (s2.l == 7 && s1.l == 9 && s1.s == 1 && s1.j == 8)
    {
        return -(3 * ppi *
                 (2 * pmag * ppmag *
                      (35 - 1260 * quick_pow(x, 2) + 6930 * quick_pow(x, 4) - 12012 * quick_pow(x, 6) +
                       6435 * quick_pow(x, 8)) *
                      (f5 - f6) +
                  quick_pow(ppmag, 2) * x *
                      (315 - 4620 * quick_pow(x, 2) + 18018 * quick_pow(x, 4) - 25740 * quick_pow(x, 6) +
                       12155 * quick_pow(x, 8)) *
                      (f5 + f6) +
                  quick_pow(pmag, 2) * x *
                      (quick_pow(ppmag, 2) *
                           (35 - 420 * quick_pow(x, 2) + 1386 * quick_pow(x, 4) - 1716 * quick_pow(x, 6) +
                            715 * quick_pow(x, 8)) *
                           f4 +
                       8 * (-35 + 315 * quick_pow(x, 2) - 693 * quick_pow(x, 4) + 429 * quick_pow(x, 6)) *
                           (f5 + f6)))) /
               (136. * sqrt(2));
    }
    else if (s2.l == 9 && s1.l == 7 && s1.s == 1 && s1.j == 8)
    {
        return -(3 * ppi *
                 (2 * pmag * ppmag *
                      (35 - 1260 * quick_pow(x, 2) + 6930 * quick_pow(x, 4) - 12012 * quick_pow(x, 6) +
                       6435 * quick_pow(x, 8)) *
                      (f5 - f6) +
                  8 * quick_pow(ppmag, 2) * x *
                      (-35 + 315 * quick_pow(x, 2) - 693 * quick_pow(x, 4) + 429 * quick_pow(x, 6)) * (f5 + f6) +
                  quick_pow(pmag, 2) * x *
                      (quick_pow(ppmag, 2) *
                           (35 - 420 * quick_pow(x, 2) + 1386 * quick_pow(x, 4) - 1716 * quick_pow(x, 6) +
                            715 * quick_pow(x, 8)) *
                           f4 +
                       (315 - 4620 * quick_pow(x, 2) + 18018 * quick_pow(x, 4) - 25740 * quick_pow(x, 6) +
                        12155 * quick_pow(x, 8)) *
                           (f5 + f6)))) /
               (136. * sqrt(2));
    }
    else if (s2.l == 9 && s1.l == 9 && s1.s == 1 && s1.j == 8)
    {
        return -0.0009191176470588235 *
               (ppi * (-413270 * pmag * ppmag * quick_pow(x, 10) * f3 +
                       206635 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 11) * f4 -
                       630 * pmag * ppmag * quick_pow(x, 2) * (85 * f3 + 4 * f5 - 4 * f6) +
                       4620 * pmag * ppmag * quick_pow(x, 4) * (85 * f3 + 3 * f5 - 3 * f6) -
                       12012 * pmag * ppmag * quick_pow(x, 6) * (85 * f3 + 2 * f5 - 2 * f6) +
                       70 * pmag * ppmag * (17 * f3 + f5 - f6) +
                       12870 * pmag * ppmag * quick_pow(x, 8) * (85 * f3 + f5 - f6) +
                       105 * quick_pow(x, 3) *
                           (748 * f1 + 748 * f2 + 863 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            44 * quick_pow(pmag, 2) * f5 - 44 * quick_pow(ppmag, 2) * f5 -
                            44 * quick_pow(pmag, 2) * f6 - 44 * quick_pow(ppmag, 2) * f6) -
                       462 * quick_pow(x, 5) *
                           (663 * f1 + 663 * f2 + 881 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            39 * quick_pow(pmag, 2) * f5 - 39 * quick_pow(ppmag, 2) * f5 -
                            39 * quick_pow(pmag, 2) * f6 - 39 * quick_pow(ppmag, 2) * f6) +
                       858 * quick_pow(x, 7) *
                           (510 * f1 + 510 * f2 + 899 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            30 * quick_pow(pmag, 2) * f5 - 30 * quick_pow(ppmag, 2) * f5 -
                            30 * quick_pow(pmag, 2) * f6 - 30 * quick_pow(ppmag, 2) * f6) -
                       715 * quick_pow(x, 9) *
                           (289 * f1 + 289 * f2 + 917 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            17 * quick_pow(pmag, 2) * f5 - 17 * quick_pow(ppmag, 2) * f5 -
                            17 * quick_pow(pmag, 2) * f6 - 17 * quick_pow(ppmag, 2) * f6) -
                       35 * x *
                           (153 * f1 + 153 * f2 + 169 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            9 * quick_pow(pmag, 2) * f5 - 9 * quick_pow(ppmag, 2) * f5 - 9 * quick_pow(pmag, 2) * f6 -
                            9 * quick_pow(ppmag, 2) * f6)));
    }
    else if (s2.l == 8 && s1.l == 8 && s1.s == 0 && s1.j == 8)
    {
        return (ppi *
                (35 - 1260 * quick_pow(x, 2) + 6930 * quick_pow(x, 4) - 12012 * quick_pow(x, 6) +
                 6435 * quick_pow(x, 8)) *
                (f1 - 3 * f2 - quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                 quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 - quick_pow(pmag, 2) * f5 -
                 quick_pow(ppmag, 2) * f5 - 2 * pmag * ppmag * x * f5 - quick_pow(pmag, 2) * f6 -
                 quick_pow(ppmag, 2) * f6 + 2 * pmag * ppmag * x * f6)) /
               64.;
    }
    else if (s2.l == 8 && s1.l == 8 && s1.s == 1 && s1.j == 8)
    {
        return (ppi * (70 * quick_pow(ppmag, 2) * quick_pow(-1 + quick_pow(x, 2), 2) *
                           (-1 + 33 * quick_pow(x, 2) - 143 * quick_pow(x, 4) + 143 * quick_pow(x, 6)) *
                           (quick_pow(pmag, 2) * f4 - f5 - f6) +
                       2 *
                           (35 - 1260 * quick_pow(x, 2) + 6930 * quick_pow(x, 4) - 12012 * quick_pow(x, 6) +
                            6435 * quick_pow(x, 8)) *
                           (f1 + f2 + quick_pow(pmag + ppmag * x, 2) * f5 + quick_pow(pmag - ppmag * x, 2) * f6) +
                       4 * ppmag * x * (1 - quick_pow(x, 2)) *
                           (-35 + 385 * quick_pow(x, 2) - 1001 * quick_pow(x, 4) + 715 * quick_pow(x, 6)) *
                           (-(pmag * (f3 - f5 + f6)) + ppmag * x * (f5 + f6)))) /
               128.;
    }
    else if (s2.l == 7 && s1.l == 7 && s1.s == 1 && s1.j == 8)
    {
        return -0.001838235294117647 *
               (ppi * (22737 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 9) * f4 +
                       35 * pmag * ppmag * (17 * f3 - f5 + f6) -
                       140 * pmag * ppmag * quick_pow(x, 2) * (119 * f3 - 9 * f5 + 9 * f6) +
                       35 * x *
                           (68 * f1 + 68 * f2 + 59 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            4 * quick_pow(pmag, 2) * f5 + 4 * quick_pow(ppmag, 2) * f5 + 4 * quick_pow(pmag, 2) * f6 +
                            4 * quick_pow(ppmag, 2) * f6) -
                       140 * quick_pow(x, 3) *
                           (153 * f1 + 153 * f2 + 143 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            9 * quick_pow(pmag, 2) * f5 + 9 * quick_pow(ppmag, 2) * f5 + 9 * quick_pow(pmag, 2) * f6 +
                            9 * quick_pow(ppmag, 2) * f6) -
                       132 * quick_pow(x, 7) *
                           (221 * f1 + 221 * f2 + 461 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            13 * quick_pow(pmag, 2) * f5 + 13 * quick_pow(ppmag, 2) * f5 +
                            13 * quick_pow(pmag, 2) * f6 + 13 * quick_pow(ppmag, 2) * f6) +
                       126 * quick_pow(x, 5) *
                           (374 * f1 + 374 * f2 + 445 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            22 * quick_pow(pmag, 2) * f5 + 22 * quick_pow(ppmag, 2) * f5 +
                            22 * quick_pow(pmag, 2) * f6 + 22 * quick_pow(ppmag, 2) * f6) +
                       630 * pmag * ppmag * quick_pow(x, 4) * (119 * f3 + 11 * (-f5 + f6)) -
                       924 * pmag * ppmag * quick_pow(x, 6) * (119 * f3 + 13 * (-f5 + f6)) +
                       429 * pmag * ppmag * quick_pow(x, 8) * (119 * f3 + 15 * (-f5 + f6))));
    }
    else if (s2.l == 7 && s1.l == 9 && s1.s == 1 && s1.j == 8)
    {
        return (-3 * ppi *
                (2 * pmag * ppmag *
                     (35 - 1260 * quick_pow(x, 2) + 6930 * quick_pow(x, 4) - 12012 * quick_pow(x, 6) +
                      6435 * quick_pow(x, 8)) *
                     (f5 - f6) +
                 quick_pow(ppmag, 2) * x *
                     (315 - 4620 * quick_pow(x, 2) + 18018 * quick_pow(x, 4) - 25740 * quick_pow(x, 6) +
                      12155 * quick_pow(x, 8)) *
                     (f5 + f6) +
                 quick_pow(pmag, 2) * x *
                     (quick_pow(ppmag, 2) *
                          (35 - 420 * quick_pow(x, 2) + 1386 * quick_pow(x, 4) - 1716 * quick_pow(x, 6) +
                           715 * quick_pow(x, 8)) *
                          f4 +
                      8 * (-35 + 315 * quick_pow(x, 2) - 693 * quick_pow(x, 4) + 429 * quick_pow(x, 6)) * (f5 + f6)))) /
               (136. * sqrt(2));
    }
    else if (s2.l == 9 && s1.l == 7 && s1.s == 1 && s1.j == 8)
    {
        return (-3 * ppi *
                (2 * pmag * ppmag *
                     (35 - 1260 * quick_pow(x, 2) + 6930 * quick_pow(x, 4) - 12012 * quick_pow(x, 6) +
                      6435 * quick_pow(x, 8)) *
                     (f5 - f6) +
                 8 * quick_pow(ppmag, 2) * x *
                     (-35 + 315 * quick_pow(x, 2) - 693 * quick_pow(x, 4) + 429 * quick_pow(x, 6)) * (f5 + f6) +
                 quick_pow(pmag, 2) * x *
                     (quick_pow(ppmag, 2) *
                          (35 - 420 * quick_pow(x, 2) + 1386 * quick_pow(x, 4) - 1716 * quick_pow(x, 6) +
                           715 * quick_pow(x, 8)) *
                          f4 +
                      (315 - 4620 * quick_pow(x, 2) + 18018 * quick_pow(x, 4) - 25740 * quick_pow(x, 6) +
                       12155 * quick_pow(x, 8)) *
                          (f5 + f6)))) /
               (136. * sqrt(2));
    }
    else if (s2.l == 9 && s1.l == 9 && s1.s == 1 && s1.j == 8)
    {
        return -0.0009191176470588235 *
               (ppi * (-413270 * pmag * ppmag * quick_pow(x, 10) * f3 +
                       206635 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 11) * f4 -
                       630 * pmag * ppmag * quick_pow(x, 2) * (85 * f3 + 4 * f5 - 4 * f6) +
                       4620 * pmag * ppmag * quick_pow(x, 4) * (85 * f3 + 3 * f5 - 3 * f6) -
                       12012 * pmag * ppmag * quick_pow(x, 6) * (85 * f3 + 2 * f5 - 2 * f6) +
                       70 * pmag * ppmag * (17 * f3 + f5 - f6) +
                       12870 * pmag * ppmag * quick_pow(x, 8) * (85 * f3 + f5 - f6) +
                       105 * quick_pow(x, 3) *
                           (748 * f1 + 748 * f2 + 863 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            44 * quick_pow(pmag, 2) * f5 - 44 * quick_pow(ppmag, 2) * f5 -
                            44 * quick_pow(pmag, 2) * f6 - 44 * quick_pow(ppmag, 2) * f6) -
                       462 * quick_pow(x, 5) *
                           (663 * f1 + 663 * f2 + 881 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            39 * quick_pow(pmag, 2) * f5 - 39 * quick_pow(ppmag, 2) * f5 -
                            39 * quick_pow(pmag, 2) * f6 - 39 * quick_pow(ppmag, 2) * f6) +
                       858 * quick_pow(x, 7) *
                           (510 * f1 + 510 * f2 + 899 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            30 * quick_pow(pmag, 2) * f5 - 30 * quick_pow(ppmag, 2) * f5 -
                            30 * quick_pow(pmag, 2) * f6 - 30 * quick_pow(ppmag, 2) * f6) -
                       715 * quick_pow(x, 9) *
                           (289 * f1 + 289 * f2 + 917 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            17 * quick_pow(pmag, 2) * f5 - 17 * quick_pow(ppmag, 2) * f5 -
                            17 * quick_pow(pmag, 2) * f6 - 17 * quick_pow(ppmag, 2) * f6) -
                       35 * x *
                           (153 * f1 + 153 * f2 + 169 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            9 * quick_pow(pmag, 2) * f5 - 9 * quick_pow(ppmag, 2) * f5 - 9 * quick_pow(pmag, 2) * f6 -
                            9 * quick_pow(ppmag, 2) * f6)));
    }
    else if (s2.l == 9 && s1.l == 9 && s1.s == 0 && s1.j == 9)
    {
        return (ppi * x *
                (315 - 4620 * quick_pow(x, 2) + 18018 * quick_pow(x, 4) - 25740 * quick_pow(x, 6) +
                 12155 * quick_pow(x, 8)) *
                (f1 - 3 * f2 - quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                 quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 - quick_pow(pmag, 2) * f5 -
                 quick_pow(ppmag, 2) * f5 - 2 * pmag * ppmag * x * f5 - quick_pow(pmag, 2) * f6 -
                 quick_pow(ppmag, 2) * f6 + 2 * pmag * ppmag * x * f6)) /
               64.;
    }
    else if (s2.l == 9 && s1.l == 9 && s1.s == 1 && s1.j == 9)
    {
        return (ppi * (9724 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 11) * f4 +
                       2431 * pmag * ppmag * quick_pow(x, 10) * (f3 + 9 * f5 - 9 * f6) -
                       6435 * pmag * ppmag * quick_pow(x, 8) * (f3 + 7 * f5 - 7 * f6) +
                       6006 * pmag * ppmag * quick_pow(x, 6) * (f3 + 5 * f5 - 5 * f6) -
                       2310 * pmag * ppmag * quick_pow(x, 4) * (f3 + 3 * f5 - 3 * f6) +
                       315 * pmag * ppmag * quick_pow(x, 2) * (f3 + f5 - f6) - 7 * pmag * ppmag * (f3 - f5 + f6) -
                       4620 * quick_pow(x, 3) *
                           (f1 + f2 - quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 + quick_pow(pmag, 2) * f5 +
                            quick_pow(ppmag, 2) * f5 + quick_pow(pmag, 2) * f6 + quick_pow(ppmag, 2) * f6) -
                       1716 * quick_pow(x, 7) *
                           (15 * f1 + 15 * f2 - 22 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            15 * quick_pow(pmag, 2) * f5 + 15 * quick_pow(ppmag, 2) * f5 +
                            15 * quick_pow(pmag, 2) * f6 + 15 * quick_pow(ppmag, 2) * f6) +
                       715 * quick_pow(x, 9) *
                           (17 * f1 + 17 * f2 - 44 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            17 * quick_pow(pmag, 2) * f5 + 17 * quick_pow(ppmag, 2) * f5 +
                            17 * quick_pow(pmag, 2) * f6 + 17 * quick_pow(ppmag, 2) * f6) +
                       462 * quick_pow(x, 5) *
                           (39 * f1 + 39 * f2 - 44 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            39 * quick_pow(pmag, 2) * f5 + 39 * quick_pow(ppmag, 2) * f5 +
                            39 * quick_pow(pmag, 2) * f6 + 39 * quick_pow(ppmag, 2) * f6) +
                       7 * x *
                           (45 * f1 + 45 * f2 - 44 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            45 * quick_pow(pmag, 2) * f5 + 45 * quick_pow(ppmag, 2) * f5 +
                            45 * quick_pow(pmag, 2) * f6 + 45 * quick_pow(ppmag, 2) * f6))) /
               64.;
    }
    else if (s2.l == 8 && s1.l == 8 && s1.s == 1 && s1.j == 9)
    {
        return (ppi *
                (19 *
                     (35 - 1260 * quick_pow(x, 2) + 6930 * quick_pow(x, 4) - 12012 * quick_pow(x, 6) +
                      6435 * quick_pow(x, 8)) *
                     f1 +
                 19 *
                     (35 - 1260 * quick_pow(x, 2) + 6930 * quick_pow(x, 4) - 12012 * quick_pow(x, 6) +
                      6435 * quick_pow(x, 8)) *
                     f2 -
                 10640 * pmag * ppmag * x * f3 + 127680 * pmag * ppmag * quick_pow(x, 3) * f3 -
                 421344 * pmag * ppmag * quick_pow(x, 5) * f3 + 521664 * pmag * ppmag * quick_pow(x, 7) * f3 -
                 217360 * pmag * ppmag * quick_pow(x, 9) * f3 + 595 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                 21455 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 +
                 132510 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 4) * f4 -
                 299838 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 6) * f4 +
                 286143 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 8) * f4 -
                 97955 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 10) * f4 +
                 35 * quick_pow(pmag, 2) * f5 + 35 * quick_pow(ppmag, 2) * f5 + 630 * pmag * ppmag * x * f5 -
                 1260 * quick_pow(pmag, 2) * quick_pow(x, 2) * f5 - 1260 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f5 -
                 9240 * pmag * ppmag * quick_pow(x, 3) * f5 + 6930 * quick_pow(pmag, 2) * quick_pow(x, 4) * f5 +
                 6930 * quick_pow(ppmag, 2) * quick_pow(x, 4) * f5 + 36036 * pmag * ppmag * quick_pow(x, 5) * f5 -
                 12012 * quick_pow(pmag, 2) * quick_pow(x, 6) * f5 -
                 12012 * quick_pow(ppmag, 2) * quick_pow(x, 6) * f5 - 51480 * pmag * ppmag * quick_pow(x, 7) * f5 +
                 6435 * quick_pow(pmag, 2) * quick_pow(x, 8) * f5 + 6435 * quick_pow(ppmag, 2) * quick_pow(x, 8) * f5 +
                 24310 * pmag * ppmag * quick_pow(x, 9) * f5 + 35 * quick_pow(pmag, 2) * f6 +
                 35 * quick_pow(ppmag, 2) * f6 - 630 * pmag * ppmag * x * f6 -
                 1260 * quick_pow(pmag, 2) * quick_pow(x, 2) * f6 - 1260 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f6 +
                 9240 * pmag * ppmag * quick_pow(x, 3) * f6 + 6930 * quick_pow(pmag, 2) * quick_pow(x, 4) * f6 +
                 6930 * quick_pow(ppmag, 2) * quick_pow(x, 4) * f6 - 36036 * pmag * ppmag * quick_pow(x, 5) * f6 -
                 12012 * quick_pow(pmag, 2) * quick_pow(x, 6) * f6 -
                 12012 * quick_pow(ppmag, 2) * quick_pow(x, 6) * f6 + 51480 * pmag * ppmag * quick_pow(x, 7) * f6 +
                 6435 * quick_pow(pmag, 2) * quick_pow(x, 8) * f6 + 6435 * quick_pow(ppmag, 2) * quick_pow(x, 8) * f6 -
                 24310 * pmag * ppmag * quick_pow(x, 9) * f6)) /
               1216.;
    }
    else if (s2.l == 8 && s1.l == 10 && s1.s == 1 && s1.j == 9)
    {
        return (-3 * sqrt(2.5) * ppi *
                (4 * pmag * ppmag * x *
                     (315 - 4620 * quick_pow(x, 2) + 18018 * quick_pow(x, 4) - 25740 * quick_pow(x, 6) +
                      12155 * quick_pow(x, 8)) *
                     (f5 - f6) +
                 quick_pow(ppmag, 2) *
                     (-63 + 3465 * quick_pow(x, 2) - 30030 * quick_pow(x, 4) + 90090 * quick_pow(x, 6) -
                      109395 * quick_pow(x, 8) + 46189 * quick_pow(x, 10)) *
                     (f5 + f6) +
                 quick_pow(pmag, 2) * (quick_pow(ppmag, 2) *
                                           (-7 + 315 * quick_pow(x, 2) - 2310 * quick_pow(x, 4) +
                                            6006 * quick_pow(x, 6) - 6435 * quick_pow(x, 8) + 2431 * quick_pow(x, 10)) *
                                           f4 +
                                       2 *
                                           (35 - 1260 * quick_pow(x, 2) + 6930 * quick_pow(x, 4) -
                                            12012 * quick_pow(x, 6) + 6435 * quick_pow(x, 8)) *
                                           (f5 + f6)))) /
               608.;
    }
    else if (s2.l == 10 && s1.l == 8 && s1.s == 1 && s1.j == 9)
    {
        return (-3 * sqrt(2.5) * ppi *
                (4 * pmag * ppmag * x *
                     (315 - 4620 * quick_pow(x, 2) + 18018 * quick_pow(x, 4) - 25740 * quick_pow(x, 6) +
                      12155 * quick_pow(x, 8)) *
                     (f5 - f6) +
                 2 * quick_pow(ppmag, 2) *
                     (35 - 1260 * quick_pow(x, 2) + 6930 * quick_pow(x, 4) - 12012 * quick_pow(x, 6) +
                      6435 * quick_pow(x, 8)) *
                     (f5 + f6) +
                 quick_pow(pmag, 2) * (quick_pow(ppmag, 2) *
                                           (-7 + 315 * quick_pow(x, 2) - 2310 * quick_pow(x, 4) +
                                            6006 * quick_pow(x, 6) - 6435 * quick_pow(x, 8) + 2431 * quick_pow(x, 10)) *
                                           f4 +
                                       (-63 + 3465 * quick_pow(x, 2) - 30030 * quick_pow(x, 4) +
                                        90090 * quick_pow(x, 6) - 109395 * quick_pow(x, 8) + 46189 * quick_pow(x, 10)) *
                                           (f5 + f6)))) /
               608.;
    }
    else if (s2.l == 10 && s1.l == 10 && s1.s == 1 && s1.j == 9)
    {
        return (ppi *
                (19 *
                     (-63 + 3465 * quick_pow(x, 2) - 30030 * quick_pow(x, 4) + 90090 * quick_pow(x, 6) -
                      109395 * quick_pow(x, 8) + 46189 * quick_pow(x, 10)) *
                     f1 +
                 19 *
                     (-63 + 3465 * quick_pow(x, 2) - 30030 * quick_pow(x, 4) + 90090 * quick_pow(x, 6) -
                      109395 * quick_pow(x, 8) + 46189 * quick_pow(x, 10)) *
                     f2 -
                 26334 * pmag * ppmag * x * f3 + 482790 * pmag * ppmag * quick_pow(x, 3) * f3 -
                 2510508 * pmag * ppmag * quick_pow(x, 5) * f3 + 5379660 * pmag * ppmag * quick_pow(x, 7) * f3 -
                 5080790 * pmag * ppmag * quick_pow(x, 9) * f3 + 1755182 * pmag * ppmag * quick_pow(x, 11) * f3 -
                 1323 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                 72702 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 -
                 677985 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 4) * f4 +
                 2390388 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 6) * f4 -
                 3906045 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 8) * f4 +
                 2999854 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 10) * f4 -
                 877591 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 12) * f4 +
                 63 * quick_pow(pmag, 2) * f5 + 63 * quick_pow(ppmag, 2) * f5 - 1260 * pmag * ppmag * x * f5 -
                 3465 * quick_pow(pmag, 2) * quick_pow(x, 2) * f5 - 3465 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f5 +
                 18480 * pmag * ppmag * quick_pow(x, 3) * f5 + 30030 * quick_pow(pmag, 2) * quick_pow(x, 4) * f5 +
                 30030 * quick_pow(ppmag, 2) * quick_pow(x, 4) * f5 - 72072 * pmag * ppmag * quick_pow(x, 5) * f5 -
                 90090 * quick_pow(pmag, 2) * quick_pow(x, 6) * f5 -
                 90090 * quick_pow(ppmag, 2) * quick_pow(x, 6) * f5 + 102960 * pmag * ppmag * quick_pow(x, 7) * f5 +
                 109395 * quick_pow(pmag, 2) * quick_pow(x, 8) * f5 +
                 109395 * quick_pow(ppmag, 2) * quick_pow(x, 8) * f5 - 48620 * pmag * ppmag * quick_pow(x, 9) * f5 -
                 46189 * quick_pow(pmag, 2) * quick_pow(x, 10) * f5 -
                 46189 * quick_pow(ppmag, 2) * quick_pow(x, 10) * f5 + 63 * quick_pow(pmag, 2) * f6 +
                 63 * quick_pow(ppmag, 2) * f6 + 1260 * pmag * ppmag * x * f6 -
                 3465 * quick_pow(pmag, 2) * quick_pow(x, 2) * f6 - 3465 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f6 -
                 18480 * pmag * ppmag * quick_pow(x, 3) * f6 + 30030 * quick_pow(pmag, 2) * quick_pow(x, 4) * f6 +
                 30030 * quick_pow(ppmag, 2) * quick_pow(x, 4) * f6 + 72072 * pmag * ppmag * quick_pow(x, 5) * f6 -
                 90090 * quick_pow(pmag, 2) * quick_pow(x, 6) * f6 -
                 90090 * quick_pow(ppmag, 2) * quick_pow(x, 6) * f6 - 102960 * pmag * ppmag * quick_pow(x, 7) * f6 +
                 109395 * quick_pow(pmag, 2) * quick_pow(x, 8) * f6 +
                 109395 * quick_pow(ppmag, 2) * quick_pow(x, 8) * f6 + 48620 * pmag * ppmag * quick_pow(x, 9) * f6 -
                 46189 * quick_pow(pmag, 2) * quick_pow(x, 10) * f6 -
                 46189 * quick_pow(ppmag, 2) * quick_pow(x, 10) * f6)) /
               2432.;
    }
    else if (s2.l == 10 && s1.l == 10 && s1.s == 0 && s1.j == 10)
    {
        return (ppi *
                (-63 + 3465 * quick_pow(x, 2) - 30030 * quick_pow(x, 4) + 90090 * quick_pow(x, 6) -
                 109395 * quick_pow(x, 8) + 46189 * quick_pow(x, 10)) *
                (f1 - 3 * f2 - quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                 quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 - quick_pow(pmag, 2) * f5 -
                 quick_pow(ppmag, 2) * f5 - 2 * pmag * ppmag * x * f5 - quick_pow(pmag, 2) * f6 -
                 quick_pow(ppmag, 2) * f6 + 2 * pmag * ppmag * x * f6)) /
               128.;
    }
    else if (s2.l == 10 && s1.l == 10 && s1.s == 1 && s1.j == 10)
    {
        return (ppi * (18 * quick_pow(ppmag, 2) * quick_pow(-1 + quick_pow(x, 2), 2) *
                           (7 - 364 * quick_pow(x, 2) + 2730 * quick_pow(x, 4) - 6188 * quick_pow(x, 6) +
                            4199 * quick_pow(x, 8)) *
                           (quick_pow(pmag, 2) * f4 - f5 - f6) +
                       2 *
                           (-63 + 3465 * quick_pow(x, 2) - 30030 * quick_pow(x, 4) + 90090 * quick_pow(x, 6) -
                            109395 * quick_pow(x, 8) + 46189 * quick_pow(x, 10)) *
                           (f1 + f2 + quick_pow(pmag + ppmag * x, 2) * f5 + quick_pow(pmag - ppmag * x, 2) * f6) +
                       4 * ppmag * x * (1 - quick_pow(x, 2)) *
                           (63 - 1092 * quick_pow(x, 2) + 4914 * quick_pow(x, 4) - 7956 * quick_pow(x, 6) +
                            4199 * quick_pow(x, 8)) *
                           (-(pmag * (f3 - f5 + f6)) + ppmag * x * (f5 + f6)))) /
               256.;
    }
    else if (s2.l == 9 && s1.l == 9 && s1.s == 1 && s1.j == 10)
    {
        return -0.000744047619047619 *
               (ppi * (209066 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 11) * f4 -
                       63 * pmag * ppmag * (21 * f3 - f5 + f6) +
                       18018 * pmag * ppmag * quick_pow(x, 6) * (63 * f3 - 5 * f5 + 5 * f6) -
                       63 * x *
                           (105 * f1 + 105 * f2 + 94 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            5 * quick_pow(pmag, 2) * f5 + 5 * quick_pow(ppmag, 2) * f5 + 5 * quick_pow(pmag, 2) * f6 +
                            5 * quick_pow(ppmag, 2) * f6) +
                       5148 * quick_pow(x, 7) *
                           (105 * f1 + 105 * f2 + 151 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            5 * quick_pow(pmag, 2) * f5 + 5 * quick_pow(ppmag, 2) * f5 + 5 * quick_pow(pmag, 2) * f6 +
                            5 * quick_pow(ppmag, 2) * f6) -
                       715 * quick_pow(x, 9) *
                           (357 * f1 + 357 * f2 + 926 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            17 * quick_pow(pmag, 2) * f5 + 17 * quick_pow(ppmag, 2) * f5 +
                            17 * quick_pow(pmag, 2) * f6 + 17 * quick_pow(ppmag, 2) * f6) +
                       210 * quick_pow(x, 3) *
                           (462 * f1 + 462 * f2 + 433 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            22 * quick_pow(pmag, 2) * f5 + 22 * quick_pow(ppmag, 2) * f5 +
                            22 * quick_pow(pmag, 2) * f6 + 22 * quick_pow(ppmag, 2) * f6) -
                       462 * quick_pow(x, 5) *
                           (819 * f1 + 819 * f2 + 886 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            39 * quick_pow(pmag, 2) * f5 + 39 * quick_pow(ppmag, 2) * f5 +
                            39 * quick_pow(pmag, 2) * f6 + 39 * quick_pow(ppmag, 2) * f6) +
                       315 * pmag * ppmag * quick_pow(x, 2) * (189 * f3 + 11 * (-f5 + f6)) -
                       2310 * pmag * ppmag * quick_pow(x, 4) * (189 * f3 + 13 * (-f5 + f6)) -
                       6435 * pmag * ppmag * quick_pow(x, 8) * (189 * f3 + 17 * (-f5 + f6)) +
                       2431 * pmag * ppmag * quick_pow(x, 10) * (189 * f3 + 19 * (-f5 + f6))));
    }
    else if (s2.l == 9 && s1.l == 11 && s1.s == 1 && s1.j == 10)
    {
        return -0.001488095238095238 *
               (sqrt(27.5) * ppi *
                (2 * pmag * ppmag *
                     (-63 + 3465 * quick_pow(x, 2) - 30030 * quick_pow(x, 4) + 90090 * quick_pow(x, 6) -
                      109395 * quick_pow(x, 8) + 46189 * quick_pow(x, 10)) *
                     (f5 - f6) +
                 quick_pow(ppmag, 2) * x *
                     (-693 + 15015 * quick_pow(x, 2) - 90090 * quick_pow(x, 4) + 218790 * quick_pow(x, 6) -
                      230945 * quick_pow(x, 8) + 88179 * quick_pow(x, 10)) *
                     (f5 + f6) +
                 quick_pow(pmag, 2) * x *
                     (quick_pow(ppmag, 2) *
                          (-63 + 1155 * quick_pow(x, 2) - 6006 * quick_pow(x, 4) + 12870 * quick_pow(x, 6) -
                           12155 * quick_pow(x, 8) + 4199 * quick_pow(x, 10)) *
                          f4 +
                      2 *
                          (315 - 4620 * quick_pow(x, 2) + 18018 * quick_pow(x, 4) - 25740 * quick_pow(x, 6) +
                           12155 * quick_pow(x, 8)) *
                          (f5 + f6))));
    }
    else if (s2.l == 11 && s1.l == 9 && s1.s == 1 && s1.j == 10)
    {
        return -0.001488095238095238 *
               (sqrt(27.5) * ppi *
                (2 * pmag * ppmag *
                     (-63 + 3465 * quick_pow(x, 2) - 30030 * quick_pow(x, 4) + 90090 * quick_pow(x, 6) -
                      109395 * quick_pow(x, 8) + 46189 * quick_pow(x, 10)) *
                     (f5 - f6) +
                 2 * quick_pow(ppmag, 2) * x *
                     (315 - 4620 * quick_pow(x, 2) + 18018 * quick_pow(x, 4) - 25740 * quick_pow(x, 6) +
                      12155 * quick_pow(x, 8)) *
                     (f5 + f6) +
                 quick_pow(pmag, 2) * x *
                     (quick_pow(ppmag, 2) *
                          (-63 + 1155 * quick_pow(x, 2) - 6006 * quick_pow(x, 4) + 12870 * quick_pow(x, 6) -
                           12155 * quick_pow(x, 8) + 4199 * quick_pow(x, 10)) *
                          f4 +
                      (-693 + 15015 * quick_pow(x, 2) - 90090 * quick_pow(x, 4) + 218790 * quick_pow(x, 6) -
                       230945 * quick_pow(x, 8) + 88179 * quick_pow(x, 10)) *
                          (f5 + f6))));
    }
    else if (s2.l == 11 && s1.l == 11 && s1.s == 1 && s1.j == 10)
    {
        return -0.0003720238095238095 *
               (ppi * (-3703518 * pmag * ppmag * quick_pow(x, 12) * f3 +
                       1851759 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 13) * f4 +
                       1386 * pmag * ppmag * quick_pow(x, 2) * (126 * f3 + 5 * f5 - 5 * f6) -
                       30030 * pmag * ppmag * quick_pow(x, 4) * (63 * f3 + 2 * f5 - 2 * f6) -
                       126 * pmag * ppmag * (21 * f3 + f5 - f6) +
                       180180 * pmag * ppmag * quick_pow(x, 6) * (42 * f3 + f5 - f6) -
                       218790 * pmag * ppmag * quick_pow(x, 8) * (63 * f3 + f5 - f6) +
                       92378 * pmag * ppmag * quick_pow(x, 10) * (126 * f3 + f5 - f6) -
                       231 * quick_pow(x, 3) *
                           (1365 * f1 + 1365 * f2 + 1528 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            65 * quick_pow(pmag, 2) * f5 - 65 * quick_pow(ppmag, 2) * f5 -
                            65 * quick_pow(pmag, 2) * f6 - 65 * quick_pow(ppmag, 2) * f6) -
                       4199 * quick_pow(x, 11) *
                           (441 * f1 + 441 * f2 + 1616 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            21 * quick_pow(pmag, 2) * f5 - 21 * quick_pow(ppmag, 2) * f5 -
                            21 * quick_pow(pmag, 2) * f6 - 21 * quick_pow(ppmag, 2) * f6) +
                       12155 * quick_pow(x, 9) *
                           (399 * f1 + 399 * f2 + 797 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            19 * quick_pow(pmag, 2) * f5 - 19 * quick_pow(ppmag, 2) * f5 -
                            19 * quick_pow(pmag, 2) * f6 - 19 * quick_pow(ppmag, 2) * f6) -
                       12870 * quick_pow(x, 7) *
                           (357 * f1 + 357 * f2 + 524 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            17 * quick_pow(pmag, 2) * f5 - 17 * quick_pow(ppmag, 2) * f5 -
                            17 * quick_pow(pmag, 2) * f6 - 17 * quick_pow(ppmag, 2) * f6) +
                       63 * x *
                           (231 * f1 + 231 * f2 + 251 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            11 * quick_pow(pmag, 2) * f5 - 11 * quick_pow(ppmag, 2) * f5 -
                            11 * quick_pow(pmag, 2) * f6 - 11 * quick_pow(ppmag, 2) * f6) +
                       15015 * quick_pow(x, 5) *
                           (126 * f1 + 126 * f2 + 155 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            6 * quick_pow(pmag, 2) * f5 - 6 * quick_pow(ppmag, 2) * f5 - 6 * quick_pow(pmag, 2) * f6 -
                            6 * quick_pow(ppmag, 2) * f6)));
    }
    else if (s2.l == 11 && s1.l == 11 && s1.s == 0 && s1.j == 11)
    {
        return (ppi * x *
                (-693 + 15015 * quick_pow(x, 2) - 90090 * quick_pow(x, 4) + 218790 * quick_pow(x, 6) -
                 230945 * quick_pow(x, 8) + 88179 * quick_pow(x, 10)) *
                (f1 - 3 * f2 - quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                 quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 - quick_pow(pmag, 2) * f5 -
                 quick_pow(ppmag, 2) * f5 - 2 * pmag * ppmag * x * f5 - quick_pow(pmag, 2) * f6 -
                 quick_pow(ppmag, 2) * f6 + 2 * pmag * ppmag * x * f6)) /
               128.;
    }
    else if (s2.l == 11 && s1.l == 11 && s1.s == 1 && s1.j == 11)
    {
        return (ppi * (146965 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 13) * f4 +
                       29393 * pmag * ppmag * quick_pow(x, 12) * (f3 + 11 * (f5 - f6)) -
                       92378 * pmag * ppmag * quick_pow(x, 10) * (f3 + 9 * f5 - 9 * f6) +
                       109395 * pmag * ppmag * quick_pow(x, 8) * (f3 + 7 * f5 - 7 * f6) -
                       60060 * pmag * ppmag * quick_pow(x, 6) * (f3 + 5 * f5 - 5 * f6) +
                       15015 * pmag * ppmag * quick_pow(x, 4) * (f3 + 3 * f5 - 3 * f6) -
                       1386 * pmag * ppmag * quick_pow(x, 2) * (f3 + f5 - f6) + 21 * pmag * ppmag * (f3 - f5 + f6) +
                       30030 * quick_pow(x, 3) *
                           (f1 + f2 - quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 + quick_pow(pmag, 2) * f5 +
                            quick_pow(ppmag, 2) * f5 + quick_pow(pmag, 2) * f6 + quick_pow(ppmag, 2) * f6) -
                       15015 * quick_pow(x, 5) *
                           (12 * f1 + 12 * f2 - 13 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            12 * quick_pow(pmag, 2) * f5 + 12 * quick_pow(ppmag, 2) * f5 +
                            12 * quick_pow(pmag, 2) * f6 + 12 * quick_pow(ppmag, 2) * f6) +
                       8398 * quick_pow(x, 11) *
                           (21 * f1 + 21 * f2 - 65 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            21 * quick_pow(pmag, 2) * f5 + 21 * quick_pow(ppmag, 2) * f5 +
                            21 * quick_pow(pmag, 2) * f6 + 21 * quick_pow(ppmag, 2) * f6) -
                       12155 * quick_pow(x, 9) *
                           (38 * f1 + 38 * f2 - 65 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            38 * quick_pow(pmag, 2) * f5 + 38 * quick_pow(ppmag, 2) * f5 +
                            38 * quick_pow(pmag, 2) * f6 + 38 * quick_pow(ppmag, 2) * f6) +
                       8580 * quick_pow(x, 7) *
                           (51 * f1 + 51 * f2 - 65 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            51 * quick_pow(pmag, 2) * f5 + 51 * quick_pow(ppmag, 2) * f5 +
                            51 * quick_pow(pmag, 2) * f6 + 51 * quick_pow(ppmag, 2) * f6) -
                       21 * x *
                           (66 * f1 + 66 * f2 - 65 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                            66 * quick_pow(pmag, 2) * f5 + 66 * quick_pow(ppmag, 2) * f5 +
                            66 * quick_pow(pmag, 2) * f6 + 66 * quick_pow(ppmag, 2) * f6))) /
               256.;
    }
    else if (s2.l == 10 && s1.l == 10 && s1.s == 1 && s1.j == 11)
    {
        return (ppi *
                (23 *
                     (-63 + 3465 * quick_pow(x, 2) - 30030 * quick_pow(x, 4) + 90090 * quick_pow(x, 6) -
                      109395 * quick_pow(x, 8) + 46189 * quick_pow(x, 10)) *
                     f1 +
                 23 *
                     (-63 + 3465 * quick_pow(x, 2) - 30030 * quick_pow(x, 4) + 90090 * quick_pow(x, 6) -
                      109395 * quick_pow(x, 8) + 46189 * quick_pow(x, 10)) *
                     f2 +
                 28980 * pmag * ppmag * x * f3 - 531300 * pmag * ppmag * quick_pow(x, 3) * f3 +
                 2762760 * pmag * ppmag * quick_pow(x, 5) * f3 - 5920200 * pmag * ppmag * quick_pow(x, 7) * f3 +
                 5591300 * pmag * ppmag * quick_pow(x, 9) * f3 - 1931540 * pmag * ppmag * quick_pow(x, 11) * f3 -
                 1323 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                 72828 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 -
                 680295 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 4) * f4 +
                 2402400 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 6) * f4 -
                 3931785 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 8) * f4 +
                 3024164 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 10) * f4 -
                 885989 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 12) * f4 -
                 63 * quick_pow(pmag, 2) * f5 - 63 * quick_pow(ppmag, 2) * f5 - 1386 * pmag * ppmag * x * f5 +
                 3465 * quick_pow(pmag, 2) * quick_pow(x, 2) * f5 + 3465 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f5 +
                 30030 * pmag * ppmag * quick_pow(x, 3) * f5 - 30030 * quick_pow(pmag, 2) * quick_pow(x, 4) * f5 -
                 30030 * quick_pow(ppmag, 2) * quick_pow(x, 4) * f5 - 180180 * pmag * ppmag * quick_pow(x, 5) * f5 +
                 90090 * quick_pow(pmag, 2) * quick_pow(x, 6) * f5 +
                 90090 * quick_pow(ppmag, 2) * quick_pow(x, 6) * f5 + 437580 * pmag * ppmag * quick_pow(x, 7) * f5 -
                 109395 * quick_pow(pmag, 2) * quick_pow(x, 8) * f5 -
                 109395 * quick_pow(ppmag, 2) * quick_pow(x, 8) * f5 - 461890 * pmag * ppmag * quick_pow(x, 9) * f5 +
                 46189 * quick_pow(pmag, 2) * quick_pow(x, 10) * f5 +
                 46189 * quick_pow(ppmag, 2) * quick_pow(x, 10) * f5 + 176358 * pmag * ppmag * quick_pow(x, 11) * f5 -
                 63 * quick_pow(pmag, 2) * f6 - 63 * quick_pow(ppmag, 2) * f6 + 1386 * pmag * ppmag * x * f6 +
                 3465 * quick_pow(pmag, 2) * quick_pow(x, 2) * f6 + 3465 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f6 -
                 30030 * pmag * ppmag * quick_pow(x, 3) * f6 - 30030 * quick_pow(pmag, 2) * quick_pow(x, 4) * f6 -
                 30030 * quick_pow(ppmag, 2) * quick_pow(x, 4) * f6 + 180180 * pmag * ppmag * quick_pow(x, 5) * f6 +
                 90090 * quick_pow(pmag, 2) * quick_pow(x, 6) * f6 +
                 90090 * quick_pow(ppmag, 2) * quick_pow(x, 6) * f6 - 437580 * pmag * ppmag * quick_pow(x, 7) * f6 -
                 109395 * quick_pow(pmag, 2) * quick_pow(x, 8) * f6 -
                 109395 * quick_pow(ppmag, 2) * quick_pow(x, 8) * f6 + 461890 * pmag * ppmag * quick_pow(x, 9) * f6 +
                 46189 * quick_pow(pmag, 2) * quick_pow(x, 10) * f6 +
                 46189 * quick_pow(ppmag, 2) * quick_pow(x, 10) * f6 - 176358 * pmag * ppmag * quick_pow(x, 11) * f6)) /
               2944.;
    }
    else if (s2.l == 10 && s1.l == 12 && s1.s == 1 && s1.j == 11)
    {
        return -0.00033967391304347825 *
               (sqrt(33) * ppi *
                (8 * pmag * ppmag * x *
                     (-693 + 15015 * quick_pow(x, 2) - 90090 * quick_pow(x, 4) + 218790 * quick_pow(x, 6) -
                      230945 * quick_pow(x, 8) + 88179 * quick_pow(x, 10)) *
                     (f5 - f6) +
                 quick_pow(ppmag, 2) *
                     (231 - 18018 * quick_pow(x, 2) + 225225 * quick_pow(x, 4) - 1021020 * quick_pow(x, 6) +
                      2078505 * quick_pow(x, 8) - 1939938 * quick_pow(x, 10) + 676039 * quick_pow(x, 12)) *
                     (f5 + f6) +
                 quick_pow(pmag, 2) *
                     (quick_pow(ppmag, 2) *
                          (21 - 1386 * quick_pow(x, 2) + 15015 * quick_pow(x, 4) - 60060 * quick_pow(x, 6) +
                           109395 * quick_pow(x, 8) - 92378 * quick_pow(x, 10) + 29393 * quick_pow(x, 12)) *
                          f4 +
                      4 *
                          (-63 + 3465 * quick_pow(x, 2) - 30030 * quick_pow(x, 4) + 90090 * quick_pow(x, 6) -
                           109395 * quick_pow(x, 8) + 46189 * quick_pow(x, 10)) *
                          (f5 + f6))));
    }
    else if (s2.l == 12 && s1.l == 10 && s1.s == 1 && s1.j == 11)
    {
        return -0.00033967391304347825 *
               (sqrt(33) * ppi *
                (8 * pmag * ppmag * x *
                     (-693 + 15015 * quick_pow(x, 2) - 90090 * quick_pow(x, 4) + 218790 * quick_pow(x, 6) -
                      230945 * quick_pow(x, 8) + 88179 * quick_pow(x, 10)) *
                     (f5 - f6) +
                 4 * quick_pow(ppmag, 2) *
                     (-63 + 3465 * quick_pow(x, 2) - 30030 * quick_pow(x, 4) + 90090 * quick_pow(x, 6) -
                      109395 * quick_pow(x, 8) + 46189 * quick_pow(x, 10)) *
                     (f5 + f6) +
                 quick_pow(pmag, 2) *
                     (quick_pow(ppmag, 2) *
                          (21 - 1386 * quick_pow(x, 2) + 15015 * quick_pow(x, 4) - 60060 * quick_pow(x, 6) +
                           109395 * quick_pow(x, 8) - 92378 * quick_pow(x, 10) + 29393 * quick_pow(x, 12)) *
                          f4 +
                      (231 - 18018 * quick_pow(x, 2) + 225225 * quick_pow(x, 4) - 1021020 * quick_pow(x, 6) +
                       2078505 * quick_pow(x, 8) - 1939938 * quick_pow(x, 10) + 676039 * quick_pow(x, 12)) *
                          (f5 + f6))));
    }
    else if (s2.l == 12 && s1.l == 12 && s1.s == 1 && s1.j == 11)
    {
        return (ppi *
                (23 *
                     (231 - 18018 * quick_pow(x, 2) + 225225 * quick_pow(x, 4) - 1021020 * quick_pow(x, 6) +
                      2078505 * quick_pow(x, 8) - 1939938 * quick_pow(x, 10) + 676039 * quick_pow(x, 12)) *
                     f1 +
                 23 *
                     (231 - 18018 * quick_pow(x, 2) + 225225 * quick_pow(x, 4) - 1021020 * quick_pow(x, 6) +
                      2078505 * quick_pow(x, 8) - 1939938 * quick_pow(x, 10) + 676039 * quick_pow(x, 12)) *
                     f2 +
                 138138 * pmag * ppmag * x * f3 - 3591588 * pmag * ppmag * quick_pow(x, 3) * f3 +
                 26936910 * pmag * ppmag * quick_pow(x, 5) * f3 - 87224280 * pmag * ppmag * quick_pow(x, 7) * f3 +
                 138105110 * pmag * ppmag * quick_pow(x, 9) * f3 - 105462084 * pmag * ppmag * quick_pow(x, 11) * f3 +
                 31097794 * pmag * ppmag * quick_pow(x, 13) * f3 +
                 5775 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                 450219 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 +
                 5924919 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 4) * f4 -
                 29984955 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 6) * f4 +
                 73695765 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 8) * f4 -
                 94456505 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 10) * f4 +
                 60814117 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 12) * f4 -
                 15548897 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 14) * f4 -
                 231 * quick_pow(pmag, 2) * f5 - 231 * quick_pow(ppmag, 2) * f5 + 5544 * pmag * ppmag * x * f5 +
                 18018 * quick_pow(pmag, 2) * quick_pow(x, 2) * f5 +
                 18018 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f5 - 120120 * pmag * ppmag * quick_pow(x, 3) * f5 -
                 225225 * quick_pow(pmag, 2) * quick_pow(x, 4) * f5 -
                 225225 * quick_pow(ppmag, 2) * quick_pow(x, 4) * f5 + 720720 * pmag * ppmag * quick_pow(x, 5) * f5 +
                 1021020 * quick_pow(pmag, 2) * quick_pow(x, 6) * f5 +
                 1021020 * quick_pow(ppmag, 2) * quick_pow(x, 6) * f5 - 1750320 * pmag * ppmag * quick_pow(x, 7) * f5 -
                 2078505 * quick_pow(pmag, 2) * quick_pow(x, 8) * f5 -
                 2078505 * quick_pow(ppmag, 2) * quick_pow(x, 8) * f5 + 1847560 * pmag * ppmag * quick_pow(x, 9) * f5 +
                 1939938 * quick_pow(pmag, 2) * quick_pow(x, 10) * f5 +
                 1939938 * quick_pow(ppmag, 2) * quick_pow(x, 10) * f5 - 705432 * pmag * ppmag * quick_pow(x, 11) * f5 -
                 676039 * quick_pow(pmag, 2) * quick_pow(x, 12) * f5 -
                 676039 * quick_pow(ppmag, 2) * quick_pow(x, 12) * f5 - 231 * quick_pow(pmag, 2) * f6 -
                 231 * quick_pow(ppmag, 2) * f6 - 5544 * pmag * ppmag * x * f6 +
                 18018 * quick_pow(pmag, 2) * quick_pow(x, 2) * f6 +
                 18018 * quick_pow(ppmag, 2) * quick_pow(x, 2) * f6 + 120120 * pmag * ppmag * quick_pow(x, 3) * f6 -
                 225225 * quick_pow(pmag, 2) * quick_pow(x, 4) * f6 -
                 225225 * quick_pow(ppmag, 2) * quick_pow(x, 4) * f6 - 720720 * pmag * ppmag * quick_pow(x, 5) * f6 +
                 1021020 * quick_pow(pmag, 2) * quick_pow(x, 6) * f6 +
                 1021020 * quick_pow(ppmag, 2) * quick_pow(x, 6) * f6 + 1750320 * pmag * ppmag * quick_pow(x, 7) * f6 -
                 2078505 * quick_pow(pmag, 2) * quick_pow(x, 8) * f6 -
                 2078505 * quick_pow(ppmag, 2) * quick_pow(x, 8) * f6 - 1847560 * pmag * ppmag * quick_pow(x, 9) * f6 +
                 1939938 * quick_pow(pmag, 2) * quick_pow(x, 10) * f6 +
                 1939938 * quick_pow(ppmag, 2) * quick_pow(x, 10) * f6 + 705432 * pmag * ppmag * quick_pow(x, 11) * f6 -
                 676039 * quick_pow(pmag, 2) * quick_pow(x, 12) * f6 -
                 676039 * quick_pow(ppmag, 2) * quick_pow(x, 12) * f6)) /
               11776.;
    }
    else if (s2.l == 12 && s1.l == 12 && s1.s == 0 && s1.j == 12)
    {
        return (ppi *
                (231 - 18018 * quick_pow(x, 2) + 225225 * quick_pow(x, 4) - 1021020 * quick_pow(x, 6) +
                 2078505 * quick_pow(x, 8) - 1939938 * quick_pow(x, 10) + 676039 * quick_pow(x, 12)) *
                (f1 - 3 * f2 - quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                 quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 2) * f4 - quick_pow(pmag, 2) * f5 -
                 quick_pow(ppmag, 2) * f5 - 2 * pmag * ppmag * x * f5 - quick_pow(pmag, 2) * f6 -
                 quick_pow(ppmag, 2) * f6 + 2 * pmag * ppmag * x * f6)) /
               512.;
    }
    else if (s2.l == 12 && s1.l == 12 && s1.s == 1 && s1.j == 12)
    {
        return (ppi * (154 * quick_pow(ppmag, 2) * quick_pow(-1 + quick_pow(x, 2), 2) *
                           (-3 + 225 * quick_pow(x, 2) - 2550 * quick_pow(x, 4) + 9690 * quick_pow(x, 6) -
                            14535 * quick_pow(x, 8) + 7429 * quick_pow(x, 10)) *
                           (quick_pow(pmag, 2) * f4 - f5 - f6) +
                       2 *
                           (231 - 18018 * quick_pow(x, 2) + 225225 * quick_pow(x, 4) - 1021020 * quick_pow(x, 6) +
                            2078505 * quick_pow(x, 8) - 1939938 * quick_pow(x, 10) + 676039 * quick_pow(x, 12)) *
                           (f1 + f2 + quick_pow(pmag + ppmag * x, 2) * f5 + quick_pow(pmag - ppmag * x, 2) * f6) +
                       4 * ppmag * x * (1 - quick_pow(x, 2)) *
                           (-231 + 5775 * quick_pow(x, 2) - 39270 * quick_pow(x, 4) + 106590 * quick_pow(x, 6) -
                            124355 * quick_pow(x, 8) + 52003 * quick_pow(x, 10)) *
                           (-(pmag * (f3 - f5 + f6)) + ppmag * x * (f5 + f6)))) /
               1024.;
    }
    else if (s2.l == 11 && s1.l == 11 && s1.s == 1 && s1.j == 12)
    {
        return -0.00015625 * (ppi * (3732911 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 13) * f4 +
                                     231 * pmag * ppmag * (25 * f3 - f5 + f6) +
                                     75075 * pmag * ppmag * quick_pow(x, 4) * (55 * f3 - 3 * f5 + 3 * f6) +
                                     231 * x *
                                         (150 * f1 + 150 * f2 + 137 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                                          6 * quick_pow(pmag, 2) * f5 + 6 * quick_pow(ppmag, 2) * f5 +
                                          6 * quick_pow(pmag, 2) * f6 + 6 * quick_pow(ppmag, 2) * f6) +
                                     15015 * quick_pow(x, 5) *
                                         (300 * f1 + 300 * f2 + 311 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                                          12 * quick_pow(pmag, 2) * f5 + 12 * quick_pow(ppmag, 2) * f5 +
                                          12 * quick_pow(pmag, 2) * f6 + 12 * quick_pow(ppmag, 2) * f6) -
                                     8398 * quick_pow(x, 11) *
                                         (525 * f1 + 525 * f2 + 1627 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                                          21 * quick_pow(pmag, 2) * f5 + 21 * quick_pow(ppmag, 2) * f5 +
                                          21 * quick_pow(pmag, 2) * f6 + 21 * quick_pow(ppmag, 2) * f6) +
                                     12155 * quick_pow(x, 9) *
                                         (950 * f1 + 950 * f2 + 1603 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                                          38 * quick_pow(pmag, 2) * f5 + 38 * quick_pow(ppmag, 2) * f5 +
                                          38 * quick_pow(pmag, 2) * f6 + 38 * quick_pow(ppmag, 2) * f6) -
                                     8580 * quick_pow(x, 7) *
                                         (1275 * f1 + 1275 * f2 + 1579 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                                          51 * quick_pow(pmag, 2) * f5 + 51 * quick_pow(ppmag, 2) * f5 +
                                          51 * quick_pow(pmag, 2) * f6 + 51 * quick_pow(ppmag, 2) * f6) -
                                     462 * quick_pow(x, 3) *
                                         (1625 * f1 + 1625 * f2 + 1531 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 +
                                          65 * quick_pow(pmag, 2) * f5 + 65 * quick_pow(ppmag, 2) * f5 +
                                          65 * quick_pow(pmag, 2) * f6 + 65 * quick_pow(ppmag, 2) * f6) -
                                     1386 * pmag * ppmag * quick_pow(x, 2) * (275 * f3 + 13 * (-f5 + f6)) -
                                     60060 * pmag * ppmag * quick_pow(x, 6) * (275 * f3 + 17 * (-f5 + f6)) +
                                     109395 * pmag * ppmag * quick_pow(x, 8) * (275 * f3 + 19 * (-f5 + f6)) -
                                     92378 * pmag * ppmag * quick_pow(x, 10) * (275 * f3 + 21 * (-f5 + f6)) +
                                     29393 * pmag * ppmag * quick_pow(x, 12) * (275 * f3 + 23 * (-f5 + f6))));
    }
    else if (s2.l == 11 && s1.l == 13 && s1.s == 1 && s1.j == 12)
    {
        return -0.0003125 *
               (sqrt(39) * ppi *
                (2 * pmag * ppmag *
                     (231 - 18018 * quick_pow(x, 2) + 225225 * quick_pow(x, 4) - 1021020 * quick_pow(x, 6) +
                      2078505 * quick_pow(x, 8) - 1939938 * quick_pow(x, 10) + 676039 * quick_pow(x, 12)) *
                     (f5 - f6) +
                 quick_pow(ppmag, 2) * x *
                     (3003 - 90090 * quick_pow(x, 2) + 765765 * quick_pow(x, 4) - 2771340 * quick_pow(x, 6) +
                      4849845 * quick_pow(x, 8) - 4056234 * quick_pow(x, 10) + 1300075 * quick_pow(x, 12)) *
                     (f5 + f6) +
                 quick_pow(pmag, 2) * x *
                     (quick_pow(ppmag, 2) *
                          (231 - 6006 * quick_pow(x, 2) + 45045 * quick_pow(x, 4) - 145860 * quick_pow(x, 6) +
                           230945 * quick_pow(x, 8) - 176358 * quick_pow(x, 10) + 52003 * quick_pow(x, 12)) *
                          f4 +
                      4 *
                          (-693 + 15015 * quick_pow(x, 2) - 90090 * quick_pow(x, 4) + 218790 * quick_pow(x, 6) -
                           230945 * quick_pow(x, 8) + 88179 * quick_pow(x, 10)) *
                          (f5 + f6))));
    }
    else if (s2.l == 13 && s1.l == 11 && s1.s == 1 && s1.j == 12)
    {
        return -0.0003125 *
               (sqrt(39) * ppi *
                (2 * pmag * ppmag *
                     (231 - 18018 * quick_pow(x, 2) + 225225 * quick_pow(x, 4) - 1021020 * quick_pow(x, 6) +
                      2078505 * quick_pow(x, 8) - 1939938 * quick_pow(x, 10) + 676039 * quick_pow(x, 12)) *
                     (f5 - f6) +
                 4 * quick_pow(ppmag, 2) * x *
                     (-693 + 15015 * quick_pow(x, 2) - 90090 * quick_pow(x, 4) + 218790 * quick_pow(x, 6) -
                      230945 * quick_pow(x, 8) + 88179 * quick_pow(x, 10)) *
                     (f5 + f6) +
                 quick_pow(pmag, 2) * x *
                     (quick_pow(ppmag, 2) *
                          (231 - 6006 * quick_pow(x, 2) + 45045 * quick_pow(x, 4) - 145860 * quick_pow(x, 6) +
                           230945 * quick_pow(x, 8) - 176358 * quick_pow(x, 10) + 52003 * quick_pow(x, 12)) *
                          f4 +
                      (3003 - 90090 * quick_pow(x, 2) + 765765 * quick_pow(x, 4) - 2771340 * quick_pow(x, 6) +
                       4849845 * quick_pow(x, 8) - 4056234 * quick_pow(x, 10) + 1300075 * quick_pow(x, 12)) *
                          (f5 + f6))));
    }
    else if (s2.l == 13 && s1.l == 13 && s1.s == 1 && s1.j == 12)
    {
        return -0.000078125 *
               (ppi * (-65003750 * pmag * ppmag * quick_pow(x, 14) * f3 +
                       32501875 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * quick_pow(x, 15) * f4 -
                       6006 * pmag * ppmag * quick_pow(x, 2) * (175 * f3 + 6 * f5 - 6 * f6) -
                       510510 * pmag * ppmag * quick_pow(x, 6) * (175 * f3 + 4 * f5 - 4 * f6) +
                       1385670 * pmag * ppmag * quick_pow(x, 8) * (175 * f3 + 3 * f5 - 3 * f6) -
                       1939938 * pmag * ppmag * quick_pow(x, 10) * (175 * f3 + 2 * f5 - 2 * f6) +
                       462 * pmag * ppmag * (25 * f3 + f5 - f6) +
                       450450 * pmag * ppmag * quick_pow(x, 4) * (35 * f3 + f5 - f6) +
                       1352078 * pmag * ppmag * quick_pow(x, 12) * (175 * f3 + f5 - f6) +
                       36465 * quick_pow(x, 7) *
                           (1900 * f1 + 1900 * f2 + 2521 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            76 * quick_pow(pmag, 2) * f5 - 76 * quick_pow(ppmag, 2) * f5 -
                            76 * quick_pow(pmag, 2) * f6 - 76 * quick_pow(ppmag, 2) * f6) +
                       88179 * quick_pow(x, 11) *
                           (1150 * f1 + 1150 * f2 + 2573 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            46 * quick_pow(pmag, 2) * f5 - 46 * quick_pow(ppmag, 2) * f5 -
                            46 * quick_pow(pmag, 2) * f6 - 46 * quick_pow(ppmag, 2) * f6) +
                       3003 * quick_pow(x, 3) *
                           (750 * f1 + 750 * f2 + 823 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            30 * quick_pow(pmag, 2) * f5 - 30 * quick_pow(ppmag, 2) * f5 -
                            30 * quick_pow(pmag, 2) * f6 - 30 * quick_pow(ppmag, 2) * f6) -
                       52003 * quick_pow(x, 13) *
                           (625 * f1 + 625 * f2 + 2599 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            25 * quick_pow(pmag, 2) * f5 - 25 * quick_pow(ppmag, 2) * f5 -
                            25 * quick_pow(pmag, 2) * f6 - 25 * quick_pow(ppmag, 2) * f6) -
                       45045 * quick_pow(x, 5) *
                           (425 * f1 + 425 * f2 + 499 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            17 * quick_pow(pmag, 2) * f5 - 17 * quick_pow(ppmag, 2) * f5 -
                            17 * quick_pow(pmag, 2) * f6 - 17 * quick_pow(ppmag, 2) * f6) -
                       231 * x *
                           (325 * f1 + 325 * f2 + 349 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            13 * quick_pow(pmag, 2) * f5 - 13 * quick_pow(ppmag, 2) * f5 -
                            13 * quick_pow(pmag, 2) * f6 - 13 * quick_pow(ppmag, 2) * f6) -
                       692835 * quick_pow(x, 9) *
                           (175 * f1 + 175 * f2 + 283 * quick_pow(pmag, 2) * quick_pow(ppmag, 2) * f4 -
                            7 * quick_pow(pmag, 2) * f5 - 7 * quick_pow(ppmag, 2) * f5 - 7 * quick_pow(pmag, 2) * f6 -
                            7 * quick_pow(ppmag, 2) * f6)));
    }
    else
    {
        error_message_print_abort("J out of aPWD method!");
    }
    return 0;
}

// regulator function used to cut-off high momentum part in the L-S equation.
double regulator_function_order(double p1, double p2, int n, NN::NN_configs configs)
{
    double temp;
    temp = exp(-(quick_pow(p1, n) + quick_pow(p2, n)) / quick_pow(configs.Lambda, n));
    return temp;
}
//------------------------------------------------------------------------------------------------------------------

// loop integral function.
double loop_function_L(double q, NN::NN_configs configs)
{
    double temp;
    double lambda = configs.Lambda_tilde;
    double mpi = configs.mass_pion_averaged;
    double w = sqrt(4.0 * mpi * mpi + q * q);
    double num = lambda * lambda * (2.0 * mpi * mpi + q * q) - 2.0 * mpi * mpi * q * q +
                 lambda * sqrt(lambda * lambda - 4.0 * mpi * mpi) * q * w;
    double den = 2.0 * mpi * mpi * (lambda * lambda + q * q);
    double fac = w / (2.0 * q);
    temp = fac * log(num / den);
    return temp;
}

double loop_function_A(double q, NN::NN_configs configs)
{
    double temp;
    double lambda = configs.Lambda_tilde;
    double mpi = configs.mass_pion_averaged;
    double num = q * (lambda - 2.0 * mpi);
    double den = q * q + 2.0 * lambda * mpi;
    double fac = 1.0 / (2.0 * q);
    temp = fac * atan(num / den);
    return temp;
}

} // end namespace interaction_aPWD

#endif // INTERACTION_aPWD_HPP
