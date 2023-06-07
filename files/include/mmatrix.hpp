#pragma once
#ifndef MMATRIX_HPP
#define MMATRIX_HPP

#include "lib_define.hpp"

// partial wave decomposition for scattering amplitude,
// for which we have S and J conservation.
namespace mmatrix
{
constexpr std::complex<double> im_unit{0, 1};
constexpr double ppi = 3.14159265358979323846;
constexpr double sqrt2 = 1.414213562373095;
constexpr double sqrt3 = 1.732050807568877;
constexpr double sqrt5 = 2.236067977499790;
constexpr double sqrt7 = 2.645751311064591;
constexpr double sqrt15 = 3.872983346207417;
using basic_math::quick_pow;

// return m11 partial wave coefficients in calculating scattering amplituide.
// t: costheta.
double coeff_m11(double t, int lp, int l, int s, int j)
{
    double costheta = t;
    double sintheta = sqrt(1.0 - t * t);
    if (lp == 0 && l == 0 && s == 0 && j == 0)
    {
        return 0;
    }
    else if (lp == 1 && l == 1 && s == 1 && j == 0)
    {
        return 0;
    }
    else if (lp == 1 && l == 1 && s == 0 && j == 1)
    {
        return 0;
    }
    else if (lp == 1 && l == 1 && s == 1 && j == 1)
    {
        return (3 * costheta) / 2.;
    }
    else if (lp == 0 && l == 0 && s == 1 && j == 1)
    {
        return 1;
    }
    else if (lp == 0 && l == 2 && s == 1 && j == 1)
    {
        return -(1 / sqrt2);
    }
    else if (lp == 2 && l == 0 && s == 1 && j == 1)
    {
        return 1 / (2. * sqrt2) - (3 * quick_pow(costheta, 2)) / (2. * sqrt2);
    }
    else if (lp == 2 && l == 2 && s == 1 && j == 1)
    {
        return -0.25 + (3 * quick_pow(costheta, 2)) / 4.;
    }
    else if (lp == 2 && l == 2 && s == 0 && j == 2)
    {
        return 0;
    }
    else if (lp == 2 && l == 2 && s == 1 && j == 2)
    {
        return -1.25 + (15 * quick_pow(costheta, 2)) / 4.;
    }
    else if (lp == 1 && l == 1 && s == 1 && j == 2)
    {
        return (3 * costheta) / 2.;
    }
    else if (lp == 1 && l == 3 && s == 1 && j == 2)
    {
        return -(sqrt(1.5) * costheta);
    }
    else if (lp == 3 && l == 1 && s == 1 && j == 2)
    {
        return (-3 * sqrt(1.5) * costheta) / 8. - (5 * sqrt(1.5) * quick_pow(costheta, 3)) / 8. +
               (15 * sqrt(1.5) * costheta * quick_pow(sintheta, 2)) / 8.;
    }
    else if (lp == 3 && l == 3 && s == 1 && j == 2)
    {
        return (3 * costheta) / 8. + (5 * quick_pow(costheta, 3)) / 8. - (15 * costheta * quick_pow(sintheta, 2)) / 8.;
    }
    else if (lp == 3 && l == 3 && s == 0 && j == 3)
    {
        return 0;
    }
    else if (lp == 3 && l == 3 && s == 1 && j == 3)
    {
        return (21 * costheta) / 16. + (35 * quick_pow(costheta, 3)) / 16. -
               (105 * costheta * quick_pow(sintheta, 2)) / 16.;
    }
    else if (lp == 2 && l == 2 && s == 1 && j == 3)
    {
        return -1 + 3 * quick_pow(costheta, 2);
    }
    else if (lp == 2 && l == 4 && s == 1 && j == 3)
    {
        return sqrt3 / 2. - (3 * sqrt3 * quick_pow(costheta, 2)) / 2.;
    }
    else if (lp == 4 && l == 2 && s == 1 && j == 3)
    {
        return (-9 * sqrt3) / 64. - (5 * sqrt3 * quick_pow(costheta, 2)) / 16. -
               (35 * sqrt3 * quick_pow(costheta, 4)) / 64. + (5 * sqrt3 * quick_pow(sintheta, 2)) / 16. +
               (105 * sqrt3 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 32. -
               (35 * sqrt3 * quick_pow(sintheta, 4)) / 64.;
    }
    else if (lp == 4 && l == 4 && s == 1 && j == 3)
    {
        return 0.2109375 + (15 * quick_pow(costheta, 2)) / 32. + (105 * quick_pow(costheta, 4)) / 128. -
               (15 * quick_pow(sintheta, 2)) / 32. - (315 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 64. +
               (105 * quick_pow(sintheta, 4)) / 128.;
    }
    else if (lp == 4 && l == 4 && s == 0 && j == 4)
    {
        return 0;
    }
    else if (lp == 4 && l == 4 && s == 1 && j == 4)
    {
        return 0.6328125 + (45 * quick_pow(costheta, 2)) / 32. + (315 * quick_pow(costheta, 4)) / 128. -
               (45 * quick_pow(sintheta, 2)) / 32. - (945 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 64. +
               (315 * quick_pow(sintheta, 4)) / 128.;
    }
    else if (lp == 3 && l == 3 && s == 1 && j == 4)
    {
        return (15 * costheta) / 16. + (25 * quick_pow(costheta, 3)) / 16. -
               (75 * costheta * quick_pow(sintheta, 2)) / 16.;
    }
    else if (lp == 3 && l == 5 && s == 1 && j == 4)
    {
        return (-3 * sqrt5 * costheta) / 8. - (5 * sqrt5 * quick_pow(costheta, 3)) / 8. +
               (15 * sqrt5 * costheta * quick_pow(sintheta, 2)) / 8.;
    }
    else if (lp == 5 && l == 3 && s == 1 && j == 4)
    {
        return (-15 * sqrt5 * costheta) / 64. - (35 * sqrt5 * quick_pow(costheta, 3)) / 128. -
               (63 * sqrt5 * quick_pow(costheta, 5)) / 128. + (105 * sqrt5 * costheta * quick_pow(sintheta, 2)) / 128. +
               (315 * sqrt5 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 64. -
               (315 * sqrt5 * costheta * quick_pow(sintheta, 4)) / 128.;
    }
    else if (lp == 5 && l == 5 && s == 1 && j == 4)
    {
        return (15 * costheta) / 32. + (35 * quick_pow(costheta, 3)) / 64. + (63 * quick_pow(costheta, 5)) / 64. -
               (105 * costheta * quick_pow(sintheta, 2)) / 64. -
               (315 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 32. +
               (315 * costheta * quick_pow(sintheta, 4)) / 64.;
    }
    else if (lp == 5 && l == 5 && s == 0 && j == 5)
    {
        return 0;
    }
    else if (lp == 5 && l == 5 && s == 1 && j == 5)
    {
        return (165 * costheta) / 128. + (385 * quick_pow(costheta, 3)) / 256. + (693 * quick_pow(costheta, 5)) / 256. -
               (1155 * costheta * quick_pow(sintheta, 2)) / 256. -
               (3465 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 128. +
               (3465 * costheta * quick_pow(sintheta, 4)) / 256.;
    }
    else if (lp == 4 && l == 4 && s == 1 && j == 5)
    {
        return 0.421875 + (15 * quick_pow(costheta, 2)) / 16. + (105 * quick_pow(costheta, 4)) / 64. -
               (15 * quick_pow(sintheta, 2)) / 16. - (315 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 32. +
               (105 * quick_pow(sintheta, 4)) / 64.;
    }
    else if (lp == 4 && l == 6 && s == 1 && j == 5)
    {
        return (-9 * sqrt(7.5)) / 64. - (5 * sqrt(7.5) * quick_pow(costheta, 2)) / 16. -
               (35 * sqrt(7.5) * quick_pow(costheta, 4)) / 64. + (5 * sqrt(7.5) * quick_pow(sintheta, 2)) / 16. +
               (105 * sqrt(7.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 32. -
               (35 * sqrt(7.5) * quick_pow(sintheta, 4)) / 64.;
    }
    else if (lp == 6 && l == 4 && s == 1 && j == 5)
    {
        return (-25 * sqrt(7.5)) / 256. - (105 * sqrt(7.5) * quick_pow(costheta, 2)) / 512. -
               (63 * sqrt(7.5) * quick_pow(costheta, 4)) / 256. - (231 * sqrt(7.5) * quick_pow(costheta, 6)) / 512. +
               (105 * sqrt(7.5) * quick_pow(sintheta, 2)) / 512. +
               (189 * sqrt(7.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 128. +
               (3465 * sqrt(7.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 512. -
               (63 * sqrt(7.5) * quick_pow(sintheta, 4)) / 256. -
               (3465 * sqrt(7.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 512. +
               (231 * sqrt(7.5) * quick_pow(sintheta, 6)) / 512.;
    }
    else if (lp == 6 && l == 6 && s == 1 && j == 5)
    {
        return 0.244140625 + (525 * quick_pow(costheta, 2)) / 1024. + (315 * quick_pow(costheta, 4)) / 512. +
               (1155 * quick_pow(costheta, 6)) / 1024. - (525 * quick_pow(sintheta, 2)) / 1024. -
               (945 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 256. -
               (17325 * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 1024. +
               (315 * quick_pow(sintheta, 4)) / 512. +
               (17325 * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 1024. -
               (1155 * quick_pow(sintheta, 6)) / 1024.;
    }
    else if (lp == 6 && l == 6 && s == 0 && j == 6)
    {
        return 0;
    }
    else if (lp == 6 && l == 6 && s == 1 && j == 6)
    {
        return 0.634765625 + (1365 * quick_pow(costheta, 2)) / 1024. + (819 * quick_pow(costheta, 4)) / 512. +
               (3003 * quick_pow(costheta, 6)) / 1024. - (1365 * quick_pow(sintheta, 2)) / 1024. -
               (2457 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 256. -
               (45045 * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 1024. +
               (819 * quick_pow(sintheta, 4)) / 512. +
               (45045 * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 1024. -
               (3003 * quick_pow(sintheta, 6)) / 1024.;
    }
    else if (lp == 5 && l == 5 && s == 1 && j == 6)
    {
        return (105 * costheta) / 128. + (245 * quick_pow(costheta, 3)) / 256. + (441 * quick_pow(costheta, 5)) / 256. -
               (735 * costheta * quick_pow(sintheta, 2)) / 256. -
               (2205 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 128. +
               (2205 * costheta * quick_pow(sintheta, 4)) / 256.;
    }
    else if (lp == 5 && l == 7 && s == 1 && j == 6)
    {
        return (-15 * sqrt(10.5) * costheta) / 64. - (35 * sqrt(10.5) * quick_pow(costheta, 3)) / 128. -
               (63 * sqrt(10.5) * quick_pow(costheta, 5)) / 128. +
               (105 * sqrt(10.5) * costheta * quick_pow(sintheta, 2)) / 128. +
               (315 * sqrt(10.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 64. -
               (315 * sqrt(10.5) * costheta * quick_pow(sintheta, 4)) / 128.;
    }
    else if (lp == 7 && l == 5 && s == 1 && j == 6)
    {
        return (-175 * sqrt(10.5) * costheta) / 1024. - (189 * sqrt(10.5) * quick_pow(costheta, 3)) / 1024. -
               (231 * sqrt(10.5) * quick_pow(costheta, 5)) / 1024. -
               (429 * sqrt(10.5) * quick_pow(costheta, 7)) / 1024. +
               (567 * sqrt(10.5) * costheta * quick_pow(sintheta, 2)) / 1024. +
               (1155 * sqrt(10.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 512. +
               (9009 * sqrt(10.5) * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 1024. -
               (1155 * sqrt(10.5) * costheta * quick_pow(sintheta, 4)) / 1024. -
               (15015 * sqrt(10.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 1024. +
               (3003 * sqrt(10.5) * costheta * quick_pow(sintheta, 6)) / 1024.;
    }
    else if (lp == 7 && l == 7 && s == 1 && j == 6)
    {
        return (525 * costheta) / 1024. + (567 * quick_pow(costheta, 3)) / 1024. +
               (693 * quick_pow(costheta, 5)) / 1024. + (1287 * quick_pow(costheta, 7)) / 1024. -
               (1701 * costheta * quick_pow(sintheta, 2)) / 1024. -
               (3465 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 512. -
               (27027 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 1024. +
               (3465 * costheta * quick_pow(sintheta, 4)) / 1024. +
               (45045 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 1024. -
               (9009 * costheta * quick_pow(sintheta, 6)) / 1024.;
    }
    else if (lp == 7 && l == 7 && s == 0 && j == 7)
    {
        return 0;
    }
    else if (lp == 7 && l == 7 && s == 1 && j == 7)
    {
        return (2625 * costheta) / 2048. + (2835 * quick_pow(costheta, 3)) / 2048. +
               (3465 * quick_pow(costheta, 5)) / 2048. + (6435 * quick_pow(costheta, 7)) / 2048. -
               (8505 * costheta * quick_pow(sintheta, 2)) / 2048. -
               (17325 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 1024. -
               (135135 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 2048. +
               (17325 * costheta * quick_pow(sintheta, 4)) / 2048. +
               (225225 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 2048. -
               (45045 * costheta * quick_pow(sintheta, 6)) / 2048.;
    }
    else if (lp == 6 && l == 6 && s == 1 && j == 7)
    {
        return 0.390625 + (105 * quick_pow(costheta, 2)) / 128. + (63 * quick_pow(costheta, 4)) / 64. +
               (231 * quick_pow(costheta, 6)) / 128. - (105 * quick_pow(sintheta, 2)) / 128. -
               (189 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 32. -
               (3465 * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 128. + (63 * quick_pow(sintheta, 4)) / 64. +
               (3465 * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 128. - (231 * quick_pow(sintheta, 6)) / 128.;
    }
    else if (lp == 6 && l == 8 && s == 1 && j == 7)
    {
        return (-25 * sqrt(3.5)) / 128. - (105 * sqrt(3.5) * quick_pow(costheta, 2)) / 256. -
               (63 * sqrt(3.5) * quick_pow(costheta, 4)) / 128. - (231 * sqrt(3.5) * quick_pow(costheta, 6)) / 256. +
               (105 * sqrt(3.5) * quick_pow(sintheta, 2)) / 256. +
               (189 * sqrt(3.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 64. +
               (3465 * sqrt(3.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 256. -
               (63 * sqrt(3.5) * quick_pow(sintheta, 4)) / 128. -
               (3465 * sqrt(3.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 256. +
               (231 * sqrt(3.5) * quick_pow(sintheta, 6)) / 256.;
    }
    else if (lp == 8 && l == 6 && s == 1 && j == 7)
    {
        return (-1225 * sqrt(3.5)) / 8192. - (315 * sqrt(3.5) * quick_pow(costheta, 2)) / 1024. -
               (693 * sqrt(3.5) * quick_pow(costheta, 4)) / 2048. - (429 * sqrt(3.5) * quick_pow(costheta, 6)) / 1024. -
               (6435 * sqrt(3.5) * quick_pow(costheta, 8)) / 8192. +
               (315 * sqrt(3.5) * quick_pow(sintheta, 2)) / 1024. +
               (2079 * sqrt(3.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 1024. +
               (6435 * sqrt(3.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 1024. +
               (45045 * sqrt(3.5) * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 2048. -
               (693 * sqrt(3.5) * quick_pow(sintheta, 4)) / 2048. -
               (6435 * sqrt(3.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 1024. -
               (225225 * sqrt(3.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 4096. +
               (429 * sqrt(3.5) * quick_pow(sintheta, 6)) / 1024. +
               (45045 * sqrt(3.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 2048. -
               (6435 * sqrt(3.5) * quick_pow(sintheta, 8)) / 8192.;
    }
    else if (lp == 8 && l == 8 && s == 1 && j == 7)
    {
        return 0.261688232421875 + (2205 * quick_pow(costheta, 2)) / 4096. + (4851 * quick_pow(costheta, 4)) / 8192. +
               (3003 * quick_pow(costheta, 6)) / 4096. + (45045 * quick_pow(costheta, 8)) / 32768. -
               (2205 * quick_pow(sintheta, 2)) / 4096. -
               (14553 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 4096. -
               (45045 * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 4096. -
               (315315 * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 8192. +
               (4851 * quick_pow(sintheta, 4)) / 8192. +
               (45045 * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 4096. +
               (1576575 * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 16384. -
               (3003 * quick_pow(sintheta, 6)) / 4096. -
               (315315 * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 8192. +
               (45045 * quick_pow(sintheta, 8)) / 32768.;
    }
    else if (lp == 8 && l == 8 && s == 0 && j == 8)
    {
        return 0;
    }
    else if (lp == 8 && l == 8 && s == 1 && j == 8)
    {
        return 0.635528564453125 + (5355 * quick_pow(costheta, 2)) / 4096. + (11781 * quick_pow(costheta, 4)) / 8192. +
               (7293 * quick_pow(costheta, 6)) / 4096. + (109395 * quick_pow(costheta, 8)) / 32768. -
               (5355 * quick_pow(sintheta, 2)) / 4096. -
               (35343 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 4096. -
               (109395 * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 4096. -
               (765765 * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 8192. +
               (11781 * quick_pow(sintheta, 4)) / 8192. +
               (109395 * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 4096. +
               (3828825 * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 16384. -
               (7293 * quick_pow(sintheta, 6)) / 4096. -
               (765765 * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 8192. +
               (109395 * quick_pow(sintheta, 8)) / 32768.;
    }
    else if (lp == 7 && l == 7 && s == 1 && j == 8)
    {
        return (1575 * costheta) / 2048. + (1701 * quick_pow(costheta, 3)) / 2048. +
               (2079 * quick_pow(costheta, 5)) / 2048. + (3861 * quick_pow(costheta, 7)) / 2048. -
               (5103 * costheta * quick_pow(sintheta, 2)) / 2048. -
               (10395 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 1024. -
               (81081 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 2048. +
               (10395 * costheta * quick_pow(sintheta, 4)) / 2048. +
               (135135 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 2048. -
               (27027 * costheta * quick_pow(sintheta, 6)) / 2048.;
    }
    else if (lp == 7 && l == 9 && s == 1 && j == 8)
    {
        return (-525 * costheta) / (512. * sqrt2) - (567 * quick_pow(costheta, 3)) / (512. * sqrt2) -
               (693 * quick_pow(costheta, 5)) / (512. * sqrt2) - (1287 * quick_pow(costheta, 7)) / (512. * sqrt2) +
               (1701 * costheta * quick_pow(sintheta, 2)) / (512. * sqrt2) +
               (3465 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / (256. * sqrt2) +
               (27027 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / (512. * sqrt2) -
               (3465 * costheta * quick_pow(sintheta, 4)) / (512. * sqrt2) -
               (45045 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / (512. * sqrt2) +
               (9009 * costheta * quick_pow(sintheta, 6)) / (512. * sqrt2);
    }
    else if (lp == 9 && l == 7 && s == 1 && j == 8)
    {
        return (-6615 * costheta) / (8192. * sqrt2) - (3465 * quick_pow(costheta, 3)) / (4096. * sqrt2) -
               (3861 * quick_pow(costheta, 5)) / (4096. * sqrt2) - (19305 * quick_pow(costheta, 7)) / (16384. * sqrt2) -
               (36465 * quick_pow(costheta, 9)) / (16384. * sqrt2) +
               (10395 * costheta * quick_pow(sintheta, 2)) / (4096. * sqrt2) +
               (19305 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / (2048. * sqrt2) +
               (405405 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / (16384. * sqrt2) +
               (328185 * quick_pow(costheta, 7) * quick_pow(sintheta, 2)) / (4096. * sqrt2) -
               (19305 * costheta * quick_pow(sintheta, 4)) / (4096. * sqrt2) -
               (675675 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / (16384. * sqrt2) -
               (2297295 * quick_pow(costheta, 5) * quick_pow(sintheta, 4)) / (8192. * sqrt2) +
               (135135 * costheta * quick_pow(sintheta, 6)) / (16384. * sqrt2) +
               (765765 * quick_pow(costheta, 3) * quick_pow(sintheta, 6)) / (4096. * sqrt2) -
               (328185 * costheta * quick_pow(sintheta, 8)) / (16384. * sqrt2);
    }
    else if (lp == 9 && l == 9 && s == 1 && j == 8)
    {
        return (2205 * costheta) / 4096. + (1155 * quick_pow(costheta, 3)) / 2048. +
               (1287 * quick_pow(costheta, 5)) / 2048. + (6435 * quick_pow(costheta, 7)) / 8192. +
               (12155 * quick_pow(costheta, 9)) / 8192. - (3465 * costheta * quick_pow(sintheta, 2)) / 2048. -
               (6435 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 1024. -
               (135135 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 8192. -
               (109395 * quick_pow(costheta, 7) * quick_pow(sintheta, 2)) / 2048. +
               (6435 * costheta * quick_pow(sintheta, 4)) / 2048. +
               (225225 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 8192. +
               (765765 * quick_pow(costheta, 5) * quick_pow(sintheta, 4)) / 4096. -
               (45045 * costheta * quick_pow(sintheta, 6)) / 8192. -
               (255255 * quick_pow(costheta, 3) * quick_pow(sintheta, 6)) / 2048. +
               (109395 * costheta * quick_pow(sintheta, 8)) / 8192.;
    }
    else if (lp == 9 && l == 9 && s == 0 && j == 9)
    {
        return 0;
    }
    else if (lp == 9 && l == 9 && s == 1 && j == 9)
    {
        return (41895 * costheta) / 32768. + (21945 * quick_pow(costheta, 3)) / 16384. +
               (24453 * quick_pow(costheta, 5)) / 16384. + (122265 * quick_pow(costheta, 7)) / 65536. +
               (230945 * quick_pow(costheta, 9)) / 65536. - (65835 * costheta * quick_pow(sintheta, 2)) / 16384. -
               (122265 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 8192. -
               (2567565 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 65536. -
               (2078505 * quick_pow(costheta, 7) * quick_pow(sintheta, 2)) / 16384. +
               (122265 * costheta * quick_pow(sintheta, 4)) / 16384. +
               (4279275 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 65536. +
               (14549535 * quick_pow(costheta, 5) * quick_pow(sintheta, 4)) / 32768. -
               (855855 * costheta * quick_pow(sintheta, 6)) / 65536. -
               (4849845 * quick_pow(costheta, 3) * quick_pow(sintheta, 6)) / 16384. +
               (2078505 * costheta * quick_pow(sintheta, 8)) / 65536.;
    }
    else if (lp == 8 && l == 8 && s == 1 && j == 9)
    {
        return 0.37384033203125 + (1575 * quick_pow(costheta, 2)) / 2048. + (3465 * quick_pow(costheta, 4)) / 4096. +
               (2145 * quick_pow(costheta, 6)) / 2048. + (32175 * quick_pow(costheta, 8)) / 16384. -
               (1575 * quick_pow(sintheta, 2)) / 2048. -
               (10395 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 2048. -
               (32175 * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 2048. -
               (225225 * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 4096. +
               (3465 * quick_pow(sintheta, 4)) / 4096. +
               (32175 * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 2048. +
               (1126125 * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 8192. -
               (2145 * quick_pow(sintheta, 6)) / 2048. -
               (225225 * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 4096. +
               (32175 * quick_pow(sintheta, 8)) / 16384.;
    }
    else if (lp == 8 && l == 10 && s == 1 && j == 9)
    {
        return (-3675 * sqrt(2.5)) / 16384. - (945 * sqrt(2.5) * quick_pow(costheta, 2)) / 2048. -
               (2079 * sqrt(2.5) * quick_pow(costheta, 4)) / 4096. -
               (1287 * sqrt(2.5) * quick_pow(costheta, 6)) / 2048. -
               (19305 * sqrt(2.5) * quick_pow(costheta, 8)) / 16384. +
               (945 * sqrt(2.5) * quick_pow(sintheta, 2)) / 2048. +
               (6237 * sqrt(2.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 2048. +
               (19305 * sqrt(2.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 2048. +
               (135135 * sqrt(2.5) * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 4096. -
               (2079 * sqrt(2.5) * quick_pow(sintheta, 4)) / 4096. -
               (19305 * sqrt(2.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 2048. -
               (675675 * sqrt(2.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 8192. +
               (1287 * sqrt(2.5) * quick_pow(sintheta, 6)) / 2048. +
               (135135 * sqrt(2.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 4096. -
               (19305 * sqrt(2.5) * quick_pow(sintheta, 8)) / 16384.;
    }
    else if (lp == 10 && l == 8 && s == 1 && j == 9)
    {
        return (-11907 * sqrt(2.5)) / 65536. - (24255 * sqrt(2.5) * quick_pow(costheta, 2)) / 65536. -
               (6435 * sqrt(2.5) * quick_pow(costheta, 4)) / 16384. -
               (57915 * sqrt(2.5) * quick_pow(costheta, 6)) / 131072. -
               (36465 * sqrt(2.5) * quick_pow(costheta, 8)) / 65536. -
               (138567 * sqrt(2.5) * quick_pow(costheta, 10)) / 131072. +
               (24255 * sqrt(2.5) * quick_pow(sintheta, 2)) / 65536. +
               (19305 * sqrt(2.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 8192. +
               (868725 * sqrt(2.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 131072. +
               (255255 * sqrt(2.5) * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 16384. +
               (6235515 * sqrt(2.5) * quick_pow(costheta, 8) * quick_pow(sintheta, 2)) / 131072. -
               (6435 * sqrt(2.5) * quick_pow(sintheta, 4)) / 16384. -
               (868725 * sqrt(2.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 131072. -
               (1276275 * sqrt(2.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 32768. -
               (14549535 * sqrt(2.5) * quick_pow(costheta, 6) * quick_pow(sintheta, 4)) / 65536. +
               (57915 * sqrt(2.5) * quick_pow(sintheta, 6)) / 131072. +
               (255255 * sqrt(2.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 16384. +
               (14549535 * sqrt(2.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 6)) / 65536. -
               (36465 * sqrt(2.5) * quick_pow(sintheta, 8)) / 65536. -
               (6235515 * sqrt(2.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 8)) / 131072. +
               (138567 * sqrt(2.5) * quick_pow(sintheta, 10)) / 131072.;
    }
    else if (lp == 10 && l == 10 && s == 1 && j == 9)
    {
        return 0.27252960205078125 + (72765 * quick_pow(costheta, 2)) / 131072. +
               (19305 * quick_pow(costheta, 4)) / 32768. + (173745 * quick_pow(costheta, 6)) / 262144. +
               (109395 * quick_pow(costheta, 8)) / 131072. + (415701 * quick_pow(costheta, 10)) / 262144. -
               (72765 * quick_pow(sintheta, 2)) / 131072. -
               (57915 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 16384. -
               (2606175 * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 262144. -
               (765765 * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 32768. -
               (18706545 * quick_pow(costheta, 8) * quick_pow(sintheta, 2)) / 262144. +
               (19305 * quick_pow(sintheta, 4)) / 32768. +
               (2606175 * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 262144. +
               (3828825 * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 65536. +
               (43648605 * quick_pow(costheta, 6) * quick_pow(sintheta, 4)) / 131072. -
               (173745 * quick_pow(sintheta, 6)) / 262144. -
               (765765 * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 32768. -
               (43648605 * quick_pow(costheta, 4) * quick_pow(sintheta, 6)) / 131072. +
               (109395 * quick_pow(sintheta, 8)) / 131072. +
               (18706545 * quick_pow(costheta, 2) * quick_pow(sintheta, 8)) / 262144. -
               (415701 * quick_pow(sintheta, 10)) / 262144.;
    }
    else if (lp == 10 && l == 10 && s == 0 && j == 10)
    {
        return 0;
    }
    else if (lp == 10 && l == 10 && s == 1 && j == 10)
    {
        return 0.6359024047851562 + (169785 * quick_pow(costheta, 2)) / 131072. +
               (45045 * quick_pow(costheta, 4)) / 32768. + (405405 * quick_pow(costheta, 6)) / 262144. +
               (255255 * quick_pow(costheta, 8)) / 131072. + (969969 * quick_pow(costheta, 10)) / 262144. -
               (169785 * quick_pow(sintheta, 2)) / 131072. -
               (135135 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 16384. -
               (6081075 * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 262144. -
               (1786785 * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 32768. -
               (43648605 * quick_pow(costheta, 8) * quick_pow(sintheta, 2)) / 262144. +
               (45045 * quick_pow(sintheta, 4)) / 32768. +
               (6081075 * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 262144. +
               (8933925 * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 65536. +
               (101846745 * quick_pow(costheta, 6) * quick_pow(sintheta, 4)) / 131072. -
               (405405 * quick_pow(sintheta, 6)) / 262144. -
               (1786785 * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 32768. -
               (101846745 * quick_pow(costheta, 4) * quick_pow(sintheta, 6)) / 131072. +
               (255255 * quick_pow(sintheta, 8)) / 131072. +
               (43648605 * quick_pow(costheta, 2) * quick_pow(sintheta, 8)) / 262144. -
               (969969 * quick_pow(sintheta, 10)) / 262144.;
    }
    else if (lp == 9 && l == 9 && s == 1 && j == 10)
    {
        return (24255 * costheta) / 32768. + (12705 * quick_pow(costheta, 3)) / 16384. +
               (14157 * quick_pow(costheta, 5)) / 16384. + (70785 * quick_pow(costheta, 7)) / 65536. +
               (133705 * quick_pow(costheta, 9)) / 65536. - (38115 * costheta * quick_pow(sintheta, 2)) / 16384. -
               (70785 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 8192. -
               (1486485 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 65536. -
               (1203345 * quick_pow(costheta, 7) * quick_pow(sintheta, 2)) / 16384. +
               (70785 * costheta * quick_pow(sintheta, 4)) / 16384. +
               (2477475 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 65536. +
               (8423415 * quick_pow(costheta, 5) * quick_pow(sintheta, 4)) / 32768. -
               (495495 * costheta * quick_pow(sintheta, 6)) / 65536. -
               (2807805 * quick_pow(costheta, 3) * quick_pow(sintheta, 6)) / 16384. +
               (1203345 * costheta * quick_pow(sintheta, 8)) / 65536.;
    }
    else if (lp == 9 && l == 11 && s == 1 && j == 10)
    {
        return (-2205 * sqrt(27.5) * costheta) / 16384. - (1155 * sqrt(27.5) * quick_pow(costheta, 3)) / 8192. -
               (1287 * sqrt(27.5) * quick_pow(costheta, 5)) / 8192. -
               (6435 * sqrt(27.5) * quick_pow(costheta, 7)) / 32768. -
               (12155 * sqrt(27.5) * quick_pow(costheta, 9)) / 32768. +
               (3465 * sqrt(27.5) * costheta * quick_pow(sintheta, 2)) / 8192. +
               (6435 * sqrt(27.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 4096. +
               (135135 * sqrt(27.5) * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 32768. +
               (109395 * sqrt(27.5) * quick_pow(costheta, 7) * quick_pow(sintheta, 2)) / 8192. -
               (6435 * sqrt(27.5) * costheta * quick_pow(sintheta, 4)) / 8192. -
               (225225 * sqrt(27.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 32768. -
               (765765 * sqrt(27.5) * quick_pow(costheta, 5) * quick_pow(sintheta, 4)) / 16384. +
               (45045 * sqrt(27.5) * costheta * quick_pow(sintheta, 6)) / 32768. +
               (255255 * sqrt(27.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 6)) / 8192. -
               (109395 * sqrt(27.5) * costheta * quick_pow(sintheta, 8)) / 32768.;
    }
    else if (lp == 11 && l == 9 && s == 1 && j == 10)
    {
        return (-14553 * sqrt(27.5) * costheta) / 131072. - (15015 * sqrt(27.5) * quick_pow(costheta, 3)) / 131072. -
               (32175 * sqrt(27.5) * quick_pow(costheta, 5)) / 262144. -
               (36465 * sqrt(27.5) * quick_pow(costheta, 7)) / 262144. -
               (46189 * sqrt(27.5) * quick_pow(costheta, 9)) / 262144. -
               (88179 * sqrt(27.5) * quick_pow(costheta, 11)) / 262144. +
               (45045 * sqrt(27.5) * costheta * quick_pow(sintheta, 2)) / 131072. +
               (160875 * sqrt(27.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 131072. +
               (765765 * sqrt(27.5) * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 262144. +
               (415701 * sqrt(27.5) * quick_pow(costheta, 7) * quick_pow(sintheta, 2)) / 65536. +
               (4849845 * sqrt(27.5) * quick_pow(costheta, 9) * quick_pow(sintheta, 2)) / 262144. -
               (160875 * sqrt(27.5) * costheta * quick_pow(sintheta, 4)) / 262144. -
               (1276275 * sqrt(27.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 262144. -
               (2909907 * sqrt(27.5) * quick_pow(costheta, 5) * quick_pow(sintheta, 4)) / 131072. -
               (14549535 * sqrt(27.5) * quick_pow(costheta, 7) * quick_pow(sintheta, 4)) / 131072. +
               (255255 * sqrt(27.5) * costheta * quick_pow(sintheta, 6)) / 262144. +
               (969969 * sqrt(27.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 6)) / 65536. +
               (20369349 * sqrt(27.5) * quick_pow(costheta, 5) * quick_pow(sintheta, 6)) / 131072. -
               (415701 * sqrt(27.5) * costheta * quick_pow(sintheta, 8)) / 262144. -
               (14549535 * sqrt(27.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 8)) / 262144. +
               (969969 * sqrt(27.5) * costheta * quick_pow(sintheta, 10)) / 262144.;
    }
    else if (lp == 11 && l == 11 && s == 1 && j == 10)
    {
        return (72765 * costheta) / 131072. + (75075 * quick_pow(costheta, 3)) / 131072. +
               (160875 * quick_pow(costheta, 5)) / 262144. + (182325 * quick_pow(costheta, 7)) / 262144. +
               (230945 * quick_pow(costheta, 9)) / 262144. + (440895 * quick_pow(costheta, 11)) / 262144. -
               (225225 * costheta * quick_pow(sintheta, 2)) / 131072. -
               (804375 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 131072. -
               (3828825 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 262144. -
               (2078505 * quick_pow(costheta, 7) * quick_pow(sintheta, 2)) / 65536. -
               (24249225 * quick_pow(costheta, 9) * quick_pow(sintheta, 2)) / 262144. +
               (804375 * costheta * quick_pow(sintheta, 4)) / 262144. +
               (6381375 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 262144. +
               (14549535 * quick_pow(costheta, 5) * quick_pow(sintheta, 4)) / 131072. +
               (72747675 * quick_pow(costheta, 7) * quick_pow(sintheta, 4)) / 131072. -
               (1276275 * costheta * quick_pow(sintheta, 6)) / 262144. -
               (4849845 * quick_pow(costheta, 3) * quick_pow(sintheta, 6)) / 65536. -
               (101846745 * quick_pow(costheta, 5) * quick_pow(sintheta, 6)) / 131072. +
               (2078505 * costheta * quick_pow(sintheta, 8)) / 262144. +
               (72747675 * quick_pow(costheta, 3) * quick_pow(sintheta, 8)) / 262144. -
               (4849845 * costheta * quick_pow(sintheta, 10)) / 262144.;
    }
    return 0;
}

double coeff_m10(double t, int lp, int l, int s, int j)
{
    double costheta = t;
    double sintheta = sqrt(1.0 - t * t);
    if (lp == 0 && l == 0 && s == 0 && j == 0)
    {
        return 0;
    }
    else if (lp == 1 && l == 1 && s == 1 && j == 0)
    {
        return -(sintheta / sqrt2);
    }
    else if (lp == 1 && l == 1 && s == 0 && j == 1)
    {
        return 0;
    }
    else if (lp == 1 && l == 1 && s == 1 && j == 1)
    {
        return 0;
    }
    else if (lp == 0 && l == 0 && s == 1 && j == 1)
    {
        return 0;
    }
    else if (lp == 0 && l == 2 && s == 1 && j == 1)
    {
        return 0;
    }
    else if (lp == 2 && l == 0 && s == 1 && j == 1)
    {
        return (-3 * costheta * sintheta) / 2.;
    }
    else if (lp == 2 && l == 2 && s == 1 && j == 1)
    {
        return (-3 * costheta * sintheta) / sqrt2;
    }
    else if (lp == 2 && l == 2 && s == 0 && j == 2)
    {
        return 0;
    }
    else if (lp == 2 && l == 2 && s == 1 && j == 2)
    {
        return 0;
    }
    else if (lp == 1 && l == 1 && s == 1 && j == 2)
    {
        return sintheta / sqrt2;
    }
    else if (lp == 1 && l == 3 && s == 1 && j == 2)
    {
        return (sqrt3 * sintheta) / 2.;
    }
    else if (lp == 3 && l == 1 && s == 1 && j == 2)
    {
        return -0.125 * (sqrt3 * sintheta) - (15 * sqrt3 * quick_pow(costheta, 2) * sintheta) / 8. +
               (5 * sqrt3 * quick_pow(sintheta, 3)) / 8.;
    }
    else if (lp == 3 && l == 3 && s == 1 && j == 2)
    {
        return (-3 * sintheta) / (8. * sqrt2) - (45 * quick_pow(costheta, 2) * sintheta) / (8. * sqrt2) +
               (15 * quick_pow(sintheta, 3)) / (8. * sqrt2);
    }
    else if (lp == 3 && l == 3 && s == 0 && j == 3)
    {
        return 0;
    }
    else if (lp == 3 && l == 3 && s == 1 && j == 3)
    {
        return 0;
    }
    else if (lp == 2 && l == 2 && s == 1 && j == 3)
    {
        return (3 * costheta * sintheta) / sqrt2;
    }
    else if (lp == 2 && l == 4 && s == 1 && j == 3)
    {
        return sqrt(6) * costheta * sintheta;
    }
    else if (lp == 4 && l == 2 && s == 1 && j == 3)
    {
        return (-5 * sqrt(1.5) * costheta * sintheta) / 8. - (35 * sqrt(1.5) * quick_pow(costheta, 3) * sintheta) / 8. +
               (35 * sqrt(1.5) * costheta * quick_pow(sintheta, 3)) / 8.;
    }
    else if (lp == 4 && l == 4 && s == 1 && j == 3)
    {
        return (-5 * costheta * sintheta) / (4. * sqrt2) - (35 * quick_pow(costheta, 3) * sintheta) / (4. * sqrt2) +
               (35 * costheta * quick_pow(sintheta, 3)) / (4. * sqrt2);
    }
    else if (lp == 4 && l == 4 && s == 0 && j == 4)
    {
        return 0;
    }
    else if (lp == 4 && l == 4 && s == 1 && j == 4)
    {
        return 0;
    }
    else if (lp == 3 && l == 3 && s == 1 && j == 4)
    {
        return (3 * sintheta) / (8. * sqrt2) + (45 * quick_pow(costheta, 2) * sintheta) / (8. * sqrt2) -
               (15 * quick_pow(sintheta, 3)) / (8. * sqrt2);
    }
    else if (lp == 3 && l == 5 && s == 1 && j == 4)
    {
        return (3 * sqrt(2.5) * sintheta) / 16. + (45 * sqrt(2.5) * quick_pow(costheta, 2) * sintheta) / 16. -
               (15 * sqrt(2.5) * quick_pow(sintheta, 3)) / 16.;
    }
    else if (lp == 5 && l == 3 && s == 1 && j == 4)
    {
        return (-3 * sqrt(2.5) * sintheta) / 32. - (63 * sqrt(2.5) * quick_pow(costheta, 2) * sintheta) / 64. -
               (315 * sqrt(2.5) * quick_pow(costheta, 4) * sintheta) / 64. +
               (21 * sqrt(2.5) * quick_pow(sintheta, 3)) / 64. +
               (315 * sqrt(2.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / 32. -
               (63 * sqrt(2.5) * quick_pow(sintheta, 5)) / 64.;
    }
    else if (lp == 5 && l == 5 && s == 1 && j == 4)
    {
        return (-15 * sintheta) / (64. * sqrt2) - (315 * quick_pow(costheta, 2) * sintheta) / (128. * sqrt2) -
               (1575 * quick_pow(costheta, 4) * sintheta) / (128. * sqrt2) +
               (105 * quick_pow(sintheta, 3)) / (128. * sqrt2) +
               (1575 * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / (64. * sqrt2) -
               (315 * quick_pow(sintheta, 5)) / (128. * sqrt2);
    }
    else if (lp == 5 && l == 5 && s == 0 && j == 5)
    {
        return 0;
    }
    else if (lp == 5 && l == 5 && s == 1 && j == 5)
    {
        return 0;
    }
    else if (lp == 4 && l == 4 && s == 1 && j == 5)
    {
        return (5 * costheta * sintheta) / (4. * sqrt2) + (35 * quick_pow(costheta, 3) * sintheta) / (4. * sqrt2) -
               (35 * costheta * quick_pow(sintheta, 3)) / (4. * sqrt2);
    }
    else if (lp == 4 && l == 6 && s == 1 && j == 5)
    {
        return (sqrt15 * costheta * sintheta) / 4. + (7 * sqrt15 * quick_pow(costheta, 3) * sintheta) / 4. -
               (7 * sqrt15 * costheta * quick_pow(sintheta, 3)) / 4.;
    }
    else if (lp == 6 && l == 4 && s == 1 && j == 5)
    {
        return (-35 * sqrt15 * costheta * sintheta) / 256. - (21 * sqrt15 * quick_pow(costheta, 3) * sintheta) / 32. -
               (693 * sqrt15 * quick_pow(costheta, 5) * sintheta) / 256. +
               (21 * sqrt15 * costheta * quick_pow(sintheta, 3)) / 32. +
               (1155 * sqrt15 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / 128. -
               (693 * sqrt15 * costheta * quick_pow(sintheta, 5)) / 256.;
    }
    else if (lp == 6 && l == 6 && s == 1 && j == 5)
    {
        return (-105 * costheta * sintheta) / (128. * sqrt2) -
               (63 * quick_pow(costheta, 3) * sintheta) / (16. * sqrt2) -
               (2079 * quick_pow(costheta, 5) * sintheta) / (128. * sqrt2) +
               (63 * costheta * quick_pow(sintheta, 3)) / (16. * sqrt2) +
               (3465 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / (64. * sqrt2) -
               (2079 * costheta * quick_pow(sintheta, 5)) / (128. * sqrt2);
    }
    else if (lp == 6 && l == 6 && s == 0 && j == 6)
    {
        return 0;
    }
    else if (lp == 6 && l == 6 && s == 1 && j == 6)
    {
        return 0;
    }
    else if (lp == 5 && l == 5 && s == 1 && j == 6)
    {
        return (15 * sintheta) / (64. * sqrt2) + (315 * quick_pow(costheta, 2) * sintheta) / (128. * sqrt2) +
               (1575 * quick_pow(costheta, 4) * sintheta) / (128. * sqrt2) -
               (105 * quick_pow(sintheta, 3)) / (128. * sqrt2) -
               (1575 * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / (64. * sqrt2) +
               (315 * quick_pow(sintheta, 5)) / (128. * sqrt2);
    }
    else if (lp == 5 && l == 7 && s == 1 && j == 6)
    {
        return (5 * sqrt(21) * sintheta) / 128. + (105 * sqrt(21) * quick_pow(costheta, 2) * sintheta) / 256. +
               (525 * sqrt(21) * quick_pow(costheta, 4) * sintheta) / 256. -
               (35 * sqrt(21) * quick_pow(sintheta, 3)) / 256. -
               (525 * sqrt(21) * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / 128. +
               (105 * sqrt(21) * quick_pow(sintheta, 5)) / 256.;
    }
    else if (lp == 7 && l == 5 && s == 1 && j == 6)
    {
        return (-25 * sqrt(21) * sintheta) / 1024. - (243 * sqrt(21) * quick_pow(costheta, 2) * sintheta) / 1024. -
               (825 * sqrt(21) * quick_pow(costheta, 4) * sintheta) / 1024. -
               (3003 * sqrt(21) * quick_pow(costheta, 6) * sintheta) / 1024. +
               (81 * sqrt(21) * quick_pow(sintheta, 3)) / 1024. +
               (825 * sqrt(21) * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / 512. +
               (15015 * sqrt(21) * quick_pow(costheta, 4) * quick_pow(sintheta, 3)) / 1024. -
               (165 * sqrt(21) * quick_pow(sintheta, 5)) / 1024. -
               (9009 * sqrt(21) * quick_pow(costheta, 2) * quick_pow(sintheta, 5)) / 1024. +
               (429 * sqrt(21) * quick_pow(sintheta, 7)) / 1024.;
    }
    else if (lp == 7 && l == 7 && s == 1 && j == 6)
    {
        return (-175 * sintheta) / (1024. * sqrt2) - (1701 * quick_pow(costheta, 2) * sintheta) / (1024. * sqrt2) -
               (5775 * quick_pow(costheta, 4) * sintheta) / (1024. * sqrt2) -
               (21021 * quick_pow(costheta, 6) * sintheta) / (1024. * sqrt2) +
               (567 * quick_pow(sintheta, 3)) / (1024. * sqrt2) +
               (5775 * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / (512. * sqrt2) +
               (105105 * quick_pow(costheta, 4) * quick_pow(sintheta, 3)) / (1024. * sqrt2) -
               (1155 * quick_pow(sintheta, 5)) / (1024. * sqrt2) -
               (63063 * quick_pow(costheta, 2) * quick_pow(sintheta, 5)) / (1024. * sqrt2) +
               (3003 * quick_pow(sintheta, 7)) / (1024. * sqrt2);
    }
    else if (lp == 7 && l == 7 && s == 0 && j == 7)
    {
        return 0;
    }
    else if (lp == 7 && l == 7 && s == 1 && j == 7)
    {
        return 0;
    }
    else if (lp == 6 && l == 6 && s == 1 && j == 7)
    {
        return (105 * costheta * sintheta) / (128. * sqrt2) + (63 * quick_pow(costheta, 3) * sintheta) / (16. * sqrt2) +
               (2079 * quick_pow(costheta, 5) * sintheta) / (128. * sqrt2) -
               (63 * costheta * quick_pow(sintheta, 3)) / (16. * sqrt2) -
               (3465 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / (64. * sqrt2) +
               (2079 * costheta * quick_pow(sintheta, 5)) / (128. * sqrt2);
    }
    else if (lp == 6 && l == 8 && s == 1 && j == 7)
    {
        return (15 * sqrt7 * costheta * sintheta) / 64. + (9 * sqrt7 * quick_pow(costheta, 3) * sintheta) / 8. +
               (297 * sqrt7 * quick_pow(costheta, 5) * sintheta) / 64. -
               (9 * sqrt7 * costheta * quick_pow(sintheta, 3)) / 8. -
               (495 * sqrt7 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / 32. +
               (297 * sqrt7 * costheta * quick_pow(sintheta, 5)) / 64.;
    }
    else if (lp == 8 && l == 6 && s == 1 && j == 7)
    {
        return (-315 * sqrt7 * costheta * sintheta) / 2048. -
               (693 * sqrt7 * quick_pow(costheta, 3) * sintheta) / 1024. -
               (3861 * sqrt7 * quick_pow(costheta, 5) * sintheta) / 2048. -
               (6435 * sqrt7 * quick_pow(costheta, 7) * sintheta) / 1024. +
               (693 * sqrt7 * costheta * quick_pow(sintheta, 3)) / 1024. +
               (6435 * sqrt7 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / 1024. +
               (45045 * sqrt7 * quick_pow(costheta, 5) * quick_pow(sintheta, 3)) / 1024. -
               (3861 * sqrt7 * costheta * quick_pow(sintheta, 5)) / 2048. -
               (45045 * sqrt7 * quick_pow(costheta, 3) * quick_pow(sintheta, 5)) / 1024. +
               (6435 * sqrt7 * costheta * quick_pow(sintheta, 7)) / 1024.;
    }
    else if (lp == 8 && l == 8 && s == 1 && j == 7)
    {
        return (-315 * costheta * sintheta) / (512. * sqrt2) -
               (693 * quick_pow(costheta, 3) * sintheta) / (256. * sqrt2) -
               (3861 * quick_pow(costheta, 5) * sintheta) / (512. * sqrt2) -
               (6435 * quick_pow(costheta, 7) * sintheta) / (256. * sqrt2) +
               (693 * costheta * quick_pow(sintheta, 3)) / (256. * sqrt2) +
               (6435 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / (256. * sqrt2) +
               (45045 * quick_pow(costheta, 5) * quick_pow(sintheta, 3)) / (256. * sqrt2) -
               (3861 * costheta * quick_pow(sintheta, 5)) / (512. * sqrt2) -
               (45045 * quick_pow(costheta, 3) * quick_pow(sintheta, 5)) / (256. * sqrt2) +
               (6435 * costheta * quick_pow(sintheta, 7)) / (256. * sqrt2);
    }
    else if (lp == 8 && l == 8 && s == 0 && j == 8)
    {
        return 0;
    }
    else if (lp == 8 && l == 8 && s == 1 && j == 8)
    {
        return 0;
    }
    else if (lp == 7 && l == 7 && s == 1 && j == 8)
    {
        return (175 * sintheta) / (1024. * sqrt2) + (1701 * quick_pow(costheta, 2) * sintheta) / (1024. * sqrt2) +
               (5775 * quick_pow(costheta, 4) * sintheta) / (1024. * sqrt2) +
               (21021 * quick_pow(costheta, 6) * sintheta) / (1024. * sqrt2) -
               (567 * quick_pow(sintheta, 3)) / (1024. * sqrt2) -
               (5775 * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / (512. * sqrt2) -
               (105105 * quick_pow(costheta, 4) * quick_pow(sintheta, 3)) / (1024. * sqrt2) +
               (1155 * quick_pow(sintheta, 5)) / (1024. * sqrt2) +
               (63063 * quick_pow(costheta, 2) * quick_pow(sintheta, 5)) / (1024. * sqrt2) -
               (3003 * quick_pow(sintheta, 7)) / (1024. * sqrt2);
    }
    else if (lp == 7 && l == 9 && s == 1 && j == 8)
    {
        return (525 * sintheta) / 4096. + (5103 * quick_pow(costheta, 2) * sintheta) / 4096. +
               (17325 * quick_pow(costheta, 4) * sintheta) / 4096. +
               (63063 * quick_pow(costheta, 6) * sintheta) / 4096. - (1701 * quick_pow(sintheta, 3)) / 4096. -
               (17325 * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / 2048. -
               (315315 * quick_pow(costheta, 4) * quick_pow(sintheta, 3)) / 4096. +
               (3465 * quick_pow(sintheta, 5)) / 4096. +
               (189189 * quick_pow(costheta, 2) * quick_pow(sintheta, 5)) / 4096. -
               (9009 * quick_pow(sintheta, 7)) / 4096.;
    }
    else if (lp == 9 && l == 7 && s == 1 && j == 8)
    {
        return (-735 * sintheta) / 8192. - (3465 * quick_pow(costheta, 2) * sintheta) / 4096. -
               (10725 * quick_pow(costheta, 4) * sintheta) / 4096. -
               (105105 * quick_pow(costheta, 6) * sintheta) / 16384. -
               (328185 * quick_pow(costheta, 8) * sintheta) / 16384. + (1155 * quick_pow(sintheta, 3)) / 4096. +
               (10725 * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / 2048. +
               (525525 * quick_pow(costheta, 4) * quick_pow(sintheta, 3)) / 16384. +
               (765765 * quick_pow(costheta, 6) * quick_pow(sintheta, 3)) / 4096. -
               (2145 * quick_pow(sintheta, 5)) / 4096. -
               (315315 * quick_pow(costheta, 2) * quick_pow(sintheta, 5)) / 16384. -
               (2297295 * quick_pow(costheta, 4) * quick_pow(sintheta, 5)) / 8192. +
               (15015 * quick_pow(sintheta, 7)) / 16384. +
               (328185 * quick_pow(costheta, 2) * quick_pow(sintheta, 7)) / 4096. -
               (36465 * quick_pow(sintheta, 9)) / 16384.;
    }
    else if (lp == 9 && l == 9 && s == 1 && j == 8)
    {
        return (-2205 * sintheta) / (16384. * sqrt2) - (10395 * quick_pow(costheta, 2) * sintheta) / (8192. * sqrt2) -
               (32175 * quick_pow(costheta, 4) * sintheta) / (8192. * sqrt2) -
               (315315 * quick_pow(costheta, 6) * sintheta) / (32768. * sqrt2) -
               (984555 * quick_pow(costheta, 8) * sintheta) / (32768. * sqrt2) +
               (3465 * quick_pow(sintheta, 3)) / (8192. * sqrt2) +
               (32175 * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / (4096. * sqrt2) +
               (1576575 * quick_pow(costheta, 4) * quick_pow(sintheta, 3)) / (32768. * sqrt2) +
               (2297295 * quick_pow(costheta, 6) * quick_pow(sintheta, 3)) / (8192. * sqrt2) -
               (6435 * quick_pow(sintheta, 5)) / (8192. * sqrt2) -
               (945945 * quick_pow(costheta, 2) * quick_pow(sintheta, 5)) / (32768. * sqrt2) -
               (6891885 * quick_pow(costheta, 4) * quick_pow(sintheta, 5)) / (16384. * sqrt2) +
               (45045 * quick_pow(sintheta, 7)) / (32768. * sqrt2) +
               (984555 * quick_pow(costheta, 2) * quick_pow(sintheta, 7)) / (8192. * sqrt2) -
               (109395 * quick_pow(sintheta, 9)) / (32768. * sqrt2);
    }
    else if (lp == 9 && l == 9 && s == 0 && j == 9)
    {
        return 0;
    }
    else if (lp == 9 && l == 9 && s == 1 && j == 9)
    {
        return 0;
    }
    else if (lp == 8 && l == 8 && s == 1 && j == 9)
    {
        return (315 * costheta * sintheta) / (512. * sqrt2) +
               (693 * quick_pow(costheta, 3) * sintheta) / (256. * sqrt2) +
               (3861 * quick_pow(costheta, 5) * sintheta) / (512. * sqrt2) +
               (6435 * quick_pow(costheta, 7) * sintheta) / (256. * sqrt2) -
               (693 * costheta * quick_pow(sintheta, 3)) / (256. * sqrt2) -
               (6435 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / (256. * sqrt2) -
               (45045 * quick_pow(costheta, 5) * quick_pow(sintheta, 3)) / (256. * sqrt2) +
               (3861 * costheta * quick_pow(sintheta, 5)) / (512. * sqrt2) +
               (45045 * quick_pow(costheta, 3) * quick_pow(sintheta, 5)) / (256. * sqrt2) -
               (6435 * costheta * quick_pow(sintheta, 7)) / (256. * sqrt2);
    }
    else if (lp == 8 && l == 10 && s == 1 && j == 9)
    {
        return (105 * sqrt5 * costheta * sintheta) / 512. + (231 * sqrt5 * quick_pow(costheta, 3) * sintheta) / 256. +
               (1287 * sqrt5 * quick_pow(costheta, 5) * sintheta) / 512. +
               (2145 * sqrt5 * quick_pow(costheta, 7) * sintheta) / 256. -
               (231 * sqrt5 * costheta * quick_pow(sintheta, 3)) / 256. -
               (2145 * sqrt5 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / 256. -
               (15015 * sqrt5 * quick_pow(costheta, 5) * quick_pow(sintheta, 3)) / 256. +
               (1287 * sqrt5 * costheta * quick_pow(sintheta, 5)) / 512. +
               (15015 * sqrt5 * quick_pow(costheta, 3) * quick_pow(sintheta, 5)) / 256. -
               (2145 * sqrt5 * costheta * quick_pow(sintheta, 7)) / 256.;
    }
    else if (lp == 10 && l == 8 && s == 1 && j == 9)
    {
        return (-4851 * sqrt5 * costheta * sintheta) / 32768. -
               (1287 * sqrt5 * quick_pow(costheta, 3) * sintheta) / 2048. -
               (104247 * sqrt5 * quick_pow(costheta, 5) * sintheta) / 65536. -
               (7293 * sqrt5 * quick_pow(costheta, 7) * sintheta) / 2048. -
               (692835 * sqrt5 * quick_pow(costheta, 9) * sintheta) / 65536. +
               (1287 * sqrt5 * costheta * quick_pow(sintheta, 3)) / 2048. +
               (173745 * sqrt5 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / 32768. +
               (51051 * sqrt5 * quick_pow(costheta, 5) * quick_pow(sintheta, 3)) / 2048. +
               (2078505 * sqrt5 * quick_pow(costheta, 7) * quick_pow(sintheta, 3)) / 16384. -
               (104247 * sqrt5 * costheta * quick_pow(sintheta, 5)) / 65536. -
               (51051 * sqrt5 * quick_pow(costheta, 3) * quick_pow(sintheta, 5)) / 2048. -
               (8729721 * sqrt5 * quick_pow(costheta, 5) * quick_pow(sintheta, 5)) / 32768. +
               (7293 * sqrt5 * costheta * quick_pow(sintheta, 7)) / 2048. +
               (2078505 * sqrt5 * quick_pow(costheta, 3) * quick_pow(sintheta, 7)) / 16384. -
               (692835 * sqrt5 * costheta * quick_pow(sintheta, 9)) / 65536.;
    }
    else if (lp == 10 && l == 10 && s == 1 && j == 9)
    {
        return (-8085 * costheta * sintheta) / (16384. * sqrt2) -
               (2145 * quick_pow(costheta, 3) * sintheta) / (1024. * sqrt2) -
               (173745 * quick_pow(costheta, 5) * sintheta) / (32768. * sqrt2) -
               (12155 * quick_pow(costheta, 7) * sintheta) / (1024. * sqrt2) -
               (1154725 * quick_pow(costheta, 9) * sintheta) / (32768. * sqrt2) +
               (2145 * costheta * quick_pow(sintheta, 3)) / (1024. * sqrt2) +
               (289575 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / (16384. * sqrt2) +
               (85085 * quick_pow(costheta, 5) * quick_pow(sintheta, 3)) / (1024. * sqrt2) +
               (3464175 * quick_pow(costheta, 7) * quick_pow(sintheta, 3)) / (8192. * sqrt2) -
               (173745 * costheta * quick_pow(sintheta, 5)) / (32768. * sqrt2) -
               (85085 * quick_pow(costheta, 3) * quick_pow(sintheta, 5)) / (1024. * sqrt2) -
               (14549535 * quick_pow(costheta, 5) * quick_pow(sintheta, 5)) / (16384. * sqrt2) +
               (12155 * costheta * quick_pow(sintheta, 7)) / (1024. * sqrt2) +
               (3464175 * quick_pow(costheta, 3) * quick_pow(sintheta, 7)) / (8192. * sqrt2) -
               (1154725 * costheta * quick_pow(sintheta, 9)) / (32768. * sqrt2);
    }
    else if (lp == 10 && l == 10 && s == 0 && j == 10)
    {
        return 0;
    }
    else if (lp == 10 && l == 10 && s == 1 && j == 10)
    {
        return 0;
    }
    else if (lp == 9 && l == 9 && s == 1 && j == 10)
    {
        return (2205 * sintheta) / (16384. * sqrt2) + (10395 * quick_pow(costheta, 2) * sintheta) / (8192. * sqrt2) +
               (32175 * quick_pow(costheta, 4) * sintheta) / (8192. * sqrt2) +
               (315315 * quick_pow(costheta, 6) * sintheta) / (32768. * sqrt2) +
               (984555 * quick_pow(costheta, 8) * sintheta) / (32768. * sqrt2) -
               (3465 * quick_pow(sintheta, 3)) / (8192. * sqrt2) -
               (32175 * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / (4096. * sqrt2) -
               (1576575 * quick_pow(costheta, 4) * quick_pow(sintheta, 3)) / (32768. * sqrt2) -
               (2297295 * quick_pow(costheta, 6) * quick_pow(sintheta, 3)) / (8192. * sqrt2) +
               (6435 * quick_pow(sintheta, 5)) / (8192. * sqrt2) +
               (945945 * quick_pow(costheta, 2) * quick_pow(sintheta, 5)) / (32768. * sqrt2) +
               (6891885 * quick_pow(costheta, 4) * quick_pow(sintheta, 5)) / (16384. * sqrt2) -
               (45045 * quick_pow(sintheta, 7)) / (32768. * sqrt2) -
               (984555 * quick_pow(costheta, 2) * quick_pow(sintheta, 7)) / (8192. * sqrt2) +
               (109395 * quick_pow(sintheta, 9)) / (32768. * sqrt2);
    }
    else if (lp == 9 && l == 11 && s == 1 && j == 10)
    {
        return (441 * sqrt(55) * sintheta) / 32768. + (2079 * sqrt(55) * quick_pow(costheta, 2) * sintheta) / 16384. +
               (6435 * sqrt(55) * quick_pow(costheta, 4) * sintheta) / 16384. +
               (63063 * sqrt(55) * quick_pow(costheta, 6) * sintheta) / 65536. +
               (196911 * sqrt(55) * quick_pow(costheta, 8) * sintheta) / 65536. -
               (693 * sqrt(55) * quick_pow(sintheta, 3)) / 16384. -
               (6435 * sqrt(55) * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / 8192. -
               (315315 * sqrt(55) * quick_pow(costheta, 4) * quick_pow(sintheta, 3)) / 65536. -
               (459459 * sqrt(55) * quick_pow(costheta, 6) * quick_pow(sintheta, 3)) / 16384. +
               (1287 * sqrt(55) * quick_pow(sintheta, 5)) / 16384. +
               (189189 * sqrt(55) * quick_pow(costheta, 2) * quick_pow(sintheta, 5)) / 65536. +
               (1378377 * sqrt(55) * quick_pow(costheta, 4) * quick_pow(sintheta, 5)) / 32768. -
               (9009 * sqrt(55) * quick_pow(sintheta, 7)) / 65536. -
               (196911 * sqrt(55) * quick_pow(costheta, 2) * quick_pow(sintheta, 7)) / 16384. +
               (21879 * sqrt(55) * quick_pow(sintheta, 9)) / 65536.;
    }
    else if (lp == 11 && l == 9 && s == 1 && j == 10)
    {
        return (-1323 * sqrt(55) * sintheta) / 131072. -
               (12285 * sqrt(55) * quick_pow(costheta, 2) * sintheta) / 131072. -
               (73125 * sqrt(55) * quick_pow(costheta, 4) * sintheta) / 262144. -
               (162435 * sqrt(55) * quick_pow(costheta, 6) * sintheta) / 262144. -
               (340119 * sqrt(55) * quick_pow(costheta, 8) * sintheta) / 262144. -
               (969969 * sqrt(55) * quick_pow(costheta, 10) * sintheta) / 262144. +
               (4095 * sqrt(55) * quick_pow(sintheta, 3)) / 131072. +
               (73125 * sqrt(55) * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / 131072. +
               (812175 * sqrt(55) * quick_pow(costheta, 4) * quick_pow(sintheta, 3)) / 262144. +
               (793611 * sqrt(55) * quick_pow(costheta, 6) * quick_pow(sintheta, 3)) / 65536. +
               (14549535 * sqrt(55) * quick_pow(costheta, 8) * quick_pow(sintheta, 3)) / 262144. -
               (14625 * sqrt(55) * quick_pow(sintheta, 5)) / 262144. -
               (487305 * sqrt(55) * quick_pow(costheta, 2) * quick_pow(sintheta, 5)) / 262144. -
               (2380833 * sqrt(55) * quick_pow(costheta, 4) * quick_pow(sintheta, 5)) / 131072. -
               (20369349 * sqrt(55) * quick_pow(costheta, 6) * quick_pow(sintheta, 5)) / 131072. +
               (23205 * sqrt(55) * quick_pow(sintheta, 7)) / 262144. +
               (340119 * sqrt(55) * quick_pow(costheta, 2) * quick_pow(sintheta, 7)) / 65536. +
               (14549535 * sqrt(55) * quick_pow(costheta, 4) * quick_pow(sintheta, 7)) / 131072. -
               (37791 * sqrt(55) * quick_pow(sintheta, 9)) / 262144. -
               (4849845 * sqrt(55) * quick_pow(costheta, 2) * quick_pow(sintheta, 9)) / 262144. +
               (88179 * sqrt(55) * quick_pow(sintheta, 11)) / 262144.;
    }
    else if (lp == 11 && l == 11 && s == 1 && j == 10)
    {
        return (-14553 * sintheta) / (131072. * sqrt2) -
               (135135 * quick_pow(costheta, 2) * sintheta) / (131072. * sqrt2) -
               (804375 * quick_pow(costheta, 4) * sintheta) / (262144. * sqrt2) -
               (1786785 * quick_pow(costheta, 6) * sintheta) / (262144. * sqrt2) -
               (3741309 * quick_pow(costheta, 8) * sintheta) / (262144. * sqrt2) -
               (10669659 * quick_pow(costheta, 10) * sintheta) / (262144. * sqrt2) +
               (45045 * quick_pow(sintheta, 3)) / (131072. * sqrt2) +
               (804375 * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / (131072. * sqrt2) +
               (8933925 * quick_pow(costheta, 4) * quick_pow(sintheta, 3)) / (262144. * sqrt2) +
               (8729721 * quick_pow(costheta, 6) * quick_pow(sintheta, 3)) / (65536. * sqrt2) +
               (160044885 * quick_pow(costheta, 8) * quick_pow(sintheta, 3)) / (262144. * sqrt2) -
               (160875 * quick_pow(sintheta, 5)) / (262144. * sqrt2) -
               (5360355 * quick_pow(costheta, 2) * quick_pow(sintheta, 5)) / (262144. * sqrt2) -
               (26189163 * quick_pow(costheta, 4) * quick_pow(sintheta, 5)) / (131072. * sqrt2) -
               (224062839 * quick_pow(costheta, 6) * quick_pow(sintheta, 5)) / (131072. * sqrt2) +
               (255255 * quick_pow(sintheta, 7)) / (262144. * sqrt2) +
               (3741309 * quick_pow(costheta, 2) * quick_pow(sintheta, 7)) / (65536. * sqrt2) +
               (160044885 * quick_pow(costheta, 4) * quick_pow(sintheta, 7)) / (131072. * sqrt2) -
               (415701 * quick_pow(sintheta, 9)) / (262144. * sqrt2) -
               (53348295 * quick_pow(costheta, 2) * quick_pow(sintheta, 9)) / (262144. * sqrt2) +
               (969969 * quick_pow(sintheta, 11)) / (262144. * sqrt2);
    }

    return 0;
}

double coeff_mpm(double t, int lp, int l, int s, int j)
{
    double costheta = t;
    double sintheta = sqrt(1.0 - t * t);
    if (lp == 0 && l == 0 && s == 0 && j == 0)
    {
        return 0;
    }
    else if (lp == 1 && l == 1 && s == 1 && j == 0)
    {
        return 0;
    }
    else if (lp == 1 && l == 1 && s == 0 && j == 1)
    {
        return 0;
    }
    else if (lp == 1 && l == 1 && s == 1 && j == 1)
    {
        return 0;
    }
    else if (lp == 0 && l == 0 && s == 1 && j == 1)
    {
        return 0;
    }
    else if (lp == 0 && l == 2 && s == 1 && j == 1)
    {
        return 0;
    }
    else if (lp == 2 && l == 0 && s == 1 && j == 1)
    {
        return (-3 * quick_pow(sintheta, 2)) / (2. * sqrt2);
    }
    else if (lp == 2 && l == 2 && s == 1 && j == 1)
    {
        return (3 * quick_pow(sintheta, 2)) / 4.;
    }
    else if (lp == 2 && l == 2 && s == 0 && j == 2)
    {
        return 0;
    }
    else if (lp == 2 && l == 2 && s == 1 && j == 2)
    {
        return (-5 * quick_pow(sintheta, 2)) / 4.;
    }
    else if (lp == 1 && l == 1 && s == 1 && j == 2)
    {
        return 0;
    }
    else if (lp == 1 && l == 3 && s == 1 && j == 2)
    {
        return 0;
    }
    else if (lp == 3 && l == 1 && s == 1 && j == 2)
    {
        return (-5 * sqrt(1.5) * costheta * quick_pow(sintheta, 2)) / 2.;
    }
    else if (lp == 3 && l == 3 && s == 1 && j == 2)
    {
        return (5 * costheta * quick_pow(sintheta, 2)) / 2.;
    }
    else if (lp == 3 && l == 3 && s == 0 && j == 3)
    {
        return 0;
    }
    else if (lp == 3 && l == 3 && s == 1 && j == 3)
    {
        return (-35 * costheta * quick_pow(sintheta, 2)) / 8.;
    }
    else if (lp == 2 && l == 2 && s == 1 && j == 3)
    {
        return quick_pow(sintheta, 2) / 2.;
    }
    else if (lp == 2 && l == 4 && s == 1 && j == 3)
    {
        return -0.25 * (sqrt3 * quick_pow(sintheta, 2));
    }
    else if (lp == 4 && l == 2 && s == 1 && j == 3)
    {
        return (-15 * sqrt3) / 64. - (5 * sqrt3 * quick_pow(costheta, 2)) / 16. +
               (35 * sqrt3 * quick_pow(costheta, 4)) / 64. + (5 * sqrt3 * quick_pow(sintheta, 2)) / 16. -
               (105 * sqrt3 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 32. +
               (35 * sqrt3 * quick_pow(sintheta, 4)) / 64.;
    }
    else if (lp == 4 && l == 4 && s == 1 && j == 3)
    {
        return (-15 * quick_pow(sintheta, 2)) / 16. + (105 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 16.;
    }
    else if (lp == 4 && l == 4 && s == 0 && j == 4)
    {
        return 0;
    }
    else if (lp == 4 && l == 4 && s == 1 && j == 4)
    {
        return (27 * quick_pow(sintheta, 2)) / 16. - (189 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 16.;
    }
    else if (lp == 3 && l == 3 && s == 1 && j == 4)
    {
        return (15 * costheta * quick_pow(sintheta, 2)) / 8.;
    }
    else if (lp == 3 && l == 5 && s == 1 && j == 4)
    {
        return (-3 * sqrt5 * costheta * quick_pow(sintheta, 2)) / 4.;
    }
    else if (lp == 5 && l == 3 && s == 1 && j == 4)
    {
        return (-21 * sqrt5 * costheta) / 64. - (21 * sqrt5 * quick_pow(costheta, 3)) / 128. +
               (63 * sqrt5 * quick_pow(costheta, 5)) / 128. + (63 * sqrt5 * costheta * quick_pow(sintheta, 2)) / 128. -
               (315 * sqrt5 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 64. +
               (315 * sqrt5 * costheta * quick_pow(sintheta, 4)) / 128.;
    }
    else if (lp == 5 && l == 5 && s == 1 && j == 4)
    {
        return (21 * costheta) / 32. + (21 * quick_pow(costheta, 3)) / 64. - (63 * quick_pow(costheta, 5)) / 64. -
               (63 * costheta * quick_pow(sintheta, 2)) / 64. +
               (315 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 32. -
               (315 * costheta * quick_pow(sintheta, 4)) / 64.;
    }
    else if (lp == 5 && l == 5 && s == 0 && j == 5)
    {
        return 0;
    }
    else if (lp == 5 && l == 5 && s == 1 && j == 5)
    {
        return (-77 * costheta) / 64. - (77 * quick_pow(costheta, 3)) / 128. + (231 * quick_pow(costheta, 5)) / 128. +
               (231 * costheta * quick_pow(sintheta, 2)) / 128. -
               (1155 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 64. +
               (1155 * costheta * quick_pow(sintheta, 4)) / 128.;
    }
    else if (lp == 4 && l == 4 && s == 1 && j == 5)
    {
        return (-3 * quick_pow(sintheta, 2)) / 4. + (21 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 4.;
    }
    else if (lp == 4 && l == 6 && s == 1 && j == 5)
    {
        return (-3 * sqrt(7.5)) / 32. - (sqrt(7.5) * quick_pow(costheta, 2)) / 8. +
               (7 * sqrt(7.5) * quick_pow(costheta, 4)) / 32. + (sqrt(7.5) * quick_pow(sintheta, 2)) / 8. -
               (21 * sqrt(7.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 16. +
               (7 * sqrt(7.5) * quick_pow(sintheta, 4)) / 32.;
    }
    else if (lp == 6 && l == 4 && s == 1 && j == 5)
    {
        return (-35 * sqrt(7.5)) / 256. - (119 * sqrt(7.5) * quick_pow(costheta, 2)) / 512. -
               (21 * sqrt(7.5) * quick_pow(costheta, 4)) / 256. + (231 * sqrt(7.5) * quick_pow(costheta, 6)) / 512. +
               (119 * sqrt(7.5) * quick_pow(sintheta, 2)) / 512. +
               (63 * sqrt(7.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 128. -
               (3465 * sqrt(7.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 512. -
               (21 * sqrt(7.5) * quick_pow(sintheta, 4)) / 256. +
               (3465 * sqrt(7.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 512. -
               (231 * sqrt(7.5) * quick_pow(sintheta, 6)) / 512.;
    }
    else if (lp == 6 && l == 6 && s == 1 && j == 5)
    {
        return 0.341796875 + (595 * quick_pow(costheta, 2)) / 1024. + (105 * quick_pow(costheta, 4)) / 512. -
               (1155 * quick_pow(costheta, 6)) / 1024. - (595 * quick_pow(sintheta, 2)) / 1024. -
               (315 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 256. +
               (17325 * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 1024. +
               (105 * quick_pow(sintheta, 4)) / 512. -
               (17325 * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 1024. +
               (1155 * quick_pow(sintheta, 6)) / 1024.;
    }
    else if (lp == 6 && l == 6 && s == 0 && j == 6)
    {
        return 0;
    }
    else if (lp == 6 && l == 6 && s == 1 && j == 6)
    {
        return -0.634765625 - (1105 * quick_pow(costheta, 2)) / 1024. - (195 * quick_pow(costheta, 4)) / 512. +
               (2145 * quick_pow(costheta, 6)) / 1024. + (1105 * quick_pow(sintheta, 2)) / 1024. +
               (585 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 256. -
               (32175 * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 1024. -
               (195 * quick_pow(sintheta, 4)) / 512. +
               (32175 * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 1024. -
               (2145 * quick_pow(sintheta, 6)) / 1024.;
    }
    else if (lp == 5 && l == 5 && s == 1 && j == 6)
    {
        return (35 * costheta) / 64. + (35 * quick_pow(costheta, 3)) / 128. - (105 * quick_pow(costheta, 5)) / 128. -
               (105 * costheta * quick_pow(sintheta, 2)) / 128. +
               (525 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 64. -
               (525 * costheta * quick_pow(sintheta, 4)) / 128.;
    }
    else if (lp == 5 && l == 7 && s == 1 && j == 6)
    {
        return (-5 * sqrt(10.5) * costheta) / 32. - (5 * sqrt(10.5) * quick_pow(costheta, 3)) / 64. +
               (15 * sqrt(10.5) * quick_pow(costheta, 5)) / 64. +
               (15 * sqrt(10.5) * costheta * quick_pow(sintheta, 2)) / 64. -
               (75 * sqrt(10.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 32. +
               (75 * sqrt(10.5) * costheta * quick_pow(sintheta, 4)) / 64.;
    }
    else if (lp == 7 && l == 5 && s == 1 && j == 6)
    {
        return (-225 * sqrt(10.5) * costheta) / 1024. - (171 * sqrt(10.5) * quick_pow(costheta, 3)) / 1024. -
               (33 * sqrt(10.5) * quick_pow(costheta, 5)) / 1024. +
               (429 * sqrt(10.5) * quick_pow(costheta, 7)) / 1024. +
               (513 * sqrt(10.5) * costheta * quick_pow(sintheta, 2)) / 1024. +
               (165 * sqrt(10.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 512. -
               (9009 * sqrt(10.5) * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 1024. -
               (165 * sqrt(10.5) * costheta * quick_pow(sintheta, 4)) / 1024. +
               (15015 * sqrt(10.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 1024. -
               (3003 * sqrt(10.5) * costheta * quick_pow(sintheta, 6)) / 1024.;
    }
    else if (lp == 7 && l == 7 && s == 1 && j == 6)
    {
        return (675 * costheta) / 1024. + (513 * quick_pow(costheta, 3)) / 1024. +
               (99 * quick_pow(costheta, 5)) / 1024. - (1287 * quick_pow(costheta, 7)) / 1024. -
               (1539 * costheta * quick_pow(sintheta, 2)) / 1024. -
               (495 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 512. +
               (27027 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 1024. +
               (495 * costheta * quick_pow(sintheta, 4)) / 1024. -
               (45045 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 1024. +
               (9009 * costheta * quick_pow(sintheta, 6)) / 1024.;
    }
    else if (lp == 7 && l == 7 && s == 0 && j == 7)
    {
        return 0;
    }
    else if (lp == 7 && l == 7 && s == 1 && j == 7)
    {
        return (-10125 * costheta) / 8192. - (7695 * quick_pow(costheta, 3)) / 8192. -
               (1485 * quick_pow(costheta, 5)) / 8192. + (19305 * quick_pow(costheta, 7)) / 8192. +
               (23085 * costheta * quick_pow(sintheta, 2)) / 8192. +
               (7425 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 4096. -
               (405405 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 8192. -
               (7425 * costheta * quick_pow(sintheta, 4)) / 8192. +
               (675675 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 8192. -
               (135135 * costheta * quick_pow(sintheta, 6)) / 8192.;
    }
    else if (lp == 6 && l == 6 && s == 1 && j == 7)
    {
        return 0.29296875 + (255 * quick_pow(costheta, 2)) / 512. + (45 * quick_pow(costheta, 4)) / 256. -
               (495 * quick_pow(costheta, 6)) / 512. - (255 * quick_pow(sintheta, 2)) / 512. -
               (135 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 128. +
               (7425 * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 512. + (45 * quick_pow(sintheta, 4)) / 256. -
               (7425 * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 512. + (495 * quick_pow(sintheta, 6)) / 512.;
    }
    else if (lp == 6 && l == 8 && s == 1 && j == 7)
    {
        return (-75 * sqrt(3.5)) / 512. - (255 * sqrt(3.5) * quick_pow(costheta, 2)) / 1024. -
               (45 * sqrt(3.5) * quick_pow(costheta, 4)) / 512. + (495 * sqrt(3.5) * quick_pow(costheta, 6)) / 1024. +
               (255 * sqrt(3.5) * quick_pow(sintheta, 2)) / 1024. +
               (135 * sqrt(3.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 256. -
               (7425 * sqrt(3.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 1024. -
               (45 * sqrt(3.5) * quick_pow(sintheta, 4)) / 512. +
               (7425 * sqrt(3.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 1024. -
               (495 * sqrt(3.5) * quick_pow(sintheta, 6)) / 1024.;
    }
    else if (lp == 8 && l == 6 && s == 1 && j == 7)
    {
        return (-1575 * sqrt(3.5)) / 8192. - (45 * sqrt(3.5) * quick_pow(costheta, 2)) / 128. -
               (495 * sqrt(3.5) * quick_pow(costheta, 4)) / 2048. +
               (6435 * sqrt(3.5) * quick_pow(costheta, 8)) / 8192. + (45 * sqrt(3.5) * quick_pow(sintheta, 2)) / 128. +
               (1485 * sqrt(3.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 1024. -
               (45045 * sqrt(3.5) * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 2048. -
               (495 * sqrt(3.5) * quick_pow(sintheta, 4)) / 2048. +
               (225225 * sqrt(3.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 4096. -
               (45045 * sqrt(3.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 2048. +
               (6435 * sqrt(3.5) * quick_pow(sintheta, 8)) / 8192.;
    }
    else if (lp == 8 && l == 8 && s == 1 && j == 7)
    {
        return 0.336456298828125 + (315 * quick_pow(costheta, 2)) / 512. + (3465 * quick_pow(costheta, 4)) / 8192. -
               (45045 * quick_pow(costheta, 8)) / 32768. - (315 * quick_pow(sintheta, 2)) / 512. -
               (10395 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 4096. +
               (315315 * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 8192. +
               (3465 * quick_pow(sintheta, 4)) / 8192. -
               (1576575 * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 16384. +
               (315315 * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 8192. -
               (45045 * quick_pow(sintheta, 8)) / 32768.;
    }
    else if (lp == 8 && l == 8 && s == 0 && j == 8)
    {
        return 0;
    }
    else if (lp == 8 && l == 8 && s == 1 && j == 8)
    {
        return -0.635528564453125 - (595 * quick_pow(costheta, 2)) / 512. - (6545 * quick_pow(costheta, 4)) / 8192. +
               (85085 * quick_pow(costheta, 8)) / 32768. + (595 * quick_pow(sintheta, 2)) / 512. +
               (19635 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 4096. -
               (595595 * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 8192. -
               (6545 * quick_pow(sintheta, 4)) / 8192. +
               (2977975 * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 16384. -
               (595595 * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 8192. +
               (85085 * quick_pow(sintheta, 8)) / 32768.;
    }
    else if (lp == 7 && l == 7 && s == 1 && j == 8)
    {
        return (4725 * costheta) / 8192. + (3591 * quick_pow(costheta, 3)) / 8192. +
               (693 * quick_pow(costheta, 5)) / 8192. - (9009 * quick_pow(costheta, 7)) / 8192. -
               (10773 * costheta * quick_pow(sintheta, 2)) / 8192. -
               (3465 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 4096. +
               (189189 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 8192. +
               (3465 * costheta * quick_pow(sintheta, 4)) / 8192. -
               (315315 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 8192. +
               (63063 * costheta * quick_pow(sintheta, 6)) / 8192.;
    }
    else if (lp == 7 && l == 9 && s == 1 && j == 8)
    {
        return (-1575 * costheta) / (2048. * sqrt2) - (1197 * quick_pow(costheta, 3)) / (2048. * sqrt2) -
               (231 * quick_pow(costheta, 5)) / (2048. * sqrt2) + (3003 * quick_pow(costheta, 7)) / (2048. * sqrt2) +
               (3591 * costheta * quick_pow(sintheta, 2)) / (2048. * sqrt2) +
               (1155 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / (1024. * sqrt2) -
               (63063 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / (2048. * sqrt2) -
               (1155 * costheta * quick_pow(sintheta, 4)) / (2048. * sqrt2) +
               (105105 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / (2048. * sqrt2) -
               (21021 * costheta * quick_pow(sintheta, 6)) / (2048. * sqrt2);
    }
    else if (lp == 9 && l == 7 && s == 1 && j == 8)
    {
        return (-8085 * costheta) / (8192. * sqrt2) - (3465 * quick_pow(costheta, 3)) / (4096. * sqrt2) -
               (2145 * quick_pow(costheta, 5)) / (4096. * sqrt2) + (2145 * quick_pow(costheta, 7)) / (16384. * sqrt2) +
               (36465 * quick_pow(costheta, 9)) / (16384. * sqrt2) +
               (10395 * costheta * quick_pow(sintheta, 2)) / (4096. * sqrt2) +
               (10725 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / (2048. * sqrt2) -
               (45045 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / (16384. * sqrt2) -
               (328185 * quick_pow(costheta, 7) * quick_pow(sintheta, 2)) / (4096. * sqrt2) -
               (10725 * costheta * quick_pow(sintheta, 4)) / (4096. * sqrt2) +
               (75075 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / (16384. * sqrt2) +
               (2297295 * quick_pow(costheta, 5) * quick_pow(sintheta, 4)) / (8192. * sqrt2) -
               (15015 * costheta * quick_pow(sintheta, 6)) / (16384. * sqrt2) -
               (765765 * quick_pow(costheta, 3) * quick_pow(sintheta, 6)) / (4096. * sqrt2) +
               (328185 * costheta * quick_pow(sintheta, 8)) / (16384. * sqrt2);
    }
    else if (lp == 9 && l == 9 && s == 1 && j == 8)
    {
        return (2695 * costheta) / 4096. + (1155 * quick_pow(costheta, 3)) / 2048. +
               (715 * quick_pow(costheta, 5)) / 2048. - (715 * quick_pow(costheta, 7)) / 8192. -
               (12155 * quick_pow(costheta, 9)) / 8192. - (3465 * costheta * quick_pow(sintheta, 2)) / 2048. -
               (3575 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 1024. +
               (15015 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 8192. +
               (109395 * quick_pow(costheta, 7) * quick_pow(sintheta, 2)) / 2048. +
               (3575 * costheta * quick_pow(sintheta, 4)) / 2048. -
               (25025 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 8192. -
               (765765 * quick_pow(costheta, 5) * quick_pow(sintheta, 4)) / 4096. +
               (5005 * costheta * quick_pow(sintheta, 6)) / 8192. +
               (255255 * quick_pow(costheta, 3) * quick_pow(sintheta, 6)) / 2048. -
               (109395 * costheta * quick_pow(sintheta, 8)) / 8192.;
    }
    else if (lp == 9 && l == 9 && s == 0 && j == 9)
    {
        return 0;
    }
    else if (lp == 9 && l == 9 && s == 1 && j == 9)
    {
        return (-10241 * costheta) / 8192. - (4389 * quick_pow(costheta, 3)) / 4096. -
               (2717 * quick_pow(costheta, 5)) / 4096. + (2717 * quick_pow(costheta, 7)) / 16384. +
               (46189 * quick_pow(costheta, 9)) / 16384. + (13167 * costheta * quick_pow(sintheta, 2)) / 4096. +
               (13585 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 2048. -
               (57057 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 16384. -
               (415701 * quick_pow(costheta, 7) * quick_pow(sintheta, 2)) / 4096. -
               (13585 * costheta * quick_pow(sintheta, 4)) / 4096. +
               (95095 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 16384. +
               (2909907 * quick_pow(costheta, 5) * quick_pow(sintheta, 4)) / 8192. -
               (19019 * costheta * quick_pow(sintheta, 6)) / 16384. -
               (969969 * quick_pow(costheta, 3) * quick_pow(sintheta, 6)) / 4096. +
               (415701 * costheta * quick_pow(sintheta, 8)) / 16384.;
    }
    else if (lp == 8 && l == 8 && s == 1 && j == 9)
    {
        return 0.299072265625 + (35 * quick_pow(costheta, 2)) / 64. + (385 * quick_pow(costheta, 4)) / 1024. -
               (5005 * quick_pow(costheta, 8)) / 4096. - (35 * quick_pow(sintheta, 2)) / 64. -
               (1155 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 512. +
               (35035 * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 1024. +
               (385 * quick_pow(sintheta, 4)) / 1024. -
               (175175 * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 2048. +
               (35035 * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 1024. -
               (5005 * quick_pow(sintheta, 8)) / 4096.;
    }
    else if (lp == 8 && l == 10 && s == 1 && j == 9)
    {
        return (-735 * sqrt(2.5)) / 4096. - (21 * sqrt(2.5) * quick_pow(costheta, 2)) / 64. -
               (231 * sqrt(2.5) * quick_pow(costheta, 4)) / 1024. +
               (3003 * sqrt(2.5) * quick_pow(costheta, 8)) / 4096. + (21 * sqrt(2.5) * quick_pow(sintheta, 2)) / 64. +
               (693 * sqrt(2.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 512. -
               (21021 * sqrt(2.5) * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 1024. -
               (231 * sqrt(2.5) * quick_pow(sintheta, 4)) / 1024. +
               (105105 * sqrt(2.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 2048. -
               (21021 * sqrt(2.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 1024. +
               (3003 * sqrt(2.5) * quick_pow(sintheta, 8)) / 4096.;
    }
    else if (lp == 10 && l == 8 && s == 1 && j == 9)
    {
        return (-14553 * sqrt(2.5)) / 65536. - (27489 * sqrt(2.5) * quick_pow(costheta, 2)) / 65536. -
               (5577 * sqrt(2.5) * quick_pow(costheta, 4)) / 16384. -
               (24453 * sqrt(2.5) * quick_pow(costheta, 6)) / 131072. +
               (7293 * sqrt(2.5) * quick_pow(costheta, 8)) / 65536. +
               (138567 * sqrt(2.5) * quick_pow(costheta, 10)) / 131072. +
               (27489 * sqrt(2.5) * quick_pow(sintheta, 2)) / 65536. +
               (16731 * sqrt(2.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 8192. +
               (366795 * sqrt(2.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 131072. -
               (51051 * sqrt(2.5) * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 16384. -
               (6235515 * sqrt(2.5) * quick_pow(costheta, 8) * quick_pow(sintheta, 2)) / 131072. -
               (5577 * sqrt(2.5) * quick_pow(sintheta, 4)) / 16384. -
               (366795 * sqrt(2.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 131072. +
               (255255 * sqrt(2.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 32768. +
               (14549535 * sqrt(2.5) * quick_pow(costheta, 6) * quick_pow(sintheta, 4)) / 65536. +
               (24453 * sqrt(2.5) * quick_pow(sintheta, 6)) / 131072. -
               (51051 * sqrt(2.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 16384. -
               (14549535 * sqrt(2.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 6)) / 65536. +
               (7293 * sqrt(2.5) * quick_pow(sintheta, 8)) / 65536. +
               (6235515 * sqrt(2.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 8)) / 131072. -
               (138567 * sqrt(2.5) * quick_pow(sintheta, 10)) / 131072.;
    }
    else if (lp == 10 && l == 10 && s == 1 && j == 9)
    {
        return 0.33309173583984375 + (82467 * quick_pow(costheta, 2)) / 131072. +
               (16731 * quick_pow(costheta, 4)) / 32768. + (73359 * quick_pow(costheta, 6)) / 262144. -
               (21879 * quick_pow(costheta, 8)) / 131072. - (415701 * quick_pow(costheta, 10)) / 262144. -
               (82467 * quick_pow(sintheta, 2)) / 131072. -
               (50193 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 16384. -
               (1100385 * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 262144. +
               (153153 * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 32768. +
               (18706545 * quick_pow(costheta, 8) * quick_pow(sintheta, 2)) / 262144. +
               (16731 * quick_pow(sintheta, 4)) / 32768. +
               (1100385 * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 262144. -
               (765765 * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 65536. -
               (43648605 * quick_pow(costheta, 6) * quick_pow(sintheta, 4)) / 131072. -
               (73359 * quick_pow(sintheta, 6)) / 262144. +
               (153153 * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 32768. +
               (43648605 * quick_pow(costheta, 4) * quick_pow(sintheta, 6)) / 131072. -
               (21879 * quick_pow(sintheta, 8)) / 131072. -
               (18706545 * quick_pow(costheta, 2) * quick_pow(sintheta, 8)) / 262144. +
               (415701 * quick_pow(sintheta, 10)) / 262144.;
    }
    else if (lp == 10 && l == 10 && s == 0 && j == 10)
    {
        return 0;
    }
    else if (lp == 10 && l == 10 && s == 1 && j == 10)
    {
        return -0.6359024047851562 - (157437 * quick_pow(costheta, 2)) / 131072. -
               (31941 * quick_pow(costheta, 4)) / 32768. - (140049 * quick_pow(costheta, 6)) / 262144. +
               (41769 * quick_pow(costheta, 8)) / 131072. + (793611 * quick_pow(costheta, 10)) / 262144. +
               (157437 * quick_pow(sintheta, 2)) / 131072. +
               (95823 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 16384. +
               (2100735 * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 262144. -
               (292383 * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 32768. -
               (35712495 * quick_pow(costheta, 8) * quick_pow(sintheta, 2)) / 262144. -
               (31941 * quick_pow(sintheta, 4)) / 32768. -
               (2100735 * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 262144. +
               (1461915 * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 65536. +
               (83329155 * quick_pow(costheta, 6) * quick_pow(sintheta, 4)) / 131072. +
               (140049 * quick_pow(sintheta, 6)) / 262144. -
               (292383 * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 32768. -
               (83329155 * quick_pow(costheta, 4) * quick_pow(sintheta, 6)) / 131072. +
               (41769 * quick_pow(sintheta, 8)) / 131072. +
               (35712495 * quick_pow(costheta, 2) * quick_pow(sintheta, 8)) / 262144. -
               (793611 * quick_pow(sintheta, 10)) / 262144.;
    }
    else if (lp == 9 && l == 9 && s == 1 && j == 10)
    {
        return (4851 * costheta) / 8192. + (2079 * quick_pow(costheta, 3)) / 4096. +
               (1287 * quick_pow(costheta, 5)) / 4096. - (1287 * quick_pow(costheta, 7)) / 16384. -
               (21879 * quick_pow(costheta, 9)) / 16384. - (6237 * costheta * quick_pow(sintheta, 2)) / 4096. -
               (6435 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 2048. +
               (27027 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 16384. +
               (196911 * quick_pow(costheta, 7) * quick_pow(sintheta, 2)) / 4096. +
               (6435 * costheta * quick_pow(sintheta, 4)) / 4096. -
               (45045 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 16384. -
               (1378377 * quick_pow(costheta, 5) * quick_pow(sintheta, 4)) / 8192. +
               (9009 * costheta * quick_pow(sintheta, 6)) / 16384. +
               (459459 * quick_pow(costheta, 3) * quick_pow(sintheta, 6)) / 4096. -
               (196911 * costheta * quick_pow(sintheta, 8)) / 16384.;
    }
    else if (lp == 9 && l == 11 && s == 1 && j == 10)
    {
        return (-441 * sqrt(27.5) * costheta) / 4096. - (189 * sqrt(27.5) * quick_pow(costheta, 3)) / 2048. -
               (117 * sqrt(27.5) * quick_pow(costheta, 5)) / 2048. +
               (117 * sqrt(27.5) * quick_pow(costheta, 7)) / 8192. +
               (1989 * sqrt(27.5) * quick_pow(costheta, 9)) / 8192. +
               (567 * sqrt(27.5) * costheta * quick_pow(sintheta, 2)) / 2048. +
               (585 * sqrt(27.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 1024. -
               (2457 * sqrt(27.5) * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 8192. -
               (17901 * sqrt(27.5) * quick_pow(costheta, 7) * quick_pow(sintheta, 2)) / 2048. -
               (585 * sqrt(27.5) * costheta * quick_pow(sintheta, 4)) / 2048. +
               (4095 * sqrt(27.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 8192. +
               (125307 * sqrt(27.5) * quick_pow(costheta, 5) * quick_pow(sintheta, 4)) / 4096. -
               (819 * sqrt(27.5) * costheta * quick_pow(sintheta, 6)) / 8192. -
               (41769 * sqrt(27.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 6)) / 2048. +
               (17901 * sqrt(27.5) * costheta * quick_pow(sintheta, 8)) / 8192.;
    }
    else if (lp == 11 && l == 9 && s == 1 && j == 10)
    {
        return (-17199 * sqrt(27.5) * costheta) / 131072. - (15561 * sqrt(27.5) * quick_pow(costheta, 3)) / 131072. -
               (23985 * sqrt(27.5) * quick_pow(costheta, 5)) / 262144. -
               (11271 * sqrt(27.5) * quick_pow(costheta, 7)) / 262144. +
               (12597 * sqrt(27.5) * quick_pow(costheta, 9)) / 262144. +
               (88179 * sqrt(27.5) * quick_pow(costheta, 11)) / 262144. +
               (46683 * sqrt(27.5) * costheta * quick_pow(sintheta, 2)) / 131072. +
               (119925 * sqrt(27.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 131072. +
               (236691 * sqrt(27.5) * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 262144. -
               (113373 * sqrt(27.5) * quick_pow(costheta, 7) * quick_pow(sintheta, 2)) / 65536. -
               (4849845 * sqrt(27.5) * quick_pow(costheta, 9) * quick_pow(sintheta, 2)) / 262144. -
               (119925 * sqrt(27.5) * costheta * quick_pow(sintheta, 4)) / 262144. -
               (394485 * sqrt(27.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 262144. +
               (793611 * sqrt(27.5) * quick_pow(costheta, 5) * quick_pow(sintheta, 4)) / 131072. +
               (14549535 * sqrt(27.5) * quick_pow(costheta, 7) * quick_pow(sintheta, 4)) / 131072. +
               (78897 * sqrt(27.5) * costheta * quick_pow(sintheta, 6)) / 262144. -
               (264537 * sqrt(27.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 6)) / 65536. -
               (20369349 * sqrt(27.5) * quick_pow(costheta, 5) * quick_pow(sintheta, 6)) / 131072. +
               (113373 * sqrt(27.5) * costheta * quick_pow(sintheta, 8)) / 262144. +
               (14549535 * sqrt(27.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 8)) / 262144. -
               (969969 * sqrt(27.5) * costheta * quick_pow(sintheta, 10)) / 262144.;
    }
    else if (lp == 11 && l == 11 && s == 1 && j == 10)
    {
        return (85995 * costheta) / 131072. + (77805 * quick_pow(costheta, 3)) / 131072. +
               (119925 * quick_pow(costheta, 5)) / 262144. + (56355 * quick_pow(costheta, 7)) / 262144. -
               (62985 * quick_pow(costheta, 9)) / 262144. - (440895 * quick_pow(costheta, 11)) / 262144. -
               (233415 * costheta * quick_pow(sintheta, 2)) / 131072. -
               (599625 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 131072. -
               (1183455 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 262144. +
               (566865 * quick_pow(costheta, 7) * quick_pow(sintheta, 2)) / 65536. +
               (24249225 * quick_pow(costheta, 9) * quick_pow(sintheta, 2)) / 262144. +
               (599625 * costheta * quick_pow(sintheta, 4)) / 262144. +
               (1972425 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 262144. -
               (3968055 * quick_pow(costheta, 5) * quick_pow(sintheta, 4)) / 131072. -
               (72747675 * quick_pow(costheta, 7) * quick_pow(sintheta, 4)) / 131072. -
               (394485 * costheta * quick_pow(sintheta, 6)) / 262144. +
               (1322685 * quick_pow(costheta, 3) * quick_pow(sintheta, 6)) / 65536. +
               (101846745 * quick_pow(costheta, 5) * quick_pow(sintheta, 6)) / 131072. -
               (566865 * costheta * quick_pow(sintheta, 8)) / 262144. -
               (72747675 * quick_pow(costheta, 3) * quick_pow(sintheta, 8)) / 262144. +
               (4849845 * costheta * quick_pow(sintheta, 10)) / 262144.;
    }

    return 0;
}

double coeff_m01(double t, int lp, int l, int s, int j)
{
    double costheta = t;
    double sintheta = sqrt(1.0 - t * t);
    if (lp == 0 && l == 0 && s == 0 && j == 0)
    {
        return 0;
    }
    else if (lp == 1 && l == 1 && s == 1 && j == 0)
    {
        return 0;
    }
    else if (lp == 1 && l == 1 && s == 0 && j == 1)
    {
        return 0;
    }
    else if (lp == 1 && l == 1 && s == 1 && j == 1)
    {
        return (3 * sintheta) / (2. * sqrt2);
    }
    else if (lp == 0 && l == 0 && s == 1 && j == 1)
    {
        return 0;
    }
    else if (lp == 0 && l == 2 && s == 1 && j == 1)
    {
        return 0;
    }
    else if (lp == 2 && l == 0 && s == 1 && j == 1)
    {
        return (-3 * costheta * sintheta) / 2.;
    }
    else if (lp == 2 && l == 2 && s == 1 && j == 1)
    {
        return (3 * costheta * sintheta) / (2. * sqrt2);
    }
    else if (lp == 2 && l == 2 && s == 0 && j == 2)
    {
        return 0;
    }
    else if (lp == 2 && l == 2 && s == 1 && j == 2)
    {
        return (5 * costheta * sintheta) / (2. * sqrt2);
    }
    else if (lp == 1 && l == 1 && s == 1 && j == 2)
    {
        return (-3 * sintheta) / (2. * sqrt2);
    }
    else if (lp == 1 && l == 3 && s == 1 && j == 2)
    {
        return (sqrt3 * sintheta) / 2.;
    }
    else if (lp == 3 && l == 1 && s == 1 && j == 2)
    {
        return -0.125 * (sqrt3 * sintheta) - (15 * sqrt3 * quick_pow(costheta, 2) * sintheta) / 8. +
               (5 * sqrt3 * quick_pow(sintheta, 3)) / 8.;
    }
    else if (lp == 3 && l == 3 && s == 1 && j == 2)
    {
        return sintheta / (4. * sqrt2) + (15 * quick_pow(costheta, 2) * sintheta) / (4. * sqrt2) -
               (5 * quick_pow(sintheta, 3)) / (4. * sqrt2);
    }
    else if (lp == 3 && l == 3 && s == 0 && j == 3)
    {
        return 0;
    }
    else if (lp == 3 && l == 3 && s == 1 && j == 3)
    {
        return (7 * sintheta) / (32. * sqrt2) + (105 * quick_pow(costheta, 2) * sintheta) / (32. * sqrt2) -
               (35 * quick_pow(sintheta, 3)) / (32. * sqrt2);
    }
    else if (lp == 2 && l == 2 && s == 1 && j == 3)
    {
        return -2 * sqrt2 * costheta * sintheta;
    }
    else if (lp == 2 && l == 4 && s == 1 && j == 3)
    {
        return sqrt(6) * costheta * sintheta;
    }
    else if (lp == 4 && l == 2 && s == 1 && j == 3)
    {
        return (-5 * sqrt(1.5) * costheta * sintheta) / 8. - (35 * sqrt(1.5) * quick_pow(costheta, 3) * sintheta) / 8. +
               (35 * sqrt(1.5) * costheta * quick_pow(sintheta, 3)) / 8.;
    }
    else if (lp == 4 && l == 4 && s == 1 && j == 3)
    {
        return (15 * costheta * sintheta) / (16. * sqrt2) + (105 * quick_pow(costheta, 3) * sintheta) / (16. * sqrt2) -
               (105 * costheta * quick_pow(sintheta, 3)) / (16. * sqrt2);
    }
    else if (lp == 4 && l == 4 && s == 0 && j == 4)
    {
        return 0;
    }
    else if (lp == 4 && l == 4 && s == 1 && j == 4)
    {
        return (9 * costheta * sintheta) / (16. * sqrt2) + (63 * quick_pow(costheta, 3) * sintheta) / (16. * sqrt2) -
               (63 * costheta * quick_pow(sintheta, 3)) / (16. * sqrt2);
    }
    else if (lp == 3 && l == 3 && s == 1 && j == 4)
    {
        return (-15 * sintheta) / (32. * sqrt2) - (225 * quick_pow(costheta, 2) * sintheta) / (32. * sqrt2) +
               (75 * quick_pow(sintheta, 3)) / (32. * sqrt2);
    }
    else if (lp == 3 && l == 5 && s == 1 && j == 4)
    {
        return (3 * sqrt(2.5) * sintheta) / 16. + (45 * sqrt(2.5) * quick_pow(costheta, 2) * sintheta) / 16. -
               (15 * sqrt(2.5) * quick_pow(sintheta, 3)) / 16.;
    }
    else if (lp == 5 && l == 3 && s == 1 && j == 4)
    {
        return (-3 * sqrt(2.5) * sintheta) / 32. - (63 * sqrt(2.5) * quick_pow(costheta, 2) * sintheta) / 64. -
               (315 * sqrt(2.5) * quick_pow(costheta, 4) * sintheta) / 64. +
               (21 * sqrt(2.5) * quick_pow(sintheta, 3)) / 64. +
               (315 * sqrt(2.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / 32. -
               (63 * sqrt(2.5) * quick_pow(sintheta, 5)) / 64.;
    }
    else if (lp == 5 && l == 5 && s == 1 && j == 4)
    {
        return (3 * sintheta) / (16. * sqrt2) + (63 * quick_pow(costheta, 2) * sintheta) / (32. * sqrt2) +
               (315 * quick_pow(costheta, 4) * sintheta) / (32. * sqrt2) -
               (21 * quick_pow(sintheta, 3)) / (32. * sqrt2) -
               (315 * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / (16. * sqrt2) +
               (63 * quick_pow(sintheta, 5)) / (32. * sqrt2);
    }
    else if (lp == 5 && l == 5 && s == 0 && j == 5)
    {
        return 0;
    }
    else if (lp == 5 && l == 5 && s == 1 && j == 5)
    {
        return (11 * sintheta) / (128. * sqrt2) + (231 * quick_pow(costheta, 2) * sintheta) / (256. * sqrt2) +
               (1155 * quick_pow(costheta, 4) * sintheta) / (256. * sqrt2) -
               (77 * quick_pow(sintheta, 3)) / (256. * sqrt2) -
               (1155 * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / (128. * sqrt2) +
               (231 * quick_pow(sintheta, 5)) / (256. * sqrt2);
    }
    else if (lp == 4 && l == 4 && s == 1 && j == 5)
    {
        return (-3 * costheta * sintheta) / (2. * sqrt2) - (21 * quick_pow(costheta, 3) * sintheta) / (2. * sqrt2) +
               (21 * costheta * quick_pow(sintheta, 3)) / (2. * sqrt2);
    }
    else if (lp == 4 && l == 6 && s == 1 && j == 5)
    {
        return (sqrt15 * costheta * sintheta) / 4. + (7 * sqrt15 * quick_pow(costheta, 3) * sintheta) / 4. -
               (7 * sqrt15 * costheta * quick_pow(sintheta, 3)) / 4.;
    }
    else if (lp == 6 && l == 4 && s == 1 && j == 5)
    {
        return (-35 * sqrt15 * costheta * sintheta) / 256. - (21 * sqrt15 * quick_pow(costheta, 3) * sintheta) / 32. -
               (693 * sqrt15 * quick_pow(costheta, 5) * sintheta) / 256. +
               (21 * sqrt15 * costheta * quick_pow(sintheta, 3)) / 32. +
               (1155 * sqrt15 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / 128. -
               (693 * sqrt15 * costheta * quick_pow(sintheta, 5)) / 256.;
    }
    else if (lp == 6 && l == 6 && s == 1 && j == 5)
    {
        return (175 * costheta * sintheta) / (256. * sqrt2) +
               (105 * quick_pow(costheta, 3) * sintheta) / (32. * sqrt2) +
               (3465 * quick_pow(costheta, 5) * sintheta) / (256. * sqrt2) -
               (105 * costheta * quick_pow(sintheta, 3)) / (32. * sqrt2) -
               (5775 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / (128. * sqrt2) +
               (3465 * costheta * quick_pow(sintheta, 5)) / (256. * sqrt2);
    }
    else if (lp == 6 && l == 6 && s == 0 && j == 6)
    {
        return 0;
    }
    else if (lp == 6 && l == 6 && s == 1 && j == 6)
    {
        return (65 * costheta * sintheta) / (256. * sqrt2) + (39 * quick_pow(costheta, 3) * sintheta) / (32. * sqrt2) +
               (1287 * quick_pow(costheta, 5) * sintheta) / (256. * sqrt2) -
               (39 * costheta * quick_pow(sintheta, 3)) / (32. * sqrt2) -
               (2145 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / (128. * sqrt2) +
               (1287 * costheta * quick_pow(sintheta, 5)) / (256. * sqrt2);
    }
    else if (lp == 5 && l == 5 && s == 1 && j == 6)
    {
        return (-35 * sintheta) / (128. * sqrt2) - (735 * quick_pow(costheta, 2) * sintheta) / (256. * sqrt2) -
               (3675 * quick_pow(costheta, 4) * sintheta) / (256. * sqrt2) +
               (245 * quick_pow(sintheta, 3)) / (256. * sqrt2) +
               (3675 * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / (128. * sqrt2) -
               (735 * quick_pow(sintheta, 5)) / (256. * sqrt2);
    }
    else if (lp == 5 && l == 7 && s == 1 && j == 6)
    {
        return (5 * sqrt(21) * sintheta) / 128. + (105 * sqrt(21) * quick_pow(costheta, 2) * sintheta) / 256. +
               (525 * sqrt(21) * quick_pow(costheta, 4) * sintheta) / 256. -
               (35 * sqrt(21) * quick_pow(sintheta, 3)) / 256. -
               (525 * sqrt(21) * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / 128. +
               (105 * sqrt(21) * quick_pow(sintheta, 5)) / 256.;
    }
    else if (lp == 7 && l == 5 && s == 1 && j == 6)
    {
        return (-25 * sqrt(21) * sintheta) / 1024. - (243 * sqrt(21) * quick_pow(costheta, 2) * sintheta) / 1024. -
               (825 * sqrt(21) * quick_pow(costheta, 4) * sintheta) / 1024. -
               (3003 * sqrt(21) * quick_pow(costheta, 6) * sintheta) / 1024. +
               (81 * sqrt(21) * quick_pow(sintheta, 3)) / 1024. +
               (825 * sqrt(21) * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / 512. +
               (15015 * sqrt(21) * quick_pow(costheta, 4) * quick_pow(sintheta, 3)) / 1024. -
               (165 * sqrt(21) * quick_pow(sintheta, 5)) / 1024. -
               (9009 * sqrt(21) * quick_pow(costheta, 2) * quick_pow(sintheta, 5)) / 1024. +
               (429 * sqrt(21) * quick_pow(sintheta, 7)) / 1024.;
    }
    else if (lp == 7 && l == 7 && s == 1 && j == 6)
    {
        return (75 * sintheta) / (512. * sqrt2) + (729 * quick_pow(costheta, 2) * sintheta) / (512. * sqrt2) +
               (2475 * quick_pow(costheta, 4) * sintheta) / (512. * sqrt2) +
               (9009 * quick_pow(costheta, 6) * sintheta) / (512. * sqrt2) -
               (243 * quick_pow(sintheta, 3)) / (512. * sqrt2) -
               (2475 * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / (256. * sqrt2) -
               (45045 * quick_pow(costheta, 4) * quick_pow(sintheta, 3)) / (512. * sqrt2) +
               (495 * quick_pow(sintheta, 5)) / (512. * sqrt2) +
               (27027 * quick_pow(costheta, 2) * quick_pow(sintheta, 5)) / (512. * sqrt2) -
               (1287 * quick_pow(sintheta, 7)) / (512. * sqrt2);
    }
    else if (lp == 7 && l == 7 && s == 0 && j == 7)
    {
        return 0;
    }
    else if (lp == 7 && l == 7 && s == 1 && j == 7)
    {
        return (375 * sintheta) / (8192. * sqrt2) + (3645 * quick_pow(costheta, 2) * sintheta) / (8192. * sqrt2) +
               (12375 * quick_pow(costheta, 4) * sintheta) / (8192. * sqrt2) +
               (45045 * quick_pow(costheta, 6) * sintheta) / (8192. * sqrt2) -
               (1215 * quick_pow(sintheta, 3)) / (8192. * sqrt2) -
               (12375 * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / (4096. * sqrt2) -
               (225225 * quick_pow(costheta, 4) * quick_pow(sintheta, 3)) / (8192. * sqrt2) +
               (2475 * quick_pow(sintheta, 5)) / (8192. * sqrt2) +
               (135135 * quick_pow(costheta, 2) * quick_pow(sintheta, 5)) / (8192. * sqrt2) -
               (6435 * quick_pow(sintheta, 7)) / (8192. * sqrt2);
    }
    else if (lp == 6 && l == 6 && s == 1 && j == 7)
    {
        return (-15 * costheta * sintheta) / (16. * sqrt2) - (9 * quick_pow(costheta, 3) * sintheta) / (2. * sqrt2) -
               (297 * quick_pow(costheta, 5) * sintheta) / (16. * sqrt2) +
               (9 * costheta * quick_pow(sintheta, 3)) / (2. * sqrt2) +
               (495 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / (8. * sqrt2) -
               (297 * costheta * quick_pow(sintheta, 5)) / (16. * sqrt2);
    }
    else if (lp == 6 && l == 8 && s == 1 && j == 7)
    {
        return (15 * sqrt7 * costheta * sintheta) / 64. + (9 * sqrt7 * quick_pow(costheta, 3) * sintheta) / 8. +
               (297 * sqrt7 * quick_pow(costheta, 5) * sintheta) / 64. -
               (9 * sqrt7 * costheta * quick_pow(sintheta, 3)) / 8. -
               (495 * sqrt7 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / 32. +
               (297 * sqrt7 * costheta * quick_pow(sintheta, 5)) / 64.;
    }
    else if (lp == 8 && l == 6 && s == 1 && j == 7)
    {
        return (-315 * sqrt7 * costheta * sintheta) / 2048. -
               (693 * sqrt7 * quick_pow(costheta, 3) * sintheta) / 1024. -
               (3861 * sqrt7 * quick_pow(costheta, 5) * sintheta) / 2048. -
               (6435 * sqrt7 * quick_pow(costheta, 7) * sintheta) / 1024. +
               (693 * sqrt7 * costheta * quick_pow(sintheta, 3)) / 1024. +
               (6435 * sqrt7 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / 1024. +
               (45045 * sqrt7 * quick_pow(costheta, 5) * quick_pow(sintheta, 3)) / 1024. -
               (3861 * sqrt7 * costheta * quick_pow(sintheta, 5)) / 2048. -
               (45045 * sqrt7 * quick_pow(costheta, 3) * quick_pow(sintheta, 5)) / 1024. +
               (6435 * sqrt7 * costheta * quick_pow(sintheta, 7)) / 1024.;
    }
    else if (lp == 8 && l == 8 && s == 1 && j == 7)
    {
        return (2205 * costheta * sintheta) / (4096. * sqrt2) +
               (4851 * quick_pow(costheta, 3) * sintheta) / (2048. * sqrt2) +
               (27027 * quick_pow(costheta, 5) * sintheta) / (4096. * sqrt2) +
               (45045 * quick_pow(costheta, 7) * sintheta) / (2048. * sqrt2) -
               (4851 * costheta * quick_pow(sintheta, 3)) / (2048. * sqrt2) -
               (45045 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / (2048. * sqrt2) -
               (315315 * quick_pow(costheta, 5) * quick_pow(sintheta, 3)) / (2048. * sqrt2) +
               (27027 * costheta * quick_pow(sintheta, 5)) / (4096. * sqrt2) +
               (315315 * quick_pow(costheta, 3) * quick_pow(sintheta, 5)) / (2048. * sqrt2) -
               (45045 * costheta * quick_pow(sintheta, 7)) / (2048. * sqrt2);
    }
    else if (lp == 8 && l == 8 && s == 0 && j == 8)
    {
        return 0;
    }
    else if (lp == 8 && l == 8 && s == 1 && j == 8)
    {
        return (595 * costheta * sintheta) / (4096. * sqrt2) +
               (1309 * quick_pow(costheta, 3) * sintheta) / (2048. * sqrt2) +
               (7293 * quick_pow(costheta, 5) * sintheta) / (4096. * sqrt2) +
               (12155 * quick_pow(costheta, 7) * sintheta) / (2048. * sqrt2) -
               (1309 * costheta * quick_pow(sintheta, 3)) / (2048. * sqrt2) -
               (12155 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / (2048. * sqrt2) -
               (85085 * quick_pow(costheta, 5) * quick_pow(sintheta, 3)) / (2048. * sqrt2) +
               (7293 * costheta * quick_pow(sintheta, 5)) / (4096. * sqrt2) +
               (85085 * quick_pow(costheta, 3) * quick_pow(sintheta, 5)) / (2048. * sqrt2) -
               (12155 * costheta * quick_pow(sintheta, 7)) / (2048. * sqrt2);
    }
    else if (lp == 7 && l == 7 && s == 1 && j == 8)
    {
        return (-1575 * sintheta) / (8192. * sqrt2) - (15309 * quick_pow(costheta, 2) * sintheta) / (8192. * sqrt2) -
               (51975 * quick_pow(costheta, 4) * sintheta) / (8192. * sqrt2) -
               (189189 * quick_pow(costheta, 6) * sintheta) / (8192. * sqrt2) +
               (5103 * quick_pow(sintheta, 3)) / (8192. * sqrt2) +
               (51975 * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / (4096. * sqrt2) +
               (945945 * quick_pow(costheta, 4) * quick_pow(sintheta, 3)) / (8192. * sqrt2) -
               (10395 * quick_pow(sintheta, 5)) / (8192. * sqrt2) -
               (567567 * quick_pow(costheta, 2) * quick_pow(sintheta, 5)) / (8192. * sqrt2) +
               (27027 * quick_pow(sintheta, 7)) / (8192. * sqrt2);
    }
    else if (lp == 7 && l == 9 && s == 1 && j == 8)
    {
        return (525 * sintheta) / 4096. + (5103 * quick_pow(costheta, 2) * sintheta) / 4096. +
               (17325 * quick_pow(costheta, 4) * sintheta) / 4096. +
               (63063 * quick_pow(costheta, 6) * sintheta) / 4096. - (1701 * quick_pow(sintheta, 3)) / 4096. -
               (17325 * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / 2048. -
               (315315 * quick_pow(costheta, 4) * quick_pow(sintheta, 3)) / 4096. +
               (3465 * quick_pow(sintheta, 5)) / 4096. +
               (189189 * quick_pow(costheta, 2) * quick_pow(sintheta, 5)) / 4096. -
               (9009 * quick_pow(sintheta, 7)) / 4096.;
    }
    else if (lp == 9 && l == 7 && s == 1 && j == 8)
    {
        return (-735 * sintheta) / 8192. - (3465 * quick_pow(costheta, 2) * sintheta) / 4096. -
               (10725 * quick_pow(costheta, 4) * sintheta) / 4096. -
               (105105 * quick_pow(costheta, 6) * sintheta) / 16384. -
               (328185 * quick_pow(costheta, 8) * sintheta) / 16384. + (1155 * quick_pow(sintheta, 3)) / 4096. +
               (10725 * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / 2048. +
               (525525 * quick_pow(costheta, 4) * quick_pow(sintheta, 3)) / 16384. +
               (765765 * quick_pow(costheta, 6) * quick_pow(sintheta, 3)) / 4096. -
               (2145 * quick_pow(sintheta, 5)) / 4096. -
               (315315 * quick_pow(costheta, 2) * quick_pow(sintheta, 5)) / 16384. -
               (2297295 * quick_pow(costheta, 4) * quick_pow(sintheta, 5)) / 8192. +
               (15015 * quick_pow(sintheta, 7)) / 16384. +
               (328185 * quick_pow(costheta, 2) * quick_pow(sintheta, 7)) / 4096. -
               (36465 * quick_pow(sintheta, 9)) / 16384.;
    }
    else if (lp == 9 && l == 9 && s == 1 && j == 8)
    {
        return (245 * sintheta) / (2048. * sqrt2) + (1155 * quick_pow(costheta, 2) * sintheta) / (1024. * sqrt2) +
               (3575 * quick_pow(costheta, 4) * sintheta) / (1024. * sqrt2) +
               (35035 * quick_pow(costheta, 6) * sintheta) / (4096. * sqrt2) +
               (109395 * quick_pow(costheta, 8) * sintheta) / (4096. * sqrt2) -
               (385 * quick_pow(sintheta, 3)) / (1024. * sqrt2) -
               (3575 * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / (512. * sqrt2) -
               (175175 * quick_pow(costheta, 4) * quick_pow(sintheta, 3)) / (4096. * sqrt2) -
               (255255 * quick_pow(costheta, 6) * quick_pow(sintheta, 3)) / (1024. * sqrt2) +
               (715 * quick_pow(sintheta, 5)) / (1024. * sqrt2) +
               (105105 * quick_pow(costheta, 2) * quick_pow(sintheta, 5)) / (4096. * sqrt2) +
               (765765 * quick_pow(costheta, 4) * quick_pow(sintheta, 5)) / (2048. * sqrt2) -
               (5005 * quick_pow(sintheta, 7)) / (4096. * sqrt2) -
               (109395 * quick_pow(costheta, 2) * quick_pow(sintheta, 7)) / (1024. * sqrt2) +
               (12155 * quick_pow(sintheta, 9)) / (4096. * sqrt2);
    }
    else if (lp == 9 && l == 9 && s == 0 && j == 9)
    {
        return 0;
    }
    else if (lp == 9 && l == 9 && s == 1 && j == 9)
    {
        return (931 * sintheta) / (32768. * sqrt2) + (4389 * quick_pow(costheta, 2) * sintheta) / (16384. * sqrt2) +
               (13585 * quick_pow(costheta, 4) * sintheta) / (16384. * sqrt2) +
               (133133 * quick_pow(costheta, 6) * sintheta) / (65536. * sqrt2) +
               (415701 * quick_pow(costheta, 8) * sintheta) / (65536. * sqrt2) -
               (1463 * quick_pow(sintheta, 3)) / (16384. * sqrt2) -
               (13585 * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / (8192. * sqrt2) -
               (665665 * quick_pow(costheta, 4) * quick_pow(sintheta, 3)) / (65536. * sqrt2) -
               (969969 * quick_pow(costheta, 6) * quick_pow(sintheta, 3)) / (16384. * sqrt2) +
               (2717 * quick_pow(sintheta, 5)) / (16384. * sqrt2) +
               (399399 * quick_pow(costheta, 2) * quick_pow(sintheta, 5)) / (65536. * sqrt2) +
               (2909907 * quick_pow(costheta, 4) * quick_pow(sintheta, 5)) / (32768. * sqrt2) -
               (19019 * quick_pow(sintheta, 7)) / (65536. * sqrt2) -
               (415701 * quick_pow(costheta, 2) * quick_pow(sintheta, 7)) / (16384. * sqrt2) +
               (46189 * quick_pow(sintheta, 9)) / (65536. * sqrt2);
    }
    else if (lp == 8 && l == 8 && s == 1 && j == 9)
    {
        return (-175 * costheta * sintheta) / (256. * sqrt2) -
               (385 * quick_pow(costheta, 3) * sintheta) / (128. * sqrt2) -
               (2145 * quick_pow(costheta, 5) * sintheta) / (256. * sqrt2) -
               (3575 * quick_pow(costheta, 7) * sintheta) / (128. * sqrt2) +
               (385 * costheta * quick_pow(sintheta, 3)) / (128. * sqrt2) +
               (3575 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / (128. * sqrt2) +
               (25025 * quick_pow(costheta, 5) * quick_pow(sintheta, 3)) / (128. * sqrt2) -
               (2145 * costheta * quick_pow(sintheta, 5)) / (256. * sqrt2) -
               (25025 * quick_pow(costheta, 3) * quick_pow(sintheta, 5)) / (128. * sqrt2) +
               (3575 * costheta * quick_pow(sintheta, 7)) / (128. * sqrt2);
    }
    else if (lp == 8 && l == 10 && s == 1 && j == 9)
    {
        return (105 * sqrt5 * costheta * sintheta) / 512. + (231 * sqrt5 * quick_pow(costheta, 3) * sintheta) / 256. +
               (1287 * sqrt5 * quick_pow(costheta, 5) * sintheta) / 512. +
               (2145 * sqrt5 * quick_pow(costheta, 7) * sintheta) / 256. -
               (231 * sqrt5 * costheta * quick_pow(sintheta, 3)) / 256. -
               (2145 * sqrt5 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / 256. -
               (15015 * sqrt5 * quick_pow(costheta, 5) * quick_pow(sintheta, 3)) / 256. +
               (1287 * sqrt5 * costheta * quick_pow(sintheta, 5)) / 512. +
               (15015 * sqrt5 * quick_pow(costheta, 3) * quick_pow(sintheta, 5)) / 256. -
               (2145 * sqrt5 * costheta * quick_pow(sintheta, 7)) / 256.;
    }
    else if (lp == 10 && l == 8 && s == 1 && j == 9)
    {
        return (-4851 * sqrt5 * costheta * sintheta) / 32768. -
               (1287 * sqrt5 * quick_pow(costheta, 3) * sintheta) / 2048. -
               (104247 * sqrt5 * quick_pow(costheta, 5) * sintheta) / 65536. -
               (7293 * sqrt5 * quick_pow(costheta, 7) * sintheta) / 2048. -
               (692835 * sqrt5 * quick_pow(costheta, 9) * sintheta) / 65536. +
               (1287 * sqrt5 * costheta * quick_pow(sintheta, 3)) / 2048. +
               (173745 * sqrt5 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / 32768. +
               (51051 * sqrt5 * quick_pow(costheta, 5) * quick_pow(sintheta, 3)) / 2048. +
               (2078505 * sqrt5 * quick_pow(costheta, 7) * quick_pow(sintheta, 3)) / 16384. -
               (104247 * sqrt5 * costheta * quick_pow(sintheta, 5)) / 65536. -
               (51051 * sqrt5 * quick_pow(costheta, 3) * quick_pow(sintheta, 5)) / 2048. -
               (8729721 * sqrt5 * quick_pow(costheta, 5) * quick_pow(sintheta, 5)) / 32768. +
               (7293 * sqrt5 * costheta * quick_pow(sintheta, 7)) / 2048. +
               (2078505 * sqrt5 * quick_pow(costheta, 3) * quick_pow(sintheta, 7)) / 16384. -
               (692835 * sqrt5 * costheta * quick_pow(sintheta, 9)) / 65536.;
    }
    else if (lp == 10 && l == 10 && s == 1 && j == 9)
    {
        return (14553 * costheta * sintheta) / (32768. * sqrt2) +
               (3861 * quick_pow(costheta, 3) * sintheta) / (2048. * sqrt2) +
               (312741 * quick_pow(costheta, 5) * sintheta) / (65536. * sqrt2) +
               (21879 * quick_pow(costheta, 7) * sintheta) / (2048. * sqrt2) +
               (2078505 * quick_pow(costheta, 9) * sintheta) / (65536. * sqrt2) -
               (3861 * costheta * quick_pow(sintheta, 3)) / (2048. * sqrt2) -
               (521235 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / (32768. * sqrt2) -
               (153153 * quick_pow(costheta, 5) * quick_pow(sintheta, 3)) / (2048. * sqrt2) -
               (6235515 * quick_pow(costheta, 7) * quick_pow(sintheta, 3)) / (16384. * sqrt2) +
               (312741 * costheta * quick_pow(sintheta, 5)) / (65536. * sqrt2) +
               (153153 * quick_pow(costheta, 3) * quick_pow(sintheta, 5)) / (2048. * sqrt2) +
               (26189163 * quick_pow(costheta, 5) * quick_pow(sintheta, 5)) / (32768. * sqrt2) -
               (21879 * costheta * quick_pow(sintheta, 7)) / (2048. * sqrt2) -
               (6235515 * quick_pow(costheta, 3) * quick_pow(sintheta, 7)) / (16384. * sqrt2) +
               (2078505 * costheta * quick_pow(sintheta, 9)) / (65536. * sqrt2);
    }
    else if (lp == 10 && l == 10 && s == 0 && j == 10)
    {
        return 0;
    }
    else if (lp == 10 && l == 10 && s == 1 && j == 10)
    {
        return (3087 * costheta * sintheta) / (32768. * sqrt2) +
               (819 * quick_pow(costheta, 3) * sintheta) / (2048. * sqrt2) +
               (66339 * quick_pow(costheta, 5) * sintheta) / (65536. * sqrt2) +
               (4641 * quick_pow(costheta, 7) * sintheta) / (2048. * sqrt2) +
               (440895 * quick_pow(costheta, 9) * sintheta) / (65536. * sqrt2) -
               (819 * costheta * quick_pow(sintheta, 3)) / (2048. * sqrt2) -
               (110565 * quick_pow(costheta, 3) * quick_pow(sintheta, 3)) / (32768. * sqrt2) -
               (32487 * quick_pow(costheta, 5) * quick_pow(sintheta, 3)) / (2048. * sqrt2) -
               (1322685 * quick_pow(costheta, 7) * quick_pow(sintheta, 3)) / (16384. * sqrt2) +
               (66339 * costheta * quick_pow(sintheta, 5)) / (65536. * sqrt2) +
               (32487 * quick_pow(costheta, 3) * quick_pow(sintheta, 5)) / (2048. * sqrt2) +
               (5555277 * quick_pow(costheta, 5) * quick_pow(sintheta, 5)) / (32768. * sqrt2) -
               (4641 * costheta * quick_pow(sintheta, 7)) / (2048. * sqrt2) -
               (1322685 * quick_pow(costheta, 3) * quick_pow(sintheta, 7)) / (16384. * sqrt2) +
               (440895 * costheta * quick_pow(sintheta, 9)) / (65536. * sqrt2);
    }
    else if (lp == 9 && l == 9 && s == 1 && j == 10)
    {
        return (-4851 * sintheta) / (32768. * sqrt2) - (22869 * quick_pow(costheta, 2) * sintheta) / (16384. * sqrt2) -
               (70785 * quick_pow(costheta, 4) * sintheta) / (16384. * sqrt2) -
               (693693 * quick_pow(costheta, 6) * sintheta) / (65536. * sqrt2) -
               (2166021 * quick_pow(costheta, 8) * sintheta) / (65536. * sqrt2) +
               (7623 * quick_pow(sintheta, 3)) / (16384. * sqrt2) +
               (70785 * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / (8192. * sqrt2) +
               (3468465 * quick_pow(costheta, 4) * quick_pow(sintheta, 3)) / (65536. * sqrt2) +
               (5054049 * quick_pow(costheta, 6) * quick_pow(sintheta, 3)) / (16384. * sqrt2) -
               (14157 * quick_pow(sintheta, 5)) / (16384. * sqrt2) -
               (2081079 * quick_pow(costheta, 2) * quick_pow(sintheta, 5)) / (65536. * sqrt2) -
               (15162147 * quick_pow(costheta, 4) * quick_pow(sintheta, 5)) / (32768. * sqrt2) +
               (99099 * quick_pow(sintheta, 7)) / (65536. * sqrt2) +
               (2166021 * quick_pow(costheta, 2) * quick_pow(sintheta, 7)) / (16384. * sqrt2) -
               (240669 * quick_pow(sintheta, 9)) / (65536. * sqrt2);
    }
    else if (lp == 9 && l == 11 && s == 1 && j == 10)
    {
        return (441 * sqrt(55) * sintheta) / 32768. + (2079 * sqrt(55) * quick_pow(costheta, 2) * sintheta) / 16384. +
               (6435 * sqrt(55) * quick_pow(costheta, 4) * sintheta) / 16384. +
               (63063 * sqrt(55) * quick_pow(costheta, 6) * sintheta) / 65536. +
               (196911 * sqrt(55) * quick_pow(costheta, 8) * sintheta) / 65536. -
               (693 * sqrt(55) * quick_pow(sintheta, 3)) / 16384. -
               (6435 * sqrt(55) * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / 8192. -
               (315315 * sqrt(55) * quick_pow(costheta, 4) * quick_pow(sintheta, 3)) / 65536. -
               (459459 * sqrt(55) * quick_pow(costheta, 6) * quick_pow(sintheta, 3)) / 16384. +
               (1287 * sqrt(55) * quick_pow(sintheta, 5)) / 16384. +
               (189189 * sqrt(55) * quick_pow(costheta, 2) * quick_pow(sintheta, 5)) / 65536. +
               (1378377 * sqrt(55) * quick_pow(costheta, 4) * quick_pow(sintheta, 5)) / 32768. -
               (9009 * sqrt(55) * quick_pow(sintheta, 7)) / 65536. -
               (196911 * sqrt(55) * quick_pow(costheta, 2) * quick_pow(sintheta, 7)) / 16384. +
               (21879 * sqrt(55) * quick_pow(sintheta, 9)) / 65536.;
    }
    else if (lp == 11 && l == 9 && s == 1 && j == 10)
    {
        return (-1323 * sqrt(55) * sintheta) / 131072. -
               (12285 * sqrt(55) * quick_pow(costheta, 2) * sintheta) / 131072. -
               (73125 * sqrt(55) * quick_pow(costheta, 4) * sintheta) / 262144. -
               (162435 * sqrt(55) * quick_pow(costheta, 6) * sintheta) / 262144. -
               (340119 * sqrt(55) * quick_pow(costheta, 8) * sintheta) / 262144. -
               (969969 * sqrt(55) * quick_pow(costheta, 10) * sintheta) / 262144. +
               (4095 * sqrt(55) * quick_pow(sintheta, 3)) / 131072. +
               (73125 * sqrt(55) * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / 131072. +
               (812175 * sqrt(55) * quick_pow(costheta, 4) * quick_pow(sintheta, 3)) / 262144. +
               (793611 * sqrt(55) * quick_pow(costheta, 6) * quick_pow(sintheta, 3)) / 65536. +
               (14549535 * sqrt(55) * quick_pow(costheta, 8) * quick_pow(sintheta, 3)) / 262144. -
               (14625 * sqrt(55) * quick_pow(sintheta, 5)) / 262144. -
               (487305 * sqrt(55) * quick_pow(costheta, 2) * quick_pow(sintheta, 5)) / 262144. -
               (2380833 * sqrt(55) * quick_pow(costheta, 4) * quick_pow(sintheta, 5)) / 131072. -
               (20369349 * sqrt(55) * quick_pow(costheta, 6) * quick_pow(sintheta, 5)) / 131072. +
               (23205 * sqrt(55) * quick_pow(sintheta, 7)) / 262144. +
               (340119 * sqrt(55) * quick_pow(costheta, 2) * quick_pow(sintheta, 7)) / 65536. +
               (14549535 * sqrt(55) * quick_pow(costheta, 4) * quick_pow(sintheta, 7)) / 131072. -
               (37791 * sqrt(55) * quick_pow(sintheta, 9)) / 262144. -
               (4849845 * sqrt(55) * quick_pow(costheta, 2) * quick_pow(sintheta, 9)) / 262144. +
               (88179 * sqrt(55) * quick_pow(sintheta, 11)) / 262144.;
    }
    else if (lp == 11 && l == 11 && s == 1 && j == 10)
    {
        return (6615 * sintheta) / (65536. * sqrt2) + (61425 * quick_pow(costheta, 2) * sintheta) / (65536. * sqrt2) +
               (365625 * quick_pow(costheta, 4) * sintheta) / (131072. * sqrt2) +
               (812175 * quick_pow(costheta, 6) * sintheta) / (131072. * sqrt2) +
               (1700595 * quick_pow(costheta, 8) * sintheta) / (131072. * sqrt2) +
               (4849845 * quick_pow(costheta, 10) * sintheta) / (131072. * sqrt2) -
               (20475 * quick_pow(sintheta, 3)) / (65536. * sqrt2) -
               (365625 * quick_pow(costheta, 2) * quick_pow(sintheta, 3)) / (65536. * sqrt2) -
               (4060875 * quick_pow(costheta, 4) * quick_pow(sintheta, 3)) / (131072. * sqrt2) -
               (3968055 * quick_pow(costheta, 6) * quick_pow(sintheta, 3)) / (32768. * sqrt2) -
               (72747675 * quick_pow(costheta, 8) * quick_pow(sintheta, 3)) / (131072. * sqrt2) +
               (73125 * quick_pow(sintheta, 5)) / (131072. * sqrt2) +
               (2436525 * quick_pow(costheta, 2) * quick_pow(sintheta, 5)) / (131072. * sqrt2) +
               (11904165 * quick_pow(costheta, 4) * quick_pow(sintheta, 5)) / (65536. * sqrt2) +
               (101846745 * quick_pow(costheta, 6) * quick_pow(sintheta, 5)) / (65536. * sqrt2) -
               (116025 * quick_pow(sintheta, 7)) / (131072. * sqrt2) -
               (1700595 * quick_pow(costheta, 2) * quick_pow(sintheta, 7)) / (32768. * sqrt2) -
               (72747675 * quick_pow(costheta, 4) * quick_pow(sintheta, 7)) / (65536. * sqrt2) +
               (188955 * quick_pow(sintheta, 9)) / (131072. * sqrt2) +
               (24249225 * quick_pow(costheta, 2) * quick_pow(sintheta, 9)) / (131072. * sqrt2) -
               (440895 * quick_pow(sintheta, 11)) / (131072. * sqrt2);
    }

    return 0;
}

double coeff_m00(double t, int lp, int l, int s, int j)
{
    double costheta = t;
    double sintheta = sqrt(1.0 - t * t);
    if (lp == 0 && l == 0 && s == 0 && j == 0)
    {
        return 0;
    }
    else if (lp == 1 && l == 1 && s == 1 && j == 0)
    {
        return costheta;
    }
    else if (lp == 1 && l == 1 && s == 0 && j == 1)
    {
        return 0;
    }
    else if (lp == 1 && l == 1 && s == 1 && j == 1)
    {
        return 0;
    }
    else if (lp == 0 && l == 0 && s == 1 && j == 1)
    {
        return 1;
    }
    else if (lp == 0 && l == 2 && s == 1 && j == 1)
    {
        return sqrt2;
    }
    else if (lp == 2 && l == 0 && s == 1 && j == 1)
    {
        return -(1 / sqrt2) + (3 * quick_pow(costheta, 2)) / sqrt2;
    }
    else if (lp == 2 && l == 2 && s == 1 && j == 1)
    {
        return -1 + 3 * quick_pow(costheta, 2);
    }
    else if (lp == 2 && l == 2 && s == 0 && j == 2)
    {
        return 0;
    }
    else if (lp == 2 && l == 2 && s == 1 && j == 2)
    {
        return 0;
    }
    else if (lp == 1 && l == 1 && s == 1 && j == 2)
    {
        return 2 * costheta;
    }
    else if (lp == 1 && l == 3 && s == 1 && j == 2)
    {
        return sqrt(6) * costheta;
    }
    else if (lp == 3 && l == 1 && s == 1 && j == 2)
    {
        return (3 * sqrt(1.5) * costheta) / 4. + (5 * sqrt(1.5) * quick_pow(costheta, 3)) / 4. -
               (15 * sqrt(1.5) * costheta * quick_pow(sintheta, 2)) / 4.;
    }
    else if (lp == 3 && l == 3 && s == 1 && j == 2)
    {
        return (9 * costheta) / 8. + (15 * quick_pow(costheta, 3)) / 8. - (45 * costheta * quick_pow(sintheta, 2)) / 8.;
    }
    else if (lp == 3 && l == 3 && s == 0 && j == 3)
    {
        return 0;
    }
    else if (lp == 3 && l == 3 && s == 1 && j == 3)
    {
        return 0;
    }
    else if (lp == 2 && l == 2 && s == 1 && j == 3)
    {
        return -1.5 + (9 * quick_pow(costheta, 2)) / 2.;
    }
    else if (lp == 2 && l == 4 && s == 1 && j == 3)
    {
        return -sqrt3 + 3 * sqrt3 * quick_pow(costheta, 2);
    }
    else if (lp == 4 && l == 2 && s == 1 && j == 3)
    {
        return (9 * sqrt3) / 32. + (5 * sqrt3 * quick_pow(costheta, 2)) / 8. +
               (35 * sqrt3 * quick_pow(costheta, 4)) / 32. - (5 * sqrt3 * quick_pow(sintheta, 2)) / 8. -
               (105 * sqrt3 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 16. +
               (35 * sqrt3 * quick_pow(sintheta, 4)) / 32.;
    }
    else if (lp == 4 && l == 4 && s == 1 && j == 3)
    {
        return 0.5625 + (5 * quick_pow(costheta, 2)) / 4. + (35 * quick_pow(costheta, 4)) / 16. -
               (5 * quick_pow(sintheta, 2)) / 4. - (105 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 8. +
               (35 * quick_pow(sintheta, 4)) / 16.;
    }
    else if (lp == 4 && l == 4 && s == 0 && j == 4)
    {
        return 0;
    }
    else if (lp == 4 && l == 4 && s == 1 && j == 4)
    {
        return 0;
    }
    else if (lp == 3 && l == 3 && s == 1 && j == 4)
    {
        return (3 * costheta) / 2. + (5 * quick_pow(costheta, 3)) / 2. - (15 * costheta * quick_pow(sintheta, 2)) / 2.;
    }
    else if (lp == 3 && l == 5 && s == 1 && j == 4)
    {
        return (3 * sqrt5 * costheta) / 4. + (5 * sqrt5 * quick_pow(costheta, 3)) / 4. -
               (15 * sqrt5 * costheta * quick_pow(sintheta, 2)) / 4.;
    }
    else if (lp == 5 && l == 3 && s == 1 && j == 4)
    {
        return (15 * sqrt5 * costheta) / 32. + (35 * sqrt5 * quick_pow(costheta, 3)) / 64. +
               (63 * sqrt5 * quick_pow(costheta, 5)) / 64. - (105 * sqrt5 * costheta * quick_pow(sintheta, 2)) / 64. -
               (315 * sqrt5 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 32. +
               (315 * sqrt5 * costheta * quick_pow(sintheta, 4)) / 64.;
    }
    else if (lp == 5 && l == 5 && s == 1 && j == 4)
    {
        return (75 * costheta) / 64. + (175 * quick_pow(costheta, 3)) / 128. + (315 * quick_pow(costheta, 5)) / 128. -
               (525 * costheta * quick_pow(sintheta, 2)) / 128. -
               (1575 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 64. +
               (1575 * costheta * quick_pow(sintheta, 4)) / 128.;
    }
    else if (lp == 5 && l == 5 && s == 0 && j == 5)
    {
        return 0;
    }
    else if (lp == 5 && l == 5 && s == 1 && j == 5)
    {
        return 0;
    }
    else if (lp == 4 && l == 4 && s == 1 && j == 5)
    {
        return 0.703125 + (25 * quick_pow(costheta, 2)) / 16. + (175 * quick_pow(costheta, 4)) / 64. -
               (25 * quick_pow(sintheta, 2)) / 16. - (525 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 32. +
               (175 * quick_pow(sintheta, 4)) / 64.;
    }
    else if (lp == 4 && l == 6 && s == 1 && j == 5)
    {
        return (9 * sqrt(7.5)) / 32. + (5 * sqrt(7.5) * quick_pow(costheta, 2)) / 8. +
               (35 * sqrt(7.5) * quick_pow(costheta, 4)) / 32. - (5 * sqrt(7.5) * quick_pow(sintheta, 2)) / 8. -
               (105 * sqrt(7.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 16. +
               (35 * sqrt(7.5) * quick_pow(sintheta, 4)) / 32.;
    }
    else if (lp == 6 && l == 4 && s == 1 && j == 5)
    {
        return (25 * sqrt(7.5)) / 128. + (105 * sqrt(7.5) * quick_pow(costheta, 2)) / 256. +
               (63 * sqrt(7.5) * quick_pow(costheta, 4)) / 128. + (231 * sqrt(7.5) * quick_pow(costheta, 6)) / 256. -
               (105 * sqrt(7.5) * quick_pow(sintheta, 2)) / 256. -
               (189 * sqrt(7.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 64. -
               (3465 * sqrt(7.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 256. +
               (63 * sqrt(7.5) * quick_pow(sintheta, 4)) / 128. +
               (3465 * sqrt(7.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 256. -
               (231 * sqrt(7.5) * quick_pow(sintheta, 6)) / 256.;
    }
    else if (lp == 6 && l == 6 && s == 1 && j == 5)
    {
        return 0.5859375 + (315 * quick_pow(costheta, 2)) / 256. + (189 * quick_pow(costheta, 4)) / 128. +
               (693 * quick_pow(costheta, 6)) / 256. - (315 * quick_pow(sintheta, 2)) / 256. -
               (567 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 64. -
               (10395 * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 256. +
               (189 * quick_pow(sintheta, 4)) / 128. +
               (10395 * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 256. - (693 * quick_pow(sintheta, 6)) / 256.;
    }
    else if (lp == 6 && l == 6 && s == 0 && j == 6)
    {
        return 0;
    }
    else if (lp == 6 && l == 6 && s == 1 && j == 6)
    {
        return 0;
    }
    else if (lp == 5 && l == 5 && s == 1 && j == 6)
    {
        return (45 * costheta) / 32. + (105 * quick_pow(costheta, 3)) / 64. + (189 * quick_pow(costheta, 5)) / 64. -
               (315 * costheta * quick_pow(sintheta, 2)) / 64. -
               (945 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 32. +
               (945 * costheta * quick_pow(sintheta, 4)) / 64.;
    }
    else if (lp == 5 && l == 7 && s == 1 && j == 6)
    {
        return (15 * sqrt(10.5) * costheta) / 32. + (35 * sqrt(10.5) * quick_pow(costheta, 3)) / 64. +
               (63 * sqrt(10.5) * quick_pow(costheta, 5)) / 64. -
               (105 * sqrt(10.5) * costheta * quick_pow(sintheta, 2)) / 64. -
               (315 * sqrt(10.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 32. +
               (315 * sqrt(10.5) * costheta * quick_pow(sintheta, 4)) / 64.;
    }
    else if (lp == 7 && l == 5 && s == 1 && j == 6)
    {
        return (175 * sqrt(10.5) * costheta) / 512. + (189 * sqrt(10.5) * quick_pow(costheta, 3)) / 512. +
               (231 * sqrt(10.5) * quick_pow(costheta, 5)) / 512. + (429 * sqrt(10.5) * quick_pow(costheta, 7)) / 512. -
               (567 * sqrt(10.5) * costheta * quick_pow(sintheta, 2)) / 512. -
               (1155 * sqrt(10.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 256. -
               (9009 * sqrt(10.5) * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 512. +
               (1155 * sqrt(10.5) * costheta * quick_pow(sintheta, 4)) / 512. +
               (15015 * sqrt(10.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 512. -
               (3003 * sqrt(10.5) * costheta * quick_pow(sintheta, 6)) / 512.;
    }
    else if (lp == 7 && l == 7 && s == 1 && j == 6)
    {
        return (1225 * costheta) / 1024. + (1323 * quick_pow(costheta, 3)) / 1024. +
               (1617 * quick_pow(costheta, 5)) / 1024. + (3003 * quick_pow(costheta, 7)) / 1024. -
               (3969 * costheta * quick_pow(sintheta, 2)) / 1024. -
               (8085 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 512. -
               (63063 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 1024. +
               (8085 * costheta * quick_pow(sintheta, 4)) / 1024. +
               (105105 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 1024. -
               (21021 * costheta * quick_pow(sintheta, 6)) / 1024.;
    }
    else if (lp == 7 && l == 7 && s == 0 && j == 7)
    {
        return 0;
    }
    else if (lp == 7 && l == 7 && s == 1 && j == 7)
    {
        return 0;
    }
    else if (lp == 6 && l == 6 && s == 1 && j == 7)
    {
        return 0.68359375 + (735 * quick_pow(costheta, 2)) / 512. + (441 * quick_pow(costheta, 4)) / 256. +
               (1617 * quick_pow(costheta, 6)) / 512. - (735 * quick_pow(sintheta, 2)) / 512. -
               (1323 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 128. -
               (24255 * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 512. +
               (441 * quick_pow(sintheta, 4)) / 256. +
               (24255 * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 512. -
               (1617 * quick_pow(sintheta, 6)) / 512.;
    }
    else if (lp == 6 && l == 8 && s == 1 && j == 7)
    {
        return (25 * sqrt(3.5)) / 64. + (105 * sqrt(3.5) * quick_pow(costheta, 2)) / 128. +
               (63 * sqrt(3.5) * quick_pow(costheta, 4)) / 64. + (231 * sqrt(3.5) * quick_pow(costheta, 6)) / 128. -
               (105 * sqrt(3.5) * quick_pow(sintheta, 2)) / 128. -
               (189 * sqrt(3.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 32. -
               (3465 * sqrt(3.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 128. +
               (63 * sqrt(3.5) * quick_pow(sintheta, 4)) / 64. +
               (3465 * sqrt(3.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 128. -
               (231 * sqrt(3.5) * quick_pow(sintheta, 6)) / 128.;
    }
    else if (lp == 8 && l == 6 && s == 1 && j == 7)
    {
        return (1225 * sqrt(3.5)) / 4096. + (315 * sqrt(3.5) * quick_pow(costheta, 2)) / 512. +
               (693 * sqrt(3.5) * quick_pow(costheta, 4)) / 1024. + (429 * sqrt(3.5) * quick_pow(costheta, 6)) / 512. +
               (6435 * sqrt(3.5) * quick_pow(costheta, 8)) / 4096. - (315 * sqrt(3.5) * quick_pow(sintheta, 2)) / 512. -
               (2079 * sqrt(3.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 512. -
               (6435 * sqrt(3.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 512. -
               (45045 * sqrt(3.5) * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 1024. +
               (693 * sqrt(3.5) * quick_pow(sintheta, 4)) / 1024. +
               (6435 * sqrt(3.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 512. +
               (225225 * sqrt(3.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 2048. -
               (429 * sqrt(3.5) * quick_pow(sintheta, 6)) / 512. -
               (45045 * sqrt(3.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 1024. +
               (6435 * sqrt(3.5) * quick_pow(sintheta, 8)) / 4096.;
    }
    else if (lp == 8 && l == 8 && s == 1 && j == 7)
    {
        return 0.59814453125 + (315 * quick_pow(costheta, 2)) / 256. + (693 * quick_pow(costheta, 4)) / 512. +
               (429 * quick_pow(costheta, 6)) / 256. + (6435 * quick_pow(costheta, 8)) / 2048. -
               (315 * quick_pow(sintheta, 2)) / 256. - (2079 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 256. -
               (6435 * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 256. -
               (45045 * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 512. +
               (693 * quick_pow(sintheta, 4)) / 512. + (6435 * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 256. +
               (225225 * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 1024. -
               (429 * quick_pow(sintheta, 6)) / 256. -
               (45045 * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 512. +
               (6435 * quick_pow(sintheta, 8)) / 2048.;
    }
    else if (lp == 8 && l == 8 && s == 0 && j == 8)
    {
        return 0;
    }
    else if (lp == 8 && l == 8 && s == 1 && j == 8)
    {
        return 0;
    }
    else if (lp == 7 && l == 7 && s == 1 && j == 8)
    {
        return (175 * costheta) / 128. + (189 * quick_pow(costheta, 3)) / 128. + (231 * quick_pow(costheta, 5)) / 128. +
               (429 * quick_pow(costheta, 7)) / 128. - (567 * costheta * quick_pow(sintheta, 2)) / 128. -
               (1155 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 64. -
               (9009 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 128. +
               (1155 * costheta * quick_pow(sintheta, 4)) / 128. +
               (15015 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 128. -
               (3003 * costheta * quick_pow(sintheta, 6)) / 128.;
    }
    else if (lp == 7 && l == 9 && s == 1 && j == 8)
    {
        return (525 * costheta) / (256. * sqrt2) + (567 * quick_pow(costheta, 3)) / (256. * sqrt2) +
               (693 * quick_pow(costheta, 5)) / (256. * sqrt2) + (1287 * quick_pow(costheta, 7)) / (256. * sqrt2) -
               (1701 * costheta * quick_pow(sintheta, 2)) / (256. * sqrt2) -
               (3465 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / (128. * sqrt2) -
               (27027 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / (256. * sqrt2) +
               (3465 * costheta * quick_pow(sintheta, 4)) / (256. * sqrt2) +
               (45045 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / (256. * sqrt2) -
               (9009 * costheta * quick_pow(sintheta, 6)) / (256. * sqrt2);
    }
    else if (lp == 9 && l == 7 && s == 1 && j == 8)
    {
        return (6615 * costheta) / (4096. * sqrt2) + (3465 * quick_pow(costheta, 3)) / (2048. * sqrt2) +
               (3861 * quick_pow(costheta, 5)) / (2048. * sqrt2) + (19305 * quick_pow(costheta, 7)) / (8192. * sqrt2) +
               (36465 * quick_pow(costheta, 9)) / (8192. * sqrt2) -
               (10395 * costheta * quick_pow(sintheta, 2)) / (2048. * sqrt2) -
               (19305 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / (1024. * sqrt2) -
               (405405 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / (8192. * sqrt2) -
               (328185 * quick_pow(costheta, 7) * quick_pow(sintheta, 2)) / (2048. * sqrt2) +
               (19305 * costheta * quick_pow(sintheta, 4)) / (2048. * sqrt2) +
               (675675 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / (8192. * sqrt2) +
               (2297295 * quick_pow(costheta, 5) * quick_pow(sintheta, 4)) / (4096. * sqrt2) -
               (135135 * costheta * quick_pow(sintheta, 6)) / (8192. * sqrt2) -
               (765765 * quick_pow(costheta, 3) * quick_pow(sintheta, 6)) / (2048. * sqrt2) +
               (328185 * costheta * quick_pow(sintheta, 8)) / (8192. * sqrt2);
    }
    else if (lp == 9 && l == 9 && s == 1 && j == 8)
    {
        return (19845 * costheta) / 16384. + (10395 * quick_pow(costheta, 3)) / 8192. +
               (11583 * quick_pow(costheta, 5)) / 8192. + (57915 * quick_pow(costheta, 7)) / 32768. +
               (109395 * quick_pow(costheta, 9)) / 32768. - (31185 * costheta * quick_pow(sintheta, 2)) / 8192. -
               (57915 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 4096. -
               (1216215 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 32768. -
               (984555 * quick_pow(costheta, 7) * quick_pow(sintheta, 2)) / 8192. +
               (57915 * costheta * quick_pow(sintheta, 4)) / 8192. +
               (2027025 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 32768. +
               (6891885 * quick_pow(costheta, 5) * quick_pow(sintheta, 4)) / 16384. -
               (405405 * costheta * quick_pow(sintheta, 6)) / 32768. -
               (2297295 * quick_pow(costheta, 3) * quick_pow(sintheta, 6)) / 8192. +
               (984555 * costheta * quick_pow(sintheta, 8)) / 32768.;
    }
    else if (lp == 9 && l == 9 && s == 0 && j == 9)
    {
        return 0;
    }
    else if (lp == 9 && l == 9 && s == 1 && j == 9)
    {
        return 0;
    }
    else if (lp == 8 && l == 8 && s == 1 && j == 9)
    {
        return 0.67291259765625 + (2835 * quick_pow(costheta, 2)) / 2048. + (6237 * quick_pow(costheta, 4)) / 4096. +
               (3861 * quick_pow(costheta, 6)) / 2048. + (57915 * quick_pow(costheta, 8)) / 16384. -
               (2835 * quick_pow(sintheta, 2)) / 2048. -
               (18711 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 2048. -
               (57915 * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 2048. -
               (405405 * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 4096. +
               (6237 * quick_pow(sintheta, 4)) / 4096. +
               (57915 * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 2048. +
               (2027025 * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 8192. -
               (3861 * quick_pow(sintheta, 6)) / 2048. -
               (405405 * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 4096. +
               (57915 * quick_pow(sintheta, 8)) / 16384.;
    }
    else if (lp == 8 && l == 10 && s == 1 && j == 9)
    {
        return (3675 * sqrt(2.5)) / 8192. + (945 * sqrt(2.5) * quick_pow(costheta, 2)) / 1024. +
               (2079 * sqrt(2.5) * quick_pow(costheta, 4)) / 2048. +
               (1287 * sqrt(2.5) * quick_pow(costheta, 6)) / 1024. +
               (19305 * sqrt(2.5) * quick_pow(costheta, 8)) / 8192. -
               (945 * sqrt(2.5) * quick_pow(sintheta, 2)) / 1024. -
               (6237 * sqrt(2.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 1024. -
               (19305 * sqrt(2.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 1024. -
               (135135 * sqrt(2.5) * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 2048. +
               (2079 * sqrt(2.5) * quick_pow(sintheta, 4)) / 2048. +
               (19305 * sqrt(2.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 1024. +
               (675675 * sqrt(2.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 4096. -
               (1287 * sqrt(2.5) * quick_pow(sintheta, 6)) / 1024. -
               (135135 * sqrt(2.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 2048. +
               (19305 * sqrt(2.5) * quick_pow(sintheta, 8)) / 8192.;
    }
    else if (lp == 10 && l == 8 && s == 1 && j == 9)
    {
        return (11907 * sqrt(2.5)) / 32768. + (24255 * sqrt(2.5) * quick_pow(costheta, 2)) / 32768. +
               (6435 * sqrt(2.5) * quick_pow(costheta, 4)) / 8192. +
               (57915 * sqrt(2.5) * quick_pow(costheta, 6)) / 65536. +
               (36465 * sqrt(2.5) * quick_pow(costheta, 8)) / 32768. +
               (138567 * sqrt(2.5) * quick_pow(costheta, 10)) / 65536. -
               (24255 * sqrt(2.5) * quick_pow(sintheta, 2)) / 32768. -
               (19305 * sqrt(2.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 4096. -
               (868725 * sqrt(2.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 65536. -
               (255255 * sqrt(2.5) * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 8192. -
               (6235515 * sqrt(2.5) * quick_pow(costheta, 8) * quick_pow(sintheta, 2)) / 65536. +
               (6435 * sqrt(2.5) * quick_pow(sintheta, 4)) / 8192. +
               (868725 * sqrt(2.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 65536. +
               (1276275 * sqrt(2.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 16384. +
               (14549535 * sqrt(2.5) * quick_pow(costheta, 6) * quick_pow(sintheta, 4)) / 32768. -
               (57915 * sqrt(2.5) * quick_pow(sintheta, 6)) / 65536. -
               (255255 * sqrt(2.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 8192. -
               (14549535 * sqrt(2.5) * quick_pow(costheta, 4) * quick_pow(sintheta, 6)) / 32768. +
               (36465 * sqrt(2.5) * quick_pow(sintheta, 8)) / 32768. +
               (6235515 * sqrt(2.5) * quick_pow(costheta, 2) * quick_pow(sintheta, 8)) / 65536. -
               (138567 * sqrt(2.5) * quick_pow(sintheta, 10)) / 65536.;
    }
    else if (lp == 10 && l == 10 && s == 1 && j == 9)
    {
        return 0.605621337890625 + (40425 * quick_pow(costheta, 2)) / 32768. +
               (10725 * quick_pow(costheta, 4)) / 8192. + (96525 * quick_pow(costheta, 6)) / 65536. +
               (60775 * quick_pow(costheta, 8)) / 32768. + (230945 * quick_pow(costheta, 10)) / 65536. -
               (40425 * quick_pow(sintheta, 2)) / 32768. -
               (32175 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 4096. -
               (1447875 * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 65536. -
               (425425 * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 8192. -
               (10392525 * quick_pow(costheta, 8) * quick_pow(sintheta, 2)) / 65536. +
               (10725 * quick_pow(sintheta, 4)) / 8192. +
               (1447875 * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 65536. +
               (2127125 * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 16384. +
               (24249225 * quick_pow(costheta, 6) * quick_pow(sintheta, 4)) / 32768. -
               (96525 * quick_pow(sintheta, 6)) / 65536. -
               (425425 * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 8192. -
               (24249225 * quick_pow(costheta, 4) * quick_pow(sintheta, 6)) / 32768. +
               (60775 * quick_pow(sintheta, 8)) / 32768. +
               (10392525 * quick_pow(costheta, 2) * quick_pow(sintheta, 8)) / 65536. -
               (230945 * quick_pow(sintheta, 10)) / 65536.;
    }
    else if (lp == 10 && l == 10 && s == 0 && j == 10)
    {
        return 0;
    }
    else if (lp == 10 && l == 10 && s == 1 && j == 10)
    {
        return 0;
    }
    else if (lp == 9 && l == 9 && s == 1 && j == 10)
    {
        return (11025 * costheta) / 8192. + (5775 * quick_pow(costheta, 3)) / 4096. +
               (6435 * quick_pow(costheta, 5)) / 4096. + (32175 * quick_pow(costheta, 7)) / 16384. +
               (60775 * quick_pow(costheta, 9)) / 16384. - (17325 * costheta * quick_pow(sintheta, 2)) / 4096. -
               (32175 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 2048. -
               (675675 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 16384. -
               (546975 * quick_pow(costheta, 7) * quick_pow(sintheta, 2)) / 4096. +
               (32175 * costheta * quick_pow(sintheta, 4)) / 4096. +
               (1126125 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 16384. +
               (3828825 * quick_pow(costheta, 5) * quick_pow(sintheta, 4)) / 8192. -
               (225225 * costheta * quick_pow(sintheta, 6)) / 16384. -
               (1276275 * quick_pow(costheta, 3) * quick_pow(sintheta, 6)) / 4096. +
               (546975 * costheta * quick_pow(sintheta, 8)) / 16384.;
    }
    else if (lp == 9 && l == 11 && s == 1 && j == 10)
    {
        return (2205 * sqrt(27.5) * costheta) / 8192. + (1155 * sqrt(27.5) * quick_pow(costheta, 3)) / 4096. +
               (1287 * sqrt(27.5) * quick_pow(costheta, 5)) / 4096. +
               (6435 * sqrt(27.5) * quick_pow(costheta, 7)) / 16384. +
               (12155 * sqrt(27.5) * quick_pow(costheta, 9)) / 16384. -
               (3465 * sqrt(27.5) * costheta * quick_pow(sintheta, 2)) / 4096. -
               (6435 * sqrt(27.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 2048. -
               (135135 * sqrt(27.5) * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 16384. -
               (109395 * sqrt(27.5) * quick_pow(costheta, 7) * quick_pow(sintheta, 2)) / 4096. +
               (6435 * sqrt(27.5) * costheta * quick_pow(sintheta, 4)) / 4096. +
               (225225 * sqrt(27.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 16384. +
               (765765 * sqrt(27.5) * quick_pow(costheta, 5) * quick_pow(sintheta, 4)) / 8192. -
               (45045 * sqrt(27.5) * costheta * quick_pow(sintheta, 6)) / 16384. -
               (255255 * sqrt(27.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 6)) / 4096. +
               (109395 * sqrt(27.5) * costheta * quick_pow(sintheta, 8)) / 16384.;
    }
    else if (lp == 11 && l == 9 && s == 1 && j == 10)
    {
        return (14553 * sqrt(27.5) * costheta) / 65536. + (15015 * sqrt(27.5) * quick_pow(costheta, 3)) / 65536. +
               (32175 * sqrt(27.5) * quick_pow(costheta, 5)) / 131072. +
               (36465 * sqrt(27.5) * quick_pow(costheta, 7)) / 131072. +
               (46189 * sqrt(27.5) * quick_pow(costheta, 9)) / 131072. +
               (88179 * sqrt(27.5) * quick_pow(costheta, 11)) / 131072. -
               (45045 * sqrt(27.5) * costheta * quick_pow(sintheta, 2)) / 65536. -
               (160875 * sqrt(27.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 65536. -
               (765765 * sqrt(27.5) * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 131072. -
               (415701 * sqrt(27.5) * quick_pow(costheta, 7) * quick_pow(sintheta, 2)) / 32768. -
               (4849845 * sqrt(27.5) * quick_pow(costheta, 9) * quick_pow(sintheta, 2)) / 131072. +
               (160875 * sqrt(27.5) * costheta * quick_pow(sintheta, 4)) / 131072. +
               (1276275 * sqrt(27.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 131072. +
               (2909907 * sqrt(27.5) * quick_pow(costheta, 5) * quick_pow(sintheta, 4)) / 65536. +
               (14549535 * sqrt(27.5) * quick_pow(costheta, 7) * quick_pow(sintheta, 4)) / 65536. -
               (255255 * sqrt(27.5) * costheta * quick_pow(sintheta, 6)) / 131072. -
               (969969 * sqrt(27.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 6)) / 32768. -
               (20369349 * sqrt(27.5) * quick_pow(costheta, 5) * quick_pow(sintheta, 6)) / 65536. +
               (415701 * sqrt(27.5) * costheta * quick_pow(sintheta, 8)) / 131072. +
               (14549535 * sqrt(27.5) * quick_pow(costheta, 3) * quick_pow(sintheta, 8)) / 131072. -
               (969969 * sqrt(27.5) * costheta * quick_pow(sintheta, 10)) / 131072.;
    }
    else if (lp == 11 && l == 11 && s == 1 && j == 10)
    {
        return (160083 * costheta) / 131072. + (165165 * quick_pow(costheta, 3)) / 131072. +
               (353925 * quick_pow(costheta, 5)) / 262144. + (401115 * quick_pow(costheta, 7)) / 262144. +
               (508079 * quick_pow(costheta, 9)) / 262144. + (969969 * quick_pow(costheta, 11)) / 262144. -
               (495495 * costheta * quick_pow(sintheta, 2)) / 131072. -
               (1769625 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 131072. -
               (8423415 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 262144. -
               (4572711 * quick_pow(costheta, 7) * quick_pow(sintheta, 2)) / 65536. -
               (53348295 * quick_pow(costheta, 9) * quick_pow(sintheta, 2)) / 262144. +
               (1769625 * costheta * quick_pow(sintheta, 4)) / 262144. +
               (14039025 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 262144. +
               (32008977 * quick_pow(costheta, 5) * quick_pow(sintheta, 4)) / 131072. +
               (160044885 * quick_pow(costheta, 7) * quick_pow(sintheta, 4)) / 131072. -
               (2807805 * costheta * quick_pow(sintheta, 6)) / 262144. -
               (10669659 * quick_pow(costheta, 3) * quick_pow(sintheta, 6)) / 65536. -
               (224062839 * quick_pow(costheta, 5) * quick_pow(sintheta, 6)) / 131072. +
               (4572711 * costheta * quick_pow(sintheta, 8)) / 262144. +
               (160044885 * quick_pow(costheta, 3) * quick_pow(sintheta, 8)) / 262144. -
               (10669659 * costheta * quick_pow(sintheta, 10)) / 262144.;
    }

    return 0;
}

double coeff_mss(double t, int lp, int l, int s, int j)
{
    double costheta = t;
    double sintheta = sqrt(1.0 - t * t);
    if (lp == 0 && l == 0 && s == 0 && j == 0)
    {
        return 1;
    }
    else if (lp == 1 && l == 1 && s == 1 && j == 0)
    {
        return 0;
    }
    else if (lp == 1 && l == 1 && s == 0 && j == 1)
    {
        return 3 * costheta;
    }
    else if (lp == 2 && l == 2 && s == 0 && j == 2)
    {
        return -2.5 + (15 * quick_pow(costheta, 2)) / 2.;
    }
    else if (lp == 3 && l == 3 && s == 0 && j == 3)
    {
        return (21 * costheta) / 8. + (35 * quick_pow(costheta, 3)) / 8. -
               (105 * costheta * quick_pow(sintheta, 2)) / 8.;
    }
    else if (lp == 4 && l == 4 && s == 0 && j == 4)
    {
        return 1.265625 + (45 * quick_pow(costheta, 2)) / 16. + (315 * quick_pow(costheta, 4)) / 64. -
               (45 * quick_pow(sintheta, 2)) / 16. - (945 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 32. +
               (315 * quick_pow(sintheta, 4)) / 64.;
    }
    else if (lp == 5 && l == 5 && s == 0 && j == 5)
    {
        return (165 * costheta) / 64. + (385 * quick_pow(costheta, 3)) / 128. + (693 * quick_pow(costheta, 5)) / 128. -
               (1155 * costheta * quick_pow(sintheta, 2)) / 128. -
               (3465 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 64. +
               (3465 * costheta * quick_pow(sintheta, 4)) / 128.;
    }
    else if (lp == 6 && l == 6 && s == 0 && j == 6)
    {
        return 1.26953125 + (1365 * quick_pow(costheta, 2)) / 512. + (819 * quick_pow(costheta, 4)) / 256. +
               (3003 * quick_pow(costheta, 6)) / 512. - (1365 * quick_pow(sintheta, 2)) / 512. -
               (2457 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 128. -
               (45045 * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 512. +
               (819 * quick_pow(sintheta, 4)) / 256. +
               (45045 * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 512. -
               (3003 * quick_pow(sintheta, 6)) / 512.;
    }
    else if (lp == 7 && l == 7 && s == 0 && j == 7)
    {
        return (2625 * costheta) / 1024. + (2835 * quick_pow(costheta, 3)) / 1024. +
               (3465 * quick_pow(costheta, 5)) / 1024. + (6435 * quick_pow(costheta, 7)) / 1024. -
               (8505 * costheta * quick_pow(sintheta, 2)) / 1024. -
               (17325 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 512. -
               (135135 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 1024. +
               (17325 * costheta * quick_pow(sintheta, 4)) / 1024. +
               (225225 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 1024. -
               (45045 * costheta * quick_pow(sintheta, 6)) / 1024.;
    }
    else if (lp == 8 && l == 8 && s == 0 && j == 8)
    {
        return 1.27105712890625 + (5355 * quick_pow(costheta, 2)) / 2048. + (11781 * quick_pow(costheta, 4)) / 4096. +
               (7293 * quick_pow(costheta, 6)) / 2048. + (109395 * quick_pow(costheta, 8)) / 16384. -
               (5355 * quick_pow(sintheta, 2)) / 2048. -
               (35343 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 2048. -
               (109395 * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 2048. -
               (765765 * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 4096. +
               (11781 * quick_pow(sintheta, 4)) / 4096. +
               (109395 * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 2048. +
               (3828825 * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 8192. -
               (7293 * quick_pow(sintheta, 6)) / 2048. -
               (765765 * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 4096. +
               (109395 * quick_pow(sintheta, 8)) / 16384.;
    }
    else if (lp == 9 && l == 9 && s == 0 && j == 9)
    {
        return (41895 * costheta) / 16384. + (21945 * quick_pow(costheta, 3)) / 8192. +
               (24453 * quick_pow(costheta, 5)) / 8192. + (122265 * quick_pow(costheta, 7)) / 32768. +
               (230945 * quick_pow(costheta, 9)) / 32768. - (65835 * costheta * quick_pow(sintheta, 2)) / 8192. -
               (122265 * quick_pow(costheta, 3) * quick_pow(sintheta, 2)) / 4096. -
               (2567565 * quick_pow(costheta, 5) * quick_pow(sintheta, 2)) / 32768. -
               (2078505 * quick_pow(costheta, 7) * quick_pow(sintheta, 2)) / 8192. +
               (122265 * costheta * quick_pow(sintheta, 4)) / 8192. +
               (4279275 * quick_pow(costheta, 3) * quick_pow(sintheta, 4)) / 32768. +
               (14549535 * quick_pow(costheta, 5) * quick_pow(sintheta, 4)) / 16384. -
               (855855 * costheta * quick_pow(sintheta, 6)) / 32768. -
               (4849845 * quick_pow(costheta, 3) * quick_pow(sintheta, 6)) / 8192. +
               (2078505 * costheta * quick_pow(sintheta, 8)) / 32768.;
    }
    else if (lp == 10 && l == 10 && s == 0 && j == 10)
    {
        return 1.2718048095703125 + (169785 * quick_pow(costheta, 2)) / 65536. +
               (45045 * quick_pow(costheta, 4)) / 16384. + (405405 * quick_pow(costheta, 6)) / 131072. +
               (255255 * quick_pow(costheta, 8)) / 65536. + (969969 * quick_pow(costheta, 10)) / 131072. -
               (169785 * quick_pow(sintheta, 2)) / 65536. -
               (135135 * quick_pow(costheta, 2) * quick_pow(sintheta, 2)) / 8192. -
               (6081075 * quick_pow(costheta, 4) * quick_pow(sintheta, 2)) / 131072. -
               (1786785 * quick_pow(costheta, 6) * quick_pow(sintheta, 2)) / 16384. -
               (43648605 * quick_pow(costheta, 8) * quick_pow(sintheta, 2)) / 131072. +
               (45045 * quick_pow(sintheta, 4)) / 16384. +
               (6081075 * quick_pow(costheta, 2) * quick_pow(sintheta, 4)) / 131072. +
               (8933925 * quick_pow(costheta, 4) * quick_pow(sintheta, 4)) / 32768. +
               (101846745 * quick_pow(costheta, 6) * quick_pow(sintheta, 4)) / 65536. -
               (405405 * quick_pow(sintheta, 6)) / 131072. -
               (1786785 * quick_pow(costheta, 2) * quick_pow(sintheta, 6)) / 16384. -
               (101846745 * quick_pow(costheta, 4) * quick_pow(sintheta, 6)) / 65536. +
               (255255 * quick_pow(sintheta, 8)) / 65536. +
               (43648605 * quick_pow(costheta, 2) * quick_pow(sintheta, 8)) / 131072. -
               (969969 * quick_pow(sintheta, 10)) / 131072.;
    }

    return 0;
}

std::complex<double> m11(double t, const NN::NN_configs &configs, std::vector<double> &phase_uncoupled,
                         std::vector<double> &phase_coupled)
{
    double t_lab_now = configs.tlab_dcs;
    double qq = sqrt(configs.mass_proton * configs.mass_proton * t_lab_now * (t_lab_now + 2.0 * configs.mass_neutron) /
                     ((configs.mass_proton + configs.mass_neutron) * (configs.mass_proton + configs.mass_neutron) +
                      2.0 * t_lab_now * configs.mass_proton));
    std::complex<double> temp{0, 0};
    // uncoupled channels
    double delta = phase_uncoupled[1];
    temp = temp + coeff_m11(t, 1, 1, 1, 0) * (exp(2.0 * im_unit * delta) - 1.0) / (2.0 * im_unit * qq); // 3P0 channel
    for (int jtemp = 1; jtemp < configs.J_uncoupled_max + 1; jtemp = jtemp + 1)
    {
        delta = phase_uncoupled[2 * jtemp + 1];
        temp = temp + coeff_m11(t, jtemp, jtemp, 1, jtemp) * (exp(2.0 * im_unit * delta) - 1.0) / (2.0 * im_unit * qq);
    }
    // coupled channels
    for (int jtemp = 1; jtemp < configs.J_coupled_max + 1; jtemp = jtemp + 1)
    {
        double deltam = phase_coupled[3 * jtemp - 3];
        double deltap = phase_coupled[3 * jtemp - 2];
        double eps = phase_coupled[3 * jtemp - 1];
        temp = temp + coeff_m11(t, jtemp - 1, jtemp - 1, 1, jtemp) *
                          (cos(2.0 * eps) * exp(2.0 * im_unit * deltam) - 1.0) / (2.0 * im_unit * qq);
        temp = temp + coeff_m11(t, jtemp + 1, jtemp + 1, 1, jtemp) *
                          (cos(2.0 * eps) * exp(2.0 * im_unit * deltap) - 1.0) / (2.0 * im_unit * qq);
        temp = temp + coeff_m11(t, jtemp - 1, jtemp + 1, 1, jtemp) * sin(2.0 * eps) / (2.0 * qq) *
                          exp(im_unit * (deltam + deltap));
        temp = temp + coeff_m11(t, jtemp + 1, jtemp - 1, 1, jtemp) * sin(2.0 * eps) / (2.0 * qq) *
                          exp(im_unit * (deltam + deltap));
    }
    return temp;
}

std::complex<double> m10(double t, const NN::NN_configs &configs, std::vector<double> &phase_uncoupled,
                         std::vector<double> &phase_coupled)
{
    double t_lab_now = configs.tlab_dcs;
    double qq = sqrt(configs.mass_proton * configs.mass_proton * t_lab_now * (t_lab_now + 2.0 * configs.mass_neutron) /
                     ((configs.mass_proton + configs.mass_neutron) * (configs.mass_proton + configs.mass_neutron) +
                      2.0 * t_lab_now * configs.mass_proton));
    std::complex<double> temp;
    // uncoupled channels
    double delta = phase_uncoupled[1];
    temp = temp + coeff_m10(t, 1, 1, 1, 0) * (exp(2.0 * im_unit * delta) - 1.0) / (2.0 * im_unit * qq); // 3P0 channel
    for (int jtemp = 1; jtemp < configs.J_uncoupled_max + 1; jtemp = jtemp + 1)
    {
        delta = phase_uncoupled[2 * jtemp + 1];
        temp = temp + coeff_m10(t, jtemp, jtemp, 1, jtemp) * (exp(2.0 * im_unit * delta) - 1.0) / (2.0 * im_unit * qq);
    }
    // coupled channels
    for (int jtemp = 1; jtemp < configs.J_coupled_max + 1; jtemp = jtemp + 1)
    {
        double deltam = phase_coupled[3 * jtemp - 3];
        double deltap = phase_coupled[3 * jtemp - 2];
        double eps = phase_coupled[3 * jtemp - 1];
        temp = temp + coeff_m10(t, jtemp - 1, jtemp - 1, 1, jtemp) *
                          (cos(2.0 * eps) * exp(2.0 * im_unit * deltam) - 1.0) / (2.0 * im_unit * qq);
        temp = temp + coeff_m10(t, jtemp + 1, jtemp + 1, 1, jtemp) *
                          (cos(2.0 * eps) * exp(2.0 * im_unit * deltap) - 1.0) / (2.0 * im_unit * qq);
        temp = temp + coeff_m10(t, jtemp - 1, jtemp + 1, 1, jtemp) * sin(2.0 * eps) / (2.0 * qq) *
                          exp(im_unit * (deltam + deltap));
        temp = temp + coeff_m10(t, jtemp + 1, jtemp - 1, 1, jtemp) * sin(2.0 * eps) / (2.0 * qq) *
                          exp(im_unit * (deltam + deltap));
    }
    return temp;
}

std::complex<double> mpm(double t, const NN::NN_configs &configs, std::vector<double> &phase_uncoupled,
                         std::vector<double> &phase_coupled)
{
    double t_lab_now = configs.tlab_dcs;
    double qq = sqrt(configs.mass_proton * configs.mass_proton * t_lab_now * (t_lab_now + 2.0 * configs.mass_neutron) /
                     ((configs.mass_proton + configs.mass_neutron) * (configs.mass_proton + configs.mass_neutron) +
                      2.0 * t_lab_now * configs.mass_proton));
    std::complex<double> temp;
    // uncoupled channels
    double delta = phase_uncoupled[1];
    temp = temp + coeff_mpm(t, 1, 1, 1, 0) * (exp(2.0 * im_unit * delta) - 1.0) / (2.0 * im_unit * qq); // 3P0 channel
    for (int jtemp = 1; jtemp < configs.J_uncoupled_max + 1; jtemp = jtemp + 1)
    {
        delta = phase_uncoupled[2 * jtemp + 1];
        temp = temp + coeff_mpm(t, jtemp, jtemp, 1, jtemp) * (exp(2.0 * im_unit * delta) - 1.0) / (2.0 * im_unit * qq);
    }
    // coupled channels
    for (int jtemp = 1; jtemp < configs.J_coupled_max + 1; jtemp = jtemp + 1)
    {
        double deltam = phase_coupled[3 * jtemp - 3];
        double deltap = phase_coupled[3 * jtemp - 2];
        double eps = phase_coupled[3 * jtemp - 1];
        temp = temp + coeff_mpm(t, jtemp - 1, jtemp - 1, 1, jtemp) *
                          (cos(2.0 * eps) * exp(2.0 * im_unit * deltam) - 1.0) / (2.0 * im_unit * qq);
        temp = temp + coeff_mpm(t, jtemp + 1, jtemp + 1, 1, jtemp) *
                          (cos(2.0 * eps) * exp(2.0 * im_unit * deltap) - 1.0) / (2.0 * im_unit * qq);
        temp = temp + coeff_mpm(t, jtemp - 1, jtemp + 1, 1, jtemp) * sin(2.0 * eps) / (2.0 * qq) *
                          exp(im_unit * (deltam + deltap));
        temp = temp + coeff_mpm(t, jtemp + 1, jtemp - 1, 1, jtemp) * sin(2.0 * eps) / (2.0 * qq) *
                          exp(im_unit * (deltam + deltap));
    }
    return temp;
}

std::complex<double> m01(double t, const NN::NN_configs &configs, std::vector<double> &phase_uncoupled,
                         std::vector<double> &phase_coupled)
{
    double t_lab_now = configs.tlab_dcs;
    double qq = sqrt(configs.mass_proton * configs.mass_proton * t_lab_now * (t_lab_now + 2.0 * configs.mass_neutron) /
                     ((configs.mass_proton + configs.mass_neutron) * (configs.mass_proton + configs.mass_neutron) +
                      2.0 * t_lab_now * configs.mass_proton));
    std::complex<double> temp;
    // uncoupled channels
    double delta = phase_uncoupled[1];
    temp = temp + coeff_m01(t, 1, 1, 1, 0) * (exp(2.0 * im_unit * delta) - 1.0) / (2.0 * im_unit * qq); // 3P0 channel
    for (int jtemp = 1; jtemp < configs.J_uncoupled_max + 1; jtemp = jtemp + 1)
    {
        delta = phase_uncoupled[2 * jtemp + 1];
        temp = temp + coeff_m01(t, jtemp, jtemp, 1, jtemp) * (exp(2.0 * im_unit * delta) - 1.0) / (2.0 * im_unit * qq);
    }
    // coupled channels
    for (int jtemp = 1; jtemp < configs.J_coupled_max + 1; jtemp = jtemp + 1)
    {
        double deltam = phase_coupled[3 * jtemp - 3];
        double deltap = phase_coupled[3 * jtemp - 2];
        double eps = phase_coupled[3 * jtemp - 1];
        temp = temp + coeff_m01(t, jtemp - 1, jtemp - 1, 1, jtemp) *
                          (cos(2.0 * eps) * exp(2.0 * im_unit * deltam) - 1.0) / (2.0 * im_unit * qq);
        temp = temp + coeff_m01(t, jtemp + 1, jtemp + 1, 1, jtemp) *
                          (cos(2.0 * eps) * exp(2.0 * im_unit * deltap) - 1.0) / (2.0 * im_unit * qq);
        temp = temp + coeff_m01(t, jtemp - 1, jtemp + 1, 1, jtemp) * sin(2.0 * eps) / (2.0 * qq) *
                          exp(im_unit * (deltam + deltap));
        temp = temp + coeff_m01(t, jtemp + 1, jtemp - 1, 1, jtemp) * sin(2.0 * eps) / (2.0 * qq) *
                          exp(im_unit * (deltam + deltap));
    }
    return temp;
}

std::complex<double> m00(double t, const NN::NN_configs &configs, std::vector<double> &phase_uncoupled,
                         std::vector<double> &phase_coupled)
{
    double t_lab_now = configs.tlab_dcs;
    double qq = sqrt(configs.mass_proton * configs.mass_proton * t_lab_now * (t_lab_now + 2.0 * configs.mass_neutron) /
                     ((configs.mass_proton + configs.mass_neutron) * (configs.mass_proton + configs.mass_neutron) +
                      2.0 * t_lab_now * configs.mass_proton));
    std::complex<double> temp;
    // uncoupled channels
    double delta = phase_uncoupled[1];
    temp = temp + coeff_m00(t, 1, 1, 1, 0) * (exp(2.0 * im_unit * delta) - 1.0) / (2.0 * im_unit * qq); // 3P0 channel
    for (int jtemp = 1; jtemp < configs.J_uncoupled_max + 1; jtemp = jtemp + 1)
    {
        delta = phase_uncoupled[2 * jtemp + 1];
        temp = temp + coeff_m00(t, jtemp, jtemp, 1, jtemp) * (exp(2.0 * im_unit * delta) - 1.0) / (2.0 * im_unit * qq);
    }
    // coupled channels
    for (int jtemp = 1; jtemp < configs.J_coupled_max + 1; jtemp = jtemp + 1)
    {
        double deltam = phase_coupled[3 * jtemp - 3];
        double deltap = phase_coupled[3 * jtemp - 2];
        double eps = phase_coupled[3 * jtemp - 1];
        temp = temp + coeff_m00(t, jtemp - 1, jtemp - 1, 1, jtemp) *
                          (cos(2.0 * eps) * exp(2.0 * im_unit * deltam) - 1.0) / (2.0 * im_unit * qq);
        temp = temp + coeff_m00(t, jtemp + 1, jtemp + 1, 1, jtemp) *
                          (cos(2.0 * eps) * exp(2.0 * im_unit * deltap) - 1.0) / (2.0 * im_unit * qq);
        temp = temp + coeff_m00(t, jtemp - 1, jtemp + 1, 1, jtemp) * sin(2.0 * eps) / (2.0 * qq) *
                          exp(im_unit * (deltam + deltap));
        temp = temp + coeff_m00(t, jtemp + 1, jtemp - 1, 1, jtemp) * sin(2.0 * eps) / (2.0 * qq) *
                          exp(im_unit * (deltam + deltap));
    }
    return temp;
}

std::complex<double> mss(double t, const NN::NN_configs &configs, std::vector<double> &phase_uncoupled,
                         std::vector<double> &phase_coupled)
{
    double t_lab_now = configs.tlab_dcs;
    double qq = sqrt(configs.mass_proton * configs.mass_proton * t_lab_now * (t_lab_now + 2.0 * configs.mass_neutron) /
                     ((configs.mass_proton + configs.mass_neutron) * (configs.mass_proton + configs.mass_neutron) +
                      2.0 * t_lab_now * configs.mass_proton));
    std::complex<double> temp;
    // only uncoupled channels
    double delta = phase_uncoupled[0];
    temp = temp + coeff_mss(t, 0, 0, 0, 0) * (exp(2.0 * im_unit * delta) - 1.0) / (2.0 * im_unit * qq); // 1S0 channel
    for (int jtemp = 1; jtemp < configs.J_uncoupled_max + 1; jtemp = jtemp + 1)
    {
        delta = phase_uncoupled[2 * jtemp];
        temp = temp + coeff_mss(t, jtemp, jtemp, 0, jtemp) * (exp(2.0 * im_unit * delta) - 1.0) / (2.0 * im_unit * qq);
    }
    return temp;
}

} // namespace mmatrix

#endif // MMATRIX_HPP
