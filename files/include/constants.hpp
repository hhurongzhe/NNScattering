#pragma once
#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

namespace constants
{

constexpr double pi = 3.14159265358979323846;
constexpr double sqrt_2 = 1.41421356237309504880;
constexpr double h_si = 6.62607004e-34;                   // m^2 kg / s
constexpr double e_si = 1.602176634e-19;                  // C
constexpr double c_si = 299792458.0;                      // m / s
constexpr double hbar_si = h_si / (2 * pi);               // J s
constexpr double hbarc_si = hbar_si * c_si;               // J m
constexpr double hbar = hbar_si / (e_si * 1e6 * 1e-21);   // MeV zs
constexpr double hbarc = hbarc_si / (e_si * 1e6 * 1e-15); // MeV fm
const double fine_struct_const = 0.00729735253;           // fine structure constant

// 1 zs = 10^-21 s

} // end namespace constants

#endif // CONSTANTS_HPP