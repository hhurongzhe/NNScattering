#pragma once
#ifndef PRECALCULATE_SFR_HPP
#define PRECALCULATE_SFR_HPP

#include "lib_define.hpp"

// this section provides precalculating corresponding values in SFR integral,
// which can efficiently improve the speed of the code.
namespace precalculate_sfr
{
constexpr double ppi = 3.14159265358979323846;

// function ImVc in two loop contributions.
double two_loop_functions_1(double mu, const NN::NN_configs &configs)
{
    double gg = pow(configs.axial_current_coupling_constant, 2);
    double gggg = gg * gg;
    double mpi = configs.mass_pion_averaged;
    double ffffff = pow(configs.pion_decay_constant, 6);
    double part_1 = 3.0 * gggg * (2.0 * mpi * mpi - mu * mu) / (ppi * mu * 4096.0 * ffffff);
    double part_2 = (mpi * mpi - 2.0 * mu * mu);
    double part_3 = 2.0 * mpi + (2.0 * mpi * mpi - mu * mu) / (2.0 * mu) * log((mu + 2.0 * mpi) / (mu - 2.0 * mpi));
    double part_4 = 4.0 * gg * mpi * (2.0 * mpi * mpi - mu * mu);
    double result = part_1 * (part_2 * part_3 + part_4);
    return result;
}

// function ImWc in two loop contributions.
double two_loop_functions_2(double mu, const NN::NN_configs &configs)
{
    double gg = pow(configs.axial_current_coupling_constant, 2);
    double gggg = gg * gg;
    double mpi = configs.mass_pion_averaged;
    double ff = pow(configs.pion_decay_constant, 2);
    double ffffff = pow(configs.pion_decay_constant, 6);
    double kk = mu * mu / 4.0 - mpi * mpi;
    double k = sqrt(kk);
    double d1d2 = configs.d1_plus_d2;
    double d3 = configs.d3;
    double d5 = configs.d5;
    double fact = 2.0 * k / (3.0 * mu * pow(8.0 * ppi, 3) * ffffff);
    double a = 0.0;
    double b = 1.0;
    double temp = 0;
    int degree = 12; //! 12 meshpoints
    double x, sqrtfart, fvalue, fpart1, fpart2, fpart3, fpart4, fpart5, fpart6, fpart7;
    // #pragma omp parallel for private(x, sqrtfart, fvalue, fpart1, fpart2, fpart3, fpart4, fpart5, fpart6, fpart7)
    // reduction(+ : temp)
    for (int i = 0; i < degree; i = i + 1)
    {
        x = 0.5 * (b - a) * gauss_legendre::zeros12[i] + 0.5 * (b + a);
        sqrtfart = sqrt(mpi * mpi + kk * x * x);
        fpart1 = gg * (mu * mu - 2.0 * mpi * mpi) + 2.0 * (1.0 - gg) * kk * x * x;
        fpart2 =
            96.0 * ppi * ppi * ff * ((2.0 * mpi * mpi - mu * mu) * d1d2 - 2.0 * kk * x * x * d3 + 4.0 * mpi * mpi * d5);
        fpart3 = (4.0 * mpi * mpi * (1.0 + 2.0 * gg) - mu * mu * (1.0 + 5.0 * gg)) * (k / mu) *
                 log((mu + 2.0 * k) / (2.0 * mpi));
        fpart4 = mu * mu / 12.0 * (5.0 + 13.0 * gg) - 2.0 * mpi * mpi * (1.0 + 2.0 * gg) - 3.0 * kk * x * x;
        fpart5 = 6.0 * k * x * sqrtfart * log((k * x + sqrtfart) / (mpi));
        fpart6 = gggg * (mu * mu - 2.0 * kk * x * x - 2.0 * mpi * mpi);
        fpart7 = 5.0 / 6.0 + mpi * mpi / (kk * x * x) -
                 pow(1.0 + mpi * mpi / (kk * x * x), 1.5) * log((k * x + sqrtfart) / mpi);
        fvalue = fpart1 * (fpart2 + fpart3 + fpart4 + fpart5 + fpart6 * fpart7);
        temp = temp + gauss_legendre::weights12[i] * fvalue;
    }
    temp = 0.5 * (b - a) * temp;
    return fact * temp;
}

// function ImVs in two loop contributions.
double two_loop_functions_3(double mu, const NN::NN_configs &configs)
{
    double gg = pow(configs.axial_current_coupling_constant, 2);
    double gggg = gg * gg;
    double gggggg = gg * gg * gg;
    double mpi = configs.mass_pion_averaged;
    double ff = pow(configs.pion_decay_constant, 2);
    double ffffff = pow(configs.pion_decay_constant, 6);
    double kk = mu * mu / 4.0 - mpi * mpi;
    double k = sqrt(kk);
    double d14d15 = configs.d14_minus_d15;
    double part0 = -gg * mu * k * k * k / (8.0 * ppi * ff * ff) * d14d15;
    double fact = 2.0 * gggggg * mu * k * k * k / (pow(8.0 * ppi, 3) * ffffff);
    double a = 0.0;
    double b = 1.0;
    double temp = 0;
    int degree = 12; //! 12 meshpoints
    double x, sqrtfart, fvalue;
    // #pragma omp parallel for private(x, sqrtfart, fvalue) reduction(+ : temp)
    for (int i = 0; i < degree; i = i + 1)
    {
        x = 0.5 * (b - a) * gauss_legendre::zeros12[i] + 0.5 * (b + a);
        sqrtfart = sqrt(mpi * mpi + kk * x * x);
        fvalue = (1.0 - x * x) * (1.0 / 6.0 - mpi * mpi / (kk * x * x) +
                                  pow(1.0 + mpi * mpi / (kk * x * x), 1.5) * log((k * x + sqrtfart) / mpi));
        temp = temp + gauss_legendre::weights12[i] * fvalue;
    }
    temp = 0.5 * (b - a) * temp;
    return part0 + fact * temp;
}

// function ImVt in two loop contributions.
double two_loop_functions_4(double mu, const NN::NN_configs &configs)
{
    double gg = pow(configs.axial_current_coupling_constant, 2);
    double gggg = gg * gg;
    double gggggg = gg * gg * gg;
    double mpi = configs.mass_pion_averaged;
    double ff = pow(configs.pion_decay_constant, 2);
    double ffffff = pow(configs.pion_decay_constant, 6);
    double kk = mu * mu / 4.0 - mpi * mpi;
    double k = sqrt(kk);
    double d14d15 = configs.d14_minus_d15;
    double part0 = -gg * mu * k * k * k / (8.0 * ppi * ff * ff) * d14d15;
    double fact = 2.0 * gggggg * mu * k * k * k / (pow(8.0 * ppi, 3) * ffffff);
    double a = 0.0;
    double b = 1.0;
    double temp = 0;
    int degree = 12; //! 12 meshpoints
    double x, sqrtfart, fvalue;
    // #pragma omp parallel for private(x, sqrtfart, fvalue) reduction(+ : temp)
    for (int i = 0; i < degree; i = i + 1)
    {
        x = 0.5 * (b - a) * gauss_legendre::zeros12[i] + 0.5 * (b + a);
        sqrtfart = sqrt(mpi * mpi + kk * x * x);
        fvalue = (1.0 - x * x) * (1.0 / 6.0 - mpi * mpi / (kk * x * x) +
                                  pow(1.0 + mpi * mpi / (kk * x * x), 1.5) * log((k * x + sqrtfart) / mpi));
        temp = temp + gauss_legendre::weights12[i] * fvalue;
    }
    temp = 0.5 * (b - a) * temp;
    return (part0 + fact * temp) / (mu * mu);
}

// function ImWs in two loop contributions.
double two_loop_functions_5(double mu, const NN::NN_configs &configs)
{
    double gg = pow(configs.axial_current_coupling_constant, 2);
    double gggg = gg * gg;
    double gggggg = gg * gg * gg;
    double mpi = configs.mass_pion_averaged;
    double ff = pow(configs.pion_decay_constant, 2);
    double ffffff = pow(configs.pion_decay_constant, 6);
    double kk = mu * mu / 4.0 - mpi * mpi;
    double k = sqrt(kk);
    double part1 = gggg * (4.0 * mpi * mpi - mu * mu) / (ppi * 4096.0 * ffffff);
    double part2 = (mpi * mpi - mu * mu / 4.0) * log((mu + 2.0 * mpi) / (mu - 2.0 * mpi)) + (1.0 + 2.0 * gg) * mu * mpi;
    return part1 * part2;
}

// function ImWt in two loop contributions.
double two_loop_functions_6(double mu, const NN::NN_configs &configs)
{
    double gg = pow(configs.axial_current_coupling_constant, 2);
    double gggg = gg * gg;
    double gggggg = gg * gg * gg;
    double mpi = configs.mass_pion_averaged;
    double ff = pow(configs.pion_decay_constant, 2);
    double ffffff = pow(configs.pion_decay_constant, 6);
    double kk = mu * mu / 4.0 - mpi * mpi;
    double k = sqrt(kk);
    double part1 = gggg * (4.0 * mpi * mpi - mu * mu) / (ppi * 4096.0 * ffffff);
    double part2 = (mpi * mpi - mu * mu / 4.0) * log((mu + 2.0 * mpi) / (mu - 2.0 * mpi)) + (1.0 + 2.0 * gg) * mu * mpi;
    return part1 * part2 / (mu * mu);
}

std::vector<double> Interaction_SFR_Container(const NN::NN_configs &configs, int index)
{
    std::vector<double> temp;
    for (int i = 0; i < 100; i = i + 1)
    {
        double mu = 0.5 * (configs.Lambda_tilde - 2 * configs.mass_pion_averaged) * gauss_legendre::zeros100[i] +
                    0.5 * (configs.Lambda_tilde + 2 * configs.mass_pion_averaged);
        double value = 0;
        if (index == 1)
        {
            value = two_loop_functions_1(mu, configs);
        }
        else if (index == 2)
        {
            value = two_loop_functions_2(mu, configs);
        }
        else if (index == 3)
        {
            value = two_loop_functions_3(mu, configs);
        }
        else if (index == 4)
        {
            value = two_loop_functions_4(mu, configs);
        }
        else if (index == 5)
        {
            value = two_loop_functions_5(mu, configs);
        }
        else if (index == 6)
        {
            value = two_loop_functions_6(mu, configs);
        }
        else
        {
            std::cerr << "error in SFR!" << std::endl;
            std::exit(-1);
        }
        temp.push_back(value);
    }
    return temp;
}

} // namespace precalculate_sfr

#endif // PRECALCULATE_SFR_HPP