#pragma once
#ifndef INTERACTION_NNLO_HPP
#define INTERACTION_NNLO_HPP

#include "interaction_aPWD.hpp"
#include "lib_define.hpp"
#include "state.hpp"

namespace interaction_NNLO
{
constexpr double ppi = 3.14159265358979323846;
constexpr double twopicubic = 248.0502134423985614038105;

// NNLO meson exchange potential.
double V_NNLO_me_np(states::LSJ_State s1, states::LSJ_State s2, const NN::NN_configs &configs)
{
    double temp = 0.0;

    double pmag = s1.momentum;
    double ppmag = s2.momentum;
    double regulator4 = interaction_aPWD::regulator_function_order(pmag, ppmag, 4, configs);
    double isospin_factor;
    if ((s1.l + s1.s) % 2 == 0)
    {
        isospin_factor = 1.0;
    }
    else
    {
        isospin_factor = -3.0;
    }

    double gaga = pow(configs.axial_current_coupling_constant, 2);
    double ffff = pow(configs.pion_decay_constant, 4);
    double mm = pow(configs.mass_pion_averaged, 2);

    for (int i = 0; i < 22; i = i + 1)
    {
        double x = gauss_legendre::zeros22[i];
        double w = gauss_legendre::weights22[i];
        double qmag = sqrt(pmag * pmag + ppmag * ppmag - 2.0 * pmag * ppmag * x);
        double f1;
        double f2;
        double f3 = 0.0;
        double f4 = 0.0;
        double f5 = 0.0;
        double f6;
        double f1_part1 = 3.0 * gaga / (16.0 * ppi * ffff);
        double f1_part2 = 2.0 * mm * (configs.C_3 - 2.0 * configs.C_1) + configs.C_3 * qmag * qmag;
        double f1_part3 = 2.0 * mm + qmag * qmag;
        double f1_part4 = interaction_aPWD::loop_function_A(qmag, configs);
        f1 = f1_part1 * f1_part2 * f1_part3 * f1_part4;
        f2 = isospin_factor * gaga * qmag * qmag / (32.0 * ppi * ffff) * configs.C_4 * (4.0 * mm + qmag * qmag) *
             interaction_aPWD::loop_function_A(qmag, configs);
        f6 = -f2 / (qmag * qmag);

        double fa = interaction_aPWD::V_kernal_SYM(s1, s2, f1, f2, f3, f4, f5, f6, x);

        temp = temp + fa * w;
    }
    return regulator4 * temp;
}

double V_NNLO_contact_np(states::LSJ_State s1, states::LSJ_State s2, const NN::NN_configs &configs) { return 0.0; }

double V_NNLO_np(states::LSJ_State s1, states::LSJ_State s2, const NN::NN_configs &configs)
{
    double muu = configs.mass_nucleon;
    double p1, p2;
    p1 = s1.momentum;
    p2 = s2.momentum;
    double e1 = sqrt(muu * muu + p1 * p1);
    double e2 = sqrt(muu * muu + p2 * p2);
    double fac = sqrt(muu / e1) * sqrt((muu / e2));
    return fac * (V_NNLO_contact_np(s1, s2, configs) + V_NNLO_me_np(s1, s2, configs)) / twopicubic;
}

} // namespace interaction_NNLO

#endif // INTERACTION_NNLO_HPP
