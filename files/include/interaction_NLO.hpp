#pragma once
#ifndef INTERACTION_NLO_HPP
#define INTERACTION_NLO_HPP

#include "interaction_aPWD.hpp"
#include "lib_define.hpp"
#include "state.hpp"

namespace interaction_NLO
{
constexpr double ppi = 3.14159265358979323846;
constexpr double twopicubic = 248.0502134423985614038105;

// NLO meson exchange potential.
double V_NLO_me_np(states::LSJ_State s1, states::LSJ_State s2, const NN::NN_configs &configs)
{
    double temp = 0.0;

    const double pmag = s1.momentum;
    const double ppmag = s2.momentum;
    const double regulator4 = interaction_aPWD::regulator_function_order(pmag, ppmag, 4, configs);
    const double ga = configs.axial_current_coupling_constant;
    const double fpi = configs.pion_decay_constant;
    const double mpi = configs.mass_pion_averaged;

    double isospin_factor;
    if ((s1.l + s1.s) % 2 == 0)
    {
        isospin_factor = 1.0;
    }
    else
    {
        isospin_factor = -3.0;
    }

    for (int i = 0; i < 22; i = i + 1)
    {
        const double x = gauss_legendre::zeros22[i];
        const double w = gauss_legendre::weights22[i];
        const double qmag = sqrt(pmag * pmag + ppmag * ppmag - 2.0 * pmag * ppmag * x);
        double f1;
        double f2;
        const double f3 = 0.0;
        const double f4 = 0.0;
        const double f5 = 0.0;
        double f6;
        double f1_fac = interaction_aPWD::loop_function_L(qmag, configs) / (384.0 * ppi * ppi * pow(fpi, 4));
        double f1_part1 = 4.0 * pow(mpi, 2) * (1.0 + 4.0 * pow(ga, 2) - 5.0 * pow(ga, 4));
        double f1_part2 = qmag * qmag * (1.0 + 10.0 * pow(ga, 2) - 23.0 * pow(ga, 4));
        double f1_part3 = -48.0 * pow(ga, 4) * pow(mpi, 4) / (4.0 * pow(mpi, 2) + qmag * qmag);
        f1 = isospin_factor * f1_fac * (f1_part1 + f1_part2 + f1_part3);
        f2 = 3.0 * pow(ga, 4) * qmag * qmag / (64.0 * ppi * ppi * pow(fpi, 4)) *
             interaction_aPWD::loop_function_L(qmag, configs);
        f6 = -f2 / (qmag * qmag);

        const double fa = interaction_aPWD::V_kernal_SYM(s1, s2, f1, f2, f3, f4, f5, f6, x);
        temp = temp + fa * w;
    }
    return regulator4 * temp;
}

// NLO contact potential.
double V_NLO_contact_np(states::LSJ_State s1, states::LSJ_State s2, const NN::NN_configs &configs)
{
    const double pmag = s1.momentum;
    const double ppmag = s2.momentum;

    const double regulator4 = interaction_aPWD::regulator_function_order(pmag, ppmag, 4, configs);
    const double regulator6 = interaction_aPWD::regulator_function_order(pmag, ppmag, 6, configs);

    if (s1 == s2 && s1.l == 0 && s1.s == 0 && s1.j == 0)
    {
        return regulator4 * configs.C_NLO_1s0 * (pmag * pmag + ppmag * ppmag); // 1S0 channel
    }
    else if (s1 == s2 && s1.l == 1 && s1.s == 1 && s1.j == 0)
    {
        return regulator4 * configs.C_NLO_3p0 * pmag * ppmag; // 3P0 channel
    }
    else if (s1 == s2 && s1.l == 1 && s1.s == 0 && s1.j == 1)
    {
        return regulator4 * configs.C_NLO_1p1 * pmag * ppmag; // 1P1 channel
    }
    else if (s1 == s2 && s1.l == 1 && s1.s == 1 && s1.j == 1)
    {
        return regulator6 * configs.C_NLO_3p1 * pmag * ppmag; // 3P1 channel, with n=3.
    }
    else if (s1 == s2 && s1.l == 0 && s1.s == 1 && s1.j == 1)
    {
        return regulator4 * configs.C_NLO_3s1 * (pmag * pmag + ppmag * ppmag); // 3S1 channel
    }
    else if (s2.l == 0 && s2.s == 1 && s2.j == 1 && s1.l == 2 && s1.s == 1 && s1.j == 1)
    {
        return regulator4 * configs.C_NLO_3s13d1 * pmag * pmag; // 3S1-3D1 channel
    }
    else if (s2.l == 2 && s2.s == 1 && s2.j == 1 && s1.l == 0 && s1.s == 1 && s1.j == 1)
    {
        return regulator4 * configs.C_NLO_3s13d1 * ppmag * ppmag; // 3D1-3S1 channel
    }
    else if (s1 == s2 && s1.l == 1 && s1.s == 1 && s1.j == 2)
    {
        return regulator4 * configs.C_NLO_3p2 * pmag * ppmag; // 3P2 channel
    }
    else
    {
        return 0.0;
    }
}

double V_NLO_np(states::LSJ_State s1, states::LSJ_State s2, const NN::NN_configs &configs)
{
    const double muu = configs.mass_nucleon;
    double p1, p2;
    p1 = s1.momentum;
    p2 = s2.momentum;
    const double e1 = sqrt(muu * muu + p1 * p1);
    const double e2 = sqrt(muu * muu + p2 * p2);
    const double fac = sqrt(muu / e1) * sqrt((muu / e2));
    return fac * (V_NLO_contact_np(s1, s2, configs) + V_NLO_me_np(s1, s2, configs)) / twopicubic;
}

} // namespace interaction_NLO

#endif // INTERACTION_NLO_HPP
