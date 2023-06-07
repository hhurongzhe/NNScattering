#pragma once
#ifndef LO_INTERACTION_HPP
#define LO_INTERACTION_HPP

#include "interaction_aPWD.hpp"
#include "lib_define.hpp"
#include "state.hpp"

namespace interaction_LO
{
constexpr double ppi = 3.14159265358979323846;
constexpr double twopicubic = 248.0502134423985614038105;

// LO one-meson exchange potential.
double V_LO_ome_np(states::LSJ_State s1, states::LSJ_State s2, const NN::NN_configs &configs)
{
    double temp = 0.0;

    double lp = -configs.mb_coupling * configs.mb_coupling;

    double pmag = s1.momentum;
    double ppmag = s2.momentum;
    const double regulator8 = interaction_aPWD::regulator_function_order(pmag, ppmag, 8, configs);

    for (int i = 0; i < 22; i = i + 1)
    {
        double x = gauss_legendre::zeros22[i];
        double w = gauss_legendre::weights22[i];
        double f1 = 0.0;
        double f2 = 0.0;
        double f3 = 0.0;
        double f4 = 0.0;
        double f5 = 0.0;

        double f6p0 = lp / (pmag * pmag + ppmag * ppmag - 2.0 * pmag * ppmag * x +
                            configs.mass_pion_neutral * configs.mass_pion_neutral);
        double f6pp = lp / (pmag * pmag + ppmag * ppmag - 2.0 * pmag * ppmag * x +
                            configs.mass_pion_charged * configs.mass_pion_charged);

        double f6_T0 = -f6p0 - 2.0 * f6pp;
        double f6_T1 = -f6p0 + 2.0 * f6pp;
        double f6;
        if ((s1.l + s1.s) % 2 == 0)
        {
            f6 = f6_T1;
        }
        else
        {
            f6 = f6_T0;
        }

        double fa = interaction_aPWD::V_kernal_SYM(s1, s2, f1, f2, f3, f4, f5, f6, x);
        temp = temp + fa * w;
    }
    return regulator8 * temp;
}

// LO contact potential.
double V_LO_contact_np(states::LSJ_State s1, states::LSJ_State s2, const NN::NN_configs &configs)
{
    const double pmag = s1.momentum;
    const double ppmag = s2.momentum;
    const double regulator6 = interaction_aPWD::regulator_function_order(pmag, ppmag, 6, configs);
    if (s1 == s2 && s1.l == 0 && s1.s == 0 && s1.j == 0)
    {
        return regulator6 * configs.C_LO_1s0; // 1S0 channel
    }
    else if (s1 == s2 && s1.l == 0 && s1.s == 1 && s1.j == 1)
    {
        return regulator6 * configs.C_LO_3s1; // 3S1 channel
    }
    else
    {
        return 0.0;
    }
}

// Leading order np interaction.
// Noting that a factor called "minimal relativity" is needed.
double V_LO_np(states::LSJ_State s1, states::LSJ_State s2, const NN::NN_configs &configs)
{
    double muu = configs.mass_nucleon;
    double p1, p2;
    p1 = s1.momentum;
    p2 = s2.momentum;
    double e1 = sqrt(muu * muu + p1 * p1);
    double e2 = sqrt(muu * muu + p2 * p2);
    double fac = sqrt(muu / e1) * sqrt((muu / e2));
    return fac * (V_LO_contact_np(s1, s2, configs) + V_LO_ome_np(s1, s2, configs)) / twopicubic;
}

} // namespace interaction_LO

#endif // LO_INTERACTION_HPP
