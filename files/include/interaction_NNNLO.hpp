#pragma once
#ifndef INTERACTION_NNNLO_HPP
#define INTERACTION_NNNLO_HPP

#include "interaction_aPWD.hpp"
#include "lib_define.hpp"
#include "state.hpp"

namespace interaction_NNNLO
{
constexpr double ppi = 3.14159265358979323846;
constexpr double twopicubic = 248.0502134423985614038105;

double spectral_intergral_CS(double q, const NN::NN_configs &configs, std::vector<double> &container_sfr_interaction)
{
    constexpr int n = 2;
    double fact = -2.0 * pow(q, 6) / (ppi);
    double a = n * configs.mass_pion_averaged;
    double b = configs.Lambda_tilde;
    double temp = 0;
    int degree = 100; //! 100 meshpoints
    double xx;
    // #pragma omp parallel for private(xx) reduction(+ : temp)
    for (int i = 0; i < degree; i = i + 1)
    {
        xx = 0.5 * (b - a) * gauss_legendre::zeros100[i] + 0.5 * (b + a);
        temp = temp + gauss_legendre::weights100[i] * container_sfr_interaction[i] / (pow(xx, 5) * (xx * xx + q * q));
    }
    temp = 0.5 * (b - a) * temp;
    return fact * temp;
}

double spectral_intergral_T(double q, const NN::NN_configs &configs, std::vector<double> &container_sfr_interaction)
{
    constexpr int n = 2;
    double fact = 2.0 * pow(q, 4) / (ppi);
    double a = n * configs.mass_pion_averaged;
    double b = configs.Lambda_tilde;
    double temp = 0;
    int degree = 100; //! 100 meshpoints
    double xx;
    // #pragma omp parallel for private(xx) reduction(+ : temp)
    for (int i = 0; i < degree; i = i + 1)
    {
        xx = 0.5 * (b - a) * gauss_legendre::zeros100[i] + 0.5 * (b + a);
        temp = temp + gauss_legendre::weights100[i] * container_sfr_interaction[i] / (pow(xx, 3) * (xx * xx + q * q));
    }
    temp = 0.5 * (b - a) * temp;
    return fact * temp;
}

double V_NNNLO_me_np(states::LSJ_State s1, states::LSJ_State s2, const NN::NN_configs &configs,
                     std::vector<double> &container_sfr_interaction_1, std::vector<double> &container_sfr_interaction_2,
                     std::vector<double> &container_sfr_interaction_3, std::vector<double> &container_sfr_interaction_4,
                     std::vector<double> &container_sfr_interaction_5, std::vector<double> &container_sfr_interaction_6)
{
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
    double gggg = gaga * gaga;
    double ff = pow(configs.pion_decay_constant, 2);
    double ffff = pow(configs.pion_decay_constant, 4);
    double mpi = configs.mass_pion_averaged;
    double mm = pow(configs.mass_pion_averaged, 2);
    double mcmc = pow(configs.mass_pion_charged, 2);
    double mmmm = pow(configs.mass_pion_averaged, 4);
    double mmmmm = pow(configs.mass_pion_averaged, 5);
    double mmmmmm = pow(configs.mass_pion_averaged, 6);
    double c1 = configs.C_1;
    double c2 = configs.C_2;
    double c3 = configs.C_3;
    double c4 = configs.C_4;
    double x, weight, qmag, ww, wwww, loopA, loopL;
    double f1 = 0, f2 = 0, f3 = 0, f4 = 0, f5 = 0, f6 = 0;
    double f1_part1 = 0, f1_part2 = 0, f1_football_diagram = 0, f1_two_loop_diagram_Vc = 0, f1_two_loop_diagram_Wc = 0,
           f1_relativistic_corrections_Vc = 0, f1_relativistic_corrections_Wc = 0, f2_football_diagram = 0,
           f2_two_loop_diagram_Vs = 0, f2_two_loop_diagram_Ws = 0, f2_relativistic_corrections_Vs = 0,
           f2_relativistic_corrections_Ws = 0, f3_relativistic_corrections_Vls = 0, f3_relativistic_corrections_Wls = 0,
           f6_football_diagram = 0, f6_two_loop_diagram_Vt = 0, f6_two_loop_diagram_Wt = 0,
           f6_relativistic_corrections_Vt = 0, f6_relativistic_corrections_Wt = 0;
    double f1_cimn_Vc = 0, f1_cimn_Wc = 0, f2_cimn_Ws = 0, f6_cimn_Wt = 0, f3_cimn_Vls = 0, f3_cimn_Wls = 0;
    double f6_pi_gamma = 0;
    double bb = 0;
    double fa = 0;
    double temp = 0.0;
#pragma omp parallel for private(                                                                                 \
        x, weight, qmag, ww, wwww, loopA, loopL, f1, f2, f3, f4, f5, f6, f1_part1, f1_part2, f1_football_diagram, \
            f1_two_loop_diagram_Vc, f1_two_loop_diagram_Wc, f1_relativistic_corrections_Vc,                       \
            f1_relativistic_corrections_Wc, f2_football_diagram, f2_two_loop_diagram_Vs, f2_two_loop_diagram_Ws,  \
            f2_relativistic_corrections_Vs, f2_relativistic_corrections_Ws, f3_relativistic_corrections_Vls,      \
            f3_relativistic_corrections_Wls, f6_football_diagram, f6_two_loop_diagram_Vt, f6_two_loop_diagram_Wt, \
            f6_relativistic_corrections_Vt, f6_relativistic_corrections_Wt, f1_cimn_Vc, f1_cimn_Wc, f2_cimn_Ws,   \
            f6_cimn_Wt, f3_cimn_Vls, f3_cimn_Wls, fa) reduction(+ : temp)
    for (int i = 0; i < 22; i = i + 1)
    {
        x = gauss_legendre::zeros22[i];
        weight = gauss_legendre::weights22[i];
        qmag = sqrt(pmag * pmag + ppmag * ppmag - 2.0 * pmag * ppmag * x);
        ww = 4.0 * mm + qmag * qmag;
        wwww = ww * ww;
        loopA = interaction_aPWD::loop_function_A(qmag, configs);
        loopL = interaction_aPWD::loop_function_L(qmag, configs);
        f1 = 0.0;
        f2 = 0.0;
        f3 = 0.0;
        f4 = 0.0;
        f5 = 0.0;
        f6 = 0.0;

        // football diagram
        f1_part1 = c2 * ww / 6.0 + c3 * (2.0 * mm + qmag * qmag) - 4.0 * c1 * mm;
        f1_part2 = c2 * c2 * wwww / 45.0;
        f1_football_diagram = (3.0 / (16.0 * ppi * ppi * ffff)) * (f1_part1 * f1_part1 + f1_part2) * loopL;
        f6_football_diagram = isospin_factor * pow(c4, 2) / (96.0 * ppi * ppi * ffff) * ww * loopL;
        f2_football_diagram = -qmag * qmag * f6_football_diagram;

        // two loop diagram
        f1_two_loop_diagram_Vc = spectral_intergral_CS(qmag, configs, container_sfr_interaction_1);
        f1_two_loop_diagram_Wc = isospin_factor * spectral_intergral_CS(qmag, configs, container_sfr_interaction_2);
        f2_two_loop_diagram_Vs = spectral_intergral_CS(qmag, configs, container_sfr_interaction_3);
        f6_two_loop_diagram_Vt = spectral_intergral_T(qmag, configs, container_sfr_interaction_4);
        f2_two_loop_diagram_Ws = isospin_factor * spectral_intergral_CS(qmag, configs, container_sfr_interaction_5);
        f6_two_loop_diagram_Wt = isospin_factor * spectral_intergral_T(qmag, configs, container_sfr_interaction_6);

        // relativistic corrections
        f1_relativistic_corrections_Vc = 3.0 * gggg / (128.0 * ppi * ffff * configs.mass_nucleon) *
                                         (mmmmm / (2.0 * ww) + (2.0 * mm + qmag * qmag) * (qmag * qmag - mm) * loopA);
        f1_relativistic_corrections_Wc =
            isospin_factor * gaga / (64.0 * ppi * ffff * configs.mass_nucleon) *
            (3.0 * gaga * mmmmm / (2.0 * ww) +
             (gaga * (3.0 * mm + 2.0 * qmag * qmag) - 2.0 * mm - qmag * qmag) * (2.0 * mm + qmag * qmag) * loopA);
        f6_relativistic_corrections_Vt =
            3.0 * gggg / (256.0 * ppi * ffff * configs.mass_nucleon) * (5.0 * mm + 2.0 * qmag * qmag) * loopA;
        f2_relativistic_corrections_Vs = -qmag * qmag * f6_relativistic_corrections_Vt;
        f6_relativistic_corrections_Wt = isospin_factor * gaga / (128.0 * ppi * ffff * configs.mass_nucleon) *
                                         (gaga * (3.0 * mm + qmag * qmag) - ww) * loopA;
        f2_relativistic_corrections_Ws = -qmag * qmag * f6_relativistic_corrections_Wt;
        f3_relativistic_corrections_Vls =
            0.5 * 3.0 * gggg / (32.0 * ppi * ffff * configs.mass_nucleon) * (2.0 * mm + qmag * qmag) * loopA;
        f3_relativistic_corrections_Wls =
            0.5 * gaga * (1.0 - gaga) / (32.0 * ppi * ffff * configs.mass_nucleon) * ww * loopA;

        // ci/MN terms
        f1_cimn_Vc = -gaga * loopL / (32.0 * ppi * ppi * configs.mass_nucleon * ffff) *
                     ((c2 - 6.0 * c3) * pow(qmag, 4) + 4.0 * (6.0 * c1 + c2 - 3.0 * c3) * qmag * qmag * mm +
                      6.0 * (c2 - 2.0 * c3) * mmmm + 24.0 * (2.0 * c1 + c3) * mmmmmm / ww);
        f1_cimn_Wc = -c4 * qmag * qmag * loopL / (192.0 * ppi * ppi * configs.mass_nucleon * ffff) *
                     (gaga * (8.0 * mm + 5.0 * qmag * qmag) + ww) * isospin_factor;
        f6_cimn_Wt = -c4 * loopL / (192.0 * ppi * ppi * configs.mass_nucleon * ffff) *
                     (gaga * (16.0 * mm + 7.0 * qmag * qmag) - ww) * isospin_factor;
        f2_cimn_Ws = -qmag * qmag * f6_cimn_Wt;
        f3_cimn_Vls = 0.5 * c2 * gaga / (8.0 * ppi * ppi * configs.mass_nucleon * ffff) * ww * loopL;
        f3_cimn_Wls = -0.5 * c4 * loopL / (48.0 * ppi * ppi * configs.mass_nucleon * ffff) *
                      (gaga * (8.0 * mm + 5.0 * qmag * qmag) + ww) * isospin_factor;

        // pi-gamma terms,
        // only for np channel!
        bb = pow(qmag / configs.mass_pion_charged, 2);
        f6_pi_gamma = -gaga / (4.0 * ff * mcmc) * (isospin_factor + 1.0) * (constants::fine_struct_const / ppi) *
                      (-(1 - bb) * (1 - bb) / (2 * bb * bb * (1 + bb)) * (log(1 + bb)) + 1 / (2 * bb));

        // total
        f1 = f1_football_diagram + f1_two_loop_diagram_Vc + f1_two_loop_diagram_Wc + f1_relativistic_corrections_Vc +
             f1_relativistic_corrections_Wc + f1_cimn_Vc + f1_cimn_Wc;
        f2 = f2_football_diagram + f2_two_loop_diagram_Vs + f2_two_loop_diagram_Ws + f2_relativistic_corrections_Vs +
             f2_relativistic_corrections_Ws + f2_cimn_Ws;
        f3 = f3_relativistic_corrections_Vls + f3_relativistic_corrections_Wls + f3_cimn_Vls + f3_cimn_Wls;
        f6 = f6_football_diagram + f6_two_loop_diagram_Vt + f6_two_loop_diagram_Wt + f6_relativistic_corrections_Vt +
             f6_relativistic_corrections_Wt + f6_cimn_Wt + f6_pi_gamma;

        fa = interaction_aPWD::V_kernal_SYM(s1, s2, f1, f2, f3, f4, f5, f6, x);
        temp = temp + fa * weight;
    }
    return regulator4 * temp;
}

double V_NNNLO_contact_np(states::LSJ_State s1, states::LSJ_State s2, const NN::NN_configs &configs)
{
    const double pmag = s1.momentum;
    const double ppmag = s2.momentum;

    double regulator4 = interaction_aPWD::regulator_function_order(pmag, ppmag, 4, configs);
    double regulator6 = interaction_aPWD::regulator_function_order(pmag, ppmag, 6, configs);
    double regulator8 = interaction_aPWD::regulator_function_order(pmag, ppmag, 8, configs);

    if (s1 == s2 && s1.l == 0 && s1.s == 0 && s1.j == 0)
    {
        return regulator4 * (configs.Dhat_NNNLO_1s0 * (pmag * pmag * pmag * pmag + ppmag * ppmag * ppmag * ppmag) +
                             configs.D_NNNLO_1s0 * pmag * pmag * ppmag * ppmag); // 1S0 channel
    }
    else if (s1 == s2 && s1.l == 1 && s1.s == 1 && s1.j == 0)
    {
        return regulator6 * configs.D_NNNLO_3p0 *
               (pmag * ppmag * ppmag * ppmag + ppmag * pmag * pmag * pmag); // 3P0 channel, n=3
    }
    else if (s1 == s2 && s1.l == 1 && s1.s == 0 && s1.j == 1)
    {
        return regulator4 * configs.D_NNNLO_1p1 *
               (pmag * ppmag * ppmag * ppmag + ppmag * pmag * pmag * pmag); // 1P1 channel
    }
    else if (s1 == s2 && s1.l == 1 && s1.s == 1 && s1.j == 1)
    {
        return regulator6 * configs.D_NNNLO_3p1 *
               (pmag * ppmag * ppmag * ppmag + ppmag * pmag * pmag * pmag); // 3P1 channel, n=3
    }
    else if (s1 == s2 && s1.l == 0 && s1.s == 1 && s1.j == 1)
    {
        return regulator4 * (configs.Dhat_NNNLO_3s1 * (pmag * pmag * pmag * pmag + ppmag * ppmag * ppmag * ppmag) +
                             configs.D_NNNLO_3s1 * pmag * pmag * ppmag * ppmag); // 3S1 channel
    }
    else if (s1 == s2 && s1.l == 2 && s1.s == 1 && s1.j == 1)
    {
        return regulator4 * configs.D_NNNLO_3d1 * pmag * pmag * ppmag * ppmag; // 3D1 channel
    }
    else if (s2.l == 0 && s2.s == 1 && s2.j == 1 && s1.l == 2 && s1.s == 1 && s1.j == 1)
    {
        return regulator4 * (configs.Dhat_NNNLO_3s13d1 * pmag * pmag * pmag * pmag +
                             configs.D_NNNLO_3s13d1 * pmag * pmag * ppmag * ppmag); // 3S1-3D1 channel
    }
    else if (s2.l == 2 && s2.s == 1 && s2.j == 1 && s1.l == 0 && s1.s == 1 && s1.j == 1)
    {
        return regulator4 * (configs.Dhat_NNNLO_3s13d1 * ppmag * ppmag * ppmag * ppmag +
                             configs.D_NNNLO_3s13d1 * pmag * pmag * ppmag * ppmag); // 3D1-3S1 channel
    }
    else if (s1 == s2 && s1.l == 2 && s1.s == 0 && s1.j == 2)
    {
        return regulator6 * configs.D_NNNLO_1d2 * pmag * pmag * ppmag * ppmag; // 1D2 channel, n=3
    }
    else if (s1 == s2 && s1.l == 2 && s1.s == 1 && s1.j == 2)
    {
        return regulator6 * configs.D_NNNLO_3d2 * pmag * pmag * ppmag * ppmag; // 3D2 channel, n=3
    }
    else if (s1 == s2 && s1.l == 1 && s1.s == 1 && s1.j == 2)
    {
        return regulator4 * configs.D_NNNLO_3p2 *
               (pmag * ppmag * ppmag * ppmag + ppmag * pmag * pmag * pmag); // 3P2 channel
    }
    else if (s2.l == 1 && s2.s == 1 && s2.j == 2 && s1.l == 3 && s1.s == 1 && s1.j == 2)
    {
        return regulator8 * configs.D_NNNLO_3p23f2 * ppmag * pmag * pmag * pmag; // 3P2-3F2 channel, n=4
    }
    else if (s2.l == 3 && s2.s == 1 && s2.j == 2 && s1.l == 1 && s1.s == 1 && s1.j == 2)
    {
        return regulator8 * configs.D_NNNLO_3p23f2 * pmag * ppmag * ppmag * ppmag; // 3F2-3P2 channel, n=4
    }
    else if (s1 == s2 && s1.l == 2 && s1.s == 1 && s1.j == 3)
    {
        return regulator4 * configs.D_NNNLO_3d3 * pmag * pmag * ppmag * ppmag; // 3D3 channel
    }
    else
    {
        return 0.0;
    }
}

double V_NNNLO_np(states::LSJ_State s1, states::LSJ_State s2, const NN::NN_configs &configs,
                  std::vector<double> &container_sfr_interaction_1, std::vector<double> &container_sfr_interaction_2,
                  std::vector<double> &container_sfr_interaction_3, std::vector<double> &container_sfr_interaction_4,
                  std::vector<double> &container_sfr_interaction_5, std::vector<double> &container_sfr_interaction_6)
{
    double muu = configs.mass_nucleon;
    double p1, p2;
    p1 = s1.momentum;
    p2 = s2.momentum;
    double e1 = sqrt(muu * muu + p1 * p1);
    double e2 = sqrt(muu * muu + p2 * p2);
    double fac = sqrt(muu / e1) * sqrt((muu / e2));
    return fac *
           (V_NNNLO_contact_np(s1, s2, configs) +
            V_NNNLO_me_np(s1, s2, configs, container_sfr_interaction_1, container_sfr_interaction_2,
                          container_sfr_interaction_3, container_sfr_interaction_4, container_sfr_interaction_5,
                          container_sfr_interaction_6)) /
           twopicubic;
}

} // namespace interaction_NNNLO

#endif // INTERACTION_NNNLO_HPP
