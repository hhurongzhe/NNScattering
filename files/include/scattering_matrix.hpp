#pragma once
#ifndef LO_SCATTERING_MATRIX_HPP
#define LO_SCATTERING_MATRIX_HPP

#include "interaction_ALL.hpp"
#include "lib_define.hpp"
#include "precalculate.hpp"

namespace scattering
{
double R_uncoupled_np(int l, int s, int j, const NN::NN_configs &configs, double t_lab_now,
                      precalculate::Interaction_Container &container_interaction,
                      std::vector<double> &container_sfr_interaction_1,
                      std::vector<double> &container_sfr_interaction_2,
                      std::vector<double> &container_sfr_interaction_3,
                      std::vector<double> &container_sfr_interaction_4,
                      std::vector<double> &container_sfr_interaction_5,
                      std::vector<double> &container_sfr_interaction_6)
{
    double mu = configs.mass_nucleon / 2.0; // reduced mass.
    double qq = sqrt(configs.mass_proton * configs.mass_proton * t_lab_now * (t_lab_now + 2.0 * configs.mass_neutron) /
                     ((configs.mass_proton + configs.mass_neutron) * (configs.mass_proton + configs.mass_neutron) +
                      2.0 * t_lab_now * configs.mass_proton));
    std::vector<double> klist = configs.momentum_mesh_points;
    klist.push_back(qq);
    int N = configs.mesh_points_LS;

    // matrix equation to be solved: AX=B, (N+1)-dimension.
    Eigen::MatrixXd AA(N + 1, N + 1);
    Eigen::MatrixXd XX;
    Eigen::VectorXd BB(N + 1);

    // initialize u_j
    std::vector<double> ulist;
    double temp = 0.0;
    for (int i = 0; i < N; i = i + 1)
    {
        double tempp = configs.momentum_mesh_weights[i] / (klist[i] * klist[i] - qq * qq);
        double ui = 2.0 * mu * klist[i] * klist[i] * tempp;
        temp = temp + tempp;
        ulist.push_back(ui);
    }
    temp = -2.0 * mu * qq * qq * temp;
    ulist.push_back(temp);

    // initialize B vector
    states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(l, s, j, j, qq);
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(l, s, j, j, klist[i]);
        double inter = interaction_ALL::V_np(
            ss1, ss2, configs, container_sfr_interaction_1, container_sfr_interaction_2, container_sfr_interaction_3,
            container_sfr_interaction_4, container_sfr_interaction_5, container_sfr_interaction_6);
        BB(i) = inter;
    }

    // initialize A matrix
    for (int i = 0; i < N + 1; i = i + 1)
    {
        for (int k = 0; k < N + 1; k = k + 1)
        {
            states::LSJ_State ss3 = states::LSJ_State::ini_LSJ_state(l, s, j, j, klist[k]);
            states::LSJ_State ss4 = states::LSJ_State::ini_LSJ_state(l, s, j, j, klist[i]);
            double inter;
            if (i != (N) && k != (N))
            {
                inter = container_interaction.get_value(configs, k, i, l, l, s, j, container_sfr_interaction_1,
                                                        container_sfr_interaction_2, container_sfr_interaction_3,
                                                        container_sfr_interaction_4, container_sfr_interaction_5,
                                                        container_sfr_interaction_6);
            }
            else
            {
                inter =
                    interaction_ALL::V_np(ss3, ss4, configs, container_sfr_interaction_1, container_sfr_interaction_2,
                                          container_sfr_interaction_3, container_sfr_interaction_4,
                                          container_sfr_interaction_5, container_sfr_interaction_6);
            }
            double a_ik = util::delta(i, k) + ulist[k] * inter;
            AA(i, k) = a_ik;
        }
    }

    XX = AA.colPivHouseholderQr().solve(BB);

    return XX(N, 0);
}

std::vector<double> R_coupled_triplet_np(
    int j, const NN::NN_configs &configs, double t_lab_now, precalculate::Interaction_Container &container_interaction,
    std::vector<double> &container_sfr_interaction_1, std::vector<double> &container_sfr_interaction_2,
    std::vector<double> &container_sfr_interaction_3, std::vector<double> &container_sfr_interaction_4,
    std::vector<double> &container_sfr_interaction_5, std::vector<double> &container_sfr_interaction_6)
{
    double mu = configs.mass_nucleon / 2.0; // reduced mass.
    double qq = sqrt(configs.mass_proton * configs.mass_proton * t_lab_now * (t_lab_now + 2.0 * configs.mass_neutron) /
                     ((configs.mass_proton + configs.mass_neutron) * (configs.mass_proton + configs.mass_neutron) +
                      2.0 * t_lab_now * configs.mass_proton));
    std::vector<double> klist = configs.momentum_mesh_points;
    klist.push_back(qq);
    const int N = configs.mesh_points_LS;

    // matrix equation to be solved: AX=B, 2(N+1)-dimension.
    Eigen::MatrixXd AA(2 * (N + 1), 2 * (N + 1));
    Eigen::MatrixXd XX;
    Eigen::MatrixXd BB(2 * (N + 1), 2);

    // initialize u_j
    std::vector<double> ulist;
    double temp = 0.0;
    for (int i = 0; i < N; i = i + 1)
    {
        double tempp = configs.momentum_mesh_weights[i] / (klist[i] * klist[i] - qq * qq);
        double ui = 2.0 * mu * klist[i] * klist[i] * tempp;
        temp = temp + tempp;
        ulist.push_back(ui);
    }
    temp = -2.0 * mu * qq * qq * temp;
    ulist.push_back(temp);

    // initialize B matrix
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, qq);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, klist[i]);
        double inter = interaction_ALL::V_np(
            ss1, ss2, configs, container_sfr_interaction_1, container_sfr_interaction_2, container_sfr_interaction_3,
            container_sfr_interaction_4, container_sfr_interaction_5, container_sfr_interaction_6);
        BB(i, 0) = inter;
    }
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, qq);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, klist[i]);
        double inter = interaction_ALL::V_np(
            ss1, ss2, configs, container_sfr_interaction_1, container_sfr_interaction_2, container_sfr_interaction_3,
            container_sfr_interaction_4, container_sfr_interaction_5, container_sfr_interaction_6);
        BB(i + N + 1, 1) = inter;
    }
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, qq);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, klist[i]);
        double inter = interaction_ALL::V_np(
            ss1, ss2, configs, container_sfr_interaction_1, container_sfr_interaction_2, container_sfr_interaction_3,
            container_sfr_interaction_4, container_sfr_interaction_5, container_sfr_interaction_6);
        BB(i, 1) = inter;
    }
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, qq);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, klist[i]);
        ss2.momentum = klist[i];
        double inter = interaction_ALL::V_np(
            ss1, ss2, configs, container_sfr_interaction_1, container_sfr_interaction_2, container_sfr_interaction_3,
            container_sfr_interaction_4, container_sfr_interaction_5, container_sfr_interaction_6);
        BB(i + N + 1, 0) = inter;
    }

    // initialize A matrix
    for (int i = 0; i < N + 1; i = i + 1)
    {
        for (int k = 0; k < N + 1; k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, klist[k]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, klist[i]);
            double inter;
            if (i != (N) && k != (N))
            {
                inter = container_interaction.get_value(configs, k, i, j + 1, j + 1, 1, j, container_sfr_interaction_1,
                                                        container_sfr_interaction_2, container_sfr_interaction_3,
                                                        container_sfr_interaction_4, container_sfr_interaction_5,
                                                        container_sfr_interaction_6);
            }
            else
            {
                inter =
                    interaction_ALL::V_np(ss1, ss2, configs, container_sfr_interaction_1, container_sfr_interaction_2,
                                          container_sfr_interaction_3, container_sfr_interaction_4,
                                          container_sfr_interaction_5, container_sfr_interaction_6);
            }
            double a_ik = util::delta(i, k) + ulist[k] * inter;
            AA(i, k) = a_ik;
        }
    }
    for (int i = N + 1; i < 2 * (N + 1); i = i + 1)
    {
        for (int k = N + 1; k < 2 * (N + 1); k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, klist[k - (N + 1)]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, klist[i - (N + 1)]);
            double inter;
            if (i != (2 * N + 1) && k != (2 * N + 1))
            {
                inter = container_interaction.get_value(configs, k - (N + 1), i - (N + 1), j - 1, j - 1, 1, j,
                                                        container_sfr_interaction_1, container_sfr_interaction_2,
                                                        container_sfr_interaction_3, container_sfr_interaction_4,
                                                        container_sfr_interaction_5, container_sfr_interaction_6);
            }
            else
            {
                inter =
                    interaction_ALL::V_np(ss1, ss2, configs, container_sfr_interaction_1, container_sfr_interaction_2,
                                          container_sfr_interaction_3, container_sfr_interaction_4,
                                          container_sfr_interaction_5, container_sfr_interaction_6);
            }
            double a_ik = util::delta(i, k) + ulist[k - (N + 1)] * inter;
            AA(i, k) = a_ik;
        }
    }
    for (int i = 0; i < (N + 1); i = i + 1)
    {
        for (int k = N + 1; k < 2 * (N + 1); k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, klist[k - (N + 1)]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, klist[i]);
            double inter;
            if (i != (N) && k != (2 * N + 1))
            {
                inter = container_interaction.get_value(configs, k - (N + 1), i, j - 1, j + 1, 1, j,
                                                        container_sfr_interaction_1, container_sfr_interaction_2,
                                                        container_sfr_interaction_3, container_sfr_interaction_4,
                                                        container_sfr_interaction_5, container_sfr_interaction_6);
            }
            else
            {
                inter =
                    interaction_ALL::V_np(ss1, ss2, configs, container_sfr_interaction_1, container_sfr_interaction_2,
                                          container_sfr_interaction_3, container_sfr_interaction_4,
                                          container_sfr_interaction_5, container_sfr_interaction_6);
            }
            double a_ik = ulist[k - (N + 1)] * inter;
            AA(i, k) = a_ik;
        }
    }
    for (int i = (N + 1); i < 2 * (N + 1); i = i + 1)
    {
        for (int k = 0; k < (N + 1); k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, klist[k]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, klist[i - (N + 1)]);
            double inter;
            if (i != (2 * N + 1) && k != (N))
            {
                inter = container_interaction.get_value(configs, k, i - (N + 1), j + 1, j - 1, 1, j,
                                                        container_sfr_interaction_1, container_sfr_interaction_2,
                                                        container_sfr_interaction_3, container_sfr_interaction_4,
                                                        container_sfr_interaction_5, container_sfr_interaction_6);
            }
            else
            {
                inter =
                    interaction_ALL::V_np(ss1, ss2, configs, container_sfr_interaction_1, container_sfr_interaction_2,
                                          container_sfr_interaction_3, container_sfr_interaction_4,
                                          container_sfr_interaction_5, container_sfr_interaction_6);
            }
            double a_ik = ulist[k] * inter;
            AA(i, k) = a_ik;
        }
    }

    XX = AA.colPivHouseholderQr().solve(BB);

    double tpp, tpm, tmp, tmm;
    tpp = XX(N, 0);
    tpm = XX(N, 1);
    tmp = XX(2 * N + 1, 0);
    tmm = XX(2 * N + 1, 1);
    std::vector<double> tem;
    tem.push_back(tpp);
    tem.push_back(tpm);
    tem.push_back(tmp);
    tem.push_back(tmm);
    return tem;
}

} // end namespace scattering

#endif // LO_SCATTERING_MATRIX_HPP
