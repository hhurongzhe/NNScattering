#pragma once
#ifndef PHASE_HPP
#define PHASE_HPP

#include "lib_define.hpp"
#include "precalculate.hpp"
#include "scattering_matrix.hpp"

namespace phases
{
constexpr std::complex<double> im_unit{0, 1};
constexpr double ppi = 3.14159265358979323846;

// phase shifts for uncoupled channels,
// arrange in order: (LSJ) = (J0J), (J1J)
// except for J=0 in which only two channels 1S0 and 3P0.
std::vector<double> get_phase_uncoupled_np(
    const NN::NN_configs &configs, double t_lab_now, precalculate::Interaction_Container &container_interaction,
    std::vector<double> &container_sfr_interaction_1, std::vector<double> &container_sfr_interaction_2,
    std::vector<double> &container_sfr_interaction_3, std::vector<double> &container_sfr_interaction_4,
    std::vector<double> &container_sfr_interaction_5, std::vector<double> &container_sfr_interaction_6)
{
    int de = configs.mesh_points_LS;
    std::vector<double> phase;

    double mu = configs.mass_nucleon / 2.0; // reduced mass.
    double qq = sqrt(configs.mass_proton * configs.mass_proton * t_lab_now * (t_lab_now + 2.0 * configs.mass_neutron) /
                     ((configs.mass_proton + configs.mass_neutron) * (configs.mass_proton + configs.mass_neutron) +
                      2.0 * t_lab_now * configs.mass_proton));
    double factor = -ppi * qq * mu;

    double tempp0, tempp1;
    double temp1s0 =
        atan(factor * scattering::R_uncoupled_np(0, 0, 0, configs, t_lab_now, container_interaction,
                                                 container_sfr_interaction_1, container_sfr_interaction_2,
                                                 container_sfr_interaction_3, container_sfr_interaction_4,
                                                 container_sfr_interaction_5, container_sfr_interaction_6));
    double temp3p0 =
        atan(factor * scattering::R_uncoupled_np(1, 1, 0, configs, t_lab_now, container_interaction,
                                                 container_sfr_interaction_1, container_sfr_interaction_2,
                                                 container_sfr_interaction_3, container_sfr_interaction_4,
                                                 container_sfr_interaction_5, container_sfr_interaction_6));
    phase.push_back(temp1s0);
    phase.push_back(temp3p0);

    const int jmax = configs.J_uncoupled_max;
    if (jmax < 0)
    {
        std::cerr << "J_max should only be non-negative integer!" << std::endl;
        std::exit(-1);
    }
    if (jmax == 0)
    {
        // only two channels for jmax=0: 1S0, 3P0.
        return phase;
    }
    else
    {
        for (int j_temp = 1; j_temp < jmax + 1; j_temp = j_temp + 1)
        {
            tempp0 =
                atan(factor * scattering::R_uncoupled_np(j_temp, 0, j_temp, configs, t_lab_now, container_interaction,
                                                         container_sfr_interaction_1, container_sfr_interaction_2,
                                                         container_sfr_interaction_3, container_sfr_interaction_4,
                                                         container_sfr_interaction_5, container_sfr_interaction_6));
            tempp1 =
                atan(factor * scattering::R_uncoupled_np(j_temp, 1, j_temp, configs, t_lab_now, container_interaction,
                                                         container_sfr_interaction_1, container_sfr_interaction_2,
                                                         container_sfr_interaction_3, container_sfr_interaction_4,
                                                         container_sfr_interaction_5, container_sfr_interaction_6));
            phase.push_back(tempp0);
            phase.push_back(tempp1);
        }
    }
    return phase;
}

std::vector<double> get_phase_bar_coupled_np(
    const NN::NN_configs &configs, double t_lab_now, precalculate::Interaction_Container &container_interaction,
    std::vector<double> &container_sfr_interaction_1, std::vector<double> &container_sfr_interaction_2,
    std::vector<double> &container_sfr_interaction_3, std::vector<double> &container_sfr_interaction_4,
    std::vector<double> &container_sfr_interaction_5, std::vector<double> &container_sfr_interaction_6)
{
    int de = configs.mesh_points_LS;
    std::vector<double> phase_bar;
    double mu = configs.mass_nucleon / 2.0; // reduced mass.
    double qq = sqrt(configs.mass_proton * configs.mass_proton * t_lab_now * (t_lab_now + 2.0 * configs.mass_neutron) /
                     ((configs.mass_proton + configs.mass_neutron) * (configs.mass_proton + configs.mass_neutron) +
                      2.0 * t_lab_now * configs.mass_proton));
    double factor = -ppi * qq * mu;
    const int jmax = configs.J_coupled_max;

    std::vector<double> coupled;
    double rpp, rpm, rmm;
    double deltap, deltam, epsilon;
    double delm, delp, epsb;
    double te1, te2;

    if (jmax < 1)
    {
        std::cerr << "error in phase shifts calculation for coupled channnels!" << std::endl;
        std::exit(-1);
    }
    else
    {
        for (int j_temp = 1; j_temp < jmax + 1; j_temp = j_temp + 1)
        {
            coupled = scattering::R_coupled_triplet_np(j_temp, configs, t_lab_now, container_interaction,
                                                       container_sfr_interaction_1, container_sfr_interaction_2,
                                                       container_sfr_interaction_3, container_sfr_interaction_4,
                                                       container_sfr_interaction_5, container_sfr_interaction_6);
            rpp = coupled[0];
            rpm = coupled[2];
            rmm = coupled[3];
            epsilon = 0.5 * atan(2.0 * rpm / (rmm - rpp));
            deltap = atan(0.5 * factor * (rmm + rpp - (rmm - rpp) / (cos(2.0 * epsilon))));
            deltam = atan(0.5 * factor * (rmm + rpp + (rmm - rpp) / (cos(2.0 * epsilon))));
            epsb = 0.5 * asin(sin(2.0 * epsilon) * sin(deltam - deltap));
            te1 = deltap + deltam;
            te2 = asin(tan(2.0 * epsb) / tan(2.0 * epsilon));
            delm = 0.5 * (te1 + te2);
            delp = 0.5 * (te1 - te2);
            phase_bar.push_back(delm);
            phase_bar.push_back(delp);
            phase_bar.push_back(epsb);
        }
    }
    return phase_bar;
}

constexpr char partial_wave_name_table[] = "SPDFGHIJKLMNOX";

std::vector<std::string> partial_wave_name_uncoupled(const NN::NN_configs &configs)
{
    std::vector<std::string> temp;
    const int jmax_uncoupled = configs.J_uncoupled_max;
    std::string name_uncoupled;
    temp.push_back("1S0");
    temp.push_back("3P0");
    if (jmax_uncoupled == 0)
    {
        return temp;
    }
    else
    {
        for (int j_temp = 1; j_temp < jmax_uncoupled + 1; j_temp = j_temp + 1)
        {
            temp.push_back(std::to_string(1) + partial_wave_name_table[std::min(j_temp, 13)] + std::to_string(j_temp));
            temp.push_back(std::to_string(3) + partial_wave_name_table[std::min(j_temp, 13)] + std::to_string(j_temp));
        }
    }
    return temp;
}

std::vector<std::string> partial_wave_name_coupled(const NN::NN_configs &configs)
{
    std::vector<std::string> temp;
    const int jmax_coupled = configs.J_coupled_max;
    std::string name_coupled;
    for (int j_temp = 1; j_temp < jmax_coupled + 1; j_temp = j_temp + 1)
    {
        temp.push_back(std::to_string(3) + partial_wave_name_table[std::min(j_temp - 1, 13)] + std::to_string(j_temp));
        temp.push_back(std::to_string(3) + partial_wave_name_table[std::min(j_temp + 1, 13)] + std::to_string(j_temp));
        temp.push_back("E" + std::to_string(j_temp));
    }
    return temp;
}

} // end namespace phases

#endif // PHASE_HPP
