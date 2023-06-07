#pragma once
#ifndef NN_CONFIGS_HPP
#define NN_CONFIGS_HPP

#include "mymath.hpp"
#include "util/inifile.hpp"

namespace NN
{

struct NN_configs
{

    // ***** interaction section *****
    double axial_current_coupling_constant;
    double pion_decay_constant;                                                            // in MeV
    double mb_coupling;
    double C_LO_1s0, C_LO_3s1;                                                             // in 10^4 GeV^(-2)
    double C_NLO_1s0, C_NLO_3p0, C_NLO_1p1, C_NLO_3p1, C_NLO_3s1, C_NLO_3s13d1, C_NLO_3p2; // in 10^4 GeV^(-4)
    double Dhat_NNNLO_1s0, D_NNNLO_1s0, D_NNNLO_3p0, D_NNNLO_1p1, D_NNNLO_3p1, Dhat_NNNLO_3s1, D_NNNLO_3s1, D_NNNLO_3d1,
        Dhat_NNNLO_3s13d1, D_NNNLO_3s13d1, D_NNNLO_1d2, D_NNNLO_3d2, D_NNNLO_3p2, D_NNNLO_3p23f2,
        D_NNNLO_3d3;                          // in 10^4 GeV^(-6)
    double C_1, C_2, C_3, C_4;                // in GeV^(-1)
    double d1_plus_d2, d3, d5, d14_minus_d15; // in GeV^(-2)
    double Lambda;                            // in MeV
    double Lambda_tilde;                      // in MeV

    // ***** kinetic parameters section ****
    double tlab_min; // all in MeV
    double tlab_max;
    double tlab_dis;
    double tlab_dcs;

    // ***** meson masses section ****
    double mass_pion_charged; // in MeV
    double mass_pion_neutral;
    double mass_pion_averaged;

    // ***** baryon masses section ****
    double mass_proton; // in MeV
    double mass_neutron;
    double mass_nucleon;

    // ***** numerical parameters section ****

    // number of points used for solving LS equation.
    int mesh_points_LS;
    std::vector<double> momentum_mesh_points;
    std::vector<double> momentum_mesh_weights;

    // constant C in Gaussian integration in Lippmann-Schwinger Equation
    // in order of 1000 MeV
    double C_Gaussian; // in MeV

    int J_uncoupled_max;
    int J_coupled_max;

    // ***** output section *****
    std::string result_dir;
    std::string result_name;

    // constructor
    NN_configs(const inifile_system::inifile &ini);

    // generate a result file name
    std::string result_file() const;
};

NN_configs::NN_configs(const inifile_system::inifile &ini)
{

    // ***** interaction section *****
    auto sec = ini.section("interaction");
    axial_current_coupling_constant = sec.get_double("axial_current_coupling_constant");
    pion_decay_constant = 0.001 * sec.get_double("pion_decay_constant");
    mb_coupling = axial_current_coupling_constant / (2.0 * pion_decay_constant);

    C_LO_1s0 = 10000.0 * sec.get_double("C_LO_1s0");
    C_LO_3s1 = 10000.0 * sec.get_double("C_LO_3s1");

    C_NLO_1s0 = 10000.0 * sec.get_double("C_NLO_1s0");
    C_NLO_3p0 = 10000.0 * sec.get_double("C_NLO_3p0");
    C_NLO_1p1 = 10000.0 * sec.get_double("C_NLO_1p1");
    C_NLO_3p1 = 10000.0 * sec.get_double("C_NLO_3p1");
    C_NLO_3s1 = 10000.0 * sec.get_double("C_NLO_3s1");
    C_NLO_3s13d1 = 10000.0 * sec.get_double("C_NLO_3s13d1");
    C_NLO_3p2 = 10000.0 * sec.get_double("C_NLO_3p2");

    Dhat_NNNLO_1s0 = 10000.0 * sec.get_double("Dhat_NNNLO_1s0");
    D_NNNLO_1s0 = 10000.0 * sec.get_double("D_NNNLO_1s0");
    D_NNNLO_3p0 = 10000.0 * sec.get_double("D_NNNLO_3p0");
    D_NNNLO_1p1 = 10000.0 * sec.get_double("D_NNNLO_1p1");
    D_NNNLO_3p1 = 10000.0 * sec.get_double("D_NNNLO_3p1");
    Dhat_NNNLO_3s1 = 10000.0 * sec.get_double("Dhat_NNNLO_3s1");
    D_NNNLO_3s1 = 10000.0 * sec.get_double("D_NNNLO_3s1");
    D_NNNLO_3d1 = 10000.0 * sec.get_double("D_NNNLO_3d1");
    Dhat_NNNLO_3s13d1 = 10000.0 * sec.get_double("Dhat_NNNLO_3s13d1");
    D_NNNLO_3s13d1 = 10000.0 * sec.get_double("D_NNNLO_3s13d1");
    D_NNNLO_1d2 = 10000.0 * sec.get_double("D_NNNLO_1d2");
    D_NNNLO_3d2 = 10000.0 * sec.get_double("D_NNNLO_3d2");
    D_NNNLO_3p2 = 10000.0 * sec.get_double("D_NNNLO_3p2");
    D_NNNLO_3p23f2 = 10000.0 * sec.get_double("D_NNNLO_3p23f2");
    D_NNNLO_3d3 = 10000.0 * sec.get_double("D_NNNLO_3d3");

    C_1 = sec.get_double("C_1");
    C_2 = sec.get_double("C_2");
    C_3 = sec.get_double("C_3");
    C_4 = sec.get_double("C_4");

    d1_plus_d2 = sec.get_double("d1_plus_d2");
    d3 = sec.get_double("d3");
    d5 = sec.get_double("d5");
    d14_minus_d15 = sec.get_double("d14_minus_d15");

    Lambda = 0.001 * sec.get_double("Lambda");
    Lambda_tilde = 0.001 * sec.get_double("Lambda_tilde");

    // ***** kinetic parameters section ****
    sec = ini.section("kinetics");
    tlab_min = 0.001 * sec.get_double("tlab_min");
    tlab_max = 0.001 * sec.get_double("tlab_max");
    tlab_dis = 0.001 * sec.get_double("tlab_dis");
    tlab_dcs = 0.001 * sec.get_double("tlab_dcs");

    // ***** meson masses section ****
    sec = ini.section("meson-masses");
    mass_pion_charged = 0.001 * sec.get_double("mass_pion_charged");
    mass_pion_neutral = 0.001 * sec.get_double("mass_pion_neutral");
    mass_pion_averaged = 0.001 * sec.get_double("mass_pion_averaged");

    // ***** baryon masses section ****
    sec = ini.section("baryon-masses");
    mass_proton = 0.001 * sec.get_double("mass_proton");
    mass_neutron = 0.001 * sec.get_double("mass_neutron");
    mass_nucleon = 0.001 * sec.get_double("mass_nucleon");

    // ***** numerical parameters section ****
    sec = ini.section("numerical-parameters");
    mesh_points_LS = sec.get_int("mesh_points_LS");
    C_Gaussian = 0.001 * sec.get_double("C_Gaussian");
    J_uncoupled_max = sec.get_int("J_uncoupled_max");
    J_coupled_max = sec.get_int("J_coupled_max");
    momentum_mesh_points = gauss_legendre::make_klist(mesh_points_LS, C_Gaussian);
    momentum_mesh_weights = gauss_legendre::make_slist(mesh_points_LS, C_Gaussian);

    // ***** output section *****
    sec = ini.section("output");
    result_dir = sec.get_string("result_dir");
    result_name = sec.get_string("result_name");
    if (result_dir.back() != '/')
    {
        result_dir += "/";
    }
};

std::string file_stem(const std::string &file)
{
    auto p1 = file.find_last_of('/') + 1;
    auto p2 = file.find_last_of('.');
    return file.substr(p1, p2 - p1);
}

std::string NN_configs::result_file() const
{
    std::ostringstream oss;
    oss << result_dir << result_name << ".trace";
    return oss.str();
}

} // namespace NN

#endif // NN_CONFIGS_HPP