#pragma once
#ifndef OBSERVABLE_HPP
#define OBSERVABLE_HPP

#include "indicators/cursor_control.hpp"
#include "indicators/progress_bar.hpp"
#include "lib_define.hpp"
#include "mmatrix.hpp"
#include "phase.hpp"
#include "precalculate.hpp"
#include "precalculate_SFR.hpp"

namespace observable
{
using namespace indicators;

constexpr std::complex<double> im_unit{0, 1};
constexpr double ppi = 3.14159265358979323846;

void differential_cross_section(const NN::NN_configs &configs)
{
    double factor = 0.3894; // transform factor from [GeV]^(-2) to mb.
    double angle_trans = 180.0 / ppi;
    double t_lab_now = configs.tlab_dcs;
    precalculate::Interaction_Container container_interaction;
    std::vector<double> container_sfr_interaction_1 = precalculate_sfr::Interaction_SFR_Container(configs, 1);
    std::vector<double> container_sfr_interaction_2 = precalculate_sfr::Interaction_SFR_Container(configs, 2);
    std::vector<double> container_sfr_interaction_3 = precalculate_sfr::Interaction_SFR_Container(configs, 3);
    std::vector<double> container_sfr_interaction_4 = precalculate_sfr::Interaction_SFR_Container(configs, 4);
    std::vector<double> container_sfr_interaction_5 = precalculate_sfr::Interaction_SFR_Container(configs, 5);
    std::vector<double> container_sfr_interaction_6 = precalculate_sfr::Interaction_SFR_Container(configs, 6);
    std::vector<double> phase1 = phases::get_phase_uncoupled_np(
        configs, t_lab_now, container_interaction, container_sfr_interaction_1, container_sfr_interaction_2,
        container_sfr_interaction_3, container_sfr_interaction_4, container_sfr_interaction_5,
        container_sfr_interaction_6);
    std::vector<double> phase2 = phases::get_phase_bar_coupled_np(
        configs, t_lab_now, container_interaction, container_sfr_interaction_1, container_sfr_interaction_2,
        container_sfr_interaction_3, container_sfr_interaction_4, container_sfr_interaction_5,
        container_sfr_interaction_6);

    std::vector<double> temp;
    double theta_min = 0;
    double theta_max = ppi;
    int interv = 180;
    std::vector<double> thetalist = util::make_table(theta_min, theta_max, interv);
    for (int i = 0; i < thetalist.size(); i = i + 1)
    {
        double t = cos(thetalist[i]);
        std::complex<double> M11 = mmatrix::m11(t, configs, phase1, phase2);
        std::complex<double> M10 = mmatrix::m10(t, configs, phase1, phase2);
        std::complex<double> M01 = mmatrix::m01(t, configs, phase1, phase2);
        std::complex<double> M00 = mmatrix::m00(t, configs, phase1, phase2);
        std::complex<double> Mpm = mmatrix::mpm(t, configs, phase1, phase2);
        std::complex<double> Mss = mmatrix::mss(t, configs, phase1, phase2);
        double X = 0.5 * (norm(M11) + norm(M10) + norm(Mpm) + norm(M01) + 0.5 * norm(M00) + 0.5 * norm(Mss));
        temp.push_back(X);
    }

    for (int i = 0; i < thetalist.size(); i = i + 1)
    {
        std::cout << angle_trans * thetalist[i] << ", ";
    }
    std::cout << "\n\n";
    for (int i = 0; i < temp.size(); i = i + 1)
    {
        std::cout << factor * temp[i] << ", ";
    }
    std::cout << "\n\n";
}

void total_cross_section(NN::NN_configs &configs)
{
    double factor = 0.3894; // transform factor from [GeV]^(-2) to mb.
    double tstart = configs.tlab_min;
    double tend = configs.tlab_max;
    int interv = int((tend - tstart) / configs.tlab_dis);
    double t_lab_now;

    precalculate::Interaction_Container container_interaction;
    std::vector<double> container_sfr_interaction_1 = precalculate_sfr::Interaction_SFR_Container(configs, 1);
    std::vector<double> container_sfr_interaction_2 = precalculate_sfr::Interaction_SFR_Container(configs, 2);
    std::vector<double> container_sfr_interaction_3 = precalculate_sfr::Interaction_SFR_Container(configs, 3);
    std::vector<double> container_sfr_interaction_4 = precalculate_sfr::Interaction_SFR_Container(configs, 4);
    std::vector<double> container_sfr_interaction_5 = precalculate_sfr::Interaction_SFR_Container(configs, 5);
    std::vector<double> container_sfr_interaction_6 = precalculate_sfr::Interaction_SFR_Container(configs, 6);

    std::vector<double> tlablist = util::make_table(tstart, tend, interv);

    // Hide cursor
    show_console_cursor(false);
    // Progress Bar settings
    ProgressBar bar{option::BarWidth{50},
                    option::Start{"["},
                    option::Fill{"■"},
                    option::Lead{"■"},
                    option::Remainder{"-"},
                    option::End{" ]"},
                    option::ForegroundColor{Color::green},
                    option::ShowPercentage{true},
                    option::ShowElapsedTime{true},
                    option::ShowRemainingTime{true},
                    option::PrefixText{"Calculating cross section "},
                    option::ProgressType{ProgressType::incremental},
                    option::MaxProgress{tlablist.size()},
                    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}};

    std::vector<double> cross_list;

    double angle_max_fixed = 1;  // max cos(theta)
    double angle_min_fixed = -1; // min cos(theta)
    int seeds_theta = 100;       // number of points in integrating over cos(theta)
    double inter = (angle_max_fixed - angle_min_fixed) / seeds_theta;

    for (int in = 0; in < tlablist.size(); in = in + 1)
    {
        // bar progress update
        bar.set_progress(in);

        t_lab_now = tlablist[in];
        configs.tlab_dcs = t_lab_now; // update tlab in configs.
        std::vector<double> phase1 = phases::get_phase_uncoupled_np(
            configs, t_lab_now, container_interaction, container_sfr_interaction_1, container_sfr_interaction_2,
            container_sfr_interaction_3, container_sfr_interaction_4, container_sfr_interaction_5,
            container_sfr_interaction_6);
        std::vector<double> phase2 = phases::get_phase_bar_coupled_np(
            configs, t_lab_now, container_interaction, container_sfr_interaction_1, container_sfr_interaction_2,
            container_sfr_interaction_3, container_sfr_interaction_4, container_sfr_interaction_5,
            container_sfr_interaction_6);

        double temp = 0;
        for (int ii = 0; ii < seeds_theta; ii = ii + 1)
        {
            double t = angle_min_fixed + (ii + 0.5) * inter; // cos(theta)
            std::complex<double> M11 = mmatrix::m11(t, configs, phase1, phase2);
            std::complex<double> M10 = mmatrix::m10(t, configs, phase1, phase2);
            std::complex<double> M01 = mmatrix::m01(t, configs, phase1, phase2);
            std::complex<double> M00 = mmatrix::m00(t, configs, phase1, phase2);
            std::complex<double> Mpm = mmatrix::mpm(t, configs, phase1, phase2);
            std::complex<double> Mss = mmatrix::mss(t, configs, phase1, phase2);
            double X = 0.5 * (norm(M11) + norm(M10) + norm(Mpm) + norm(M01) + 0.5 * norm(M00) + 0.5 * norm(Mss));
            temp = temp + 2.0 * ppi * factor * inter * X;
        }
        temp = temp * (2.0 / (angle_max_fixed - angle_min_fixed)); // normalize
        cross_list.push_back(temp);
    }
    bar.mark_as_completed();
    // Show cursor
    show_console_cursor(true);

    std::cout.precision(6);
    for (int i = 0; i < tlablist.size(); i = i + 1)
    {
        std::cout << tlablist[i] << ", ";
    }
    std::cout << "\n\n";
    for (int i = 0; i < cross_list.size(); i = i + 1)
    {
        std::cout << cross_list[i] << ", ";
    }
    std::cout << "\n\n";
}

void spin_observables(const NN::NN_configs &configs)
{
    double factor = 0.3894; // transform factor from [GeV]^(-2) to mb.
    double angle_trans = 180.0 / ppi;
    double t_lab_now = configs.tlab_dcs;
    precalculate::Interaction_Container container_interaction;
    std::vector<double> container_sfr_interaction_1 = precalculate_sfr::Interaction_SFR_Container(configs, 1);
    std::vector<double> container_sfr_interaction_2 = precalculate_sfr::Interaction_SFR_Container(configs, 2);
    std::vector<double> container_sfr_interaction_3 = precalculate_sfr::Interaction_SFR_Container(configs, 3);
    std::vector<double> container_sfr_interaction_4 = precalculate_sfr::Interaction_SFR_Container(configs, 4);
    std::vector<double> container_sfr_interaction_5 = precalculate_sfr::Interaction_SFR_Container(configs, 5);
    std::vector<double> container_sfr_interaction_6 = precalculate_sfr::Interaction_SFR_Container(configs, 6);
    std::vector<double> phase1 = phases::get_phase_uncoupled_np(
        configs, t_lab_now, container_interaction, container_sfr_interaction_1, container_sfr_interaction_2,
        container_sfr_interaction_3, container_sfr_interaction_4, container_sfr_interaction_5,
        container_sfr_interaction_6);
    std::vector<double> phase2 = phases::get_phase_bar_coupled_np(
        configs, t_lab_now, container_interaction, container_sfr_interaction_1, container_sfr_interaction_2,
        container_sfr_interaction_3, container_sfr_interaction_4, container_sfr_interaction_5,
        container_sfr_interaction_6);

    std::vector<double> DSG, D, P, A, R, Rp, Axx, Azz, Axz;
    double theta_min = 0;
    double theta_max = ppi;
    int interv = 180;
    std::vector<double> thetalist = util::make_table(theta_min, theta_max, interv);
    for (int i = 0; i < thetalist.size(); i = i + 1)
    {
        double the = thetalist[i];
        double t = cos(the);
        double tt = sin(the);
        std::complex<double> M11 = mmatrix::m11(t, configs, phase1, phase2);
        std::complex<double> M10 = mmatrix::m10(t, configs, phase1, phase2);
        std::complex<double> M01 = mmatrix::m01(t, configs, phase1, phase2);
        std::complex<double> M00 = mmatrix::m00(t, configs, phase1, phase2);
        std::complex<double> Mpm = mmatrix::mpm(t, configs, phase1, phase2);
        std::complex<double> Mss = mmatrix::mss(t, configs, phase1, phase2);
        double I0 = 0.5 * (norm(M11) + norm(M10) + norm(Mpm) + norm(M01) + 0.5 * norm(M00) + 0.5 * norm(Mss)); // I0
        DSG.push_back(I0);
        double I01D = 0.25 * norm(M11 + Mpm - Mss) + 0.25 * norm(M11 - Mpm - M00) + 0.5 * norm(M10 + M01); // I0(1-D)
        D.push_back(1 - I01D / I0);
        double I0P = sqrt(2) / 4.0 * real(im_unit * (M10 - M01) * conj(M11 - Mpm + M00));                  // I0P
        P.push_back(I0P / I0);
        double I0A = -0.5 *
                     real((M00 + sqrt(2) * (t + 1) / tt * M10) * conj(M11 + Mpm + Mss) -
                          sqrt(2) / tt * (M10 + M01) * conj(M11 + Mpm)) *
                     sin(0.5 * the); // I0A
        A.push_back(I0A / I0);
        double I0R = 0.5 *
                     real((M00 + sqrt(2) * (t - 1) / tt * M10) * conj(M11 + Mpm + Mss) +
                          sqrt(2) / tt * (M10 + M01) * conj(Mss)) *
                     cos(0.5 * the); // I0R
        R.push_back(I0R / I0);
        double I0Rp = 0.5 *
                      real((M00 + sqrt(2) * (t + 1) / tt * M10) * conj(M11 + Mpm + Mss) -
                           sqrt(2) / tt * (M10 + M01) * conj(Mss)) *
                      sin(0.5 * the); // I0Rp
        Rp.push_back(I0Rp / I0);
        double I0Axx = 0.25 * norm(M00) - 0.25 * norm(Mss) - 0.5 * norm(M01) + 0.5 * norm(M10) + real(M11 * conj(Mpm));
        Axx.push_back(I0Axx / I0);
        double I0Azz =
            0.5 * norm(M11) - 0.25 * norm(M00) - 0.25 * norm(Mss) + 0.5 * norm(M01) - 0.5 * norm(M10) + 0.5 * norm(Mpm);
        Azz.push_back(I0Azz / I0);
        double I0Axz = 0.25 * tan(the) * (norm(M11 - Mpm) - norm(M00)) - 0.5 / tan(the) * (norm(M01) - norm(M10));
        Axz.push_back(I0Axz / I0);
    }

    std::cout << "theta:\n";
    for (int i = 0; i < thetalist.size(); i = i + 1)
    {
        std::cout << angle_trans * thetalist[i] << ", ";
    }
    std::cout << "\n\n";

    std::cout << "DSG:\n";
    for (int i = 0; i < DSG.size(); i = i + 1)
    {
        std::cout << factor * DSG[i] << ", ";
    }
    std::cout << "\n\n";

    std::cout << "D:\n";
    for (int i = 0; i < D.size(); i = i + 1)
    {
        std::cout << D[i] << ", ";
    }
    std::cout << "\n\n";

    std::cout << "P:\n";
    for (int i = 0; i < P.size(); i = i + 1)
    {
        std::cout << P[i] << ", ";
    }
    std::cout << "\n\n";

    std::cout << "A:\n";
    for (int i = 0; i < A.size(); i = i + 1)
    {
        std::cout << A[i] << ", ";
    }
    std::cout << "\n\n";

    std::cout << "R:\n";
    for (int i = 0; i < R.size(); i = i + 1)
    {
        std::cout << R[i] << ", ";
    }
    std::cout << "\n\n";

    std::cout << "R':\n";
    for (int i = 0; i < Rp.size(); i = i + 1)
    {
        std::cout << Rp[i] << ", ";
    }
    std::cout << "\n\n";

    std::cout << "Axx:\n";
    for (int i = 0; i < Axx.size(); i = i + 1)
    {
        std::cout << Axx[i] << ", ";
    }
    std::cout << "\n\n";

    std::cout << "Azz:\n";
    for (int i = 0; i < Azz.size(); i = i + 1)
    {
        std::cout << Azz[i] << ", ";
    }
    std::cout << "\n\n";

    std::cout << "Axz:\n";
    for (int i = 0; i < Axz.size(); i = i + 1)
    {
        std::cout << Axz[i] << ", ";
    }
    std::cout << "\n\n";
}

} // namespace observable

#endif // OBSERVABLE_HPP
