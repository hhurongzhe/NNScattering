#pragma once
#ifndef EVALUATE_HPP
#define EVALUATE_HPP

#include "indicators/cursor_control.hpp"
#include "indicators/progress_bar.hpp"
#include "lib_define.hpp"
#include "phase.hpp"
#include "precalculate.hpp"
#include "precalculate_SFR.hpp"

namespace evaluate
{
using namespace indicators;

void evaluate_phase_np(const NN::NN_configs &configs)
{
    double tstart = configs.tlab_min;
    double tend = configs.tlab_max;
    int interv = int((tend - tstart) / configs.tlab_dis);
    int de = configs.mesh_points_LS;
    double angle_trans = 180.0 / constants::pi;
    double t_lab_now;

    // interaction container, which is used to store precalculated interaction matrix elements.
    precalculate::Interaction_Container container_interaction;
    std::vector<double> container_sfr_interaction_1 = precalculate_sfr::Interaction_SFR_Container(configs, 1);
    std::vector<double> container_sfr_interaction_2 = precalculate_sfr::Interaction_SFR_Container(configs, 2);
    std::vector<double> container_sfr_interaction_3 = precalculate_sfr::Interaction_SFR_Container(configs, 3);
    std::vector<double> container_sfr_interaction_4 = precalculate_sfr::Interaction_SFR_Container(configs, 4);
    std::vector<double> container_sfr_interaction_5 = precalculate_sfr::Interaction_SFR_Container(configs, 5);
    std::vector<double> container_sfr_interaction_6 = precalculate_sfr::Interaction_SFR_Container(configs, 6);

    std::vector<double> tlablist = util::make_table(tstart, tend, interv);

    //* if you just need to evaluate some specific points:
    // std::vector<double> tlablist = {1., 5., 10., 25., 50., 100., 150., 200., 250., 300., 350., 400., 450., 500.};
    // for (int t = 0; t < tlablist.size(); t = t + 1)
    // {
    //     tlablist[t] = tlablist[t] * 0.001;
    // }

    //* store calculated phase shifts in matrix by Eigen:
    Eigen::MatrixXd phase_uncoupled(tlablist.size(), 2 * configs.J_uncoupled_max + 2);
    Eigen::MatrixXd phase_coupled(tlablist.size(), 3 * configs.J_uncoupled_max);

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
                    option::PrefixText{"Calculating phase shifts "},
                    option::ProgressType{ProgressType::incremental},
                    option::MaxProgress{tlablist.size()},
                    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}};

    // note that in, ii and  jj are dummy variables.
    for (int in = 0; in < tlablist.size(); in = in + 1)
    {
        t_lab_now = tlablist[in];

        // bar progress update
        bar.set_progress(in);

        std::vector<double> phase1 = phases::get_phase_uncoupled_np(
            configs, t_lab_now, container_interaction, container_sfr_interaction_1, container_sfr_interaction_2,
            container_sfr_interaction_3, container_sfr_interaction_4, container_sfr_interaction_5,
            container_sfr_interaction_6);
        std::vector<double> phase2 = phases::get_phase_bar_coupled_np(
            configs, t_lab_now, container_interaction, container_sfr_interaction_1, container_sfr_interaction_2,
            container_sfr_interaction_3, container_sfr_interaction_4, container_sfr_interaction_5,
            container_sfr_interaction_6);
        for (int ii = 0; ii < 2 * configs.J_uncoupled_max + 2; ii = ii + 1)
        {
            phase_uncoupled(in, ii) = angle_trans * phase1[ii];
        }
        for (int ii = 0; ii < 3 * configs.J_coupled_max; ii = ii + 1)
        {
            phase_coupled(in, ii) = angle_trans * phase2[ii];
        }
    }

    bar.mark_as_completed();
    // Show cursor
    show_console_cursor(true);

    // write results in the output file.
    std::ofstream fp(configs.result_file());
    std::cout.precision(6);
    auto name_uncoupled = phases::partial_wave_name_uncoupled(configs);
    auto name_coupled = phases::partial_wave_name_coupled(configs);
    // print partial wave names on the first line in output file:
    fp << std::fixed << "Tlab";
    for (int ii = 0; ii < name_uncoupled.size(); ii = ii + 1)
    {
        fp << "\t" << std::setw(15) << name_uncoupled[ii];
    }
    for (int ii = 0; ii < name_coupled.size(); ii = ii + 1)
    {
        fp << "\t" << std::setw(15) << name_coupled[ii];
    }
    fp << "\n";
    for (int ii = 0; ii < tlablist.size(); ii = ii + 1)
    {
        fp << std::fixed << 1000.0 * tlablist[ii];
        for (int jj = 0; jj < 2 * configs.J_uncoupled_max + 2; jj = jj + 1)
        {
            fp << "\t" << std::setw(15) << phase_uncoupled(ii, jj);
        }
        for (int jj = 0; jj < 3 * configs.J_coupled_max; jj = jj + 1)
        {
            fp << "\t" << std::setw(15) << phase_coupled(ii, jj);
        }
        fp << '\n';
    }
    fp.close();
}

} // namespace evaluate

#endif