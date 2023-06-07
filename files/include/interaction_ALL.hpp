#pragma once
#ifndef ALL_INTERACTION_HPP
#define ALL_INTERACTION_HPP

#include "interaction_LO.hpp"
#include "interaction_NLO.hpp"
#include "interaction_NNLO.hpp"
#include "interaction_NNNLO.hpp"
#include "lib_define.hpp"

namespace interaction_ALL
{

double V_np(states::LSJ_State s1, states::LSJ_State s2, const NN::NN_configs &configs,
            std::vector<double> &container_sfr_interaction_1, std::vector<double> &container_sfr_interaction_2,
            std::vector<double> &container_sfr_interaction_3, std::vector<double> &container_sfr_interaction_4,
            std::vector<double> &container_sfr_interaction_5, std::vector<double> &container_sfr_interaction_6)
{
    double temp = 0.0;
    temp = temp + interaction_LO::V_LO_np(s1, s2, configs);
    temp = temp + interaction_NLO::V_NLO_np(s1, s2, configs);
    temp = temp + interaction_NNLO::V_NNLO_np(s1, s2, configs);
    temp =
        temp + interaction_NNNLO::V_NNNLO_np(s1, s2, configs, container_sfr_interaction_1, container_sfr_interaction_2,
                                             container_sfr_interaction_3, container_sfr_interaction_4,
                                             container_sfr_interaction_5, container_sfr_interaction_6);
    return temp;
}

} // namespace interaction_ALL

#endif // INTERACTION_ALL_HPP
