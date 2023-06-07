#pragma once
#ifndef PRECALCULATE_HPP
#define PRECALCULATE_HPP

#include "interaction_ALL.hpp"
#include "lib_define.hpp"

namespace precalculate
{

struct Interaction_Container
{
  private:
    std::unordered_map<unsigned long long int, double> container;
    unsigned long long int get_key(int i, int j, int L1, int L2, int S, int J)
    {
        // compute key using bitwise XOR.(more efficient than addition or multiplication in general)
        unsigned long long int key = ((unsigned long long int)(i)) ^ (((unsigned long long int)(j)) << 12) ^
                                     (((unsigned long long int)(L1)) << 24) ^ (((unsigned long long int)(L2)) << 36) ^
                                     (((unsigned long long int)(S)) << 48) ^ (((unsigned long long int)(J)) << 60);
        return key;
    }

  public:
    Interaction_Container() {}
    double get_value(const NN::NN_configs &configs, int i, int j, int L1, int L2, int S, int J,
                     std::vector<double> &container_sfr_interaction_1, std::vector<double> &container_sfr_interaction_2,
                     std::vector<double> &container_sfr_interaction_3, std::vector<double> &container_sfr_interaction_4,
                     std::vector<double> &container_sfr_interaction_5, std::vector<double> &container_sfr_interaction_6)
    {
        auto key = get_key(i, j, L1, L2, S, J);
        auto it = container.find(key);
        if (it != container.end())
        {
            return it->second;
        }
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(L1, S, J, J, configs.momentum_mesh_points[i]);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(L2, S, J, J, configs.momentum_mesh_points[j]);
        double value = interaction_ALL::V_np(
            ss1, ss2, configs, container_sfr_interaction_1, container_sfr_interaction_2, container_sfr_interaction_3,
            container_sfr_interaction_4, container_sfr_interaction_5, container_sfr_interaction_6);
        container[key] = value;
        return value;
    }
};

} // namespace precalculate

#endif // PRECALCULATE_HPP