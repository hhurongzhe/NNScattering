#pragma once
#ifndef LO_STATE_HPP
#define LO_STATE_HPP

#include "lib_define.hpp"
#include "mymath.hpp"

namespace states
{

struct LSJ_State
{
    int l; // total orbital momentum of a LSJ basis
    int s; // total spin of a LSJ basis
    int j; // j=l+s
    int m; // j in z axis
    // int p;           //parity factor \eta of a LSJ basis, +1 or -1, i.e P = \eta (-1)^J
    double momentum; // momentum of a LSJ basis
    friend bool operator==(const LSJ_State &s1, const LSJ_State &s2)
    {
        return s1.j == s2.j && s1.l == s2.l && s1.s == s2.s;
    }

    // initialize a LSJ state
    static states::LSJ_State ini_LSJ_state(int l, int s, int j, int m, double pp)
    {
        states::LSJ_State state;
        state.l = l;
        state.s = s;
        state.j = j;
        state.m = m;
        state.momentum = pp;
        return state;
    }
};

} // end namespace states

#endif // LO_STATE_HPP