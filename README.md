# NNScattering(will be updated this year)

Nucleon-Nucleon Scattering Calculator, written in C++.

--Phase shifts are calculated by solving LS equation with matrix inversion method in momentum space. Two-body states are defined by SYM-LSJ convention, under which matrix elements are calculated.

--Other two-body scattering observables can be calculated, too. See "observable.hpp"

--The evaluation of interaction matrix elements is by aPWD method.

--"Eigen.zip" should be unzipped before compiling.

--As an example, phase shifts of np phase shifts up to 300 MeV are calculated with chiral N3LO interaction (EMN500new).

--All parameters can be found or changed in "NN.ini".

--Originary bulit by xmake, see "xmake.lua".
