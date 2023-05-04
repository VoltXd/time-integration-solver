#pragma once 

#include "OrdinaryDifferentialEquationSolver.hpp"


/// @brief Solver for 1st order differential equations using Heun Method.
/// Solving y' = a*y + u, knowing y(0), and where u is the input and y the output.
class HeunSolver : public OrdinaryDifferentialEquationSolver
{
    public:
    HeunSolver();
    
    double iterate(double u) override;
    void iterateAll(const std::vector<double>& uVector) override;
};