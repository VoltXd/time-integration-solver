#pragma once

#include "FirstOrderDifferentialEquationSolver.hpp"

/// @brief Solver for 1st order differential equations using Euler Method.
/// Solving y' = a*y + u, knowing y(0), and where u is the input and y the output.
class EulerSolver : public FirstOrderDifferentialEquationSolver
{
    public:
    EulerSolver();
    
    double iterate(double u) override;
    void iterateAll(const std::vector<double>& uVector) override;
};