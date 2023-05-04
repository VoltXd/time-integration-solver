#pragma once

#include "FirstOrderDifferentialEquationSolver.hpp"

/// @brief Solver for 1st order differential equations using 4-Stage Runge-Kutta Method.
/// Solving y' = a*y + u, knowing y(0), and where u is the input and y the output.
class RungeKutta4Solver : public FirstOrderDifferentialEquationSolver
{
    public:
    RungeKutta4Solver();
    
    double iterate(double u) override;
    void iterateAll(const std::vector<double>& uVector) override;
};