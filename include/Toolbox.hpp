#pragma once

#include <cmath>
#include <vector>

/// @brief Compute the step response of a first order linear ordinary differential equation with constant coefficients (y' = a*y + u)
/// @param time The time when we want the response
/// @param a The equation's constant coefficient
/// @param y0 The initial value of the output
/// @param constantInput The input value
/// @return y(t), the solution of the equation
inline double firstOrderLODEwCCStepResponse(double time, double a, double y0, double constantInput)
{
    return (y0 + constantInput / a) * exp(a * time) - constantInput / a;
}

/// @brief Compute the step response of a second order linear ordinary differential equation with constant coefficients (y'' = a1*y' + a0*y + u).
/// @brief Only for pseudo-oscillating systems
/// @param time The time when we want the response
/// @param a_i The equation's constant coefficient
/// @param y0 The initial value of the output
/// @param constantInput The input value
/// @return y(t), the solution of the equation
inline double secondOrderLODEwCCStepResponse(double time, double a0, double a1, double constantInput)
{
    double A = constantInput / a0;
    double B = - a1 * constantInput / a0 / sqrt(-a1*a1 - 4*a0);
    double omega = sqrt(-a1*a1 - 4*a0) * 0.5;
    return (A * cos(omega * time) + B * sin(omega * time)) * exp(0.5 * a1 * time) - constantInput / a0; 
}

double linearODEwCC(double time, double input, const std::vector<double>& outputDerivativesVector, const std::vector<double>& coefficientVector);
