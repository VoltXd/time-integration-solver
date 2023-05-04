#pragma once

#include <cmath>

/// @brief Compute the step response of a first order linear ordinary differential equation with constant coefficients (y' = a*y + u)
/// @param time The time when we want the response
/// @param a The equation's constant coefficient
/// @param y0 The initial value of the output
/// @param constantInput The input value
/// @return y(t), the solution of the equation
inline double FirstOrderLODEwCCStepResponse(double time, double a, double y0, double constantInput)
{
    return (y0 + constantInput / a) * exp(a * time) - constantInput / a;
}
