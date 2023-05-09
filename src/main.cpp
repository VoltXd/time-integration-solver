#include <iostream>
#include <fstream>
#include <vector>

#include "EulerSolver.hpp"
#include "RungeKutta4Solver.hpp"
#include "Timer.hpp"
#include "HeunSolver.hpp"
#include "Toolbox.hpp" 

int main(int argc, char* argv[])
{
    double RC = 10;

    // Input vector
    unsigned long numberOfSamples = 1024;
    std::vector<double> inputVector(numberOfSamples, RC * RC);
    std::vector<double> coefficientsVector { -RC * RC, -RC };
    std::vector<double> initialConditionsVector { 0, 0 };

    // Euler Solver
    std::cout << "Euler Method\n";
    EulerSolver eulerSolver;
    eulerSolver.initialise(linearODEwCC, coefficientsVector, initialConditionsVector, 1e-3, numberOfSamples);
    {
        Timer t;
        eulerSolver.iterateAll_CUDA(inputVector);
    }
    eulerSolver.filePrintOutputSolutionAndError("results/Euler.csv", RC * RC);

    // Heun Solver 
    std::cout << "Heun Method\n";
    HeunSolver heunSolver;
    heunSolver.initialise(linearODEwCC, coefficientsVector, initialConditionsVector, 1e-3, numberOfSamples);
    {
        Timer t;
        heunSolver.iterateAll_CUDA(inputVector);
    }
    heunSolver.filePrintOutputSolutionAndError("results/Heun.csv", RC * RC);

    // 4-stage Runge-Kutta Solver 
    std::cout << "RK4 Method\n";
    RungeKutta4Solver rk4Solver;
    rk4Solver.initialise(linearODEwCC, coefficientsVector, initialConditionsVector, 1e-3, numberOfSamples);
    {
        Timer t;
        rk4Solver.iterateAll_CUDA(inputVector);
    }
    rk4Solver.filePrintOutputSolutionAndError("results/RK4.csv", RC * RC);

    return 0;
}