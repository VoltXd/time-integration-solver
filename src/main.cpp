#include <iostream>
#include <fstream>
#include <vector>

#include "EulerSolver.hpp"
#include "RungeKutta4Solver.hpp"
#include "HeunSolver.hpp"
#include "Timer.hpp"

int main(int argc, char* argv[])
{
    double RC = 10;

    // Input vector
    unsigned long numberOfSamples = 1024;
    std::vector<double> inputVector(numberOfSamples, RC);

    // Euler Solver
    std::cout << "Euler Method\n";
    EulerSolver eulerSolver;
    eulerSolver.initialise(-RC, 0, 1e-3, numberOfSamples);
    {
        Timer t;
        eulerSolver.iterateAll(inputVector);
    }
    eulerSolver.filePrintOutputSolutionAndError("Euler.csv", RC);

    // 4-stage Runge-Kutta Solver 
    std::cout << "RK4 Method\n";
    RungeKutta4Solver rk4Solver;
    rk4Solver.initialise(-RC, 0, 1e-3, numberOfSamples);
    {
        Timer t;
        rk4Solver.iterateAll(inputVector);
    }
    rk4Solver.filePrintOutputSolutionAndError("RK4.csv", RC);

    // Heun Solver 
    std::cout << "Heun Method\n";
    HeunSolver heunSolver;
    heunSolver.initialise(-RC, 0, 1e-3, numberOfSamples);
    {
        Timer t;
        heunSolver.iterateAll(inputVector);
    }
    heunSolver.filePrintOutputSolutionAndError("Heun.csv", RC);

    return 0;
}