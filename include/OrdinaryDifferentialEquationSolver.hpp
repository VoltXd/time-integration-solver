#pragma once

#include <iostream>
#include <vector>

#include "Toolbox.hpp"

class OrdinaryDifferentialEquationSolver
{
    public:
    void initialise(ODE* highestDerivative, 
                    const std::vector<double>& coefficientsVector, 
                    const std::vector<double>& initialConditionsVector, 
                    double samplePeriod, 
                    unsigned long numberOfSamples);
    virtual double iterate(double u) = 0;
    virtual void iterateAll(const std::vector<double>& uVector) = 0;

    void filePrintOutputSolutionAndError(const std::string& filename, double stepValue);

    inline std::vector<double> getYVector() const { return m_yVector; } 
    inline std::vector<double> getTimeVector() const { return m_timeVector; } 

    protected:
    ODE* fptr_highestDerivative;

    std::vector<double> m_timeVector;  
    std::vector<double> m_yVector;

    /// @brief (a_0 a_1 ... a_N-1)^T. Example => d^N(y)/dt^N = a_N-1 d^N-1(y)/dt^N-1 + ... + a_0 y + u
    std::vector<double> m_coefficientsVector;
    /// @brief (y(0) y'(0) ... d^N-1(y)/dt^N-1(0))^T
    std::vector<double> m_outputDerivativesVector;

    bool m_isSolverReady;
    double m_samplePeriod;
    unsigned long m_currentSampleIndex;



    friend inline std::ostream& operator<<(std::ostream& stream, const OrdinaryDifferentialEquationSolver& solver);
};

inline std::ostream& operator<<(std::ostream& stream, const OrdinaryDifferentialEquationSolver& solver)
{
    unsigned long numberOfSamples = solver.m_timeVector.size();
    for (unsigned long i = 0; i < numberOfSamples; i++)
    {
        stream << solver.m_timeVector.at(i) << ';' << solver.m_yVector.at(i) << '\n';
    }

    return stream;
}