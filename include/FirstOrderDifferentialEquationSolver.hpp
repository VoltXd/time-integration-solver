#pragma once

#include <iostream>
#include <vector>

#include "Toolbox.hpp"

class FirstOrderDifferentialEquationSolver
{
    public:
    void initialise(double a, double y0, double samplePeriod, unsigned long numberOfSamples);
    virtual double iterate(double u) = 0;
    virtual void iterateAll(const std::vector<double>& uVector) = 0;
    void unitStepError(double stepValue);

    void filePrintOutputSolutionAndError(const std::string& filename, double stepValue);

    protected:
    std::vector<double> m_timeVector;  
    std::vector<double> m_yVector;
    
    double m_a;
    double m_samplePeriod;

    bool m_isSolverReady;

    unsigned long m_currentSampleIndex;

    friend inline std::ostream& operator<<(std::ostream& stream, const FirstOrderDifferentialEquationSolver& solver);
};

inline std::ostream& operator<<(std::ostream& stream, const FirstOrderDifferentialEquationSolver& solver)
{
    unsigned long numberOfSamples = solver.m_timeVector.size();
    for (unsigned long i = 0; i < numberOfSamples; i++)
    {
        stream << solver.m_timeVector.at(i) << ';' << solver.m_yVector.at(i) << '\n';
    }

    return stream;
}