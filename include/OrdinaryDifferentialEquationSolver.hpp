#pragma once

#include <iostream>
#include <vector>

#include "Toolbox.hpp"

class OrdinaryDifferentialEquationSolver
{
    public:
    void initialise(double (*highestDerivative)(double time, double input, const std::vector<double>& outputDerivativesVector, 
                    const std::vector<double>& coefficientsVector), 
                    const std::vector<double>& coefficientsVector, 
                    const std::vector<double>& initialConditionsVector, 
                    double samplePeriod, 
                    unsigned long numberOfSamples);
    virtual double iterate(double u) = 0;
    virtual void iterateAll(const std::vector<double>& uVector) = 0;
    void unitStepError(double stepValue);

    void filePrintOutputSolutionAndError(const std::string& filename, double stepValue);

    inline std::vector<double> getYVector() const { return m_yVector; } 
    inline std::vector<double> getTimeVector() const { return m_timeVector; } 

    protected:
    double (*fptr_highestDerivative)(double time, double input, const std::vector<double>& outputDerivativesVector, const std::vector<double>& coefficientsVector);

    std::vector<double> m_timeVector;  
    std::vector<double> m_yVector;

    std::vector<double> m_coefficientsVector;
    std::vector<double> m_outputDerivativesVector;

    bool m_isSolverReady;
    unsigned long m_equationOrder;
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