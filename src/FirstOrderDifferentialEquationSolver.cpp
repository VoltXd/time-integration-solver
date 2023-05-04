#include "FirstOrderDifferentialEquationSolver.hpp"

#include <fstream>

void FirstOrderDifferentialEquationSolver::initialise(double a, double y0, double samplePeriod, unsigned long numberOfSamples)
{
    // Ensure that we have a correct number of samples
    if (numberOfSamples == 0)
        return;

    m_isSolverReady = true;

    // Set the sample rate and the output coefficient.
    m_a = a;
    m_samplePeriod = samplePeriod;

    // Allocate the number of samples in the time & output vector, compute the time vector and set the initial condition.
    m_timeVector.resize(numberOfSamples);
    m_timeVector.at(0) = 0;
    for (int i = 1; i < numberOfSamples; i++)
        m_timeVector.at(i) = m_timeVector.at(i - 1) + samplePeriod;

    m_yVector.resize(numberOfSamples);
    m_yVector.at(0) = y0;

    m_currentSampleIndex = 1;
}

/// @brief Write in a file the time vector, the computed output, the exact solution and the error (in CSV format).
/// @param filename File where the datas are wrote.
/// @param stepValue Value 
/// @warning Relevant only when using step responses of the same magnitude
void FirstOrderDifferentialEquationSolver::filePrintOutputSolutionAndError(const std::string& filename, double constantInput)
{
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cout << "Cannot open " << filename << '\n';
        return;
    }

    unsigned long numberOfSamples = m_timeVector.size();
    if (numberOfSamples == 0)
        return;
    
    double time;
    double y0 = m_yVector.at(0);

    double solverY;
    double exactY;


    for (unsigned long i = 0; i < numberOfSamples; i++)
    {
        time = m_timeVector.at(i);
        solverY = m_yVector.at(i);
        exactY = FirstOrderLODEwCCStepResponse(time, m_a, y0, constantInput);
        file << time << ';' << solverY << ';' << exactY << ';' << exactY - solverY << '\n';
    }

}



void FirstOrderDifferentialEquationSolver::unitStepError(double stepValue)
{
    unsigned long vectorsSize = m_timeVector.size();

    std::vector<double> inputVector(vectorsSize, stepValue);
    iterateAll(inputVector);

    for (int i = 0; i < vectorsSize; i++)
        m_yVector.at(i) -= FirstOrderLODEwCCStepResponse(m_timeVector.at(i), m_a, m_yVector.at(0), stepValue);
}