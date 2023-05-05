#include "OrdinaryDifferentialEquationSolver.hpp"

#include <fstream>

void OrdinaryDifferentialEquationSolver::initialise(double (*highestDerivative)(double time, double input, const std::vector<double>& outputDerivativesVector, const std::vector<double>& coefficientsVector),
                                                    const std::vector<double>& coefficientsVector,
                                                    const std::vector<double>& initialConditionsVector, 
                                                    double samplePeriod, 
                                                    unsigned long numberOfSamples)
{
    // Ensure that we have a correct number of samples
    if (numberOfSamples == 0)
    {
        std::cout << "Solver init. error: numberOfSamples == 0\n\n";
        return;
    }

    // Ensure that there as much coefficients as initial conditions
    if (coefficientsVector.size() != initialConditionsVector.size())
    {
        std::cout << "Solver init. error: coefficient and initial conditions vectors have different size !\n\n";
        return;
    }

    // Ensure the function pointer passed in the arguments is not a null pointer
    if (highestDerivative == nullptr)
    {
        std::cout << "Solver init. error: function pointer passed for highest derivative function is a null pointer !\n\n";
        return;
    } 

    m_isSolverReady = true;

    // Retrieve the highest derivative function pointer
    fptr_highestDerivative = highestDerivative;

    // Set the sample rate and copy the coefficient vector.
    m_coefficientsVector = coefficientsVector;
    m_samplePeriod = samplePeriod;

    // Allocate the number of samples in the time & output vector and compute the time vector
    m_timeVector.resize(numberOfSamples);
    m_yVector.resize(numberOfSamples);

    m_timeVector.at(0) = 0;
    for (int i = 1; i < numberOfSamples; i++)
        m_timeVector.at(i) = m_timeVector.at(i - 1) + samplePeriod;

    
    // Set the initial conditions.
    m_yVector.at(0) = initialConditionsVector.at(0);
    m_outputDerivativesVector = initialConditionsVector;

    // Init. iteration index
    m_currentSampleIndex = 1;
}


/// @brief Write in a file the time vector, the computed output, the exact solution and the error (in CSV format).
/// @param filename File where the datas are wrote.
/// @param stepValue Value 
/// @warning Relevant only when using step responses of the same magnitude
void OrdinaryDifferentialEquationSolver::filePrintOutputSolutionAndError(const std::string& filename, double constantInput)
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
        exactY = secondOrderLODEwCCStepResponse(time, m_coefficientsVector.at(0), m_coefficientsVector.at(1), constantInput);
        // exactY = firstOrderLODEwCCStepResponse(time, m_coefficientsVector.at(0), 0, constantInput);
        file << time << ';' << solverY << ';' << exactY << ';' << exactY - solverY << '\n';
    }

}
