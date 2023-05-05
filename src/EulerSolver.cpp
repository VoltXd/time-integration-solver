#include "EulerSolver.hpp"
#include "Toolbox.hpp"

#include <iostream>
#include <limits>

EulerSolver::EulerSolver()
{
    fptr_highestDerivative = nullptr;
    m_isSolverReady = false;
}

double EulerSolver::iterate(double u)
{
    // Iterations are possible only if the solver is properly set up.
    if (!m_isSolverReady)
    {
        std::cout << "The solver is not ready...\n";
        std::cout << "Number of samples allocated: " << m_yVector.size() << "\n\n";
        return std::numeric_limits<double>::infinity();
    }

    // Cannot iterate if the sample index is out of bounds
    if (m_currentSampleIndex >= m_yVector.size())
    {
        std::cout << "You are trying to run the solver but every iterations are computed...\n";
        std::cout << "Current index: " << m_currentSampleIndex << '\n';
        std::cout << "Number of samples allocated: " << m_yVector.size() << "\n\n";
        return std::numeric_limits<double>::infinity();
    }

    // Compute the highest derivative value
    double highestDerivativeValue = fptr_highestDerivative(m_timeVector.at(m_currentSampleIndex - 1), u, m_outputDerivativesVector, m_coefficientsVector);
    m_outputDerivativesVector.back() += m_samplePeriod * highestDerivativeValue;
    for (int i = m_outputDerivativesVector.size() - 2; 0 <= i; i--)
        m_outputDerivativesVector.at(i) += m_samplePeriod * m_outputDerivativesVector.at(i + 1);  

    m_yVector.at(m_currentSampleIndex) = m_outputDerivativesVector.at(0);

    return m_yVector.at(m_currentSampleIndex++);
}

void EulerSolver::iterateAll(const std::vector<double>& uVector)
{
    // Iterations are possible only if the solver is properly set up.
    if (!m_isSolverReady)
    {
        std::cout << "The solver is not ready...\n";
        std::cout << "Number of samples allocated: " << m_yVector.size() << "\n\n";
        return;
    }

    // Cannot iterate if the sample index is out of bounds
    if (m_currentSampleIndex >= m_yVector.size())
    {
        std::cout << "You are trying to run the solver but every iterations are computed...\n";
        std::cout << "Current index: " << m_currentSampleIndex << '\n';
        std::cout << "Number of samples allocated: " << m_yVector.size() << "\n\n";
        return;
    }

    unsigned long numberOfSamples = m_yVector.size();
    double highestDerivativeValue;
    for (; m_currentSampleIndex < numberOfSamples; m_currentSampleIndex++)
    {
        highestDerivativeValue = fptr_highestDerivative(m_timeVector.at(m_currentSampleIndex - 1), uVector.at(m_currentSampleIndex - 1), m_outputDerivativesVector, m_coefficientsVector);
        m_outputDerivativesVector.back() += m_samplePeriod * highestDerivativeValue;
        for (int i = m_outputDerivativesVector.size() - 2; 0 <= i; i--)
        {
            m_outputDerivativesVector.at(i) += m_samplePeriod * m_outputDerivativesVector.at(i + 1);  
        }
        m_yVector.at(m_currentSampleIndex) = m_outputDerivativesVector.at(0);
    }
}