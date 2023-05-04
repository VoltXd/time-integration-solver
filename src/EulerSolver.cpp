#include "EulerSolver.hpp"
#include "Toolbox.hpp"

#include <iostream>
#include <limits>

EulerSolver::EulerSolver()
{
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

    double yPrime = m_a * m_yVector.at(m_currentSampleIndex - 1) + u;
    m_yVector.at(m_currentSampleIndex) = m_samplePeriod * yPrime + m_yVector.at(m_currentSampleIndex - 1);

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
    double yPrime;
    for (; m_currentSampleIndex < numberOfSamples; m_currentSampleIndex++)
    {
        yPrime = m_a * m_yVector.at(m_currentSampleIndex - 1) + uVector.at(m_currentSampleIndex - 1);
        m_yVector.at(m_currentSampleIndex) = m_samplePeriod * yPrime + m_yVector.at(m_currentSampleIndex - 1);
    }
}