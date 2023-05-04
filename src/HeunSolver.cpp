#include "HeunSolver.hpp"

HeunSolver::HeunSolver()
{
    fptr_highestDerivative = nullptr;
    m_isSolverReady = false;
}

double HeunSolver::iterate(double u)
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

    double thirdSamplePeriod = m_samplePeriod / 3.0;
    double previousY = m_yVector.at(m_currentSampleIndex - 1);

    // Compute y_n = y_n-1 + Ts (k1 + 3*k3) / 4
    // k1 = f(x_n-1             , y_n-1)
    // k2 = f(x_n-1 +   Ts/3    , y_n-1 +   k1*Ts / 3)
    // k3 = f(x_n-1 + 2*Ts/3    , y_n-1 + 2*k2*Ts / 3)
    // f(x, y) = a y + u
    double k1 = m_a * previousY + u; // Same as euler method
    double k2 = m_a * (previousY + k1 * thirdSamplePeriod) + u;
    double k3 = m_a * (previousY + 2.0 * k2 * thirdSamplePeriod) + u;
    
    m_yVector.at(m_currentSampleIndex) = previousY + m_samplePeriod * (k1 + 3.0 * k3) * 0.25;

    return m_yVector.at(m_currentSampleIndex++);
}

void HeunSolver::iterateAll(const std::vector<double>& uVector)
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

    // Iterations variables
    double previousY, previousU;
    double k1, k2, k3, k4;
    double thirdSamplePeriod = m_samplePeriod / 3.0;

    unsigned long numberOfSamples = m_timeVector.size();

    // Iterations
    for (; m_currentSampleIndex < numberOfSamples; m_currentSampleIndex++)
    {
        previousY = m_yVector.at(m_currentSampleIndex - 1);
        previousU = uVector.at(m_currentSampleIndex - 1);

        
        // Compute y_n = y_n-1 + Ts (k1 + 3*k3) / 4
        // k1 = f(x_n-1             , y_n-1)
        // k2 = f(x_n-1 +   Ts/3    , y_n-1 +   k1*Ts / 3)
        // k3 = f(x_n-1 + 2*Ts/3    , y_n-1 + 2*k2*Ts / 3)
        // f(x, y) = a y + u
        k1 = m_a * previousY + previousU; // Same as euler method
        k2 = m_a * (previousY + k1 * thirdSamplePeriod) + previousU;
        k3 = m_a * (previousY + 2.0 * k2 * thirdSamplePeriod) + previousU;
        
        m_yVector.at(m_currentSampleIndex) = previousY + m_samplePeriod * (k1 + 3.0 * k3) * 0.25;
    }
}