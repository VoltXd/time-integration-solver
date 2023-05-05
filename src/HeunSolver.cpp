#include "HeunSolver.hpp"

HeunSolver::HeunSolver()
{
    fptr_highestDerivative = nullptr;
    m_isSolverReady = false;
}

// TODO: Correct this function or delete it
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

    std::vector<double> outputVector_k2(m_outputDerivativesVector.size());
    std::vector<double> outputVector_k3(m_outputDerivativesVector.size());


    // Compute y_n = y_n-1 + Ts (k1 + 3*k3) / 4
    // k1 = f(x_n-1             , y_n-1)
    // k2 = f(x_n-1 +   Ts/3    , y_n-1 +   k1*Ts / 3)
    // k3 = f(x_n-1 + 2*Ts/3    , y_n-1 + 2*k2*Ts / 3)
    // f(x, y) = a y + u
    double k1 = fptr_highestDerivative(m_timeVector.at(m_currentSampleIndex - 1), u, m_outputDerivativesVector, m_coefficientsVector); // Same as euler method

    // Fill the k2 temporaries output vectors
    for (int i = 0; i < m_outputDerivativesVector.size(); i++)
        outputVector_k2.at(i) = m_outputDerivativesVector.at(i) + k1 * thirdSamplePeriod;
    double k2 = fptr_highestDerivative(m_timeVector.at(m_currentSampleIndex - 1), u, outputVector_k2, m_coefficientsVector);

    // Fill the k2 temporaries output vectors
    for (int i = 0; i < m_outputDerivativesVector.size(); i++)
        outputVector_k3.at(i) = m_outputDerivativesVector.at(i) + 2.0 * k2 * thirdSamplePeriod;
    double k3 = fptr_highestDerivative(m_timeVector.at(m_currentSampleIndex - 1), u, outputVector_k3, m_coefficientsVector);

    m_outputDerivativesVector.back() += m_samplePeriod * (k1 + 3.0 * k3) * 0.25;
    for (int i = m_outputDerivativesVector.size() - 2; 0 <= i; i--)
    {
        k1 = m_outputDerivativesVector.at(i + 1);
        k2 = m_outputDerivativesVector.at(i + 1) + k1 * thirdSamplePeriod;
        k3 = m_outputDerivativesVector.at(i + 1) + 2.0 * k2 * thirdSamplePeriod;
        m_outputDerivativesVector.at(i) += m_samplePeriod * (k1 + 3.0 * k3) * 0.25;
    }    
    m_yVector.at(m_currentSampleIndex) = m_outputDerivativesVector.at(0);

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
    double previousU;
    double k1, k2, k3;
    double thirdSamplePeriod = m_samplePeriod / 3.0;

    unsigned long numberOfSamples = m_timeVector.size();
   
    std::vector<double> outputVector_k2(m_outputDerivativesVector.size());
    std::vector<double> outputVector_k3(m_outputDerivativesVector.size());

    // Iterations
    for (; m_currentSampleIndex < numberOfSamples; m_currentSampleIndex++)
    {
        previousU = uVector.at(m_currentSampleIndex - 1);

        outputVector_k2 = m_outputDerivativesVector;
        outputVector_k3 = m_outputDerivativesVector;

        // Compute y_n = y_n-1 + Ts (k1 + 3*k3) / 4
        // k1 = f(x_n-1             , y_n-1)
        // k2 = f(x_n-1 +   Ts/3    , y_n-1 +   k1*Ts / 3)
        // k3 = f(x_n-1 + 2*Ts/3    , y_n-1 + 2*k2*Ts / 3)
        // f(x, y) = a y + u
        k1 = fptr_highestDerivative(m_timeVector.at(m_currentSampleIndex - 1), previousU, m_outputDerivativesVector, m_coefficientsVector); // Same as euler method

        outputVector_k2.back() += k1 * thirdSamplePeriod;
        k2 = fptr_highestDerivative(m_timeVector.at(m_currentSampleIndex - 1), previousU, outputVector_k2, m_coefficientsVector);

        outputVector_k3.back() += 2.0 * k2 * thirdSamplePeriod;
        k3 = fptr_highestDerivative(m_timeVector.at(m_currentSampleIndex - 1), previousU, outputVector_k3, m_coefficientsVector);

        m_outputDerivativesVector.back() += m_samplePeriod * (k1 + 3.0 * k3) * 0.25;
        for (int i = m_outputDerivativesVector.size() - 2; 0 <= i; i--)
        {
            k1 = m_outputDerivativesVector.at(i + 1);
            k2 = m_outputDerivativesVector.at(i + 1) + k1 * thirdSamplePeriod;
            k3 = m_outputDerivativesVector.at(i + 1) + 2.0 * k2 * thirdSamplePeriod;
            m_outputDerivativesVector.at(i) += m_samplePeriod * (k1 + 3.0 * k3) * 0.25;
        }    
        m_yVector.at(m_currentSampleIndex) = m_outputDerivativesVector.at(0);
    }
}