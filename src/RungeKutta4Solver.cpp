#include "RungeKutta4Solver.hpp"

RungeKutta4Solver::RungeKutta4Solver()
{
    fptr_highestDerivative = nullptr;
    m_isSolverReady = false;
}

// TODO: Correct this function or delete it
double RungeKutta4Solver::iterate(double u)
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

    double halfSamplePeriod = m_samplePeriod * 0.5;
    double previousY = m_yVector.at(m_currentSampleIndex - 1);

    // // Compute y_n = y_n-1 + Ts (k1 + 2k2 + 2k3 + k4) / 6
    // // k1 = f(x_n-1         , y_n-1)
    // // k2 = f(x_n-1 + Ts/2  , y_n-1 + k1 Ts / 2)
    // // k3 = f(x_n-1 + Ts/2  , y_n-1 + k2 Ts / 2)
    // // k4 = f(x_n-1 + Ts    , y_n-1 + k3 Ts)
    // double k1 = m_a * previousY + u; // Same as euler method
    // double k2 = m_a * (previousY + k1 * halfSamplePeriod) + u;
    // double k3 = m_a * (previousY + k2 * halfSamplePeriod) + u;
    // double k4 = m_a * (previousY + k3 * m_samplePeriod) + u;
    
    // m_yVector.at(m_currentSampleIndex) = previousY + m_samplePeriod * (k1 + 2.0 * (k2 + k3) + k4) / 6.0;

    return m_yVector.at(m_currentSampleIndex++);
}

void RungeKutta4Solver::iterateAll(const std::vector<double>& uVector)
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
    double halfSamplePeriod = m_samplePeriod * 0.5;

    unsigned long numberOfSamples = m_timeVector.size();

    std::vector<double> outputVector_k2(m_outputDerivativesVector.size());
    std::vector<double> outputVector_k3(m_outputDerivativesVector.size());
    std::vector<double> outputVector_k4(m_outputDerivativesVector.size());


    // Iterations
    for (; m_currentSampleIndex < numberOfSamples; m_currentSampleIndex++)
    {
        previousU = uVector.at(m_currentSampleIndex - 1);

        outputVector_k2 = m_outputDerivativesVector;
        outputVector_k3 = m_outputDerivativesVector;
        outputVector_k4 = m_outputDerivativesVector;

        // Compute y_n = y_n-1 + Ts (k1 + 2k2 + 2k3 + k4) / 6, with:
        // k1 = f(x_n-1         , y_n-1)
        // k2 = f(x_n-1 + Ts/2  , y_n-1 + k1 Ts / 2)
        // k3 = f(x_n-1 + Ts/2  , y_n-1 + k2 Ts / 2)
        // k4 = f(x_n-1 + Ts    , y_n-1 + k3 Ts)
        // f(x, y) = a y + u
        k1 = fptr_highestDerivative(m_timeVector.at(m_currentSampleIndex - 1), previousU, m_outputDerivativesVector, m_coefficientsVector); // Same as euler method

        outputVector_k2.back() += k1 * halfSamplePeriod;
        k2 = fptr_highestDerivative(m_timeVector.at(m_currentSampleIndex - 1), previousU, outputVector_k2, m_coefficientsVector);

        outputVector_k3.back() += k2 * halfSamplePeriod;
        k3 = fptr_highestDerivative(m_timeVector.at(m_currentSampleIndex - 1), previousU, outputVector_k3, m_coefficientsVector);

        outputVector_k4.back() += k3 * m_samplePeriod;
        k4 = fptr_highestDerivative(m_timeVector.at(m_currentSampleIndex - 1), previousU, outputVector_k4, m_coefficientsVector);

        m_outputDerivativesVector.back() += m_samplePeriod * (k1 + 2.0 * (k2 + k3) + k4) / 6.0;
        for (int i = m_outputDerivativesVector.size() - 2; 0 <= i; i--)
        {
            k1 = m_outputDerivativesVector.at(i + 1);
            k2 = m_outputDerivativesVector.at(i + 1) + k1 * halfSamplePeriod;
            k3 = m_outputDerivativesVector.at(i + 1) + k2 * halfSamplePeriod;
            k4 = m_outputDerivativesVector.at(i + 1) + k3 * m_samplePeriod;
            m_outputDerivativesVector.at(i) += m_samplePeriod * (k1 + 2.0 * (k2 + k3) + k4) / 6.0;
        }    
        m_yVector.at(m_currentSampleIndex) = m_outputDerivativesVector.at(0);
    }
}