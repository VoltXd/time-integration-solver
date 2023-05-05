#include "Toolbox.hpp"

double linearODEwCC(double time, double input, const std::vector<double>& outputDerivativesVector, const std::vector<double>& coefficientVector)
{
    double result = input;
    for (int i = 0; i < outputDerivativesVector.size(); i++)
        result += coefficientVector.at(i) * outputDerivativesVector.at(i);
    return result;
}
