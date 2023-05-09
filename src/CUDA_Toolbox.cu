#include "CUDA_Toolbox.cuh"

__device__ double device_linearODEwCC(double input, const double* outputDerivativesVector, const double* coefficientVector, unsigned long equationOrder)
{
    double result = input;
    for (int i = 0; i < equationOrder; i++)
        result += coefficientVector[i] * outputDerivativesVector[i];
    return result;
}

