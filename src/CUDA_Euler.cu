#include "CUDA_Euler.cuh"

#include <stdio.h>
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "CUDA_Toolbox.cuh"

__global__ void kernel_eulerSolver(unsigned long numberOfSamples, double* outputDerivativesVector, double* coefficientsVector, double* uVector, double* yVector, const unsigned long equationOrder, double samplePeriod)
{
    double highestDerivativeValue;
    unsigned long currentSampleIndex = 1;
    for (; currentSampleIndex < numberOfSamples; currentSampleIndex++)
    {
        highestDerivativeValue = device_linearODEwCC(uVector[currentSampleIndex - 1], outputDerivativesVector, coefficientsVector, equationOrder);
        outputDerivativesVector[equationOrder - 1] += samplePeriod * highestDerivativeValue;
        for (int i = equationOrder - 2; 0 <= i; i--)
        {
            outputDerivativesVector[i] += samplePeriod * outputDerivativesVector[i + 1];  
        }
        yVector[currentSampleIndex] = outputDerivativesVector[0];
    }
}

void euler_iterateAll_CUDA(std::vector<double>& outputDerivativesVector, const std::vector<double>& coefficientsVector, const std::vector<double>& uVector, std::vector<double>& yVector, double samplePeriod)
{
    // I won't use the fPtr for the highest derivative value in CUDA
    // I'll use a hardcoded function (linear ODE w/ CC) for the moment

    // The init. and the result part should be wrapped in two functions, so that this function have only 3 calls (not sure the structure is the same).

    unsigned long numberOfSamples = yVector.size();
    unsigned long equationOrder = coefficientsVector.size();

    // Allocate GPU memory
    double* cuda_outputDerivativesVector;
    double* cuda_coefficientsVector;
    double* cuda_uVector;
    double* cuda_yVector;

    cudaMalloc(&cuda_outputDerivativesVector, outputDerivativesVector.size() * sizeof(double));
    cudaMalloc(&cuda_coefficientsVector, coefficientsVector.size() * sizeof(double));
    cudaMalloc(&cuda_uVector, uVector.size() * sizeof(double));
    cudaMalloc(&cuda_yVector, yVector.size() * sizeof(double));

    // Copy datas in the GPU's memory
    cudaMemcpy(cuda_outputDerivativesVector, outputDerivativesVector.data(), outputDerivativesVector.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_coefficientsVector, coefficientsVector.data(), coefficientsVector.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_uVector, uVector.data(), uVector.size() * sizeof(double), cudaMemcpyHostToDevice);

    // Run the kernel
    kernel_eulerSolver<<<1, 1>>>(numberOfSamples, cuda_outputDerivativesVector, cuda_coefficientsVector, cuda_uVector, cuda_yVector, equationOrder, samplePeriod);

    // Copy the result in the CPU memory
    cudaMemcpy(yVector.data(), cuda_yVector, yVector.size() * sizeof(double), cudaMemcpyDeviceToHost);

    // FREE GPU MEMORY!!!!! (IMPORTANT)
    cudaFree(cuda_outputDerivativesVector);
    cudaFree(cuda_coefficientsVector);
    cudaFree(cuda_uVector);
    cudaFree(cuda_yVector);
}