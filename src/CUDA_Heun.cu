#include "CUDA_Heun.cuh"

#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "CUDA_Toolbox.cuh"

__global__ void kernel_heunSolver(unsigned long numberOfSamples, double* outputDerivativesVector, double* coefficientsVector, double* uVector, double* yVector, unsigned long equationOrder, double samplePeriod)
{
    // Iterations variables
    double previousU;
    double k1, k2, k3;
    double thirdSamplePeriod = samplePeriod / 3.0;

    double* outputVector_k2 = (double*)malloc(equationOrder * sizeof(double)); // y vector used to compute k2
    double* outputVector_k3 = (double*)malloc(equationOrder * sizeof(double)); // y vector used to compute k3

    // Iterations
    for (unsigned long currentSampleIndex = 1; currentSampleIndex < numberOfSamples; currentSampleIndex++)
    {
        previousU = uVector[currentSampleIndex - 1];

        memcpy(outputVector_k2, outputDerivativesVector, sizeof(double) * equationOrder);
        memcpy(outputVector_k3, outputDerivativesVector, sizeof(double) * equationOrder);

        // Compute y_n = y_n-1 + Ts (k1 + 3*k3) / 4
        // k1 = f(x_n-1             , y_n-1)
        // k2 = f(x_n-1 +   Ts/3    , y_n-1 +   k1*Ts / 3)
        // k3 = f(x_n-1 + 2*Ts/3    , y_n-1 + 2*k2*Ts / 3)
        // f(x, y) = a y + u
        k1 = device_linearODEwCC(previousU, outputDerivativesVector, coefficientsVector, equationOrder); // Same as euler method

        outputVector_k2[equationOrder - 1] += k1 * thirdSamplePeriod;
        k2 = device_linearODEwCC(previousU, outputVector_k2, coefficientsVector, equationOrder);

        outputVector_k3[equationOrder - 1] += 2.0 * k2 * thirdSamplePeriod;
        k3 = device_linearODEwCC(previousU, outputVector_k3, coefficientsVector, equationOrder);

        outputDerivativesVector[equationOrder - 1] += samplePeriod * (k1 + 3.0 * k3) * 0.25;
        for (int i = equationOrder - 2; 0 <= i; i--)
        {
            k1 = outputDerivativesVector[i + 1];
            k2 = outputDerivativesVector[i + 1] + k1 * thirdSamplePeriod;
            k3 = outputDerivativesVector[i + 1] + 2.0 * k2 * thirdSamplePeriod;
            outputDerivativesVector[i] += samplePeriod * (k1 + 3.0 * k3) * 0.25;
        }    
        yVector[currentSampleIndex] = outputDerivativesVector[0];
    }

    free(outputVector_k2);
    free(outputVector_k3);
}

void heun_iterateAll_CUDA(std::vector<double>& outputDerivativesVector, const std::vector<double>& coefficientsVector, const std::vector<double>& uVector, std::vector<double>& yVector, double samplePeriod)
{
    // I won't use the fPtr for the highest derivative value in CUDA
    // I'll use a hardcoded function (linear ODE w/ CC) for the moment

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
    kernel_heunSolver<<<1, 1>>>(numberOfSamples, cuda_outputDerivativesVector, cuda_coefficientsVector, cuda_uVector, cuda_yVector, equationOrder, samplePeriod);

    // Copy the result in the CPU memory
    cudaMemcpy(yVector.data(), cuda_yVector, yVector.size() * sizeof(double), cudaMemcpyDeviceToHost);

    // FREE GPU MEMORY!!!!! (IMPORTANT)
    cudaFree(cuda_outputDerivativesVector);
    cudaFree(cuda_coefficientsVector);
    cudaFree(cuda_uVector);
    cudaFree(cuda_yVector);
}