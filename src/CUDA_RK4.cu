#include "CUDA_RK4.cuh"

#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "CUDA_Toolbox.cuh"

__global__ void kernel_rk4Solver(unsigned long numberOfSamples, double* outputDerivativesVector, double* coefficientsVector, double* uVector, double* yVector, unsigned long equationOrder, double samplePeriod)
{
    // Iterations variables
    double previousY, previousU;
    double k1, k2, k3, k4;
    double halfSamplePeriod = samplePeriod * 0.5;

    double* outputVector_k2 = (double*)malloc(equationOrder * sizeof(double));
    double* outputVector_k3 = (double*)malloc(equationOrder * sizeof(double));
    double* outputVector_k4 = (double*)malloc(equationOrder * sizeof(double));


    // Iterations
    for (unsigned long currentSampleIndex = 1; currentSampleIndex < numberOfSamples; currentSampleIndex++)
    {
        previousU = uVector[currentSampleIndex - 1];

        memcpy(outputVector_k2, outputDerivativesVector, equationOrder * sizeof(double));
        memcpy(outputVector_k3, outputDerivativesVector, equationOrder * sizeof(double));
        memcpy(outputVector_k4, outputDerivativesVector, equationOrder * sizeof(double));

        // Compute y_n = y_n-1 + Ts (k1 + 2k2 + 2k3 + k4) / 6, with:
        // k1 = f(x_n-1         , y_n-1)
        // k2 = f(x_n-1 + Ts/2  , y_n-1 + k1 Ts / 2)
        // k3 = f(x_n-1 + Ts/2  , y_n-1 + k2 Ts / 2)
        // k4 = f(x_n-1 + Ts    , y_n-1 + k3 Ts)
        // f(x, y) = a y + u
        k1 = device_linearODEwCC(previousU, outputDerivativesVector, coefficientsVector, equationOrder); // Same as euler method

        outputVector_k2[equationOrder - 1] += k1 * halfSamplePeriod;
        k2 = device_linearODEwCC(previousU, outputVector_k2, coefficientsVector, equationOrder);

        outputVector_k3[equationOrder - 1] += k2 * halfSamplePeriod;
        k3 = device_linearODEwCC(previousU, outputVector_k3, coefficientsVector, equationOrder);

        outputVector_k4[equationOrder - 1] += k3 * samplePeriod;
        k4 = device_linearODEwCC(previousU, outputVector_k4, coefficientsVector, equationOrder);

        outputDerivativesVector[equationOrder - 1] += samplePeriod * (k1 + 2.0 * (k2 + k3) + k4) / 6.0;
        for (int i = equationOrder - 2; 0 <= i; i--)
        {
            k1 = outputDerivativesVector[i + 1];
            k2 = outputDerivativesVector[i + 1] + k1 * halfSamplePeriod;
            k3 = outputDerivativesVector[i + 1] + k2 * halfSamplePeriod;
            k4 = outputDerivativesVector[i + 1] + k3 * samplePeriod;
            outputDerivativesVector[i] += samplePeriod * (k1 + 2.0 * (k2 + k3) + k4) / 6.0;
        }    
        yVector[currentSampleIndex] = outputDerivativesVector[0];
    }

    free(outputVector_k2);
    free(outputVector_k3);
    free(outputVector_k4);
}

void rk4_iterateAll_CUDA(std::vector<double>& outputDerivativesVector, const std::vector<double>& coefficientsVector, const std::vector<double>& uVector, std::vector<double>& yVector, double samplePeriod)
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
    kernel_rk4Solver<<<1, 1>>>(numberOfSamples, cuda_outputDerivativesVector, cuda_coefficientsVector, cuda_uVector, cuda_yVector, equationOrder, samplePeriod);

    // Copy the result in the CPU memory
    cudaMemcpy(yVector.data(), cuda_yVector, yVector.size() * sizeof(double), cudaMemcpyDeviceToHost);

    // FREE GPU MEMORY!!!!! (IMPORTANT)
    cudaFree(cuda_outputDerivativesVector);
    cudaFree(cuda_coefficientsVector);
    cudaFree(cuda_uVector);
    cudaFree(cuda_yVector);
}