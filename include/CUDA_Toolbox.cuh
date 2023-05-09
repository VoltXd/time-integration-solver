#pragma once

__device__ double device_linearODEwCC(double input, const double* outputDerivativesVector, const double* coefficientVector, unsigned long equationOrder);