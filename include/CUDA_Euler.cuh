#pragma once

#include <vector>

void euler_iterateAll_CUDA(std::vector<double>& outputDerivativesVector, const std::vector<double>& coefficientsVector, const std::vector<double>& uVector, std::vector<double>& yVector, double samplePeriod);