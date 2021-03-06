#pragma once

#include <precision.hpp>

void copyGPUConstants(
    int nG, int nX, int nN, int nZ,
    real oodz, real oodx, real oodz2,
    real aspectRatio, real wavelength,
    mode xSinDerivativeFactor, mode xCosDerivativeFactor,
    real Ra, real Pr, real RaXi, real tau
  );
