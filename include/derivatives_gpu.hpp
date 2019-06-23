#pragma once

#include <variable_gpu.hpp>
#include <constants.hpp>

void computeLinearTemperatureDerivativeGPU(VariableGPU &dTmpdt, const VariableGPU &tmp, const Constants &c);
void computeLinearVorticityDerivativeGPU(VariableGPU &dOmgdt, const VariableGPU &omg, const VariableGPU &tmp, const Constants &c);
void addAdvectionApproximationGPU(
    VariableGPU &dTmpdt, const VariableGPU &tmp,
    VariableGPU &dOmgdt, const VariableGPU &omg,
    VariableGPU &dXidt, const VariableGPU &xi,
    const VariableGPU &psi,
    const Constants &c);
