#pragma once

#include <variable.hpp>
#include <constants.hpp>

void computeLinearTemperatureDerivativeGPU(Variable &dTmpdt, const Variable &tmp, const Constants &c);
void computeLinearVorticityDerivativeGPU(Variable &dOmgdt, const Variable &omg, const Variable &tmp, const Constants &c);
void addAdvectionApproximationGPU(
    Variable &dTmpdt, const Variable &tmp,
    Variable &dOmgdt, const Variable &omg,
    Variable &dXidt, const Variable &xi,
    const Variable &psi,
    const Constants &c);
