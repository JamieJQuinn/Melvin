#pragma once

#include <variable.hpp>

void computeLinearTemperatureDerivative(Variable &dTmpdt, const Variable &tmp, const Constants &c);
void computeLinearVorticityDerivative(Variable &dOmgdt, const Variable &omg, const Variable &tmp, const Constants &c);
void computeLinearXiDerivative(Variable &dXidt, const Variable &xi, Variable &dOmgdt, const Constants &c);
void computeLinearDerivatives(
    Variable &dTmpdt, const Variable &tmp,
    Variable &dOmgdt, const Variable &omg,
    Variable &dXidt, const Variable &xi,
    const Constants &c);
void addAdvectionApproximation(
    Variable &dTmpdt, const Variable &tmp,
    Variable &dOmgdt, const Variable &omg,
    Variable &dXidt, const Variable &xi,
    const Variable &psi,
    const Constants &c);
void computeNonlinearTemperatureDerivative(
    Variable &dTmpdt, const Variable &tmp,
    const Variable &psi,
    const Constants &c);
void computeNonlinearXiDerivative(
    Variable &dXidt, const Variable &xi,
    const Variable &psi,
    const Constants &c);
void computeNonlinearVorticityDerivative(
    Variable &dOmgdt, const Variable &omg,
    const Variable &psi,
    const Constants &c);
void computeNonlinearDerivatives(
    Variable &dTmpdt, const Variable &tmp,
    Variable &dOmgdt, const Variable &omg,
    Variable &dXidt, const Variable &xi,
    const Variable &psi,
    const Constants &c);
