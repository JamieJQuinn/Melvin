#pragma once

#include <precision.hpp>

void require_within_error(const real x, const real y, const real margin=1e-3);
void require_within_error(const mode x, const mode y, const real margin=1e-3);
void check_equal(const real x, const real y);
void require_equal(const real x, const real y);
void check_equal(const mode x, const mode y);
void require_equal(const mode x, const mode y);
