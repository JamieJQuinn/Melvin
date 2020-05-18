#pragma once
#include <cfloat>
#include <complex>

#ifdef CUDA
#include <cuComplex.h>
#endif

#define PRECISION_DOUBLE
// NOTE SINGLE PRECISION NOT SUPPORTED YET
//#define PRECISION_SINGLE

#if defined PRECISION_DOUBLE
typedef double real;
const real EPSILON = FLT_EPSILON;
#ifdef CUDA
typedef cuDoubleComplex gpu_mode;
#endif
#elif defined PRECISION_SINGLE
typedef float real;
const real EPSILON = DBL_EPSILON;
#ifdef CUDA
typedef cuComplex gpu_mode;
#endif
#endif

typedef std::complex<real> mode;
