#pragma once

#include <precision.hpp>
#ifdef CUDA
#include <cuComplex.h>
#endif

__device__ __host__ cuDoubleComplex  operator*(const cuDoubleComplex a, const cuDoubleComplex b);
__device__ __host__ cuDoubleComplex  operator+(const cuDoubleComplex a, const cuDoubleComplex b);
__device__ __host__ cuDoubleComplex  operator-(const cuDoubleComplex a, const cuDoubleComplex b);
__device__ __host__ cuDoubleComplex  operator/(const cuDoubleComplex a, const cuDoubleComplex b);

__device__ __host__ cuComplex  operator*(const cuComplex a, const cuComplex b);
__device__ __host__ cuComplex  operator+(const cuComplex a, const cuComplex b);
__device__ __host__ cuComplex  operator-(const cuComplex a, const cuComplex b);
__device__ __host__ cuComplex  operator/(const cuComplex a, const cuComplex b);

__device__ __host__ gpu_mode& operator+=(gpu_mode& a, const gpu_mode b);
__device__ __host__ gpu_mode& operator-=(gpu_mode& a, const gpu_mode b);


// Generic complex number creation
__device__ __host__ gpu_mode makeComplex(real a, real b);

// Float-complex multiplication
__device__ __host__ gpu_mode  operator*(const real a, const gpu_mode b);

__device__ __host__ gpu_mode operator*(const gpu_mode a, const real b);
