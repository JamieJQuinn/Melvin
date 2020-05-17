#include <complex_gpu.hpp>
#include <precision.hpp>

__device__ __host__ cuDoubleComplex  operator*(const cuDoubleComplex a, const cuDoubleComplex b) { return cuCmul(a,b); }
__device__ __host__ cuDoubleComplex  operator+(const cuDoubleComplex a, const cuDoubleComplex b) { return cuCadd(a,b); }
__device__ __host__ cuDoubleComplex  operator-(const cuDoubleComplex a, const cuDoubleComplex b) { return cuCsub(a,b); }
__device__ __host__ cuDoubleComplex  operator/(const cuDoubleComplex a, const cuDoubleComplex b) { return cuCdiv(a,b); }

__device__ __host__ cuComplex  operator*(const cuComplex a, const cuComplex b) { return cuCmulf(a,b); }
__device__ __host__ cuComplex  operator+(const cuComplex a, const cuComplex b) { return cuCaddf(a,b); }
__device__ __host__ cuComplex  operator-(const cuComplex a, const cuComplex b) { return cuCsubf(a,b); }
__device__ __host__ cuComplex  operator/(const cuComplex a, const cuComplex b) { return cuCdivf(a,b); }

__device__ __host__ gpu_mode& operator+=(gpu_mode& a, const gpu_mode b) {
  a = a+b;
  return a;
}

__device__ __host__ gpu_mode& operator-=(gpu_mode& a, const gpu_mode b) {
  a = a-b;
  return a;
}

__device__ __host__ gpu_mode makeComplex(real a, real b) {
#if defined PRECISION_DOUBLE
  return make_cuDoubleComplex(a, b);
#elif defined PRECISION_SINGLE
  return make_cuComplex(a, b);
#endif
}

__device__ __host__ gpu_mode  operator*(const real a, const gpu_mode b) {
#if defined PRECISION_DOUBLE
  return cuCmul(makeComplex(a, 0.0f), b);
#elif defined PRECISION_SINGLE
  return cuCmulf(makeComplex(a, 0.0f), b);
#endif
}

__device__ __host__ gpu_mode  operator*(const gpu_mode a, const real b) { return b*a; }


