#pragma once

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=false)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %d %s %s %d\n", code, cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}
