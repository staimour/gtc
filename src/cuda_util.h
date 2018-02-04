/* This file is part of GTC version 3 
! GTC version 3 is released under the 3-Clause BSD license:

! Copyright (c) 2002,2010,2016, GTC Team (team leader: Zhihong Lin, zhihongl@uci.edu)
! All rights reserved.

! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the following conditions are met:

! 1. Redistributions of source code must retain the above copyright notice, 
!    this list of conditions and the following disclaimer.

! 2. Redistributions in binary form must reproduce the above copyright notice, 
!    this list of conditions and the following disclaimer in the documentation 
!    and/or other materials provided with the distribution.

! 3. Neither the name of the GTC Team nor the names of its contributors may be 
!    used to endorse or promote products derived from this software without 
!    specific prior written permission. */

#ifndef CUDA_UTIL_H
#define CUDA_UTIL_H

#define CUDA_SAFE_CALL(call) {                                    \
  cudaError err = call;                                                    \
  if( cudaSuccess != err) {                                                \
  fprintf(stderr, "Cuda error in file '%s' in line %i : %d, %s.\n",        \
          __FILE__, __LINE__, err, cudaGetErrorString( err) );           \
  fflush(stderr); \
  exit(EXIT_FAILURE);                                                  \
  } }

template <class T>
class gpu_array
{
public:
  gpu_array()
  {
    size_ = 0;
    ptr_ = NULL;
  }
  gpu_array(size_t size)
  {
    size_ = size;
    CUDA_SAFE_CALL( cudaMalloc((void**)&ptr_, size*sizeof(T)) );
  }
  gpu_array(T* data, size_t size)
  {
    size_ = size;
    CUDA_SAFE_CALL( cudaMalloc((void**)&ptr_, size*sizeof(T)) );
    CUDA_SAFE_CALL( cudaMemcpy(ptr_, data, size*sizeof(T), cudaMemcpyHostToDevice) );
  }
  ~gpu_array()
  {
    if (ptr_ != NULL) {
      CUDA_SAFE_CALL( cudaFree(ptr_) );
    }
  }
  void malloc(size_t size)
  {
    size_ = size;
    CUDA_SAFE_CALL( cudaMalloc((void**)&ptr_, size*sizeof(T)) );
  }
  void free()
  {
    if (ptr_ != NULL) {
      CUDA_SAFE_CALL( cudaFree(ptr_) );
      ptr_ = NULL;
    }
  }
  void D2H(T* data)
  {
    CUDA_SAFE_CALL( cudaMemcpy(data, ptr_, size_*sizeof(T), cudaMemcpyDeviceToHost) );
  }
  void D2H(T* data, size_t size)
  {
    CUDA_SAFE_CALL( cudaMemcpy(data, ptr_, size*sizeof(T), cudaMemcpyDeviceToHost) );
  }
  void H2D(T *data)
  {
    CUDA_SAFE_CALL( cudaMemcpy(ptr_, data, size_*sizeof(T), cudaMemcpyHostToDevice) );
  }
  void H2D(T *data, size_t size)
  {
    CUDA_SAFE_CALL( cudaMemcpy(ptr_, data, size*sizeof(T), cudaMemcpyHostToDevice) );
  }
  void D2D(T *data, size_t size)
  {
    CUDA_SAFE_CALL( cudaMemcpy(ptr_, data, size*sizeof(T), cudaMemcpyDeviceToDevice) );
  }
  void D2HAsync(T* data, cudaStream_t stream) 
  {
    CUDA_SAFE_CALL( cudaMemcpyAsync(data, ptr_, size_*sizeof(T), cudaMemcpyDeviceToHost, stream) );
  }
  void D2HAsync(T* data, size_t size, cudaStream_t stream) 
  {
    CUDA_SAFE_CALL( cudaMemcpyAsync(data, ptr_, size*sizeof(T), cudaMemcpyDeviceToHost, stream) );
  }
  void H2DAsync(T *data, cudaStream_t stream) 
  {
    CUDA_SAFE_CALL( cudaMemcpyAsync(ptr_, data, size_*sizeof(T), cudaMemcpyHostToDevice, stream) );
  }
  void H2DAsync(T *data, size_t size, cudaStream_t stream) 
  {
    CUDA_SAFE_CALL( cudaMemcpyAsync(ptr_, data, size*sizeof(T), cudaMemcpyHostToDevice, stream) );
  }
  T* ptr() { return ptr_; }
  size_t size() { return size_; }
private:
  T *ptr_;
  size_t size_;
};

#include <sys/time.h>

class CTimer
{
protected:
  timeval start,end;
public:
  void Start() {gettimeofday(&start, NULL);}
  double GetET() {
    gettimeofday(&end,NULL);
    double et=(end.tv_sec+end.tv_usec*0.000001)-(start.tv_sec+start.tv_usec*0.000001);
    return et;
  }
};


#endif
