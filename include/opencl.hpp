#ifndef __OPENCL_UTIL__
#define __OPENCL_UTIL__

#ifdef GPU_ENABLED

#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_TARGET_OPENCL_VERSION 200
#include <CL/opencl.hpp>

#endif

#endif