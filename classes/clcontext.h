#ifndef H_CLCONTEXT
#define H_CLCONTEXT

#include <stdio.h>
#include <sys/utsname.h>

// avoid warnings about deprecated functions - mac os only supports opencl 1.2
#define CL_USE_DEPRECATED_OPENCL_1_1_APIS
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#define CL_USE_DEPRECATED_OPENCL_2_0_APIS

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#include "cllogger.h"

class clbuffer;
class clkernel;

class clcontext
{
  public:
    cl_int ret;
    cl_device_id device_id;
    cl_context context;
    cl_command_queue command_queue;
    cllogger *logger;
    double profiling_time;
    long profiling_count;
    bool profiling;

    std::vector<clbuffer *> v_bufferList;
    std::vector<clkernel *> v_kernelList;

    clcontext(cllogger *loggern, int device_nr = 0, bool profiling_enabled = false);
    ~clcontext();
    void finish();
    void full_info();
    void cl_device_info(cl_device_id device);
    unsigned long used_memory();
};


#endif
