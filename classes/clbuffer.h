#ifndef H_BUFFER
#define H_BUFFER

#include "stdio.h"
#include "math.h"
#include "clcontext.h"
#include "cllogger.h"


#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

class clbuffer
{
  public:
    clcontext *context;
    cl_mem buffer;
    float *data;
    size_t datasize;
    cl_int ret;
    cllogger *logger;
    std::string s_name;
    //buffer size - not euqal to system size eg. par.sx
    int sx, sy, sz;


    clbuffer(clcontext *context_, std::string name_, int sx_, int sy_, int sz_);
    ~clbuffer();
    void set(int x, int y, int z, int c, float v);
    float get(int x, int y, int z, int c);
    void debug();
    void ram2device();
    void device2ram();

    void read_raw(std::string s_filename);
    void write_raw(std::string s_filename);
    void test_reading();

  private:
    int getindex(int x, int y, int z);
};

#endif
