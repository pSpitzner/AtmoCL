#ifndef H_KERNEL
#define H_KERNEL

#include <stdio.h>
#include <vector>
#include <sys/stat.h>       // check if file is readable
#include <algorithm>
#include "clcontext.h"
#include "clbuffer.h"
#include "cllogger.h"
#include "parameters.h"

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

class clkernel
{
  public:
    cl_program program;
    cl_kernel kernel;
    cl_int ret;
    parameters par;
    std::string s_name;
    std::string s_sourceCode;
    std::vector<std::string> v_argCl;
    std::vector<clbuffer *> v_argBuffer;
    int binding_count;
    bool bindings_checked;

    clcontext *context;
    cllogger *logger;

    clkernel(clcontext *context, parameters parn, std::string filePath = "", std::string headerPathOne = "./classes/parameters.h", std::string headerPathTwo = "./kernels/k_header.cl");
    ~clkernel();
    std::string read_file(std::string inputPath);
    void read_argument_list(std::string inputPath);
    void replace_string(std::string &oldString, std::string from, std::string to);
    size_t get_pos_of_argument_from_src(std::string s_kernelArg);
    void bind(const int pos, clbuffer *b);
    void bind(std::string s_kernelArg, unsigned int myInt);
    void bind(std::string s_kernelArg, clbuffer *b);
    void bind_custom(std::string s_kernelArg, void *s, size_t custom_size);
    void check_bindings();

       // void bind(const int pos, float dt);
    void set_par(parameters parn);
    void step(int kx, int ky, int kz);

    //new
    void bind();
};

#endif
