#ifndef H_EXPORT
#define H_EXPORT

#include <stdio.h>
#include <stdlib.h>
#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision
#include <sstream>      // std::stringstream, std::stringbuf
#include <vector>
#include "gd.h"

#ifdef __APPLE__
  #include <OpenCL/opencl.h>
#else
  #include <CL/cl.h>
#endif

#include "cllogger.h"
#include "clbuffer.h"
#include "clcontext.h"
#include "clkernel.h"
#include "parameters.h"

class clexport {
  public:
    std::string s_name, s_kname;
    cllogger *logger;
    clcontext *context;
    parameters par;

    // cpu based averaging for time series and vertical profiles
    std::string s_ts_path;
    std::string s_img_path;
    std::string s_vp_path;
    std::ofstream ts_stream;
    std::ofstream vp_stream;

    // what files to write
    bool img; // mode[0]
    bool ts;  // mode[1]
    bool vp;  // mode[2]

    // when to write files in seconds of simulation time (call export every mis step!)
    int t_img;
    int t_ts;
    int t_vp;
    float img_every;
    float ts_every;
    float vp_every;


    // cut plane and offset
    int dim, ref;

    // kernel 2d size
    int kx, ky, kz;

    // images to be exported for desired variables
    std::vector<clbuffer *> be_export_images;

    // set once, when kernel is compiled. for different fields to be exported, create new exporter.
    clbuffer *b_source_scalars[3];
    clbuffer *b_source_momenta;

    clkernel *ke_render;

    clexport(clcontext *context_, std::string name_, std::string kname_, parameters par_, clbuffer *mom, clbuffer *sc0, clbuffer *sc1, clbuffer *sc2, int d_, int r_, bool img_, bool ts_, bool vp_);
    ~clexport();

    void write_files(int it);

};

#endif
