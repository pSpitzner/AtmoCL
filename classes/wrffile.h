#ifndef H_WRF
#define H_WRF

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include "math.h"
#include "string.h"
#include "gd.h"

#include "cllogger.h"
#include "clbuffer.h"
#include "clkernel.h"
#include "clcontext.h"
#include "parameters.h"

class wrffile {
  public:
    parameters par;
    wrfparameters wrf;
    clcontext *context;
    cllogger *logger;
    cl_int ret;

    // wrf texture variables
    float *T2,*PSFC,*U10,*V10,*Q2,*TSK,*HSFC,*RRTOT,*HFX,*LH,*TK,*PRESSURE,*UMET,*VMET,*W,*Q,*QC;

    // wrf buffer matching wrf domainsize
    clbuffer *b_wrf_source_vc[3];
    clbuffer *b_wrf_flux;
    clbuffer *b_hz;

    clkernel *k_interpolate;


    std::string file_name;



    wrffile(clcontext *contextn, cllogger *loggern, parameters parn, std::string file_namen, clbuffer *b_target_scalars[3]);
    ~wrffile();
    int index(int x, int y, int z);
    void load(int wrfindex=1);
};

#endif
