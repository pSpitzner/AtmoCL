#ifndef H_ASYSTEM
#define H_ASYSTEM

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include "string.h"
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
#include "clexport.h"
#include "parameters.h"

class asystem {
  public:
    parameters par;
    clcontext *context;
    cllogger *logger;
    cl_int ret;
    int debugStepCount;
    long frame_index;
    long fast_index;


    clbuffer *b_temp[4];

    clbuffer *bf_momenta_fc_a;
    clbuffer *bf_momenta_fc_b;

    // scalar field requires 2 (for ice 3) buffers for7 variables +4 for ice and snow
    clbuffer *bf_scalars_vc_a[3];
    clbuffer *bf_scalars_vc_b[3];

    clbuffer *bs_momenta_fc[3];

    clbuffer *bsRhs_momenta_fc[3];

    // 3 stages, 2 scalar buffers + 1 for ice
    // bs_sclars_vc_stageindex_scalarfieldindex
    clbuffer *bs_scalars_vc[3][3];

    // 3 stages, 2 scalar buffers + 1 for ice
    // bs_sclars_fc_direction_stageindex_scalarfieldindex
    clbuffer *bs_scalars_fc_x[3][3];
    clbuffer *bs_scalars_fc_y[3][3];
    clbuffer *bs_scalars_fc_z[3][3];

    clkernel *k_init_scalars;
    clkernel *k_init_momenta;

    clkernel *k_damping;
    clkernel *k_nesting;
    clkernel *k_perturb;
    clkernel *k_clone[4];

    clkernel *ks_ext_forcings;
    clkernel *ks_f2c;
    clkernel *ks_c2f[3];
    clkernel *ks_adv_momenta;
    clkernel *ks_adv_scalars[3][3];
    clkernel *ks_prep_fast[3][4];

    // 4 fields: density, particlenumbers, ice, momenta
    clkernel *kf_copy[4];


    // 3 stages, 4 fields: density, particlenumbers, ice, momenta
    clkernel *ks_copy[3][4];

    clkernel *kf_pressure;
    clkernel *kf_microphys;
    clkernel *kf_step_momenta[3];

    // again 3 stages, 3 fields
    clkernel *kf_step_scalars[3][3];

    // clexport *exporter[3];
    std::vector<clexport *> v_exporter;

    asystem(clcontext *contextn, cllogger *loggern, int timescheme = 1);
    ~asystem();

    void read_single_profile(clbuffer *b, FILE *f, int c, float dz, float offset = 0.0);
    void set_par(int timescheme);
    void init_from_kernel();
    void init_from_file(std::string s_filePath);

    void equilibrate();
    void perturb();
    void mis_step();
    void mis_step(int damping, int kx, int ky, int kz);
    void slow_stage(int s, int kx, int ky, int kz);
    void fast_stage(int s, int kx, int ky, int kz);

    void write_files(int step);

    void write_state(std::string s_filename);
    void read_state(std::string s_filename);
};

#endif
