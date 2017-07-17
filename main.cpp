#include <stdio.h>
#include <stdlib.h>
#include <iostream>     // std::cout, std::endl
#include <iomanip>      // std::setw
#include <sstream>      // std::stringstream
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
#include "asystem.h"


int main(int argc, char * argv[]) {

  cl_device_type DeviceType = CL_DEVICE_TYPE_CPU;
  bool profiling_enabled = false;
  // std::string s_profilePath = "./profiles/data125";
  std::string s_profilePath = "./profiles/dycoms_2";
  std::string s_statePath = "";
  int timescheme = 1;
  int device = 0;

  for (int i = 1; i < argc; i++) {
    if (std::string(argv[i]) == "--device" || std::string(argv[i]) == "-dev") {
      device = (int)atof(argv[i+1]);
    } else if (std::string(argv[i]) == "-profiling") {
      profiling_enabled = true;
    } else if (std::string(argv[i]) == "-p") {
      s_profilePath = (argv[i+1]);
    } else if (std::string(argv[i]) == "-b") {
      s_statePath = (argv[i+1]);
    } else if (std::string(argv[i]) == "-rk3") {
      timescheme = 0;
    } else if (std::string(argv[i]) == "-mis") {
      timescheme = 1;
    }
  }

  clcontext *context;
  cllogger *logger;
  asystem *sys;


  logger  = new cllogger(1, profiling_enabled);
  context = new clcontext(logger, device, profiling_enabled);
  sys     = new asystem(context, logger, timescheme);
  long long ns;

  // context->info();

  std::stringstream FileName;

  context->used_memory();

  sys->init_from_wrf("test");
  // return 0;
  /*
  if (s_statePath.length() == 0) {
    logger->log(0, "No state files given\n");
      if (s_profilePath.length() == 0) {
        logger->log(0, "No profile provided\n");
        return -1;
      } else {
        sys->init_from_file(s_profilePath);
      }
    logger->log(0, "Equilibrating...");
    sys->equilibrate();
    sys->write_state("./snapshots/equil");
  } else {
    logger->log(0, "Reading state files from %s\n", s_statePath.c_str());
    sys->read_state(s_statePath);
  }
  */

  // sys->write_files(0);
  // return 0;

  logger->start_new_timer();
  logger->log(2, "Estimating remaining time...\n");


  int iterations = 50000;
  for (int i = 0; i < iterations; i++) {
    logger->start_new_timer();
    sys->mis_step();

    long long ns = logger->end_last_timer();
    double seconds = double(ns)/1.0E9;
    logger->log_estimated_time(seconds, i, iterations);
    logger->log(2,"  -  %d/%d",i,iterations);
    if (profiling_enabled) logger->avg_plog();
  }

  double seconds_taken = double(logger->end_last_timer())/1.0E9;
  int hours_taken = int(seconds_taken/3600.0);
  seconds_taken = fmod(seconds_taken,3600.0);
  int minutes_taken = int(seconds_taken/60.0);
  seconds_taken = fmod(seconds_taken,60.0);
  logger->log(0, "\ndone after %02dh %02dm %02ds\n",hours_taken,minutes_taken,int(seconds_taken));


  if (profiling_enabled) {
    context->full_info();
    logger->log(0, "\tcl_time: %f\n", context->profiling_time * 1.0e-9);
  }

  delete sys;
  delete context;
  delete logger;
  return 0;
}

