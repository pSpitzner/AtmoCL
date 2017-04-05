#ifndef H_ERROR
#define H_ERROR

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision
#include <sstream>      // std::stringstream, std::stringbuf
#include <fstream>
#include <ctime>
#include <chrono>
#include <cmath>
#include <vector>
#include <memory>       // std::unique_ptr
#include "assert.h"

#ifdef __APPLE__
  #include <OpenCL/opencl.h>
#else
  #include <CL/cl.h>
#endif

class cllogger {
  public:

    // loglevel:
    // 0 screen and log file
    // 1 log file only
    // 2 screen only

    cllogger(int logleveln = 0, bool profiling=false, std::string path = "./output/logs/",std::string filenamen = "");
    ~cllogger();
    int loglevel;
    std::string Log_path;
    std::string filename;
    std::ofstream Debug_stream;
    std::ofstream Profiling_stream;
    std::ofstream ProfAvg_stream;
    bool profiling;
    FILE *LogFile;

    // Error and Warning logs
    void log(int msglvl, cl_int err, const char* comment, ...);
    void log(int msglvl, const char* comment, ...);
    std::string clGetErrorString(cl_int err);

    // export of additional log files
    void open_debugstream(std::string s_debug, int step);
    void close_debugstream();

    // profiling
    template<typename ... Args>
    std::string string_printf(const std::string& format, Args ... args);
    void plog(long long ns, std::string kernel_name = "", std::string note = "");
    void avg_plog();

    std::vector<std::string> v_kernel_name;
    std::vector<std::string> v_note;
    std::vector<long long>   v_time;

    // want the following to survive outside of avg_plog() so we can clear v_kernel_name at the end
    // to avoid repeating the iteration through the whole log history
    std::vector<std::string> v_name_slow;
    std::vector<std::string> v_name_fast;
    std::vector<std::string> v_name_other;
    std::vector<long long> v_time_slow;
    std::vector<long long> v_time_fast;
    std::vector<long long> v_time_other;
    std::vector<long long> v_totl_slow;
    std::vector<long long> v_totl_fast;
    std::vector<long long> v_totl_other;
    std::vector<long long> v_entr_slow;
    std::vector<long long> v_entr_fast;
    std::vector<long long> v_entr_other;
    long long total_time;

    // some timing
    int t_measured;
    double t_avg;
    typedef std::chrono::high_resolution_clock clock_cr;
    std::vector<clock_cr::time_point> timers;
    void start_new_timer();
    long long end_last_timer();
    void log_estimated_time(double t_last_step, int i, int iterations);
};

#endif
