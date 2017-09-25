#include "cllogger.h"

cllogger::cllogger(int logleveln, bool profilingn, std::string s_outputn) {
  t_measured = 0;
  t_avg = 0.0;
  loglevel = logleveln;
  profiling = profilingn;
  // 0: only log lvl 0 infos to screen
  // 1: log everything to file, lvl 0 infos to screen
  if (loglevel == 0) {
    return;
  }

  char mbstr[100];
  char datestr[100];
  std::time_t now = std::time(NULL);
  std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d_%H-%M-%S", std::localtime(&now));
  std::strftime(datestr, sizeof(datestr), "%Y-%m-%d", std::localtime(&now));

  if (s_outputn == "") { // set default
    s_output = "./output/" + std::string(datestr) + "/";
  } else {
    std::string::size_type found = s_outputn.find_last_of("/");
    if (found == s_outputn.length()-1) s_output = s_outputn;
    else s_output = s_outputn + "/";
  }
  s_logPath = s_output + "logs/";

  printf("Cleaning old files...\n");
  // clean old files if present. deleting directories may be slow on network drive, delete files first
  int ret;
  ret = system(("find "+s_output+" -name '*.png' -type f -delete 2>/dev/null").c_str());
  ret = system(("find "+s_output+" -name '*.vp' -type f -delete 2>/dev/null").c_str());
  ret = system(("rm -r "+s_output+"img/ 2>/dev/null").c_str());
  ret = system(("rm -r "+s_output+"timeseries/ 2>/dev/null").c_str());
  ret = system(("rm -r "+s_output+"verticalprofiles/ 2>/dev/null").c_str());




  ret = system(("mkdir -p " + s_logPath).c_str());
  std::string filename = s_logPath + "AtmoCL.log";

  if (profiling) {
    Profiling_stream.open(s_logPath + "Profiling.log", std::ofstream::out);
    int mysize = 5000*100;
    v_kernel_name.reserve(mysize);
    v_note.reserve(mysize);
    v_time.reserve(mysize);
  }

  LogFile = fopen(filename.c_str(), "w");
  // fseek(LogFile, 0, SEEK_SET); // overwrite
  // ofs.open(path+filename.c_str(), std::ofstream::out);

  log(0, "------------------ %s ------------------\n", mbstr);
  log(0, "Writing files to %s\n\n", s_output.c_str());

  // profiling
  if (profiling) Profiling_stream << "------------------ " << std::string(mbstr) << " ------------------\n";
}

cllogger::~cllogger() {
  if (LogFile != NULL) {
    // ofs.close();
    std::time_t now = std::time(NULL);
    char mbstr[100];
    std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d_%H-%M-%S", std::localtime(&now));

    log(1,"------------------ %s ------------------\n", mbstr);
    fflush(LogFile);
    fclose(LogFile);
  }

  if (profiling) {
    Profiling_stream.close();
    Profiling_stream.clear();
  }
}

std::string cllogger::clGetErrorString(cl_int err)
{
  switch (err) {
    case CL_SUCCESS:                            return "Success!                             ";
    case CL_DEVICE_NOT_FOUND:                   return "ERROR: Device not found              ";
    case CL_DEVICE_NOT_AVAILABLE:               return "ERROR: Device not available          ";
    case CL_COMPILER_NOT_AVAILABLE:             return "ERROR: Compiler not available        ";
    case CL_MEM_OBJECT_ALLOCATION_FAILURE:      return "ERROR: Memory obj alloc failed       ";
    case CL_OUT_OF_RESOURCES:                   return "ERROR: Out of resource               ";
    case CL_OUT_OF_HOST_MEMORY:                 return "ERROR: Out of host memory            ";
    case CL_PROFILING_INFO_NOT_AVAILABLE:       return "ERROR: Profiling info n.a.           ";
    case CL_MEM_COPY_OVERLAP:                   return "ERROR: Memory copy overlap           ";
    case CL_IMAGE_FORMAT_MISMATCH:              return "ERROR: Image format mismatch         ";
    case CL_IMAGE_FORMAT_NOT_SUPPORTED:         return "ERROR: Image format not suppted      ";
    case CL_BUILD_PROGRAM_FAILURE:              return "ERROR: Program build failure         ";
    case CL_MAP_FAILURE:                        return "ERROR: Map failure                   ";
    case CL_INVALID_VALUE:                      return "ERROR: Invalid value                 ";
    case CL_INVALID_DEVICE_TYPE:                return "ERROR: Invalid device type           ";
    case CL_INVALID_PLATFORM:                   return "ERROR: Invalid platform              ";
    case CL_INVALID_DEVICE:                     return "ERROR: Invalid device                ";
    case CL_INVALID_CONTEXT:                    return "ERROR: Invalid context               ";
    case CL_INVALID_QUEUE_PROPERTIES:           return "ERROR: Invalid queue properts        ";
    case CL_INVALID_COMMAND_QUEUE:              return "ERROR: Invalid command queue         ";
    case CL_INVALID_HOST_PTR:                   return "ERROR: Invalid host pointer          ";
    case CL_INVALID_MEM_OBJECT:                 return "ERROR: Invalid memory object         ";
    case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:    return "ERROR: Invalid img format description";
    case CL_INVALID_IMAGE_SIZE:                 return "ERROR: Invalid img size              ";
    case CL_INVALID_SAMPLER:                    return "ERROR: Invalid sampler               ";
    case CL_INVALID_BINARY:                     return "ERROR: Invalid binary                ";
    case CL_INVALID_BUILD_OPTIONS:              return "ERROR: Invalid build options         ";
    case CL_INVALID_PROGRAM:                    return "ERROR: Invalid program               ";
    case CL_INVALID_PROGRAM_EXECUTABLE:         return "ERROR: Invalid program executable    ";
    case CL_INVALID_KERNEL_NAME:                return "ERROR: Invalid kernel name           ";
    case CL_INVALID_KERNEL_DEFINITION:          return "ERROR: Invalid kernel definition     ";
    case CL_INVALID_KERNEL:                     return "ERROR: Invalid kernel                ";
    case CL_INVALID_ARG_INDEX:                  return "ERROR: Invalid argument index        ";
    case CL_INVALID_ARG_VALUE:                  return "ERROR: Invalid argument value        ";
    case CL_INVALID_ARG_SIZE:                   return "ERROR: Invalid argument size         ";
    case CL_INVALID_KERNEL_ARGS:                return "ERROR: Invalid kernel args           ";
    case CL_INVALID_WORK_DIMENSION:             return "ERROR: Invalid work dimension        ";
    case CL_INVALID_WORK_GROUP_SIZE:            return "ERROR: Invalid work group size       ";
    case CL_INVALID_WORK_ITEM_SIZE:             return "ERROR: Invalid work item size        ";
    case CL_INVALID_GLOBAL_OFFSET:              return "ERROR: Invalid global offset         ";
    case CL_INVALID_EVENT_WAIT_LIST:            return "ERROR: Invalid event wait list       ";
    case CL_INVALID_EVENT:                      return "ERROR: Invalid event                 ";
    case CL_INVALID_OPERATION:                  return "ERROR: Invalid operation             ";
    case CL_INVALID_GL_OBJECT:                  return "ERROR: Invalid OpenGL object         ";
    case CL_INVALID_BUFFER_SIZE:                return "ERROR: Invalid buffer size           ";
    default:                                    return "ERROR: Unknown                       ";
  }
}



void cllogger::log(int msglvl, const char* comment, ...)
{

  if (msglvl == 0 || msglvl == 2) {
    va_list args;
    va_start(args, comment);
    vfprintf(stderr, comment, args);
    va_end(args);
    fflush(stderr);
    fflush(stdout);
    if (msglvl == 2) {
      // dont write to logfile
      return;
    }
  }
  if (LogFile != NULL) {
    va_list args;
    va_start(args, comment);
    // fprintf(LogFile, "%d: ", msglvl);
    vfprintf(LogFile, comment, args);
    va_end(args);
    fflush(LogFile);
  }


}

void cllogger::log(int msglvl, cl_int err, const char* comment, ... )
{

  if (msglvl == 0 || msglvl == 2) {
    va_list args;
    va_start(args, comment);
    fprintf(stderr, "%s\t", clGetErrorString(err).c_str());
    vfprintf(stderr, comment, args);
    va_end(args);
    fflush(stderr);
    fflush(stdout);
    if (msglvl == 2) {
      return;
    }
  }
  if (LogFile != NULL) {
    va_list args;
    va_start(args, comment);
    // fprintf(LogFile, "%d: %s \t", msglvl, clGetErrorString(err).c_str());
    fprintf(LogFile, "%s\t", clGetErrorString(err).c_str());
    vfprintf(LogFile, comment, args);
    va_end(args);
    fflush(LogFile);
    // if Error is encountered, always print to screen
    if (err != CL_SUCCESS) {
      va_list args;
      va_start(args, comment);
      fprintf(stderr, "%s\t", clGetErrorString(err).c_str());
      vfprintf(stderr, comment, args);
      va_end(args);
      fflush(LogFile);
    }
  }
}

void cllogger::open_debugstream(std::string s_debug, int step) {
  std::stringstream Temp;
  Temp << s_logPath << "Debug_" << s_debug << "_" << std::setw(5) << std::setfill('0') << step << ".log";
  Debug_stream.open(Temp.str(), std::ofstream::out); //needs c++11 oterwise try .c_str()
}

void cllogger::close_debugstream() {
  assert(Debug_stream.is_open());
  Debug_stream.flush();
  Debug_stream.close();
  Debug_stream.clear();
}

void cllogger::start_new_timer() {
  clock_cr::time_point now = clock_cr::now();
  timers.push_back(now);
}

long long cllogger::end_last_timer() {
  assert(timers.size() > 0);
  clock_cr::time_point start = timers.back();
  clock_cr::time_point end = clock_cr::now();
  timers.pop_back();
  std::chrono::duration<long long, std::nano> time_span = end-start;
  // std::chrono::duration<double, std::milli> time_span = end-start;
  // return time_span.count()/1000.0;
  return time_span.count();
}

// printf into string
template<typename ... Args>
std::string cllogger::string_printf(const std::string& format, Args ... args) {
  size_t size = 1 + std::snprintf(nullptr, 0, format.c_str(), args ...);
  std::unique_ptr<char[]> buf(new char[size]);
  std::snprintf(buf.get(), size, format.c_str(), args ...);
  // return std::string(buf.get(), buf.get() + size);
  // -1 avoids null-termination, which prevents appanding more strings with + operator etc.
  return std::string(buf.get(), buf.get() + size - 1);
}

// template<typename ... Args>
// void cllogger::plog(long long ns, const std::string& note, Args ... args) {
void cllogger::plog(long long ns, std::string kernel_name, std::string note) {
  // pushing into vector for grouped export
  v_kernel_name.push_back(kernel_name);
  v_note.push_back(note);
  v_time.push_back(ns);

  // writing to file
  long long s, ms, mus;
  s   = ns/1e9;
  ns  = fmod(ns,1e9);
  ms  = ns/1e6;
  ns  = fmod(ns,1e6);
  mus = ns/1e3;
  ns  = fmod(ns,1e3);

  std::string out;
  if (s==0) {
    if (ms==0) {
      if (mus==0) {
        out=string_printf("%15lld",ns);
      } else out=string_printf("%11lld %03lld",mus,ns);
    } else out=string_printf("%7lld %03lld %03lld",ms,mus,ns);
  } else out=string_printf("%3lld %03lld %03lld %03lld",s,ms,mus,ns);

  Profiling_stream << out;
  Profiling_stream << " ns    " << std::setw(25) << std::left << kernel_name << note << std::endl;

  Profiling_stream.flush();
}

void cllogger::avg_plog() {

  if (v_kernel_name.size() == 0) return;

  long long slow_time = 0;
  long long fast_time = 0;
  long long other_time = 0;
  long long fast_entries = 0;
  long long slow_entries = 0;
  long long other_entries = 0;

  // read history and sort into categories
  for (int i=0; i < v_kernel_name.size(); i++) {
    std::string myname = v_kernel_name.at(i);
    long long mytime = v_time.at(i);
    total_time += mytime;
    if (myname.find("kf_") != std::string::npos) {
      int pos = -1;
      for (int j = 0; j < v_name_fast.size(); j++) {
        if (myname.compare(v_name_fast.at(j)) == 0) {
          pos = j;
          break;
        }
      }
      if (pos == -1) {
        v_name_fast.push_back(myname);
        v_time_fast.push_back(mytime);
        v_totl_fast.push_back(mytime);
        v_entr_fast.push_back(1);
      } else {
        v_time_fast.at(pos) = (v_entr_fast.at(pos)*v_time_fast.at(pos) + mytime)/(v_entr_fast.at(pos) + 1);
        v_totl_fast.at(pos) += mytime;
        v_entr_fast.at(pos) += 1;
      }
    }
    else if (myname.find("ks_") != std::string::npos) {
      int pos = -1;
      for (int j = 0; j < v_name_slow.size(); j++) {
        if (myname.compare(v_name_slow.at(j)) == 0) {
          pos = j;
          break;
        }
      }
      if (pos == -1) {
        v_name_slow.push_back(myname); // names
        v_time_slow.push_back(mytime); // average time per name
        v_totl_slow.push_back(mytime); // total time per name
        v_entr_slow.push_back(1);
      } else {
        v_time_slow.at(pos) = (v_entr_slow.at(pos)*v_time_slow.at(pos) + mytime)/(v_entr_slow.at(pos) + 1);
        v_totl_slow.at(pos) += mytime;
        v_entr_slow.at(pos) += 1;
      }
    }
    else {
      int pos = -1;
      for (int j = 0; j < v_name_other.size(); j++) {
        if (v_name_other.at(j).compare(myname) == 0) {
          pos = j;
          break;
        }
      }
      if (pos == -1) {
        v_name_other.push_back(myname);
        v_time_other.push_back(mytime);
        v_totl_other.push_back(mytime);
        v_entr_other.push_back(1);
        fflush(stdout);
      } else {
        v_time_other.at(pos) = (v_entr_other.at(pos)*v_time_other.at(pos) + mytime)/(v_entr_other.at(pos) + 1);
        v_totl_other.at(pos) += mytime;
        v_entr_other.at(pos) += 1;
      }
    }
  }

  if (v_time_fast.size() > 0) {
    for (int i = 0; i < v_totl_fast.size(); i++) {
      fast_time += v_totl_fast.at(i);
    }

    ProfAvg_stream.open(s_logPath + "Profiling_Average_fast.log", std::ofstream::out | std::ofstream::trunc);
    ProfAvg_stream << std::setw(25) << std::left << "#process";
    ProfAvg_stream << std::setw(20) << std::right << "total_time";
    ProfAvg_stream << std::setw(20) << std::right << "average_time";
    ProfAvg_stream << std::setw(15) << std::right << "repeat";
    ProfAvg_stream << std::setw(10) << std::right << "%_fast";
    ProfAvg_stream << std::setw(10) << std::right << "%_total";
    ProfAvg_stream << std::endl;
    for (int i = 0; i < v_name_fast.size(); i++) {
      ProfAvg_stream << std::setw(25) << std::left << v_name_fast.at(i);
      ProfAvg_stream << std::setw(20) << std::right << v_totl_fast.at(i);
      ProfAvg_stream << std::setw(20) << std::right << v_time_fast.at(i);
      ProfAvg_stream << std::setw(15) << std::right << v_entr_fast.at(i);
      ProfAvg_stream << std::fixed << std::setprecision(1) << std::setw(10) << std::right << 100.0*double(v_totl_fast.at(i))/double(fast_time);
      ProfAvg_stream << std::fixed << std::setprecision(1) << std::setw(10) << std::right << 100.0*double(v_totl_fast.at(i))/double(total_time);
      ProfAvg_stream << std::endl;
    }
    ProfAvg_stream.close();
    ProfAvg_stream.clear();
  }

  if (v_time_slow.size() > 0) {
    for (int i = 0; i < v_totl_slow.size(); i++) {
      slow_time += v_totl_slow.at(i);
    }

    ProfAvg_stream.open(s_logPath + "Profiling_Average_slow.log", std::ofstream::out | std::ofstream::trunc);
    ProfAvg_stream << std::setw(25) << std::left << "#process";
    ProfAvg_stream << std::setw(20) << std::right << "total_time";
    ProfAvg_stream << std::setw(20) << std::right << "average_time";
    ProfAvg_stream << std::setw(15) << std::right << "repeat";
    ProfAvg_stream << std::setw(10) << std::right << "%_slow";
    ProfAvg_stream << std::setw(10) << std::right << "%_total";
    ProfAvg_stream << std::endl;
    for (int i = 0; i < v_name_slow.size(); i++) {
      ProfAvg_stream << std::setw(25) << std::left << v_name_slow.at(i);
      ProfAvg_stream << std::setw(20) << std::right << v_totl_slow.at(i);
      ProfAvg_stream << std::setw(20) << std::right << v_time_slow.at(i);
      ProfAvg_stream << std::setw(15) << std::right << v_entr_slow.at(i);
      ProfAvg_stream << std::fixed << std::setprecision(1) << std::setw(10) << std::right << 100.0*double(v_totl_slow.at(i))/double(slow_time);
      ProfAvg_stream << std::fixed << std::setprecision(1) << std::setw(10) << std::right << 100.0*double(v_totl_slow.at(i))/double(total_time);
      ProfAvg_stream << std::endl;
    }
    ProfAvg_stream.close();
    ProfAvg_stream.clear();
  }

  if (v_time_other.size() > 0) {
    for (int i = 0; i < v_totl_other.size(); i++) {
      other_time += v_totl_other.at(i);
    }

    ProfAvg_stream.open(s_logPath + "Profiling_Average_other.log", std::ofstream::out | std::ofstream::trunc);
    ProfAvg_stream << std::setw(25) << std::left << "#process";
    ProfAvg_stream << std::setw(20) << std::right << "total_time";
    ProfAvg_stream << std::setw(20) << std::right << "average_time";
    ProfAvg_stream << std::setw(15) << std::right << "repeat";
    ProfAvg_stream << std::setw(10) << std::right << "%_other";
    ProfAvg_stream << std::setw(10) << std::right << "%_total";
    ProfAvg_stream << std::endl;
    for (int i = 0; i < v_name_other.size(); i++) {
      ProfAvg_stream << std::setw(25) << std::left << v_name_other.at(i);
      ProfAvg_stream << std::setw(20) << std::right << v_totl_other.at(i);
      ProfAvg_stream << std::setw(20) << std::right << v_time_other.at(i);
      ProfAvg_stream << std::setw(15) << std::right << v_entr_other.at(i);
      ProfAvg_stream << std::fixed << std::setprecision(1) << std::setw(10) << std::right << 100.0*double(v_totl_other.at(i))/double(other_time);
      ProfAvg_stream << std::fixed << std::setprecision(1) << std::setw(10) << std::right << 100.0*double(v_totl_other.at(i))/double(total_time);
      ProfAvg_stream << std::endl;
    }
    ProfAvg_stream.close();
    ProfAvg_stream.clear();
  }

  // this might break the whole thing, (vectors not copying correctly)
  v_kernel_name.clear();
  v_note.clear();
  v_time.clear();

}

void cllogger::log_estimated_time(double t_last_step, int i, int iterations) {
  t_avg = (i*t_avg + t_last_step)/double(i+1);
  double seconds_left = double(iterations-i)*t_avg;
  double seconds_used = t_last_step;
  int hours_left = int(seconds_left/3600.0);
  int hours_used = int(seconds_used/3600.0);
  seconds_left = fmod(seconds_left,3600.0);
  seconds_used = fmod(seconds_used,3600.0);
  int minutes_left = int(seconds_left/60.0);
  int minutes_used = int(seconds_used/60.0);
  seconds_left = fmod(seconds_left,60.0);
  seconds_used = fmod(seconds_used,60.0);
  log(2, "\rEstimated time  -  last step: %02dh %02dm %02ds  -  remaining: %02dh %02dm %02ds",hours_used,minutes_used,int(seconds_used),hours_left,minutes_left,int(seconds_left));
}
