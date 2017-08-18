#include "clkernel.h"


clkernel::clkernel(clcontext *context, parameters parn, std::string filePath, std::string headerPathOne, std::string headerPathTwo):
  context(context), logger(context->logger) {

  s_name = filePath;
  replace_string(s_name, "./kernels/exporter/", "");
  replace_string(s_name, "./kernels/", "");
  replace_string(s_name, ".cl", "");

  context->v_kernelList.push_back(this);
  s_sourceCode = read_file(headerPathOne) + read_file(headerPathTwo) + read_file(filePath);
  // if (s_name.compare("ks_adv_scalars") == 0) printf("\n\n%s\n\n", s_sourceCode.c_str());

  bool print_log = false;

  const char* cStrSource[1] = { s_sourceCode.data() };

  program = clCreateProgramWithSource(context->context, 1, cStrSource, NULL, &ret);
  logger->log(1, ret, "%s\tCreating kernel program\n", s_name.c_str());
  ret = clBuildProgram(program, 1, &(context->device_id), NULL, NULL, NULL);
  if (ret != CL_SUCCESS) print_log = true;
  // logger->log(1, ret, "%s\tBuilding kernel program\n", s_name.c_str());
  //   if (print_log) logger->log(1, "source code:\n\n%s\n", s_sourceCode.c_str());

  char buildlog[320768];
  ret = clGetProgramBuildInfo(program,  context->device_id, CL_PROGRAM_BUILD_LOG, sizeof(buildlog), buildlog, NULL);
  std::string logString = buildlog;

  if (logString.size()!=0) {
    replace_string(logString, "\n", "\n\t");
    logString.insert(0,"\t");
  }


  if (print_log) logger->log(0, ret, "%s\tBuild log for kernel:\n%s\n", s_name.c_str(), logString.c_str());
  else logger->log(1, ret, "%s\tBuild log for kernel:\n%s\n", s_name.c_str(), logString.c_str());

  /* Create kernel */
  std::string kernel_function = s_name + "_kernel_main";
  kernel = clCreateKernel(program, kernel_function.c_str(), &ret);
  logger->log(1, ret, "%s\tCreating kernel from program\n",s_name.c_str());

  set_par(parn);
  binding_count = 0;
  bindings_checked = false;
}

clkernel::~clkernel() {
  ret = clReleaseKernel(kernel);
  ret = clReleaseProgram(program);
}

void clkernel::replace_string(std::string &oldString, std::string from, std::string to) {
  // some string manipulation
  // printf("replacing from: %s \n", from.c_str());
  std::string newString;
  newString.reserve(oldString.length());  // avoids a few memory allocations
  std::string::size_type lastPos = 0;
  std::string::size_type findPos;
  while (std::string::npos != ( findPos = oldString.find(from, lastPos))) {
    newString.append(oldString, lastPos, findPos - lastPos);
    newString += to;
    lastPos = findPos + from.length();
  }
  newString += oldString.substr(lastPos);
  oldString.swap(newString);
  // if (oldString.compare(newString) != 0) {
  //   printf("changed to: %s\n",oldString.c_str());
  // }
}

std::string clkernel::read_file(std::string importPath) {
  std::ifstream inputStream (importPath);
  std::stringstream streamBuffer;

  if (inputStream.is_open()) {
    streamBuffer << inputStream.rdbuf();
    inputStream.close();
  } else {
    logger->log(0, "Unable to read: %s\n", importPath.c_str());
  }
  std::string myString = streamBuffer.str();
  return myString;
}

void clkernel::read_argument_list(std::string importPath) {
  std::ifstream is_file (importPath);
  std::string s_line;

  if (v_argCl.size() != 0 || v_argBuffer.size() != 0) {
    logger->log(0, "List of Arguments of %s is not empty.\n Are you sure you want to read them?");
  }

  // printf("reading: %s with\n",s_name.c_str());
  if (is_file.is_open()) {
    while (std::getline(is_file, s_line)) {
      std::stringstream ss_line(s_line);
      std::string s_readWrite, s_argType, s_argName, s_buffer;
      if (ss_line >> s_readWrite >> s_argType >> s_argName >> s_buffer) {
        v_argCl.push_back(s_readWrite + " " + s_argType + " " + s_argName);
        for (int i = 0; i<context->v_bufferList.size(); i++) {
          clbuffer *b = context->v_bufferList[i];
          if (s_buffer.compare(b->s_name) == 0) v_argBuffer.push_back(b);
        }
        // printf("%s\t%s\n",v_argCl.back().c_str(), v_argBuffer.back()->s_name.c_str());
      }
    }
  } else {
    logger->log(0, "Unable to read: %s\n", importPath.c_str());
  }

  if (v_argCl.size() != v_argBuffer.size()) {
    logger->log(0, "Number of buffers did not match Number of arguments for %s!");
  }
}

void clkernel::set_par(parameters parn) {
  par = parn;
  //  cl_mem par_buffer = clCreateBuffer(context->context, CL_MEM_READ_ONLY, sizeof(parameters), NULL, &ret);
  //  ret = clEnqueueWriteBuffer(context->command_queue, par_buffer, CL_TRUE, 0, sizeof(parameters), &par, 0, NULL, NULL);
  ret = clSetKernelArg(kernel, 0, sizeof(parameters), (void *) &par);
  //  logger->log(1, ret, "Setting parameters: %d %f", par.sx, par.dt );
}

void clkernel::check_bindings() {
  std::string::size_type start = s_sourceCode.find("kernel_main(", 0);
  std::string::size_type end = s_sourceCode.find(")", start);
  std::string s_argList = s_sourceCode.substr(start,end-start);
  int needed_bindings = std::count(s_argList.begin(), s_argList.end(), ',');
  if ( binding_count != needed_bindings) logger->log(0, "ERROR: Missing Bindings for %s (%d/%d)\n", s_name.c_str(), binding_count, needed_bindings);
  bindings_checked = true;
}

size_t clkernel::get_pos_of_argument_from_src(std::string s_kernelArg) {
  // find the index of kernel argument in kernel source
  std::string::size_type start = s_sourceCode.find("kernel_main(", 0);
  std::string::size_type end = s_sourceCode.find(s_kernelArg, start);
  std::string s_argList = s_sourceCode.substr(start,end-start);

  return std::count(s_argList.begin(), s_argList.end(), ',');
}

void clkernel::bind_custom(std::string s_kernelArg, void *s, size_t custom_size) {
  // use this with care. kernels need to know what s looks like. void pointer are usually a bad idea
  // provides direct access to the opencl api so to say
  size_t pos = get_pos_of_argument_from_src(s_kernelArg);
  ret = clSetKernelArg(kernel, pos, custom_size, s);
  binding_count++;
  logger->log(1, ret, "binding custom data type to kernel %s at argument %d to %s\n",  s_name.c_str(), pos, s_kernelArg.c_str());
}

void clkernel::bind(const int pos, clbuffer *b) {
  ret = clSetKernelArg(kernel, pos, sizeof(b->buffer), (void *) &b->buffer);
  binding_count++;
  // logger->log(1, ret, "binding to %d: buffer %s with size %d\n", pos, b->s_name.c_str(), sizeof(b->buffer));
}

// overload to set integers
void clkernel::bind(std::string s_kernelArg, unsigned int myInt) {
  size_t pos = get_pos_of_argument_from_src(s_kernelArg);

  ret = clSetKernelArg(kernel, pos, sizeof(unsigned int), (void *) &myInt);
  binding_count++;
  logger->log(1, ret, "binding int to kernel %s at argument %d - int %d to %s\n",  s_name.c_str(), pos, myInt, s_kernelArg.c_str());
}

void clkernel::bind() {
  if (v_argBuffer.size()>0) {
    // 1 reserved for parameters
    for (int i = 1; i<=v_argBuffer.size(); i++) {
      clbuffer *b = v_argBuffer.at(i-1);
      ret = clSetKernelArg(kernel, i, sizeof(b->buffer), (void *) &b->buffer);
      logger->log(1, ret, "binding buffer %s to %d of %s\n", b->s_name.c_str(), i, s_name.c_str());
    }
  }
}

void clkernel::bind(std::string s_kernelArg, clbuffer *b) {
  size_t pos = get_pos_of_argument_from_src(s_kernelArg);
  ret = clSetKernelArg(kernel, pos, sizeof(b->buffer), (void *) &b->buffer);
  binding_count++;
  logger->log(1, ret, "binding buffer to kernel %s at argument %d - buffer %s to %s\n",  s_name.c_str(), pos, b->s_name.c_str(), s_kernelArg.c_str());
}

void clkernel::step(int kx, int ky, int kz) {
  if (!bindings_checked) check_bindings();
  // specify dimensions, number of kernels
  //  std::cout << par.sx << "" << par.sy << std::endl;
  size_t global_item_size[3] = {size_t(kx), size_t(ky), size_t(kz)};
  // size_t local_item_size[3] = {1, 1, 1};


  if (context->profiling) {
    // logger->log(1, "kernel step (%d %d %d) %s ...\n", kx, ky, kz, s_name.c_str());
    cl_ulong ev_start_time, ev_end_time;
    cl_event profiling_event;
    ret = clEnqueueNDRangeKernel(context->command_queue, kernel, 3, NULL, global_item_size, NULL, 0, NULL, &profiling_event);
    ret = clWaitForEvents(1, &profiling_event);
    ev_start_time = (cl_ulong)0;
    ev_end_time = (cl_ulong)0;
    size_t return_bytes;
    ret = clGetEventProfilingInfo(profiling_event, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &ev_start_time, &return_bytes);
    ret = clGetEventProfilingInfo(profiling_event, CL_PROFILING_COMMAND_END,   sizeof(cl_ulong), &ev_end_time,  &return_bytes);
    context->profiling_time += (double)(ev_end_time - ev_start_time);
    context->profiling_count += 1;
    clReleaseEvent(profiling_event);
    logger->plog((long long)(ev_end_time - ev_start_time), s_name, "cl     kernel_step");
    // logger->plog((long long)(ev_end_time - ev_start_time), "%s", s_name.c_str());
    // logger->string_printf("foo %d", 3);
  } else {
    logger->log(1, "kernel step (%d %d %d) %s ...\n", global_item_size[0], global_item_size[1], global_item_size[2], s_name.c_str());
    cl_event dummy_event;
    ret = clEnqueueNDRangeKernel(context->command_queue, kernel, 3, NULL, global_item_size, NULL, 0, NULL, &dummy_event);
    // clEnqueueNDRangeKernel is non-blocking... which can cause confusion
    // ret = clWaitForEvents(1, &dummy_event);
    // dont wait unless necessary -> let the cpu do stuff meanwhile!
    clReleaseEvent(dummy_event);
  }
  // logger->log(1, ret, "kernel step (%d %d %d) %s\n", kx, ky, kz, s_name.c_str());
}
