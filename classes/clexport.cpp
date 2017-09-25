#include "clexport.h"

clexport::clexport(clcontext *context_, std::string name_, std::string kname_, parameters par_, clbuffer *mom, clbuffer *sc0, clbuffer *sc1, clbuffer *sc2, int d_, int r_, bool img_, bool ts_, bool vp_) {
  context = context_;
  par = par_;
  logger = context->logger;
  s_name = name_;
  s_kname = kname_;
  b_source_scalars[0] = sc0;
  b_source_scalars[1] = sc1;
  b_source_scalars[2] = sc2;
  b_source_momenta = mom;
  dim = d_;
  ref = r_;

  kx = par.sx;
  ky = par.sy;
  kz = par.sz;

  if      (dim == 0) kx = 1;
  else if (dim == 1) ky = 1;
  else if (dim == 2) kz = 1;
  else logger->log(0, "ERROR: no cut plane given for %s\n", s_name.c_str());


  img = img_;
  ts  = ts_;
  vp  = vp_;

  img_every = 60.0;
  ts_every  = 60.0;
  vp_every  = 10.0*120.0;
  t_img = int(std::ceil(img_every/par.dT));
  t_ts  = int(std::ceil( ts_every/par.dT));
  t_vp  = int(std::ceil( vp_every/par.dT));

  // overwrite for exporting every step
  // t_img = 1;
  // t_ts = 1;
  // t_vp = 1;

  logger->log(1, "%s size: %03d   %03d   %03d\n",  s_name.c_str(), kx, ky, kz);
  logger->log(1, "%s mode: img %d   ts %d   vp %d\n",s_name.c_str(), img, ts, vp);
  logger->log(1, "%s writing files every: %d (%2.0fs)   %d (%2.0fs)   %d (%2.0fs)\n",s_name.c_str(), t_img, img_every, t_ts, ts_every, t_vp, vp_every);

  int ret; // just to shut up warnings from the compiler
  s_img_path = logger->s_output+"img/" + s_name;
  if (img) ret = system(("mkdir -p " + s_img_path).c_str());

  s_ts_path = logger->s_output+"timeseries/";
  if (ts) {
    ret = system(("mkdir -p " + s_ts_path).c_str());
    s_ts_path += s_name + ".ts";
    ts_stream.open(s_ts_path, std::ofstream::out);
    ts_stream << "#new entry every " << ts_every << " seconds" << std::endl;
    ts_stream << "#iteration\t" << s_name << std::endl;
    ts_stream << std::scientific;
    ts_stream << std::setprecision(10);
  }

  s_vp_path = logger->s_output+"verticalprofiles/" + s_name;
  if (vp) ret = system(("mkdir -p " + s_vp_path).c_str());

  be_export_images.push_back(new clbuffer(context, "be_" + s_name, kx, ky, kz));

  ke_render = new clkernel(context, par, s_kname, "./classes/parameters.h", "./kernels/k_header.cl");

  ke_render->bind("b_source_scalars_0", b_source_scalars[0]);
  ke_render->bind("b_source_scalars_1", b_source_scalars[1]);
  ke_render->bind("b_source_scalars_2", b_source_scalars[2]);
  ke_render->bind("b_source_momenta",   b_source_momenta);
  ke_render->bind("b_target",           be_export_images[0]);
  ke_render->bind("ref", ref);
  ke_render->bind("dim", dim);
}

clexport::~clexport() {
  delete ke_render;
  if (ts) {
    ts_stream.close();
    ts_stream.clear();
  }
  // buffers will be deleted from context
}

void clexport::write_files(int it) {

  if (it%t_img!=0 && it%t_vp!=0 && it%t_ts!=0) return;

  // logger->log(1, "%s writing files\n", s_name.c_str());

  ke_render->step(kx, ky, kz);
  if (context->profiling) logger->start_new_timer();

  for (int e = 0; e < be_export_images.size(); e++) {
    be_export_images[e]->device2ram();

    if (vp && it%t_vp==0) {
      std::stringstream ss_vpName;
      ss_vpName.str("");
      ss_vpName << s_vp_path;
      ss_vpName << "/" << std::setw(5) << std::setfill('0') << it/t_vp << ".vp";
      vp_stream.open(ss_vpName.str(), std::ofstream::out);
      vp_stream << std::fixed << std::setprecision(2);
      vp_stream << "#new file every " << vp_every << " seconds" << std::endl;
      vp_stream << "#height(m)\t" << s_name << std::endl;

      // also write averaged profiling information in this increment
      if (context->profiling) logger->avg_plog();
    }

    float r, g, b, a;
    float ts_var = 0.0;
    float vp_var = 0.0;
    int ki, kj;
    if      (dim == 0) {ki = ky; kj = kz;} // YZ
    else if (dim == 1) {ki = kx; kj = kz;} // XZ
    else if (dim == 2) {ki = kx; kj = ky;} // XY
    gdImagePtr im;
    im = gdImageCreateTrueColor(ki, kj);

    for (int j = 0; j < kj; j++) {
      for (int i = 0; i < ki; i++) {
        r = be_export_images[e]->get(i,j,0,0);
        g = be_export_images[e]->get(i,j,0,1);
        b = be_export_images[e]->get(i,j,0,2);
        a = be_export_images[e]->get(i,j,0,3);

        if (img && it%t_img==0) {
          gdImageSetPixel(im, i, kj-(j+1), gdImageColorExact(im, int(r), int(g), int(b)));
        }
        if (ts && it%t_ts==0) {
          // running average to avoid numerical problems
          ts_var = (float(i+j*ki)*ts_var + a)/float(i+j*ki+1);
        }
        if (vp && it%t_vp==0) {
          // want discrete values for z and average over x or y. z is always j for dim 0 and 1!
          vp_var = (float(i)*vp_var + a)/float(i+1);
        }
      }

      if (vp && it%t_vp==0) {
        vp_stream << std::fixed << std::setprecision(2);
        vp_stream << j*par.dz << "\t";
        vp_stream << std::scientific  << std::setprecision(5);
        vp_stream << vp_var << std::endl;
        vp_stream.flush();
        vp_var = 0.0;
      }
    }

    if (img && it%t_img==0) {
      std::stringstream ss_imageName;
      ss_imageName.str("");
      ss_imageName << s_img_path;
      ss_imageName << "/" << std::setw(5) << std::setfill('0') << it/t_img << ".png";
      FILE *out = fopen(ss_imageName.str().c_str(), "w");
      gdImagePng(im, out);
      fclose(out);
      gdImageDestroy(im);
    }

    if (ts && it%t_ts==0) {
      ts_stream << std::fixed << std::setprecision(0) << std::setw(5) << std::setfill('0') << it/t_ts;
      ts_stream << std::scientific << std::setprecision(5) << "\t" << ts_var << std::endl;
      ts_stream.flush();
    }

    if (vp && it%t_vp==0) {
      vp_stream.close();
      vp_stream.clear();
    }
  }

  if (context->profiling) {
    long long ns = logger->end_last_timer();
    logger->plog(ns, s_name, "host   writing_files");
  }
}
