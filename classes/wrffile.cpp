#include "wrffile.h"

wrffile::wrffile(clcontext *contextn, cllogger *loggern, parameters parn, std::string file_namen, clbuffer *b_target_scalars[3], clbuffer *b_target_momenta) {
  context = contextn;
  logger = loggern;
  par = parn;
  file_name = file_namen;

  logger->log(0,"WRF Import from %s\n", file_name.c_str());

  // wrf integer domain size
  wrf.sx = 136;
  wrf.sy = 91;
  wrf.sz = 26;

  // wrf level height [m]
  wrf.hz[0]  =   125;
  wrf.hz[1]  =   250;
  wrf.hz[2]  =   375;
  wrf.hz[3]  =   500;
  wrf.hz[4]  =   750;
  wrf.hz[5]  =  1000;
  wrf.hz[6]  =  1250;
  wrf.hz[7]  =  1500;
  wrf.hz[8]  =  1750;
  wrf.hz[9]  =  2000;
  wrf.hz[10] =  2500;
  wrf.hz[11] =  3000;
  wrf.hz[12] =  3500;
  wrf.hz[13] =  4000;
  wrf.hz[14] =  4500;
  wrf.hz[15] =  5000;
  wrf.hz[16] =  6000;
  wrf.hz[17] =  7000;
  wrf.hz[18] =  8000;
  wrf.hz[19] =  9000;
  wrf.hz[20] = 10000;
  wrf.hz[21] = 11000;
  wrf.hz[22] = 12000;
  wrf.hz[23] = 14000;
  wrf.hz[24] = 16000;
  wrf.hz[25] = 18000;

  wrf.zoffset = 3;

  // wrf angular cell size
  wrf.dlat = 0.02227027;
  wrf.dlon = 0.02227027;

  // wrf lower left corner lat/long postion
  wrf.lon0 = 11.4946;
  wrf.lat0 = 49.9988;

  // angular to linear [m]
  wrf.dx = wrf.dlon*78847.0;
  wrf.dy = wrf.dlat*111132.0;

  // target domaincenter lat/long
  wrf.domaincenterx = 12.927716;
  wrf.domaincentery = 51.525345;

  // target domainsize [m]
  wrf.domainsizex = par.sx*par.dx;
  wrf.domainsizey = par.sy*par.dy;

  // target domainsize [cells]
  wrf.dsx = wrf.domainsizex/wrf.dx;
  wrf.dsy = wrf.domainsizey/wrf.dy;

  // target ll [cells] in wrf cs
  wrf.dx0 = (wrf.domaincenterx-wrf.lon0)/wrf.dlon-wrf.dsx/2.0;
  wrf.dy0 = (wrf.domaincentery-wrf.lat0)/wrf.dlat-wrf.dsy/2.0;

  // wrf source texture
  T2       = new float[wrf.sx*wrf.sy];
  PSFC     = new float[wrf.sx*wrf.sy];
  U10      = new float[wrf.sx*wrf.sy];
  V10      = new float[wrf.sx*wrf.sy];
  Q2       = new float[wrf.sx*wrf.sy];
  TSK      = new float[wrf.sx*wrf.sy];
  HSFC     = new float[wrf.sx*wrf.sy];
  RRTOT    = new float[wrf.sx*wrf.sy];
  HFX      = new float[wrf.sx*wrf.sy];
  LH       = new float[wrf.sx*wrf.sy];

  TK       = new float[wrf.sx*wrf.sy*wrf.sz];
  PRESSURE = new float[wrf.sx*wrf.sy*wrf.sz];
  UMET     = new float[wrf.sx*wrf.sy*wrf.sz];
  VMET     = new float[wrf.sx*wrf.sy*wrf.sz];
  W        = new float[wrf.sx*wrf.sy*wrf.sz];
  Q        = new float[wrf.sx*wrf.sy*wrf.sz];
  QC       = new float[wrf.sx*wrf.sy*wrf.sz];

  b_wrf_source_scalars_vc[0]  = new clbuffer(context, "b_wrf_source_scalars_vc_0", wrf.sx, wrf.sy, wrf.sz);
  b_wrf_source_scalars_vc[1]  = new clbuffer(context, "b_wrf_source_scalars_vc_1", wrf.sx, wrf.sy, wrf.sz);
  b_wrf_source_velocities_vc  = new clbuffer(context, "b_wrf_source_velocities_vc", wrf.sx, wrf.sy, wrf.sz);
  b_wrf_flux  = new clbuffer(context, "b_wrf_flux", wrf.sx, wrf.sy, 1);


  k_interpolate_scalars = new clkernel(context, par, "./kernels/k_wrf_interpolate_scalars.cl");
  k_interpolate_scalars->bind_custom("wrfparameters wrf", &wrf, sizeof(wrf)); // whitespaces need to match kernel
  k_interpolate_scalars->bind("b_wrf_source_scalars_0", b_wrf_source_scalars_vc[0]);
  k_interpolate_scalars->bind("b_wrf_source_scalars_1", b_wrf_source_scalars_vc[1]);
  k_interpolate_scalars->bind("b_target_scalars_0", b_target_scalars[0]);
  k_interpolate_scalars->bind("b_target_scalars_1", b_target_scalars[1]);
  k_interpolate_scalars->bind("b_target_scalars_2", b_target_scalars[2]);

  k_interpolate_momenta = new clkernel(context, par, "./kernels/k_wrf_interpolate_momenta.cl");
  k_interpolate_momenta->bind_custom("wrfparameters wrf", &wrf, sizeof(wrf));

  k_interpolate_momenta->bind("b_source_velocites", b_wrf_source_velocities_vc);
  k_interpolate_momenta->bind("b_scalars_vc", b_target_scalars[0]);
  k_interpolate_momenta->bind("b_target_momenta", b_target_momenta);

  logger->log(0,"Domain in wrf: %g %g %g %g\n",wrf.dx0,wrf.dy0,wrf.dx0+wrf.dsx,wrf.dy0+wrf.dsy);
  logger->log(0,"Cellsize in wrf: %g %g\n",wrf.dx,wrf.dy);
}

inline int wrffile::index(int x,int y,int z) {
  return int(x+wrf.sx*y+wrf.sx*wrf.sy*z);
}

void wrffile::load(int wrfindex) {
  char  cname[255];
  sprintf(cname,"%s_%02d.dat",file_name.c_str(),wrfindex);
  float u,v,w,p,q,qc,T,sig,rhod,rhov,rho,rhoc,Rml,Cpml;

  printf("loading: %s\n",cname);
  FILE *f=fopen(cname,"rb");
  fread(T2,sizeof(float),wrf.sx*wrf.sy,f);
  fread(PSFC,sizeof(float),wrf.sx*wrf.sy,f);
  fread(U10,sizeof(float),wrf.sx*wrf.sy,f);
  fread(V10,sizeof(float),wrf.sx*wrf.sy,f);
  fread(Q2,sizeof(float),wrf.sx*wrf.sy,f);
  fread(TSK,sizeof(float),wrf.sx*wrf.sy,f);
  fread(HSFC,sizeof(float),wrf.sx*wrf.sy,f);
  fread(RRTOT,sizeof(float),wrf.sx*wrf.sy,f);
  fread(HFX,sizeof(float),wrf.sx*wrf.sy,f);
  fread(LH,sizeof(float),wrf.sx*wrf.sy,f);
  fread(TK,sizeof(float),wrf.sx*wrf.sy*wrf.sz,f);
  fread(PRESSURE,sizeof(float),wrf.sx*wrf.sy*wrf.sz,f);
  fread(UMET,sizeof(float),wrf.sx*wrf.sy*wrf.sz,f);
  fread(VMET,sizeof(float),wrf.sx*wrf.sy*wrf.sz,f);
  fread(W,sizeof(float),wrf.sx*wrf.sy*wrf.sz,f);
  fread(Q,sizeof(float),wrf.sx*wrf.sy*wrf.sz,f);
  fread(QC,sizeof(float),wrf.sx*wrf.sy*wrf.sz,f);
  fclose(f);

  for (int x=0; x<wrf.sx; x++) {
    for (int y=0; y<wrf.sy; y++) {
      for (int z=0; z<wrf.sz; z++) {
        p =PRESSURE[index(x,y,z)]*100.;
        q =Q[index(x,y,z)];
        qc=QC[index(x,y,z)];
        T =TK[index(x,y,z)];
        u =UMET[index(x,y,z)];
        v =VMET[index(x,y,z)];
        w =W[index(x,y,z)];

        rhod =p/(par.rd+q*par.rv)/T;
        rhov =q*rhod;
        rhoc =qc*rhod;
        rho  =rhov+rhod+rhoc;
        Rml  =(rhod*par.rd+rhov*par.rv)/rho;
        Cpml =(rhod*par.cpd+rhov*par.cpv+rhoc*par.cpl)/rho;
        sig  =Cpml*log(T)-Rml*log(p);

        b_wrf_source_scalars_vc[0]->set(x,y,z,0,sig*rho);  // sig
        b_wrf_source_scalars_vc[0]->set(x,y,z,1,rho);      // rho
        b_wrf_source_scalars_vc[0]->set(x,y,z,2,rhov);     // rhov
        b_wrf_source_scalars_vc[0]->set(x,y,z,3,rhoc);     // rhoc

        b_wrf_source_scalars_vc[1]->set(x,y,z,0,0.0);            // rhor
        b_wrf_source_scalars_vc[1]->set(x,y,z,1,1000e6*rho);     // nd
        b_wrf_source_scalars_vc[1]->set(x,y,z,2,800e6*1e5*rhoc); // nc
        b_wrf_source_scalars_vc[1]->set(x,y,z,3,0.0);            // nr

        b_wrf_source_velocities_vc->set(x,y,z,0,u);    // v
        b_wrf_source_velocities_vc->set(x,y,z,1,v);    // u
        b_wrf_source_velocities_vc->set(x,y,z,2,w);    // w
        b_wrf_source_velocities_vc->set(x,y,z,3,0.0);  // unused
      }
    }
  }

  // fluxes for k_ext_forcings
  float hfxm=0.; //direct
  float lhm=0.;  //latent
  for (int x=0; x<wrf.sx; x++)
  {
    for (int y=0; y<wrf.sy; y++)
    {
      b_wrf_flux->set(x,y,0,0,HFX[index(x,y,0)]);
      b_wrf_flux->set(x,y,0,1, LH[index(x,y,0)]);
      b_wrf_flux->set(x,y,0,2,0.0);
      b_wrf_flux->set(x,y,0,3,0.0);
      hfxm+=HFX[index(x,y,0)];
      lhm += LH[index(x,y,0)];
    }
  }

  logger->log(0, "total flux hfxm: %g lhm: %g\n",hfxm/wrf.sx/wrf.sy,lhm/wrf.sx/wrf.sy);

  b_wrf_source_scalars_vc[0]->ram2device();
  b_wrf_source_scalars_vc[1]->ram2device();
  b_wrf_source_velocities_vc->ram2device();
  b_wrf_flux->ram2device();

  k_interpolate_scalars->step(1,1,1);
  k_interpolate_momenta->step(1,1,1);
}

