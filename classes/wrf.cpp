#include "wrf.h"

wrf::wrf(clcontext *contextn, cllogger *loggern, parameters parn, std::string file_namen, clbuffer *b_target_scalars[3]) {
  context = contextn;
  logger = loggern;
  par = parn;
  file_name = file_namen;

  logger->log(0,"WRF Import from %s\n", file_name.c_str());

  // wrf integer domain size
  sx = 136;
  sy = 91;
  sz = 26;

  // wrf level height [m]
  hz = new float[sz];
  // hz[0]  =     0;
  hz[0]  =   125;
  hz[1]  =   250;
  hz[2]  =   375;
  hz[3]  =   500;
  hz[4]  =   750;
  hz[5]  =  1000;
  hz[6]  =  1250;
  hz[7]  =  1500;
  hz[8]  =  1750;
  hz[9]  =  2000;
  hz[10] =  2500;
  hz[11] =  3000;
  hz[12] =  3500;
  hz[13] =  4000;
  hz[14] =  4500;
  hz[15] =  5000;
  hz[16] =  6000;
  hz[17] =  7000;
  hz[18] =  8000;
  hz[19] =  9000;
  hz[20] = 10000;
  hz[21] = 11000;
  hz[22] = 12000;
  hz[23] = 14000;
  hz[24] = 16000;
  hz[25] = 18000;

  // wrf angular cell size
  dlat = 0.02227027;
  dlon = 0.02227027;

  // wrf lower left corner lat/long postion
  lon0 = 11.4946;
  lat0 = 49.9988;

  // angular to linear [m]
  dx = dlon*78847.0;
  dy = dlat*111132.0;

  // target domaincenter lat/long
  domaincenterx = 12.927716;
  domaincentery = 51.525345;

  // target domainsize [m]
  domainsizex = par.sx*par.dx;
  domainsizey = par.sy*par.dy;

  // target domainsize [cells]
  dsx = domainsizex/dx;
  dsy = domainsizey/dy;

  // target ll [cells] in wrf cs
  dx0 = (domaincenterx-lon0)/dlon-dsx/2.;
  dy0 = (domaincentery-lat0)/dlat-dsy/2.;

  // wrf source texture
  T2       = new float[sx*sy];
  PSFC     = new float[sx*sy];
  U10      = new float[sx*sy];
  V10      = new float[sx*sy];
  Q2       = new float[sx*sy];
  TSK      = new float[sx*sy];
  HSFC     = new float[sx*sy];
  RRTOT    = new float[sx*sy];
  HFX      = new float[sx*sy];
  LH       = new float[sx*sy];

  TK       = new float[sx*sy*sz];
  PRESSURE = new float[sx*sy*sz];
  UMET     = new float[sx*sy*sz];
  VMET     = new float[sx*sy*sz];
  W        = new float[sx*sy*sz];
  Q        = new float[sx*sy*sz];
  QC       = new float[sx*sy*sz];

  b_wrf_source_vc[0]  = new clbuffer(context, "b_wrf_source_vc_0", sx, sy, sz);
  b_wrf_source_vc[1]  = new clbuffer(context, "b_wrf_source_vc_1", sx, sy, sz);
  b_wrf_source_vc[2]  = new clbuffer(context, "b_wrf_source_vc_2", sx, sy, sz);
  b_wrf_flux  = new clbuffer(context, "b_wrf_flux", sx, sy, 1);

  k_interpolate = new clkernel(context, par, "./kernels/k_initterpolate.cl");

  logger->log(0,"Domain in wrf: %g %g %g %g\n",dx0,dy0,dx0+dsx,dy0+dsy);
  logger->log(0,"Cellsize in wrf: %g %g\n",dx,dy);
}

inline int wrf::index(int x,int y,int z) {
  return int(x+sx*y+sx*sy*z);
}

void wrf::load(int wrfindex) {
  char  cname[255];
  sprintf(cname,"%s_%02d.dat",file_name.c_str(),wrfindex);
  float u,v,w,p,q,qc,T,sig,rhod,rhov,rho,rhoc,Rml,Cpml;

  printf("loading: %s\n",cname);
  FILE *f=fopen(cname,"rb");
  fread(T2,sizeof(float),sx*sy,f);
  fread(PSFC,sizeof(float),sx*sy,f);
  fread(U10,sizeof(float),sx*sy,f);
  fread(V10,sizeof(float),sx*sy,f);
  fread(Q2,sizeof(float),sx*sy,f);
  fread(TSK,sizeof(float),sx*sy,f);
  fread(HSFC,sizeof(float),sx*sy,f);
  fread(RRTOT,sizeof(float),sx*sy,f);
  fread(HFX,sizeof(float),sx*sy,f);
  fread(LH,sizeof(float),sx*sy,f);
  fread(TK,sizeof(float),sx*sy*sz,f);
  fread(PRESSURE,sizeof(float),sx*sy*sz,f);
  fread(UMET,sizeof(float),sx*sy*sz,f);
  fread(VMET,sizeof(float),sx*sy*sz,f);
  fread(W,sizeof(float),sx*sy*sz,f);
  fread(Q,sizeof(float),sx*sy*sz,f);
  fread(QC,sizeof(float),sx*sy*sz,f);
  fclose(f);

  for (int x=0; x<sx; x++)
  {
    for (int y=0; y<sy; y++)
    {
      for (int z=0; z<sz; z++)
      {
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

        b_wrf_source_vc[0]->set(x,y,z,0,u);   // v
        b_wrf_source_vc[0]->set(x,y,z,1,v);   // u
        b_wrf_source_vc[0]->set(x,y,z,2,w);   // w
        b_wrf_source_vc[0]->set(x,y,z,3,0.0);  // unused

        b_wrf_source_vc[1]->set(x,y,z,0,sig*rho);  // sig
        b_wrf_source_vc[1]->set(x,y,z,1,rho);  // rho
        b_wrf_source_vc[1]->set(x,y,z,2,rhov); // rhov
        b_wrf_source_vc[1]->set(x,y,z,3,rhoc); // rhoc

        b_wrf_source_vc[2]->set(x,y,z,0,0.0);  // rhoc // rhor?
        b_wrf_source_vc[2]->set(x,y,z,1,800e6*1e5*rhoc);  // nc
        b_wrf_source_vc[2]->set(x,y,z,2,0.0);  // rhos
        b_wrf_source_vc[2]->set(x,y,z,3,1000e6*rho); // Nv

        // template
        // (*cRhs).s0 += dsig;
        // (*cRhs).s1 += drho;
        // (*cRhs).s2 += drho_v;
        // (*cRhs).s3 += drho_c;
        // (*cRhs).s4 += drho_r;
        // (*cRhs).s5 += dn_d;
        // (*cRhs).s6 += dn_c;
        // (*cRhs).s7 += dn_r;
        // (*cRhs_ice).s0 += drho_i;
        // (*cRhs_ice).s1 += drho_s;
        // (*cRhs_ice).s2 += dn_i;
        // (*cRhs_ice).s3 += dn_s;
      }
    }
  }

  // fluxes for k_ext_forcings
  float hfxm=0.; //direct
  float lhm=0.;  //latent
  for (int x=0; x<sx; x++)
  {
    for (int y=0; y<sy; y++)
    {
      b_wrf_flux->set(x,y,0,0,HFX[index(x,y,0)]);
      b_wrf_flux->set(x,y,0,1, LH[index(x,y,0)]);
      b_wrf_flux->set(x,y,0,2,0.0);
      b_wrf_flux->set(x,y,0,3,0.0);
      hfxm+=HFX[index(x,y,0)];
      lhm += LH[index(x,y,0)];
    }
  }

  logger->log(0, "total flux hfxm: %g lhm: %g\n",hfxm/sx/sy,lhm/sx/sy);

  b_wrf_source_vc[0]->ram2device();
  b_wrf_source_vc[1]->ram2device();
  b_wrf_source_vc[2]->ram2device();
  b_wrf_flux->ram2device();
}
