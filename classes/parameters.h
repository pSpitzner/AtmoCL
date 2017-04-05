#ifndef Tropos_parameters_h
#define Tropos_parameters_h

typedef struct {
  int sx, sy, sz; //unsigned long = size_t
  float dx, dy, dz;
  float dT;           // large timestep
  float dt;           // maximal small timestep
  float dti;          // temporary small timestep
  int nout;           // export picture every nout steps
  int ns;             // timesplitting
  // read from profile now
  float ui, vi, wi;   // initial velocities
  // float ri;           // initial density
  // float ti;           // initial potential temperature
  // float pi;           // initial pressure
  // float nu;           // divergence damping constant
  float gr;           // gravity
  float cpd;          // specific heat for dry air
  float cpv;         // vapour
  float cpl;         // liquid
  float cpi;         // snow

  float rd;           // Gas constant for dry air
  float rv;

  float rr;           // rho_0 seifert 2005
  float pr;           // reference surface pressure
  float sr;           // rho*sigma offset to avoid numerical errors
  float lr;   // reference latent heat of evaporation at reference temperature
  float lre0;  // reference latent heat of evap at 0 kelvin
  float lrs0;  // reference latent heat of sublimation at 0 kelvin
  float lvf;  // latent heat of freezing
  float tr;   // reference temperature

  float svr;
  float ri;
  float zi;

  int timescheme;

  float b21;          // constants for MIS
  float b31, b32;
  float b41, b42, b43;
  float a32;
  float a42, a43;
  float g32;
  float g42, g43;
  float b[3][3];
  float a[3][3];
  float g[3][3];
  float dtis[3];
  int nsi[3];


} parameters;




#endif











// have 50 lines long header to find kernel debug messages easily
