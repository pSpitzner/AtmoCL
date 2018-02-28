constant sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP | CLK_FILTER_NEAREST;

// header
typedef struct {
  int x, y, z;
  int xr, yr, zr;
  int xl, yl, zl;
  int xrr, yrr, zrr;
  int xll, yll, zll;
  int s_xl, s_xr, s_yl, s_yr, s_zl, s_zr;
  int s_yll, s_yrr, s_xll, s_xrr, s_zll, s_zrr;
  float ulf, urf, vlf, vrf, wlf, wrf;
} position;

typedef struct {
  float rho;     // rho total
  float rho_v;   // rho vapour
  float rho_c;   // rho liquid cloud
  float rho_r;   // rho liquid rain
  float rho_i;   // rho liquid ice
  float rho_s;   // rho liquid snow
  float rho_l;   // rho liquid total
  float rho_d;   // rho dry
  float rho_f;   // rho frozen total, basicly ice and snow
  float sig;     // rho*sigma
  float cpml;    // mixed specific heat
  float rml;     // mixed gas constant
  float lnP;     // log. total pressure
  float lnT;     // log temperature
  float T;
  float P;
  float pv;      // vapour pressure
  float sv;      // saturation vapour pressure of water
  float svi;     // saturation vapour pressure of ice
  float lv;      // latent heat
  float n_d;     // nr. dirt
  float n_c;     // nr. cloud droplets
  float n_r;     // nr. rain droplets
  float n_i;     // nr. ice droplets
  float n_s;     // nr. snow droplets
  float sat;     // saturation
} state;

__constant bool include_ice=true;
// __constant bool include_ice=false;

__private float rand(unsigned int seed);
__private position get_pos_bc(parameters *par);
__private state init_state(parameters *par, float8 *c);
__private state init_state_with_ice(parameters *par, float8 *c, float4 *cice);

void central_dif(parameters *par, position *pos, float4 *fyn, float4 *x_c, float4 *x_cr, float4 *y_c, float4 *y_cr, float4 *z_c, float4 *z_cr);
__private float4 read_f4(int x, int y, int z, __read_only image3d_t buf_0);
__private float8 read_f8(int x, int y, int z, __read_only image3d_t buf_0, __read_only image3d_t buf_1);
// __private float4 read_f12(int x, int y, int z, __read_only image3d_t buf_0, __read_only image3d_t buf_1, __read_only image3d_t buf_2);
// __private float4 read_f16(int x, int y, int z, __read_only image3d_t buf_0, __read_only image3d_t buf_1, __read_only image3d_t buf_2, __read_only image3d_t buf_3);

void write_f4(int x, int y, int z, float4 *f, __write_only image3d_t buf_0);
void write_f8(int x, int y, int z, float8 *f, __write_only image3d_t buf_0, __write_only image3d_t buf_1);
// void write_f12(int x, int y, int z, float12 f, __write_only image3d_t buf_0, __write_only image3d_t buf_1, __write_only image3d_t buf_2);
// void write_f16(int x, int y, int z, float16 f, __write_only image3d_t buf_0, __write_only image3d_t buf_1, __write_only image3d_t buf_2, __write_only image3d_t buf_3);

__private float rand(unsigned int seed) {
  unsigned int rand_a = 16807;
  unsigned int rand_m = 2147483647;
  unsigned int rand_c = 1043;
  unsigned int rand_z = seed;
  for (int i = 0; i < 10; i++) rand_z = (rand_a*rand_z+rand_c)%rand_m;
  return (float)(rand_z)/(float)(rand_m);
}

__private position get_pos_bc(parameters *par) {
  position pos;
  pos.x = get_global_id(0);
  pos.y = get_global_id(1);
  pos.z = get_global_id(2);

  size_t sx = get_global_size(0);
  size_t sy = get_global_size(1);
  size_t sz = get_global_size(2);

  // periodic bc
  pos.xr  = ((pos.x  + 1 == sx) ?      0 : pos.x + 1 );
  pos.xl  = ((pos.x  - 1 < 0)   ? sx - 1 : pos.x - 1 );
  pos.yr  = ((pos.y  + 1 == sy) ?      0 : pos.y + 1 );
  pos.yl  = ((pos.y  - 1 < 0)   ? sy - 1 : pos.y - 1 );
  pos.zr  = ((pos.z  + 1 == sz) ?      0 : pos.z + 1 );
  pos.zl  = ((pos.z  - 1 < 0)   ? sz - 1 : pos.z - 1 );

  pos.xrr = ((pos.xr + 1 == sx) ?      0 : pos.xr + 1 );
  pos.xll = ((pos.xl - 1 < 0)   ? sx - 1 : pos.xl - 1 );
  pos.yrr = ((pos.yr + 1 == sy) ?      0 : pos.yr + 1 );
  pos.yll = ((pos.yl - 1 < 0)   ? sy - 1 : pos.yl - 1 );
  pos.zrr = ((pos.zr + 1 == sz) ?      0 : pos.zr + 1 );
  pos.zll = ((pos.zl - 1 < 0)   ? sz - 1 : pos.zl - 1 );

  // fixing select bc
  // 1 normal cell, -1 boundary cell
  pos.s_xl = 1;
  pos.s_xr = 1;
  pos.s_yl = 1;
  pos.s_yr = 1;
  pos.s_zl = 1;
  pos.s_zr = 1;

  pos.s_xll = 1;
  pos.s_xrr = 1;
  pos.s_yll = 1;
  pos.s_yrr = 1;
  pos.s_zll = 1;
  pos.s_zrr = 1;

  if (pos.z==sz-1) {pos.zr  = pos.z; pos.zrr = pos.zl; pos.s_zr = -1; pos.s_zrr=-1;}
  if (pos.z==sz-2) {pos.zrr = pos.zr; pos.s_zrr=-1;}
  if (pos.z==   0) {pos.zl  = pos.z; pos.zll = pos.zr; pos.s_zl = -1; pos.s_zll=-1;}
  if (pos.z==   1) {pos.zll = pos.zl; pos.s_zll=-1;}

  return pos;
}

__private state init_state(parameters *par, float8 *c) {
  // c does not support default arguments, work around
  float4 *cice = 0;
  return init_state_with_ice(par, c, cice);
}
__private state init_state_with_ice(parameters *par, float8 *c, float4 *cice) {
  state st;

  st.sig   = (*c).s0;
  st.rho   = (*c).s1;
  st.rho_v = (*c).s2;
  st.rho_c = (*c).s3;
  st.rho_r = max(0.0f, (*c).s4);
  // clamp to avoid negative ns
  st.n_d   = max(0.0f, (*c).s5);
  st.n_c   = max(0.0f, (*c).s6);
  st.n_r   = max(0.0f, (*c).s7);

  //ice
  st.rho_i = (*cice).s0;
  st.rho_s = (*cice).s1;
  st.n_i   = max(0.0f, (*cice).s2);
  st.n_s   = max(0.0f, (*cice).s3);

  st.rho_l = st.rho_c+st.rho_r;
  st.rho_f = st.rho_i+st.rho_s;
  st.rho_d = st.rho-st.rho_v-st.rho_l-st.rho_f;
  st.cpml  = st.rho_d*par->cpd+st.rho_v*par->cpv+st.rho_l*par->cpl+st.rho_f*par->cpi;
  st.rml   = st.rho_d*par->rd+st.rho_v*par->rv;
  st.lnP   = st.sig/(st.cpml-st.rml)+1.0f/(1.0f-st.rml/st.cpml)*log(st.rml);
  st.lnT   = st.sig/(st.cpml-st.rml)+log(st.rml)/(st.cpml/st.rml-1.0f);
  st.T     = exp(st.lnT);
  st.P     = exp(st.lnP);
  st.pv    = st.T*par->rv*st.rho_v;
  st.sv    = par->svr*pow(st.T/par->tr, (par->cpv-par->cpl)/par->rv)*exp(par->lre0/par->rv*(1.0f/par->tr-1.0f/st.T));
  st.svi   = par->svr*pow(st.T/par->tr, (par->cpv-par->cpi)/par->rv)*exp(par->lrs0/par->rv*(1.0f/par->tr-1.0f/st.T));
  st.lv    = par->lre0+(par->cpv-par->cpl)*st.T;
  st.sat   = 100.0f*(st.pv/st.sv-1.0f); // relative luftfecuhte über 100 => sat=1 entspricht 101 rel f.

  return st;
}

__private float4 read_f4(int x, int y, int z, __read_only image3d_t buf_0) {
  int4 pos = (int4)(x, y, z, 0);
  float4 f_0 = read_imagef(buf_0, pos);
  return f_0;
}
__private float8 read_f8(int x, int y, int z, __read_only image3d_t buf_0, __read_only image3d_t buf_1) {
  int4 pos = (int4)(x, y, z, 0);
  float4 f_0 = read_imagef(buf_0, pos);
  float4 f_1 = read_imagef(buf_1, pos);
  return (float8)(f_0, f_1);
}

void write_f4(int x, int y, int z, float4 *f, __write_only image3d_t buf_0) {
  int4 pos = (int4)(x, y, z, 0);
  write_imagef(buf_0, pos, (float4)((*f).s0, (*f).s1, (*f).s2, (*f).s3));
}
void write_f8(int x, int y, int z, float8 *f, __write_only image3d_t buf_0, __write_only image3d_t buf_1) {
  int4 pos = (int4)(x, y, z, 0);
  write_imagef(buf_0, pos, (float4)((*f).s0, (*f).s1, (*f).s2, (*f).s3));
  write_imagef(buf_1, pos, (float4)((*f).s4, (*f).s5, (*f).s6, (*f).s7));
}

void central_dif(parameters *par, position *pos,
                 float4 *fyn,
                 float4 *x_c, float4 *x_cr,
                 float4 *y_c, float4 *y_cr,
                 float4 *z_c, float4 *z_cr)
{
  float4 Fx, Fy, Fz;
  // Fx = (x_c*pos.ulf - x_cr*pos.urf);
  // Fx /= par->dx;

  // Fy = (y_c*pos.vlf - y_cr*pos.vrf);
  // Fy /= par->dy;

  // Fz = (z_c*pos.wlf - z_cr*pos.wrf);
  // Fz /= par->dz;

  Fx = ((*x_c)*pos->ulf/par->dx - (*x_cr)*pos->urf/par->dx);
  Fy = ((*y_c)*pos->vlf/par->dy - (*y_cr)*pos->vrf/par->dy);
  Fz = ((*z_c)*pos->wlf/par->dz - (*z_cr)*pos->wrf/par->dz);


  *fyn = Fx + Fy + Fz;
}



















// have some decent number of lines in header to find kernel debug messages easily
