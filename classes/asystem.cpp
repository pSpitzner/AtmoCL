#include "asystem.h"

void asystem::set_par(int timescheme_) {
  // dycoms
  // par.sx  = 128;     // system size
  // par.sy  = 128;
  // par.sz  = 32;
  // par.dx  = 50.0;
  // par.dy  = 50.0;
  // par.dz  = 50.0;

  // isdac
  par.sx  = 128;     // system size
  par.sy  = 1;
  par.sz  = 64;
  par.dx  = 50.0;
  par.dy  = 50.0;
  par.dz  = 20.0;

  par.ui = 0.0;
  par.vi = 0.0;
  par.wi = 0.0;

  // time steps
  par.dT = 1.0;                   // large timestep
  par.ns = 2*2*std::ceil(par.dT*350.0/sqrt(float(par.dx*par.dx + par.dy*par.dy + par.dz*par.dz)));
  // par.ns  = 12;                // max nr small per large
  logger->log(0, "ns: %d\n", par.ns);
  par.nout = 10;


  par.timescheme = timescheme_;   // mis = 1, rk3 = 0
  if (par.timescheme == 0) logger->log(0, "using rk3 time integration scheme\n");
  else if (par.timescheme == 1) logger->log(0, "using mis time integration scheme\n");

  par.gr  = 9.81;

  par.cpd = 1004.0;  // specific heat for dry air
  par.cpv = 1885.0;  // vapour
  par.cpl = 4186.0;  // liquid
  par.cpi = 2110.0;  // ice

  par.rd  = 287.0;   // Gas constant for dry air
  par.rv  = 461.0;

  par.rr  = 1.225;   // air density at surface conditions, rho_0 in seifert 2005, needed for microphys
  par.pr  = 1.0e5;   // reference pressure
  // par.sr  = 2400;    // rho*sigma offset to avoid numerical errors, not needed when using microphys -> sigma profile is more pronounced
  par.tr  = 273.15;  // reference temperature

  par.lr  = 2.5e6;   // reference latent heat of evaporation at reference temperature
  par.lre0 = par.lr+(par.cpv-par.cpl)*(-par.tr);     // reference latent heat of evap at 0 kelvin
  par.lvf = 3.33e5;   // latent heat of freezing

  par.lrs0 = par.lr+par.lvf+(par.cpv-par.cpi)*(-par.tr); // latent heat of sublimation at 0 kelvin


  par.svr = 610.7;   // ref saturation vapour pressure at ref temp

  par.ri = 1.13;      // dycoms density
  par.zi = 840.0;     // dycoms, inversion height

  // set all coefficients 0.0
  for (int i = 0; i<3; i++) {
    for (int j = 0; j<3; j++) {
      par.b[i][j] = 0.0;
      par.a[i][j] = 0.0;
      par.g[i][j] = 0.0;
    }
  }

  // constants for mis2
  // par.b[0][0] =  0.126848494553;
  // par.b[1][0] = -0.784838278826;
  // par.b[1][1] =  1.37442675268;
  // par.b[2][0] = -0.0456727081749;
  // par.b[2][1] = -0.00875082271190;
  // par.b[2][2] =  0.524775788629;

  // par.a[1][1] =  0.536946566710;
  // par.a[2][1] =  0.480892968551;
  // par.a[2][2] =  0.500561163566;

  // par.g[1][1] =  0.652465126004;
  // par.g[2][1] = -0.0732769849457;
  // par.g[2][2] =  0.144902430420;

  if (par.timescheme == 1) {
    // constants for mis3c
    par.b[0][0] =  0.397525189225;      //b21
    par.b[1][0] = -0.227036463644;      //b31
    par.b[1][1] =  0.624528794618;      //b32
    par.b[2][0] = -0.00295238076840;    //b41
    par.b[2][1] = -0.270971764284;      //b42
    par.b[2][2] =  0.671323159437;      //b43

    par.a[1][1] =  0.589557277145;      //a32
    par.a[2][1] =  0.544036601551;      //a42
    par.a[2][2] =  0.565511042564;      //a43

    par.g[1][1] =  0.142798786398;      //g32
    par.g[2][1] = -0.0428918957402;     //g42
    par.g[2][2] =  0.0202720980282;     //g43
  }

  else if (par.timescheme == 0) {
    // constants for legacy rk
    par.b[0][0] =  1.0/3.0;
    par.b[1][1] =  1.0/2.0;
    par.b[2][2] =  1.0/1.0;
  }


  // set small time steps
  for (int s = 0; s<3; s++) {
    float di = par.b[s][0] + par.b[s][1] + par.b[s][2];
    int nsi = std::ceil(di*par.ns);
    par.dtis[s] = di*par.dT/float(nsi);
    logger->log(2, "nsi: %d\tdti(%d): %f\n", nsi, s, par.dtis[s]);
  }

  par.nsi[0] = std::ceil((par.b[0][0] + par.b[0][1] + par.b[0][2])*par.ns);
  par.nsi[1] = std::ceil((par.b[1][0] + par.b[1][1] + par.b[1][2])*par.ns);
  par.nsi[2] = std::ceil((par.b[2][0] + par.b[2][1] + par.b[2][2])*par.ns);
}

asystem::asystem(clcontext *contextn, cllogger *loggern, int timescheme) {
  context = contextn;
  logger = loggern;
  debugStepCount = 0;
  frame_index = 0;
  fast_index = 0;

  set_par(timescheme);

  // ----------------------------------------------------------------- //
  // buffer //
  // ----------------------------------------------------------------- //

  b_temp[0]           = new clbuffer(context, "b_temp_0", par.sx, par.sy, par.sz);
  b_temp[1]           = new clbuffer(context, "b_temp_1", par.sx, par.sy, par.sz);
  b_temp[2]           = new clbuffer(context, "b_temp_2", par.sx, par.sy, par.sz);
  b_temp[3]           = new clbuffer(context, "b_temp_3", par.sx, par.sy, par.sz);

  bf_momenta_fc_a     = new clbuffer(context, "bf_momenta_fc_a", par.sx, par.sy, par.sz);
  bf_momenta_fc_b     = new clbuffer(context, "bf_momenta_fc_b", par.sx, par.sy, par.sz);

  bf_scalars_vc_a[0]  = new clbuffer(context, "bf_scalars_vc_a_0", par.sx, par.sy, par.sz);
  bf_scalars_vc_b[0]  = new clbuffer(context, "bf_scalars_vc_b_0", par.sx, par.sy, par.sz);
  bf_scalars_vc_a[1]  = new clbuffer(context, "bf_scalars_vc_a_1", par.sx, par.sy, par.sz);
  bf_scalars_vc_b[1]  = new clbuffer(context, "bf_scalars_vc_b_1", par.sx, par.sy, par.sz);
  bf_scalars_vc_a[2]  = new clbuffer(context, "bf_scalars_vc_a_2", par.sx, par.sy, par.sz); // ice
  bf_scalars_vc_b[2]  = new clbuffer(context, "bf_scalars_vc_b_2", par.sx, par.sy, par.sz);

  bs_momenta_fc[0]    = new clbuffer(context, "bs_momenta_fc_0", par.sx, par.sy, par.sz);
  bs_momenta_fc[1]    = new clbuffer(context, "bs_momenta_fc_1", par.sx, par.sy, par.sz);
  bs_momenta_fc[2]    = new clbuffer(context, "bs_momenta_fc_2", par.sx, par.sy, par.sz);

  bsRhs_momenta_fc[0] = new clbuffer(context, "bsRhs_momenta_fc_0", par.sx, par.sy, par.sz);
  bsRhs_momenta_fc[1] = new clbuffer(context, "bsRhs_momenta_fc_1", par.sx, par.sy, par.sz);
  bsRhs_momenta_fc[2] = new clbuffer(context, "bsRhs_momenta_fc_2", par.sx, par.sy, par.sz);

  bs_scalars_vc[0][0] = new clbuffer(context, "bs_scalars_vc_0_0", par.sx, par.sy, par.sz);
  bs_scalars_vc[1][0] = new clbuffer(context, "bs_scalars_vc_1_0", par.sx, par.sy, par.sz);
  bs_scalars_vc[2][0] = new clbuffer(context, "bs_scalars_vc_2_0", par.sx, par.sy, par.sz);

  bs_scalars_vc[0][1] = new clbuffer(context, "bs_scalars_vc_0_1", par.sx, par.sy, par.sz);
  bs_scalars_vc[1][1] = new clbuffer(context, "bs_scalars_vc_1_1", par.sx, par.sy, par.sz);
  bs_scalars_vc[2][1] = new clbuffer(context, "bs_scalars_vc_2_1", par.sx, par.sy, par.sz);

  bs_scalars_vc[0][2] = new clbuffer(context, "bs_scalars_vc_0_2", par.sx, par.sy, par.sz); // ice
  bs_scalars_vc[1][2] = new clbuffer(context, "bs_scalars_vc_1_2", par.sx, par.sy, par.sz);
  bs_scalars_vc[2][2] = new clbuffer(context, "bs_scalars_vc_2_2", par.sx, par.sy, par.sz);

  bs_scalars_fc_x[0][0]  = new clbuffer(context, "bs_scalars_fc_x_0_0", par.sx, par.sy, par.sz);
  bs_scalars_fc_x[1][0]  = new clbuffer(context, "bs_scalars_fc_x_1_0", par.sx, par.sy, par.sz);
  bs_scalars_fc_x[2][0]  = new clbuffer(context, "bs_scalars_fc_x_2_0", par.sx, par.sy, par.sz);
  bs_scalars_fc_y[0][0]  = new clbuffer(context, "bs_scalars_fc_y_0_0", par.sx, par.sy, par.sz);
  bs_scalars_fc_y[1][0]  = new clbuffer(context, "bs_scalars_fc_y_1_0", par.sx, par.sy, par.sz);
  bs_scalars_fc_y[2][0]  = new clbuffer(context, "bs_scalars_fc_y_2_0", par.sx, par.sy, par.sz);
  bs_scalars_fc_z[0][0]  = new clbuffer(context, "bs_scalars_fc_z_0_0", par.sx, par.sy, par.sz);
  bs_scalars_fc_z[1][0]  = new clbuffer(context, "bs_scalars_fc_z_1_0", par.sx, par.sy, par.sz);
  bs_scalars_fc_z[2][0]  = new clbuffer(context, "bs_scalars_fc_z_2_0", par.sx, par.sy, par.sz);

  bs_scalars_fc_x[0][1]  = new clbuffer(context, "bs_scalars_fc_x_0_1", par.sx, par.sy, par.sz);
  bs_scalars_fc_x[1][1]  = new clbuffer(context, "bs_scalars_fc_x_1_1", par.sx, par.sy, par.sz);
  bs_scalars_fc_x[2][1]  = new clbuffer(context, "bs_scalars_fc_x_2_1", par.sx, par.sy, par.sz);
  bs_scalars_fc_y[0][1]  = new clbuffer(context, "bs_scalars_fc_y_0_1", par.sx, par.sy, par.sz);
  bs_scalars_fc_y[1][1]  = new clbuffer(context, "bs_scalars_fc_y_1_1", par.sx, par.sy, par.sz);
  bs_scalars_fc_y[2][1]  = new clbuffer(context, "bs_scalars_fc_y_2_1", par.sx, par.sy, par.sz);
  bs_scalars_fc_z[0][1]  = new clbuffer(context, "bs_scalars_fc_z_0_1", par.sx, par.sy, par.sz);
  bs_scalars_fc_z[1][1]  = new clbuffer(context, "bs_scalars_fc_z_1_1", par.sx, par.sy, par.sz);
  bs_scalars_fc_z[2][1]  = new clbuffer(context, "bs_scalars_fc_z_2_1", par.sx, par.sy, par.sz);

  bs_scalars_fc_x[0][2]  = new clbuffer(context, "bs_scalars_fc_x_0_2", par.sx, par.sy, par.sz); // ice
  bs_scalars_fc_x[1][2]  = new clbuffer(context, "bs_scalars_fc_x_1_2", par.sx, par.sy, par.sz);
  bs_scalars_fc_x[2][2]  = new clbuffer(context, "bs_scalars_fc_x_2_2", par.sx, par.sy, par.sz);
  bs_scalars_fc_y[0][2]  = new clbuffer(context, "bs_scalars_fc_y_0_2", par.sx, par.sy, par.sz);
  bs_scalars_fc_y[1][2]  = new clbuffer(context, "bs_scalars_fc_y_1_2", par.sx, par.sy, par.sz);
  bs_scalars_fc_y[2][2]  = new clbuffer(context, "bs_scalars_fc_y_2_2", par.sx, par.sy, par.sz);
  bs_scalars_fc_z[0][2]  = new clbuffer(context, "bs_scalars_fc_z_0_2", par.sx, par.sy, par.sz);
  bs_scalars_fc_z[1][2]  = new clbuffer(context, "bs_scalars_fc_z_1_2", par.sx, par.sy, par.sz);
  bs_scalars_fc_z[2][2]  = new clbuffer(context, "bs_scalars_fc_z_2_2", par.sx, par.sy, par.sz);

  // ----------------------------------------------------------------- //
  // kernel //
  // ----------------------------------------------------------------- //

  k_init_scalars          = new clkernel(context, par, "./kernels/k_init_scalars_dycoms.cl");
  k_init_momenta          = new clkernel(context, par, "./kernels/k_init_momenta.cl");
  ks_ext_forcings         = new clkernel(context, par, "./kernels/k_ext_forcings_isdac.cl");
  ks_copy[0][0]           = new clkernel(context, par, "./kernels/k_copy_one.cl");
  ks_copy[0][1]           = new clkernel(context, par, "./kernels/k_copy_one.cl");
  ks_copy[0][2]           = new clkernel(context, par, "./kernels/k_copy_one.cl");
  ks_copy[0][3]           = new clkernel(context, par, "./kernels/k_copy_one.cl");
  ks_copy[1][0]           = new clkernel(context, par, "./kernels/k_copy_one.cl");
  ks_copy[1][1]           = new clkernel(context, par, "./kernels/k_copy_one.cl");
  ks_copy[1][2]           = new clkernel(context, par, "./kernels/k_copy_one.cl");
  ks_copy[1][3]           = new clkernel(context, par, "./kernels/k_copy_one.cl");
  ks_copy[2][0]           = new clkernel(context, par, "./kernels/k_copy_one.cl");
  ks_copy[2][1]           = new clkernel(context, par, "./kernels/k_copy_one.cl");
  ks_copy[2][2]           = new clkernel(context, par, "./kernels/k_copy_one.cl");
  ks_copy[2][3]           = new clkernel(context, par, "./kernels/k_copy_one.cl");
  ks_f2c                  = new clkernel(context, par, "./kernels/ks_facetocell.cl");
  ks_c2f[0]               = new clkernel(context, par, "./kernels/ks_celltoface.cl");
  ks_c2f[1]               = new clkernel(context, par, "./kernels/ks_celltoface.cl");
  ks_c2f[2]               = new clkernel(context, par, "./kernels/ks_celltoface.cl");
  ks_adv_momenta          = new clkernel(context, par, "./kernels/ks_adv_momenta.cl");
  ks_adv_scalars[0][0]    = new clkernel(context, par, "./kernels/ks_adv_scalars.cl");
  ks_adv_scalars[0][1]    = new clkernel(context, par, "./kernels/ks_adv_scalars.cl");
  ks_adv_scalars[0][2]    = new clkernel(context, par, "./kernels/ks_adv_scalars.cl");
  ks_adv_scalars[1][0]    = new clkernel(context, par, "./kernels/ks_adv_scalars.cl");
  ks_adv_scalars[1][1]    = new clkernel(context, par, "./kernels/ks_adv_scalars.cl");
  ks_adv_scalars[1][2]    = new clkernel(context, par, "./kernels/ks_adv_scalars.cl");
  ks_adv_scalars[2][0]    = new clkernel(context, par, "./kernels/ks_adv_scalars.cl");
  ks_adv_scalars[2][1]    = new clkernel(context, par, "./kernels/ks_adv_scalars.cl");
  ks_adv_scalars[2][2]    = new clkernel(context, par, "./kernels/ks_adv_scalars.cl");
  ks_prep_fast[0][0]      = new clkernel(context, par, "./kernels/ks_prep_fast_0.cl");
  ks_prep_fast[0][1]      = new clkernel(context, par, "./kernels/ks_prep_fast_0.cl");
  ks_prep_fast[0][2]      = new clkernel(context, par, "./kernels/ks_prep_fast_0.cl");
  ks_prep_fast[0][3]      = new clkernel(context, par, "./kernels/ks_prep_fast_0.cl");
  ks_prep_fast[1][0]      = new clkernel(context, par, "./kernels/ks_prep_fast_1.cl");
  ks_prep_fast[1][1]      = new clkernel(context, par, "./kernels/ks_prep_fast_1.cl");
  ks_prep_fast[1][2]      = new clkernel(context, par, "./kernels/ks_prep_fast_1.cl");
  ks_prep_fast[1][3]      = new clkernel(context, par, "./kernels/ks_prep_fast_1.cl");
  ks_prep_fast[2][0]      = new clkernel(context, par, "./kernels/ks_prep_fast_2.cl");
  ks_prep_fast[2][1]      = new clkernel(context, par, "./kernels/ks_prep_fast_2.cl");
  ks_prep_fast[2][2]      = new clkernel(context, par, "./kernels/ks_prep_fast_2.cl");
  ks_prep_fast[2][3]      = new clkernel(context, par, "./kernels/ks_prep_fast_2.cl");
  kf_pressure             = new clkernel(context, par, "./kernels/kf_pressure.cl");
  kf_microphys            = new clkernel(context, par, "./kernels/kf_microphys.cl");
  kf_step_momenta[0]      = new clkernel(context, par, "./kernels/kf_step_momenta_0.cl");
  kf_step_momenta[1]      = new clkernel(context, par, "./kernels/kf_step_momenta_1.cl");
  kf_step_momenta[2]      = new clkernel(context, par, "./kernels/kf_step_momenta_2.cl");
  kf_step_scalars[0][0]   = new clkernel(context, par, "./kernels/kf_step_scalars_0.cl");
  kf_step_scalars[0][1]   = new clkernel(context, par, "./kernels/kf_step_scalars_0.cl");
  kf_step_scalars[0][2]   = new clkernel(context, par, "./kernels/kf_step_scalars_0.cl");
  kf_step_scalars[1][0]   = new clkernel(context, par, "./kernels/kf_step_scalars_1.cl");
  kf_step_scalars[1][1]   = new clkernel(context, par, "./kernels/kf_step_scalars_1.cl");
  kf_step_scalars[1][2]   = new clkernel(context, par, "./kernels/kf_step_scalars_1.cl");
  kf_step_scalars[2][0]   = new clkernel(context, par, "./kernels/kf_step_scalars_2.cl");
  kf_step_scalars[2][1]   = new clkernel(context, par, "./kernels/kf_step_scalars_2.cl");
  kf_step_scalars[2][2]   = new clkernel(context, par, "./kernels/kf_step_scalars_2.cl");
  kf_copy[0]              = new clkernel(context, par, "./kernels/k_copy_one.cl");
  kf_copy[1]              = new clkernel(context, par, "./kernels/k_copy_one.cl");
  kf_copy[2]              = new clkernel(context, par, "./kernels/k_copy_one.cl");
  kf_copy[3]              = new clkernel(context, par, "./kernels/k_copy_one.cl");
  k_damping               = new clkernel(context, par, "./kernels/k_damping.cl");
  k_nesting               = new clkernel(context, par, "./kernels/k_nesting_isdac.cl");
  k_clone[0]              = new clkernel(context, par, "./kernels/k_clone.cl");
  k_clone[1]              = new clkernel(context, par, "./kernels/k_clone.cl");
  k_clone[2]              = new clkernel(context, par, "./kernels/k_clone.cl");
  k_clone[3]              = new clkernel(context, par, "./kernels/k_clone.cl");

  // ----------------------------------------------------------------- //
  // bindings //
  // ----------------------------------------------------------------- //

  // init
  k_init_scalars->bind("bf_scalars_vc_a_0", bf_scalars_vc_a[0]);
  k_init_scalars->bind("bf_scalars_vc_a_1", bf_scalars_vc_a[1]);
  k_init_scalars->bind("bf_scalars_vc_a_2", bf_scalars_vc_a[2]); // ice
  k_init_scalars->bind("bf_scalars_vc_b_0", bf_scalars_vc_b[0]);
  k_init_scalars->bind("bf_scalars_vc_b_1", bf_scalars_vc_b[1]);
  k_init_scalars->bind("bf_scalars_vc_b_2", bf_scalars_vc_b[2]); // ice

  k_init_momenta->bind("bf_scalars_vc_a_0", bf_scalars_vc_a[0]);
  k_init_momenta->bind("bf_momenta_fc_a", bf_momenta_fc_a);

  k_damping->bind("b_source_momenta",   bf_momenta_fc_a);
  k_damping->bind("b_target_momenta",   bf_momenta_fc_b);

  k_nesting->bind("b_source_scalars_0", bf_scalars_vc_a[0]);
  k_nesting->bind("b_source_scalars_1", bf_scalars_vc_a[1]);
  k_nesting->bind("b_source_scalars_2", bf_scalars_vc_a[2]);
  k_nesting->bind("b_target_scalars_0", bf_scalars_vc_b[0]);
  k_nesting->bind("b_target_scalars_1", bf_scalars_vc_b[1]);
  k_nesting->bind("b_target_scalars_2", bf_scalars_vc_b[2]);
  k_nesting->bind("b_source_momenta",   bf_momenta_fc_a);
  k_nesting->bind("b_target_momenta",   bf_momenta_fc_b);

  // slow
  ks_ext_forcings->bind("b_source_scalars_0", bf_scalars_vc_a[0]);
  ks_ext_forcings->bind("b_source_scalars_1", bf_scalars_vc_a[1]);
  ks_ext_forcings->bind("b_source_scalars_2", bf_scalars_vc_a[2]); // ice
  ks_ext_forcings->bind("b_source_momenta",   bf_momenta_fc_a);
  ks_ext_forcings->bind("b_target_scalars_0", bf_scalars_vc_b[0]);
  ks_ext_forcings->bind("b_target_scalars_1", bf_scalars_vc_b[1]);
  ks_ext_forcings->bind("b_target_scalars_2", bf_scalars_vc_b[2]); // ice
  ks_ext_forcings->bind("b_target_momenta",   bf_momenta_fc_b);

  // copy scalars [stages][fields]
  for (int s = 0; s < 3; s++) {
    for (int f = 0; f < 3; f++) {
      ks_copy[s][f]->bind("b_source", bf_scalars_vc_a[f]);
      ks_copy[s][f]->bind("b_target", bs_scalars_vc[s][f]);
    }
  }

  // copy momenta [stages][fields]
  for (int s = 0; s < 3; s++) {
    ks_copy[s][3]->bind("b_source", bf_momenta_fc_a);
    ks_copy[s][3]->bind("b_target", bs_momenta_fc[s]);
  }

  ks_f2c->bind("bf_momenta_fc_a",  bf_momenta_fc_a);
  ks_f2c->bind("b_m_vc_uv",        b_temp[0]);
  ks_f2c->bind("b_m_vc_w" ,        b_temp[1]);

  ks_adv_momenta->bind("bf_scalars_vc_a_0", bf_scalars_vc_a[0]);
  ks_adv_momenta->bind("bf_momenta_fc_a",   bf_momenta_fc_a);
  ks_adv_momenta->bind("b_m_vc_uv",         b_temp[0]);
  ks_adv_momenta->bind("b_m_vc_w",          b_temp[1]);
  ks_adv_momenta->bind("bRhs_m_vc_uv",      b_temp[2]);
  ks_adv_momenta->bind("bRhs_m_vc_w",       b_temp[3]);

  ks_c2f[0]->bind("bRhs_m_vc_uv",       b_temp[2]);
  ks_c2f[0]->bind("bRhs_m_vc_w",        b_temp[3]);
  ks_c2f[0]->bind("bsRhs_momenta_fc_s", bsRhs_momenta_fc[0]);

  ks_c2f[1]->bind("bRhs_m_vc_uv",       b_temp[2]);
  ks_c2f[1]->bind("bRhs_m_vc_w",        b_temp[3]);
  ks_c2f[1]->bind("bsRhs_momenta_fc_s", bsRhs_momenta_fc[1]);

  ks_c2f[2]->bind("bRhs_m_vc_uv",       b_temp[2]);
  ks_c2f[2]->bind("bRhs_m_vc_w",        b_temp[3]);
  ks_c2f[2]->bind("bsRhs_momenta_fc_s", bsRhs_momenta_fc[2]);

  // slow advection [stages][fields]
  for (int s = 0; s < 3; s++) {
    for (int f = 0; f < 3; f++) {
      ks_adv_scalars[s][f]->bind("bf_density_vc_a",   bf_scalars_vc_a[0]);
      ks_adv_scalars[s][f]->bind("bf_scalars_vc_a",   bf_scalars_vc_a[f]);
      ks_adv_scalars[s][f]->bind("bf_momenta_fc_a",   bf_momenta_fc_a);
      ks_adv_scalars[s][f]->bind("bs_scalars_fc_x_s", bs_scalars_fc_x[s][f]);
      ks_adv_scalars[s][f]->bind("bs_scalars_fc_y_s", bs_scalars_fc_y[s][f]);
      ks_adv_scalars[s][f]->bind("bs_scalars_fc_z_s", bs_scalars_fc_z[s][f]);
    }
  }

  // prep fast stage 0
  for (int f = 0; f < 3; f++) {
    ks_prep_fast[0][f]->bind("bs_source_0", bs_scalars_vc[0][f]);
    ks_prep_fast[0][f]->bind("bf_target_a", bf_scalars_vc_a[f]);
  }
  ks_prep_fast[0][3]->bind("bs_source_0", bs_momenta_fc[0]);
  ks_prep_fast[0][3]->bind("bf_target_a", bf_momenta_fc_a);

  // prep fast stage 1
  for (int f = 0; f < 3; f++) {
    ks_prep_fast[1][f]->bind("bs_source_0", bs_scalars_vc[0][f]);
    ks_prep_fast[1][f]->bind("bs_source_1", bs_scalars_vc[1][f]);
    ks_prep_fast[1][f]->bind("bf_target_a", bf_scalars_vc_a[f]);
  }
  ks_prep_fast[1][3]->bind("bs_source_0", bs_momenta_fc[0]);
  ks_prep_fast[1][3]->bind("bs_source_1", bs_momenta_fc[1]);
  ks_prep_fast[1][3]->bind("bf_target_a", bf_momenta_fc_a);

  // prep fast stage 2
  for (int f = 0; f < 3; f++) {
    ks_prep_fast[2][f]->bind("bs_source_0", bs_scalars_vc[0][f]);
    ks_prep_fast[2][f]->bind("bs_source_1", bs_scalars_vc[1][f]);
    ks_prep_fast[2][f]->bind("bs_source_2", bs_scalars_vc[2][f]);
    ks_prep_fast[2][f]->bind("bf_target_a", bf_scalars_vc_a[f]);
  }
  ks_prep_fast[2][3]->bind("bs_source_0", bs_momenta_fc[0]);
  ks_prep_fast[2][3]->bind("bs_source_1", bs_momenta_fc[1]);
  ks_prep_fast[2][3]->bind("bs_source_2", bs_momenta_fc[2]);
  ks_prep_fast[2][3]->bind("bf_target_a", bf_momenta_fc_a);

  // fast
  kf_pressure->bind("bf_scalars_vc_a_0", bf_scalars_vc_a[0]);
  kf_pressure->bind("bf_scalars_vc_a_1", bf_scalars_vc_a[1]);
  kf_pressure->bind("bf_scalars_vc_a_2", bf_scalars_vc_a[2]); // ice
  kf_pressure->bind("bRhs_p_fc",         b_temp[0]);

  kf_step_momenta[0]->bind("bf_momenta_fc_a",    bf_momenta_fc_a);
  kf_step_momenta[0]->bind("bsRhs_momenta_fc_0", bsRhs_momenta_fc[0]);
  kf_step_momenta[0]->bind("bRhs_p_fc",          b_temp[0]);
  kf_step_momenta[0]->bind("bf_momenta_fc_b",    bf_momenta_fc_b);

  kf_microphys->bind("bf_scalars_vc_a_0", bf_scalars_vc_a[0]);
  kf_microphys->bind("bf_scalars_vc_a_1", bf_scalars_vc_a[1]);
  kf_microphys->bind("bf_scalars_vc_a_2", bf_scalars_vc_a[2]); // ice
  kf_microphys->bind("bf_momenta_fc_b",   bf_momenta_fc_b);
  kf_microphys->bind("bRhs_mp_vc_0",      b_temp[0]);
  kf_microphys->bind("bRhs_mp_vc_1",      b_temp[1]);
  kf_microphys->bind("bRhs_mp_vc_2",      b_temp[2]); // ice
  kf_microphys->bind("phys",              (unsigned int)(1023)); // enable all per default

  for (int f = 0; f < 3; f++) {
    kf_step_scalars[0][f]->bind("bf_scalars_vc_a",     bf_scalars_vc_a[f]);
    kf_step_scalars[0][f]->bind("bf_momenta_fc_b",     bf_momenta_fc_b);
    kf_step_scalars[0][f]->bind("bs_scalars_fc_x_0",   bs_scalars_fc_x[0][f]);
    kf_step_scalars[0][f]->bind("bs_scalars_fc_y_0",   bs_scalars_fc_y[0][f]);
    kf_step_scalars[0][f]->bind("bs_scalars_fc_z_0",   bs_scalars_fc_z[0][f]);
    kf_step_scalars[0][f]->bind("bRhs_mp_vc",          b_temp[f]);
    kf_step_scalars[0][f]->bind("bf_scalars_vc_b",     bf_scalars_vc_b[f]);
  }

  kf_copy[0]->bind("b_source", bf_scalars_vc_b[0]);
  kf_copy[0]->bind("b_target", bf_scalars_vc_a[0]);
  kf_copy[1]->bind("b_source", bf_scalars_vc_b[1]);
  kf_copy[1]->bind("b_target", bf_scalars_vc_a[1]);
  kf_copy[2]->bind("b_source", bf_scalars_vc_b[2]);
  kf_copy[2]->bind("b_target", bf_scalars_vc_a[2]);
  kf_copy[3]->bind("b_source", bf_momenta_fc_b);
  kf_copy[3]->bind("b_target", bf_momenta_fc_a);

  k_clone[0]->bind("b_source", bf_scalars_vc_a[0]);
  k_clone[0]->bind("b_target", bf_scalars_vc_b[0]);
  k_clone[1]->bind("b_source", bf_scalars_vc_a[1]);
  k_clone[1]->bind("b_target", bf_scalars_vc_b[1]);
  k_clone[2]->bind("b_source", bf_scalars_vc_a[2]);
  k_clone[2]->bind("b_target", bf_scalars_vc_b[2]);
  k_clone[3]->bind("b_source", bf_momenta_fc_a);
  k_clone[3]->bind("b_target", bf_momenta_fc_b);

  kf_step_momenta[1]->bind("bf_momenta_fc_a",    bf_momenta_fc_a);
  kf_step_momenta[1]->bind("bsRhs_momenta_fc_0", bsRhs_momenta_fc[0]);
  kf_step_momenta[1]->bind("bsRhs_momenta_fc_1", bsRhs_momenta_fc[1]);
  kf_step_momenta[1]->bind("bs_momenta_fc_0",    bs_momenta_fc[0]);
  kf_step_momenta[1]->bind("bs_momenta_fc_1",    bs_momenta_fc[1]);
  kf_step_momenta[1]->bind("bRhs_p_fc",          b_temp[0]);
  kf_step_momenta[1]->bind("bf_momenta_fc_b",    bf_momenta_fc_b);

  for (int f = 0; f < 3; f++) {
    kf_step_scalars[1][f]->bind("bf_scalars_vc_a",     bf_scalars_vc_a[f]);
    kf_step_scalars[1][f]->bind("bf_momenta_fc_b",     bf_momenta_fc_b);
    kf_step_scalars[1][f]->bind("bs_scalars_fc_x_0",   bs_scalars_fc_x[0][f]);
    kf_step_scalars[1][f]->bind("bs_scalars_fc_x_1",   bs_scalars_fc_x[1][f]);
    kf_step_scalars[1][f]->bind("bs_scalars_fc_y_0",   bs_scalars_fc_y[0][f]);
    kf_step_scalars[1][f]->bind("bs_scalars_fc_y_1",   bs_scalars_fc_y[1][f]);
    kf_step_scalars[1][f]->bind("bs_scalars_fc_z_0",   bs_scalars_fc_z[0][f]);
    kf_step_scalars[1][f]->bind("bs_scalars_fc_z_1",   bs_scalars_fc_z[1][f]);
    kf_step_scalars[1][f]->bind("bs_scalars_vc_0",     bs_scalars_vc[0][f]);
    kf_step_scalars[1][f]->bind("bs_scalars_vc_1",     bs_scalars_vc[1][f]);
    kf_step_scalars[1][f]->bind("bRhs_mp_vc",          b_temp[f]);
    kf_step_scalars[1][f]->bind("bf_scalars_vc_b",     bf_scalars_vc_b[f]);
  }

  kf_step_momenta[2]->bind("bf_momenta_fc_a",    bf_momenta_fc_a);
  kf_step_momenta[2]->bind("bsRhs_momenta_fc_0", bsRhs_momenta_fc[0]);
  kf_step_momenta[2]->bind("bsRhs_momenta_fc_1", bsRhs_momenta_fc[1]);
  kf_step_momenta[2]->bind("bsRhs_momenta_fc_2", bsRhs_momenta_fc[2]);
  kf_step_momenta[2]->bind("bs_momenta_fc_0",    bs_momenta_fc[0]);
  kf_step_momenta[2]->bind("bs_momenta_fc_1",    bs_momenta_fc[1]);
  kf_step_momenta[2]->bind("bs_momenta_fc_2",    bs_momenta_fc[2]);
  kf_step_momenta[2]->bind("bRhs_p_fc",          b_temp[0]);
  kf_step_momenta[2]->bind("bf_momenta_fc_b",    bf_momenta_fc_b);

  for (int f = 0; f < 3; f++) {
    kf_step_scalars[2][f]->bind("bf_scalars_vc_a",     bf_scalars_vc_a[f]);
    kf_step_scalars[2][f]->bind("bf_momenta_fc_b",     bf_momenta_fc_b);
    kf_step_scalars[2][f]->bind("bs_scalars_fc_x_0",   bs_scalars_fc_x[0][f]);
    kf_step_scalars[2][f]->bind("bs_scalars_fc_x_1",   bs_scalars_fc_x[1][f]);
    kf_step_scalars[2][f]->bind("bs_scalars_fc_x_2",   bs_scalars_fc_x[2][f]);
    kf_step_scalars[2][f]->bind("bs_scalars_fc_y_0",   bs_scalars_fc_y[0][f]);
    kf_step_scalars[2][f]->bind("bs_scalars_fc_y_1",   bs_scalars_fc_y[1][f]);
    kf_step_scalars[2][f]->bind("bs_scalars_fc_y_2",   bs_scalars_fc_y[2][f]);
    kf_step_scalars[2][f]->bind("bs_scalars_fc_z_0",   bs_scalars_fc_z[0][f]);
    kf_step_scalars[2][f]->bind("bs_scalars_fc_z_1",   bs_scalars_fc_z[1][f]);
    kf_step_scalars[2][f]->bind("bs_scalars_fc_z_2",   bs_scalars_fc_z[2][f]);
    kf_step_scalars[2][f]->bind("bs_scalars_vc_0",     bs_scalars_vc[0][f]);
    kf_step_scalars[2][f]->bind("bs_scalars_vc_1",     bs_scalars_vc[1][f]);
    kf_step_scalars[2][f]->bind("bs_scalars_vc_2",     bs_scalars_vc[2][f]);
    kf_step_scalars[2][f]->bind("bRhs_mp_vc",          b_temp[f]);
    kf_step_scalars[2][f]->bind("bf_scalars_vc_b",     bf_scalars_vc_b[f]);
  }

  // ----------------------------------------------------------------- //
  // exporter //
  // ----------------------------------------------------------------- //

  // exporter[0] = new clexport(context, "YZ", "./kernels/exporter/ke_render.cl", par, bf_momenta_fc_a, bf_scalars_vc_a[0], bf_scalars_vc_a[1], 0, par.sx/2);
  // exporter[1] = new clexport(context, "XZ", "./kernels/exporter/ke_render.cl", par, bf_momenta_fc_a, bf_scalars_vc_a[0], bf_scalars_vc_a[1], 1, par.sy/2);
  // exporter[2] = new clexport(context, "XY", "./kernels/exporter/ke_render.cl", par, bf_momenta_fc_a, bf_scalars_vc_a[0], bf_scalars_vc_a[1], 2, par.sz/2);

  v_exporter.push_back(new clexport(context, "XZ_w",     "./kernels/exporter/ke_w.cl",     par, bf_momenta_fc_a, bf_scalars_vc_a[0], bf_scalars_vc_a[1], bf_scalars_vc_a[2], 1, par.sy/2,  1,1,0  ));

  v_exporter.push_back(new clexport(context, "XZ_rho_c", "./kernels/exporter/ke_rho_c.cl", par, bf_momenta_fc_a, bf_scalars_vc_a[0], bf_scalars_vc_a[1], bf_scalars_vc_a[2], 1, par.sy/2,  1,1,0  ));
  v_exporter.push_back(new clexport(context, "XZ_rho_r", "./kernels/exporter/ke_rho_r.cl", par, bf_momenta_fc_a, bf_scalars_vc_a[0], bf_scalars_vc_a[1], bf_scalars_vc_a[2], 1, par.sy/2,  1,1,0  ));
  v_exporter.push_back(new clexport(context, "XZ_rho_i", "./kernels/exporter/ke_rho_i.cl", par, bf_momenta_fc_a, bf_scalars_vc_a[0], bf_scalars_vc_a[1], bf_scalars_vc_a[2], 1, par.sy/2,  1,1,1  ));
  v_exporter.push_back(new clexport(context, "XZ_rho_s", "./kernels/exporter/ke_rho_s.cl", par, bf_momenta_fc_a, bf_scalars_vc_a[0], bf_scalars_vc_a[1], bf_scalars_vc_a[2], 1, par.sy/2,  1,1,1  ));

  v_exporter.push_back(new clexport(context, "XZ_n_c",   "./kernels/exporter/ke_n_c.cl",   par, bf_momenta_fc_a, bf_scalars_vc_a[0], bf_scalars_vc_a[1], bf_scalars_vc_a[2], 1, par.sy/2,  1,1,0  ));
  v_exporter.push_back(new clexport(context, "XZ_n_r",   "./kernels/exporter/ke_n_r.cl",   par, bf_momenta_fc_a, bf_scalars_vc_a[0], bf_scalars_vc_a[1], bf_scalars_vc_a[2], 1, par.sy/2,  1,1,0  ));
  v_exporter.push_back(new clexport(context, "XZ_n_d",   "./kernels/exporter/ke_n_d.cl",   par, bf_momenta_fc_a, bf_scalars_vc_a[0], bf_scalars_vc_a[1], bf_scalars_vc_a[2], 1, par.sy/2,  1,1,0  ));

  v_exporter.push_back(new clexport(context, "XZ_n_i",   "./kernels/exporter/ke_n_i.cl",   par, bf_momenta_fc_a, bf_scalars_vc_a[0], bf_scalars_vc_a[1], bf_scalars_vc_a[2], 1, par.sy/2,  1,1,1  ));
  v_exporter.push_back(new clexport(context, "XZ_n_s",   "./kernels/exporter/ke_n_s.cl",   par, bf_momenta_fc_a, bf_scalars_vc_a[0], bf_scalars_vc_a[1], bf_scalars_vc_a[2], 1, par.sy/2,  1,1,1  ));

  v_exporter.push_back(new clexport(context, "XY_lwp",   "./kernels/exporter/ke_int_lwp.cl",   par, bf_momenta_fc_a, bf_scalars_vc_a[0], bf_scalars_vc_a[1], bf_scalars_vc_a[2], 2, par.sz/2,  1,1,0  ));


  // vertical profiles
  v_exporter.push_back(new clexport(context, "VP_ql",    "./kernels/exporter/ke_int_ql.cl",par, bf_momenta_fc_a, bf_scalars_vc_a[0], bf_scalars_vc_a[1], bf_scalars_vc_a[2], 1, par.sy/2,  0,0,1  ));
  v_exporter.push_back(new clexport(context, "VP_qt",    "./kernels/exporter/ke_int_qt.cl",par, bf_momenta_fc_a, bf_scalars_vc_a[0], bf_scalars_vc_a[1], bf_scalars_vc_a[2], 1, par.sy/2,  0,0,1  ));
  v_exporter.push_back(new clexport(context, "VP_qi",    "./kernels/exporter/ke_int_qi.cl",par, bf_momenta_fc_a, bf_scalars_vc_a[0], bf_scalars_vc_a[1], bf_scalars_vc_a[2], 1, par.sy/2,  0,0,1  ));
  v_exporter.push_back(new clexport(context, "VP_thetal","./kernels/exporter/ke_int_thetal.cl",par, bf_momenta_fc_a, bf_scalars_vc_a[0], bf_scalars_vc_a[1], bf_scalars_vc_a[2], 1, par.sy/2,  0,0,1  ));
  v_exporter.push_back(new clexport(context, "VP_temperature","./kernels/exporter/ke_int_temperature.cl",par, bf_momenta_fc_a, bf_scalars_vc_a[0], bf_scalars_vc_a[1], bf_scalars_vc_a[2], 1, par.sy/2,  0,0,1  ));
  v_exporter.push_back(new clexport(context, "VP_verlocity_variance","./kernels/exporter/ke_int_velocity_variance.cl",par, bf_momenta_fc_a, bf_scalars_vc_a[0], bf_scalars_vc_a[1], bf_scalars_vc_a[2], 1, par.sy/2,  0,0,1  ));

}

asystem::~asystem() {
}

void asystem::read_single_profile(clbuffer *b, FILE *f, int c, float dz, float offset) {
  char buf[255];
  char *back;
  int numlvl,dummy,k;
  float *prof, *height;

  back=fgets(buf,255,f);
  sscanf(buf,"%d", &numlvl);
  prof  =new float[numlvl];
  height=new float[numlvl];
  back=fgets(buf,255,f);
  for (k=0; k<numlvl; k++) {
    back=fgets(buf,255,f);
    sscanf(buf,"%f %f",&height[k],&prof[k]);
  }
  for (int z=0; z<par.sz; z++) {
    k=0;
    while (height[k]<(z+0.5)*dz)k++;
    // printf("%d %d %d %f %f\n", c, z, k, height[k], prof[k] );
    for (int x=0; x<par.sx; x++)
      for (int y=0; y<par.sy; y++)
        b->set(x, y, z, c, prof[k]+offset);
  }
  delete prof;
  delete height;
}

void asystem::write_files(int step) {
  for (int e = 0; e<v_exporter.size(); e++) {
    v_exporter.at(e)->write_files(step);
  }
}

void asystem::write_state(std::string s_filename) {
  logger->log(0, "Writing states to %s\n", s_filename.c_str());
  int ret = system(("mkdir -p " + s_filename).c_str());
  bf_scalars_vc_a[0]->write_raw(s_filename+"/s0.dat");
  bf_scalars_vc_a[1]->write_raw(s_filename+"/s1.dat");
  bf_scalars_vc_a[2]->write_raw(s_filename+"/s2.dat");
  bf_momenta_fc_a   ->write_raw(s_filename+"/m0.dat");
}

void asystem::read_state(std::string s_filename) {
  logger->log(0, "Reading states from %s\n", s_filename.c_str());
  bf_scalars_vc_a[0]->read_raw(s_filename+"/s0.dat");
  bf_scalars_vc_a[1]->read_raw(s_filename+"/s1.dat");
  bf_scalars_vc_a[2]->read_raw(s_filename+"/s2.dat");
  bf_momenta_fc_a   ->read_raw(s_filename+"/m0.dat");
  bf_scalars_vc_a[0]->ram2device();
  bf_scalars_vc_a[1]->ram2device();
  bf_scalars_vc_a[2]->ram2device();
  bf_momenta_fc_a   ->ram2device();

  write_files(frame_index);
  frame_index+=1;
}

void asystem::init_from_kernel() {
  logger->log(0, "Init from kernel\n");
  k_init_scalars->step(par.sx, par.sy, par.sz);
  k_init_momenta->step(par.sx, par.sy, par.sz);
}

void asystem::init_from_file(std::string s_filePath) {
  logger->log(0, "Init from file %s\n", s_filePath.c_str());

  char buf[255];
  FILE *f;
  f=fopen(s_filePath.c_str(),"r");

  while (fgets(buf,255,f)!=NULL) {
    if (strstr(buf,"ThMProf")!=NULL)  read_single_profile(bf_scalars_vc_b[0], f, 0, par.dz);
    if (strstr(buf,"RhoProf")!=NULL)  read_single_profile(bf_scalars_vc_b[0], f, 1, par.dz);
    if (strstr(buf,"RhoVProf")!=NULL) read_single_profile(bf_scalars_vc_b[0], f, 2, par.dz);
    if (strstr(buf,"RhoLProf")!=NULL) read_single_profile(bf_scalars_vc_b[0], f, 3, par.dz);
    if (strstr(buf,"uProf")!=NULL)    read_single_profile(bf_momenta_fc_a, f, 0, par.dz);
    if (strstr(buf,"vProf")!=NULL)    read_single_profile(bf_momenta_fc_a, f, 1, par.dz);
  }
  fclose(f);
  bf_scalars_vc_b[0]->ram2device();
  bf_momenta_fc_a->ram2device();

  k_init_scalars->step(par.sx, par.sy, par.sz);
  k_init_momenta->step(par.sx, par.sy, par.sz);
}

void asystem::equilibrate() {
  // custom mis step to create initial profile with hydrostaic equilibrium
  int kx = 1;
  int ky = 1;
  int kz = par.sz;

  kf_microphys->bind("phys", (unsigned int)(2)); // only enable condensation

  int equil_steps = 5000;

  for (int f = 0; f < equil_steps; f++) {
    if (f%10==0) logger->log(2,"\rEquilibrating  -  %d/%d",f,equil_steps);
    if (par.timescheme == 0) {
      for (int s=0; s<3; s++) {
        slow_stage(0, kx, ky, kz);
        for (int i=0; i<par.nsi[s]; i++) fast_stage(0, kx, ky, kz);
      }
    } else if (par.timescheme == 1) {
      for (int s=0; s<3; s++) {
        slow_stage(s, kx, ky, kz);
        for (int i=0; i<par.nsi[s]; i++) {
          fast_stage(s, kx, ky, kz);
          k_nesting->bind("damping_strength", frame_index);
          k_nesting->step(kx, ky, kz);
          kf_copy[0]->step(kx, ky, kz);
          kf_copy[1]->step(kx, ky, kz);
          kf_copy[2]->step(kx, ky, kz);
          kf_copy[3]->step(kx, ky, kz);
        }
        // strip to full
        k_clone[0]->step(par.sx, par.sy, par.sz);
        k_clone[1]->step(par.sx, par.sy, par.sz);
        k_clone[2]->step(par.sx, par.sy, par.sz);
        k_clone[3]->step(par.sx, par.sy, par.sz);
        kf_copy[0]->step(par.sx, par.sy, par.sz);
        kf_copy[1]->step(par.sx, par.sy, par.sz);
        kf_copy[2]->step(par.sx, par.sy, par.sz);
        kf_copy[3]->step(par.sx, par.sy, par.sz);

        // write_files(frame_index*3+s);
      }
    }
    write_files(frame_index);
    frame_index += 1;
  }
  kf_microphys->bind("phys", (unsigned int)(1023));
  logger->log(2,"\rEquilibrating  -  done\n");
}

void asystem::mis_step() {
  // overload for default arguments
  mis_step(0, par.sx, par.sy, par.sz);
}

void asystem::mis_step(int damping, int kx, int ky, int kz) {

  ks_ext_forcings->bind("frame_index", frame_index);
  ks_ext_forcings->step(kx, ky, kz);
  kf_copy[0]->step(kx, ky, kz);
  kf_copy[1]->step(kx, ky, kz);
  kf_copy[2]->step(kx, ky, kz);
  kf_copy[3]->step(kx, ky, kz);

  if (par.timescheme == 0) {
    for (int s=0; s<3; s++) {
      slow_stage(0, kx, ky, kz);
      for (int i=0; i<par.nsi[s]; i++) fast_stage(0, kx, ky, kz);
    }
  } else if (par.timescheme == 1) {
    for (int s=0; s<3; s++) {
      slow_stage(s, kx, ky, kz);
      for (int i=0; i<par.nsi[s]; i++) {
        fast_stage(s, kx, ky, kz);
        // k_damping->bind("damping_strength", (unsigned int)frame_index);
        // k_damping->step(kx, ky, kz);
        // kf_copy[0]->step(kx, ky, kz);
        // kf_copy[1]->step(kx, ky, kz);
        // kf_copy[2]->step(kx, ky, kz);
        // kf_copy[3]->step(kx, ky, kz);
      }
    }
  }

  write_files(frame_index);
  if (frame_index*par.dT == 3600.0) {
    int hour = (int)(frame_index*par.dT/3600.0);
    write_state("./snapshots/"+std::to_string(hour)+"h");
  }
  frame_index += 1;
}

void asystem::slow_stage(int s, int kx, int ky, int kz) {

  // stage init
  ks_copy[s][0]->step(kx, ky, kz);
  ks_copy[s][1]->step(kx, ky, kz);
  ks_copy[s][2]->step(kx, ky, kz);
  ks_copy[s][3]->step(kx, ky, kz);

  // advection of momenta
  ks_f2c->step(kx, ky, kz);
  ks_adv_momenta->step(kx, ky, kz);
  ks_c2f[s]->step(kx, ky, kz);

  // advection of scalars
  ks_adv_scalars[s][0]->step(kx, ky, kz);
  ks_adv_scalars[s][1]->step(kx, ky, kz);
  ks_adv_scalars[s][2]->step(kx, ky, kz);
  ks_prep_fast[s][0]->step(kx, ky, kz);
  ks_prep_fast[s][1]->step(kx, ky, kz);
  ks_prep_fast[s][2]->step(kx, ky, kz);
  ks_prep_fast[s][3]->step(kx, ky, kz);
}

void asystem::fast_stage(int s, int kx, int ky, int kz) {
  kf_pressure->step(kx, ky, kz);
  kf_step_momenta[s]->step(kx, ky, kz);
  kf_microphys->step(kx, ky, kz);
  kf_step_scalars[s][0]->step(kx, ky, kz);
  kf_step_scalars[s][1]->step(kx, ky, kz);
  kf_step_scalars[s][2]->step(kx, ky, kz);
  kf_copy[0]->step(kx, ky, kz);
  kf_copy[1]->step(kx, ky, kz);
  kf_copy[2]->step(kx, ky, kz);
  kf_copy[3]->step(kx, ky, kz);
}


