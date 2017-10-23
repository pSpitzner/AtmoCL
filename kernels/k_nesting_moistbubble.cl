float rhovs(float T, parameters par) {
  // density of vapour at saturation
  float sv = par.svr*pow(T/par.tr, (par.cpv-par.cpl)/par.rv)*exp(par.lre0/par.rv*(1.0f/par.tr-1.0f/T));
  return sv/par.rv/T;
}

__kernel void k_nesting_moistbubble_kernel_main(__private parameters par,
                                                __private uint damping_strength,
                                                __read_only image3d_t b_source_scalars_0,
                                                __read_only image3d_t b_source_scalars_1,
                                                __read_only image3d_t b_source_scalars_2,
                                                __read_only image3d_t b_source_momenta,
                                                __write_only image3d_t b_target_scalars_0,
                                                __write_only image3d_t b_target_scalars_1,
                                                __write_only image3d_t b_target_scalars_2,
                                                __write_only image3d_t b_target_momenta)
{
  position pos = get_pos_bc(par, get_global_id(0), get_global_id(1), get_global_id(2));

  // just create hydrostatic equilibrium according to profile. bubble needs to be added later to 'break' it

  float8 c;
  float4 cice, mom;
  state st;
  float theta_e, theta_e_prof;
  float q_t, q_t_prof;
  float q_v, q_v_prof;
  float q_l, q_l_prof;
  float z;
  float pd; // partial pressure of dry air

  c    = read_f8(pos.x, pos.y, pos.z, b_source_scalars_0, b_source_scalars_1);
  cice = (float4)(0.0f);
  mom  = read_f4(pos.x, pos.y, pos.z, b_source_momenta);
  st = init_state_with_ice(par, c, cice);

  q_t = (st.rho_l+st.rho_v)/st.rho;
  q_v = (st.rho_v)/st.rho;
  q_l = q_t - q_v;

  pd = st.rho_d*par.rd*st.T; // ?
  theta_e = st.T*pow(pd/par.pr,-par.rd/(par.cpd+par.cpl*q_t))*exp(par.lre0*q_v/(par.cpd+par.cpl*q_t));

  float theta, theta_l;
  theta = exp((st.sig+st.rml*log(par.pr))/st.cpml);
  theta_l = theta - (theta/st.T*par.lr/st.cpml)*st.rho_l/st.rho;

  theta_e_prof = 320.0f;
  q_t_prof = 0.02f;
  q_v_prof = rhovs(st.T, par)/st.rho;
  q_l_prof = q_t_prof - q_v_prof;

  // if (pos.x == par.sx/2 && pos.z == 0) printf("%d %g %g (%g) |%g (%g) %g (%g) | %g | %g\n", pos.z, theta, theta_e, theta_e_prof, q_t, q_t_prof, q_v, q_v_prof, st.P, offset);

  // decrease damping over time
  float ds = (1.0f - 0.3f*exp(-(float)(damping_strength)/20.0f));
  float dm = (1.0f - 0.3f*exp(-(float)(damping_strength)/2000.0f));

  // forcing
  // sig
  c.s0 += (theta_e_prof-theta_e)*1e-0f*ds;
  c.s2 += (q_v_prof-q_v)*1e-2f*ds;
  c.s3 += (q_l_prof-q_l)*1e-2f*ds;

  // fix density at ground according to reference pressure
  if (pos.z == 0) c.s1 += (par.pr - st.P)*1.0e-7f*ds;

  cice = (float4)(0.0f);
  mom *= dm;

  write_f8(pos.x, pos.y, pos.z, c,    b_target_scalars_0, b_target_scalars_1);
  write_f4(pos.x, pos.y, pos.z, cice, b_target_scalars_2);
  write_f4(pos.x, pos.y, pos.z, mom,  b_target_momenta);
}

