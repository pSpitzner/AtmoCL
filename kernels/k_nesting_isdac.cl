float rhovs(float T, parameters par) {
  // density of vapour at saturation
  float sv = par.svr*pow(T/par.tr, (par.cpv-par.cpl)/par.rv)*exp(par.lre0/par.rv*(1.0f/par.tr-1.0f/T));
  return sv/par.rv/T;
}

__kernel void k_nesting_isdac_kernel_main(__private parameters par,
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
  position pos = get_pos_bc(&par);




  float8 c;
  float4 cice, mom;
  state st;
  float theta, theta_l, theta_l_prof;
  float q_t, q_t_prof;
  float z;

  c    = read_f8(pos.x, pos.y, pos.z, b_source_scalars_0, b_source_scalars_1);
  cice = read_f4(pos.x, pos.y, pos.z, b_source_scalars_2);
  mom  = read_f4(pos.x, pos.y, pos.z, b_source_momenta);
  st = init_state_with_ice(&par, &c, &cice);

  q_t = (st.rho_l+st.rho_v)/st.rho;
  theta = exp((st.sig+st.rml*log(par.pr))/st.cpml);
  theta_l = theta - (theta/st.T*par.lr/st.cpml)*st.rho_l/st.rho;

  z = ((float)(pos.z)+0.5f)*par.dz;

  // theta_l profile according to isdac
  if (z < 400.0f) {
    theta_l_prof = 265.0f + 0.004f*(z-400.0f);
  } else if (z >= 400.0f && z < 825.0f) {
    theta_l_prof = 265.0f;
  } else if (z >= 825.0f && z < 2045.0f) {
    theta_l_prof = 266.0f + pow(z-825.0f, 0.3f);
  } else {
    theta_l_prof = 271.0f + pow(z-2000.0f, 0.33f);
  }

  // q_t_prof profile according to isdac
  if (z < 400.0f) {
    q_t_prof = 1.5f - 0.00075f*(z-400.0f);
  } else if (z >= 400.0f && z < 825.0f) {
    q_t_prof = 1.5f;
  } else if (z >= 825.0f && z < 2045.0f) {
    q_t_prof = 1.2f;
  } else {
    q_t_prof = 0.5f - 0.000075f*(z-2045.0f);
  }
  // correct units
  q_t_prof *= 1.0e-3f;

  // decrease damping over time
  float ds = (1.0f - 0.3f*exp(-(float)(damping_strength)/20.0f));
  float dm = (1.0f - 0.3f*exp(-(float)(damping_strength)/2000.0f));

  // forcing
  // sig
  c.s0 += (theta_l_prof-theta_l)*1e-1f*ds;
  c.s2 += (q_t_prof-q_t)*1e-1f*ds;

  // this section needs work, temperature profile has a cut near z=10m with negative gradient, that should not occur.
  // fix density at ground according to A2 of isdac
  if (pos.z == 0) c.s1 += (102000.0f - st.P)*1.0e-7f*ds;

  cice = (float4)(0.0f);
  mom *= dm;

  write_f8(pos.x, pos.y, pos.z, &c,    b_target_scalars_0, b_target_scalars_1);
  write_f4(pos.x, pos.y, pos.z, &cice, b_target_scalars_2);
  write_f4(pos.x, pos.y, pos.z, &mom,  b_target_momenta);
}

