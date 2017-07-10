float rhovs(float T, parameters par) {
  // density of vapour at saturation
  float sv = par.svr*pow(T/par.tr, (par.cpv-par.cpl)/par.rv)*exp(par.lre0/par.rv*(1.0f/par.tr-1.0f/T));
  return sv/par.rv/T;
}

__kernel void k_nesting_cell_kernel_main(__private parameters par,
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




  float8 c;
  float4 cice, mom;
  state st;
  float theta, theta_prof;
  float z, z_tr;
  float h, h_prof;
  float pv_prof;

  c    = read_f8(pos.x, pos.y, pos.z, b_source_scalars_0, b_source_scalars_1);
  cice = read_f4(pos.x, pos.y, pos.z, b_source_scalars_2);
  mom  = read_f4(pos.x, pos.y, pos.z, b_source_momenta);
  st = init_state_with_ice(par, c, cice);

  theta = exp((st.sig+st.rml*log(par.pr))/st.cpml);
  h = st.pv/st.sv;

  z = ((float)(pos.z)+0.5f)*par.dz;
  z_tr = 12000.0f;

  // theta profile according to weisman_klemp
  if (z <= z_tr) {
    theta_prof = (300.0f - (343.0f - 300.0f)*pow(z/z_tr, 5.0f/4.0f));
  } else {
    theta_prof = 343.0f*exp(par.gr/par.cpd/213.0f*(z-z_tr));
  }

  // relative humidity profile according to weisman_klemp
  // H ~ st.pv/satvap
  if (z <= z_tr) {
    h_prof = 1.0f - 3.0f/4.0f*pow(z/z_tr, 5.0f/4.0f);
  } else {
    h_prof = 0.25f;
  }
  pv_prof = h_prof*st.sv;


  // decrease damping over time
  float ds = (1.0f - 0.3f*exp(-(float)(damping_strength)/0.1f));
  float dm = (1.0f - 0.3f*exp(-(float)(damping_strength)/2000.0f));


  // if(pos.x == par.sx/2 && pos.y == par.sy/2) {
  // printf("%f | %f %f (%f %f) | %f %f \n", pos.z*par.dz, st.rho_v, h_prof, h, st.pv, st.sv, theta_prof, theta);
  // }

  // forcing
  // sig
  c.s0 += (theta_prof-theta)*1e-2f*ds;
  c.s2 += (pv_prof-st.pv)*1e-6f*ds;

  // this section needs work, temperature profile has a cut near z=10m with negative gradient, that should not occur.
  // fix density at ground according to A2 of isdac
  // if (pos.z == 0) c.s1 += (102000.0f - st.P)*1.0e-7f*ds;

  cice = (float4)(0.0f);
  mom *= dm;

  write_f8(pos.x, pos.y, pos.z, c,    b_target_scalars_0, b_target_scalars_1);
  write_f4(pos.x, pos.y, pos.z, cice, b_target_scalars_2);
  write_f4(pos.x, pos.y, pos.z, mom,  b_target_momenta);
}

