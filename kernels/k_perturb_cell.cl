float rhovs(float T, parameters par) {
  // density of vapour at saturation
  float sv = par.svr*pow(T/par.tr, (par.cpv-par.cpl)/par.rv)*exp(par.lre0/par.rv*(1.0f/par.tr-1.0f/T));
  return sv/par.rv/T;
}

__kernel void k_perturb_cell_kernel_main(__private parameters par,
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

  float8 c    = read_f8(pos.x, pos.y, pos.z, b_source_scalars_0, b_source_scalars_1);
  float4 cice = read_f4(pos.x, pos.y, pos.z, b_source_scalars_2);
  float4 mom  = read_f4(pos.x, pos.y, pos.z, b_source_momenta);
  state st = init_state_with_ice(par, c, cice);

  //cosine squared, r in meters
  float r_h = 12000.0f;
  float r_v = 1400.0f;

  float4 arg = (float4)(0.0f);
  arg.s0 = ((pos.x+0.5f)*par.dx-par.sx*0.5f*par.dx)/r_h;
  arg.s1 = ((pos.y+0.5f)*par.dy-par.sy*0.5f*par.dy)/r_h;
  arg.s2 = ((pos.z+0.5f)*par.dz-r_v               )/r_v;

  float theta_vp;
  float offset, len;
  float T;
  offset = 0.0f;
  len = length(arg);
  if (len <= 1.0f) {
    offset = 2.0f*pow(cospi(len*0.5f), 2.0f);
    T = st.T + offset;

    float sig = (log(T)-log(st.rml)/(st.cpml/st.rml-1.0f))*(st.cpml-st.rml);
    // if (len <= 0.2f) printf("%f %f | %f %f\n", st.sig, sig, st.T, T);
    c.s0 = sig;
  }
  state stn = init_state_with_ice(par, c, cice);

  // if (pos.x == 64 && pos.y == 128) printf("%f %f %f %f %f\n", pos.z*par.dz, st.T, stn.T, st.sig, stn.sig);

  write_f8(pos.x, pos.y, pos.z, c,    b_target_scalars_0, b_target_scalars_1);
  write_f4(pos.x, pos.y, pos.z, cice, b_target_scalars_2);
  write_f4(pos.x, pos.y, pos.z, mom,  b_target_momenta);
}

