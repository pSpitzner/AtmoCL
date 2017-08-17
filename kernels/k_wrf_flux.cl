__kernel void k_wrf_flux_kernel_main(__private parameters par,
                                     __private uint frame_index,
                                     __read_only image3d_t b_source_scalars_0,
                                     __read_only image3d_t b_source_scalars_1,
                                     __read_only image3d_t b_source_scalars_2,
                                     __read_only image3d_t b_source_momenta,
                                     __read_only image3d_t bwrf_flux_new,
                                     __read_only image3d_t bwrf_flux_old,
                                     __write_only image3d_t b_target_scalars_0)
{
  position pos = get_pos_bc(par, get_global_id(0), get_global_id(1), get_global_id(2));


  float8    c = read_f8(pos.x, pos.y, pos.z, b_source_scalars_0, b_source_scalars_1);
  float4 cice = read_f4(pos.x, pos.y, pos.z, b_source_scalars_2);
  state    st = init_state_with_ice(par, c, cice);

  float4 wrf_new = read_f4(pos.x, pos.y, pos.z, bwrf_flux_new);
  float4 wrf_old = read_f4(pos.x, pos.y, pos.z, bwrf_flux_old);

  float wrfdt = 15.0f*60.0f;
  float weight = fmod((float)(frame_index)*par.dT/wrfdt,wrfdt);
  float4 flux = (1.0f-weight)*wrf_old + (weight)*wrf_new;

  // should always be true, only cal this kernel for ground layer
  if (pos.s_zl==-1) {
    // Watt per m**2
    heat_flux  = flux.s0;
    moist_flux = flux.s1;
  }
  printf("%d %d %d %f\n", pos.x, pos.y, pos.z, heat_flux);

  float4 cRhs = (float4)(0.0f);
  // d(rhosig) nach drho_v  // trockene dichte weg, wasserdampf hin
  cRhs.s0 += (st.lnT*(par.cpv-par.cpd)-st.lnP*(par.rv-par.rd))*moist_flux;
  cRhs.s2 += moist_flux;
  cRhs.s0 += heat_flux;

  c += cRhs*par.dT;

  write_imagef(b_target_scalars_0, (int4)(pos.x, pos.y, pos.z, 0), c.s0123);
}
