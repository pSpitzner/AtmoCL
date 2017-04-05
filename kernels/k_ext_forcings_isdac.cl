__kernel void k_ext_forcings_isdac_kernel_main(__private parameters par,
                                               __private uint frame_index,
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

  float4 TempScalars_0 = read_imagef(b_source_scalars_0, (int4)(pos.x, pos.y, pos.z, 0));
  float4 TempScalars_1 = read_imagef(b_source_scalars_1, (int4)(pos.x, pos.y, pos.z, 0));
  float4 TempMomenta = read_imagef(b_source_momenta, (int4)(pos.x, pos.y, pos.z, 0));

  float8 c;
  state st;

  c = read_f8(pos.x, pos.y, pos.z, b_source_scalars_0, b_source_scalars_1);
  if (include_ice) {
    float4 cice = read_f4(pos.x, pos.y, pos.z, b_source_scalars_2);
    st = init_state_with_ice(par, c, cice);
  } else {
    st = init_state(par, c);
  }


  float theta_p = st.T*pow(par.pr/exp(st.lnP),par.rd/par.cpv);
  float moist_flux = 0.0f;
  float heat_flux = 0.0f;

  if (pos.s_zl==-1) {
    // Watt per m**2
    moist_flux = 0.0f;
    heat_flux  = 0.0f;
  }

  float8 cRhs = (float8)(0.0f);

  float qtop, qbot;
  float rc1, rc2;
  int h;
  float hf;

  // simple radiation according to isadc
  if (pos.s_zl==1 && pos.s_zr==1) {
    qtop=0.0f; qbot=0.0f;
    for (h=pos.z; h<par.sz; h++) {
      // add up rho_l(h) = rho_c(h) + rho_d(h) for all h above this cell
      // treat ice here - or dont, use different forcings!
      qtop +=max(0.0f, read_f4(pos.x, pos.y, h, b_source_scalars_0).s3 +
                 read_f4(pos.x, pos.y, h, b_source_scalars_1).s0)*par.dz;
    }

    for (h=pos.z; h>0; h--) {
      qbot +=max(0.0f, read_f4(pos.x, pos.y, h, b_source_scalars_0).s3 +
                 read_f4(pos.x, pos.y, h, b_source_scalars_1).s0)*par.dz;
    }
    rc1 = 72.0f*exp(-170.0f*qtop)+15.0f*exp(-170.0f*qbot);

    qtop -=max(0.0f, st.rho_l)*par.dz;
    qbot +=max(0.0f, st.rho_l)*par.dz;

    rc2 =72.0f*exp(-170.0f*qtop)+15.0f*exp(-170.0f*qbot);

    cRhs.s0 += (rc1-rc2)/st.T/par.dz;
  }

  // d(rhosig) nach drho_v  // trockene dichte weg, wasserdampf hin
  // cRhs.s0 += (st.lnT*(par.cpv-par.cpd)-st.lnP*(par.rv-par.rd))*moist_flux;
  // cRhs.s2 += moist_flux;

  // cRhs.s0 += heat_flux;

  TempScalars_0.s0 += cRhs.s0*par.dT;
  TempScalars_0.s2 += cRhs.s2*par.dT;

  // add random disturbance
  if (frame_index == 1) {
    float rhosig = c.s0;
    if (pos.z <= 200.0f/par.dz) {
       unsigned long seed = get_global_id(0) + get_global_id(1)*par.sx + get_global_id(2)*par.sx*par.sy;
      TempScalars_0.s0 -=(rhosig/1.0e7f)*(rand(seed)-0.5f);
    }
  }

  // asystem expects this kernel to do a copy for all scalars...
  float4 TempScalars_2 = read_imagef(b_source_scalars_2, (int4)(pos.x, pos.y, pos.z, 0));

  write_imagef(b_target_scalars_0, (int4)(pos.x, pos.y, pos.z, 0), TempScalars_0);
  write_imagef(b_target_scalars_1, (int4)(pos.x, pos.y, pos.z, 0), TempScalars_1);
  write_imagef(b_target_scalars_2, (int4)(pos.x, pos.y, pos.z, 0), TempScalars_2);
  write_imagef(b_target_momenta,   (int4)(pos.x, pos.y, pos.z, 0), TempMomenta);
}
