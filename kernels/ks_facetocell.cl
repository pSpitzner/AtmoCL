__kernel void ks_facetocell_kernel_main(__private parameters par,
                                        __read_only image3d_t bf_momenta_fc_a,
                                        __write_only image3d_t b_m_vc_uv,
                                        __write_only image3d_t b_m_vc_w)
{
  position pos = get_pos_bc(par, get_global_id(0), get_global_id(1), get_global_id(2));

  float ul,vl,ur,vr,wl,wr;

  ul = read_imagef(bf_momenta_fc_a, (int4)(pos.x , pos.y , pos.z , 0)).x;
  ur = read_imagef(bf_momenta_fc_a, (int4)(pos.xr, pos.y , pos.z , 0)).x;
  vl = read_imagef(bf_momenta_fc_a, (int4)(pos.x , pos.y , pos.z , 0)).y;
  vr = read_imagef(bf_momenta_fc_a, (int4)(pos.x , pos.yr, pos.z , 0)).y;
  wl = read_imagef(bf_momenta_fc_a, (int4)(pos.x , pos.y , pos.z , 0)).z;
  wr = read_imagef(bf_momenta_fc_a, (int4)(pos.x , pos.y , pos.zr, 0)).z;

  // fixed bc
  if (pos.z == 0) {
    wl = 0.0f;
  } else if (pos.z==par.sz-1) {
    wr = 0.0f;
  }

  // if (pos.x == 50 && pos.y ==43 && pos.z == 19) printf("uls: %f %f %f %f %f %f\n", ul, ur, vl, vr, wl, wr );
  write_imagef(b_m_vc_uv, (int4)(pos.x, pos.y, pos.z, 0), (float4)(ul, ur, vl, vr));
  write_imagef(b_m_vc_w , (int4)(pos.x, pos.y, pos.z, 0), (float4)(wl, wr, 0.0f, 0.0f));
}
