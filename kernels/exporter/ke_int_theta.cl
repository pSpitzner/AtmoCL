__kernel void ke_int_theta_kernel_main(__private parameters par,
                                       __private uint ref,
                                       __private uint dim,
                                       __read_only image3d_t b_source_scalars_0,
                                       __read_only image3d_t b_source_scalars_1,
                                       __read_only image3d_t b_source_scalars_2,
                                       __read_only image3d_t b_source_momenta,
                                       __write_only image3d_t b_target)
{
  position pos = get_pos_bc(par, get_global_id(0), get_global_id(1), get_global_id(2));

  float8 c     = (float8)(0.0f);
  float4 c_ice = (float4)(0.0f);
  state st;

  float theta = 0.0f;
  float theta_temp = 0.0f;

  // YZ
  if      (dim == 0) {
    for (int x = 0; x < par.sx; x++) {
      c     = read_f8(x, pos.y, pos.z, b_source_scalars_0, b_source_scalars_1);
      c_ice = read_f4(x, pos.y, pos.z, b_source_scalars_2);
      st = init_state_with_ice(par, c, c_ice);
      theta_temp = exp((st.sig+st.rml*log(par.pr))/st.cpml);
      theta = ((float)(x)*theta + theta_temp)/(float)(x+1);
    }
  }
  // XZ
  else if (dim == 1) {
    for (int y = 0; y < par.sy; y++) {
      c     = read_f8(pos.x, y, pos.z, b_source_scalars_0, b_source_scalars_1);
      c_ice = read_f4(pos.x, y, pos.z, b_source_scalars_2);
      st = init_state_with_ice(par, c, c_ice);
      theta_temp = exp((st.sig+st.rml*log(par.pr))/st.cpml);
      theta = ((float)(y)*theta + theta_temp)/(float)(y+1);
    }
  }

  // abusing unused alpha channel to get the actual variable in full precision into host memory for further processing
  float4 rgba = (float4)(0.0f, 0.0f, 0.0f, theta);
  write_f4(pos.x, pos.y, pos.z, rgba, b_target);
}
