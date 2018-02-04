__kernel void ke_int_ql_kernel_main(__private parameters par,
                                    __private uint ref,
                                    __private uint dim,
                                    __read_only image3d_t b_source_scalars_0,
                                    __read_only image3d_t b_source_scalars_1,
                                    __read_only image3d_t b_source_scalars_2,
                                    __read_only image3d_t b_source_momenta,
                                    __write_only image3d_t b_target)
{
  position pos = get_pos_bc(&par);

  // vertical profile for liquid water mixing ratio
  float8 c     = (float8)(0.0f);
  float4 cice = (float4)(0.0f);
  state st;

  float ql = 0.0f;
  float q_temp = 0.0f;

  // YZ
  if      (dim == 0) {
    for (int x = 0; x < par.sx; x++) {
      c     = read_f8(x, pos.y, pos.z, b_source_scalars_0, b_source_scalars_1);
      cice = read_f4(x, pos.y, pos.z, b_source_scalars_2);
      st = init_state_with_ice(&par, &c, &cice);
      q_temp = (st.rho_d != 0.0f ? st.rho_l/st.rho_d : 0.0f);
      ql = ((float)(x)*ql + q_temp)/(float)(x+1);
    }
  }
  // XZ
  else if (dim == 1) {
    for (int y = 0; y < par.sy; y++) {
      c     = read_f8(pos.x, y, pos.z, b_source_scalars_0, b_source_scalars_1);
      cice = read_f4(pos.x, y, pos.z, b_source_scalars_2);
      st = init_state_with_ice(&par, &c, &cice);
      q_temp = (st.rho_d != 0.0f ? st.rho_l/st.rho_d : 0.0f);
      ql = ((float)(y)*ql + q_temp)/(float)(y+1);
    }
  }

  // abusing unused alpha channel to get the actual variable in full precision into host memory for further processing
  float4 rgba = (float4)(0.0f, 0.0f, 0.0f, ql);
  write_f4(pos.x, pos.y, pos.z, &rgba, b_target);
}
