__private float4 map_rgba(float var, float offset, float factor) {
  float r = 0.0f;
  float g = 0.0f;
  float b = 0.0f;
  float a = var;

  r = (var+offset)*factor*255.0f;
  if (r<0.0f) b = - r;

  if (r > 255.0f) r = 255.0f;
  if (g > 255.0f) g = 255.0f;
  if (b > 255.0f) b = 255.0f;
  if (r < 0.0f) r = 0.0f;
  if (g < 0.0f) g = 0.0f;
  if (b < 0.0f) b = 0.0f;
  float4 rgba = (float4)(r, g, b, a);
  return rgba;
}

__kernel void ke_theta_e_kernel_main(__private parameters par,
                                     __private uint ref,
                                     __private uint dim,
                                     __read_only image3d_t b_source_scalars_0,
                                     __read_only image3d_t b_source_scalars_1,
                                     __read_only image3d_t b_source_scalars_2,
                                     __read_only image3d_t b_source_momenta,
                                     __write_only image3d_t b_target)
{
  position pos = get_pos_bc(&par);

  float8 c     = (float8)(0.0f);
  float4 c_ice = (float4)(0.0f);
  state st;

  float theta_e = 0.0f;
  float pd, r_t, r_v;

  // YZ
  if      (dim == 0) {
    c     = read_f8(pos.x+ref, pos.y, pos.z, b_source_scalars_0, b_source_scalars_1);
    c_ice = read_f4(pos.x+ref, pos.y, pos.z, b_source_scalars_2);
    st = init_state_with_ice(&par, &c, &c_ice);
    pd = st.rho_d*par.rd*st.T;
    r_t = (st.rho_l+st.rho_v)/st.rho_d;
    r_v = (st.rho_v)/st.rho_d;
    theta_e = st.T*pow(pd/par.pr,-par.rd/(par.cpd+par.cpl*r_t))*exp(par.lre0*r_v/(par.cpd+par.cpl*r_t)/st.T);
  }
  // XZ
  else if (dim == 1) {
    c     = read_f8(pos.x, pos.y+ref, pos.z, b_source_scalars_0, b_source_scalars_1);
    c_ice = read_f4(pos.x, pos.y+ref, pos.z, b_source_scalars_2);
    st = init_state_with_ice(&par, &c, &c_ice);
    pd = st.rho_d*par.rd*st.T;
    r_t = (st.rho_l+st.rho_v)/st.rho_d;
    r_v = (st.rho_v)/st.rho_d;
    theta_e = st.T*pow(pd/par.pr,-par.rd/(par.cpd+par.cpl*r_t))*exp(par.lre0*r_v/(par.cpd+par.cpl*r_t)/st.T);
  }

  // if (pos.x == par.sx/2) printf("%d %d %d %g\n", pos.x, pos.y, pos.z, theta_e);
  // printf("%d %f\n", pos.z, theta_e);
  float4 rgba;
  rgba = map_rgba(theta_e, -320.0f, 5.0e-1f);
  write_f4(pos.x, pos.y, pos.z, &rgba, b_target);
}
