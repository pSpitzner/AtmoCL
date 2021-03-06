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

__kernel void ke_rho_r_kernel_main(__private parameters par,
                                   __private uint ref,
                                   __private uint dim,
                                   __read_only image3d_t b_source_scalars_0,
                                   __read_only image3d_t b_source_scalars_1,
                                   __read_only image3d_t b_source_scalars_2,
                                   __read_only image3d_t b_source_momenta,
                                   __write_only image3d_t b_target)
{
  position pos = get_pos_bc(&par);

  float8 c;
  if      (dim == 0) c = read_f8(pos.x+ref, pos.y,     pos.z,     b_source_scalars_0, b_source_scalars_1);
  else if (dim == 1) c = read_f8(pos.x,     pos.y+ref, pos.z,     b_source_scalars_0, b_source_scalars_1);
  else if (dim == 2) c = read_f8(pos.x,     pos.y,     pos.z+ref, b_source_scalars_0, b_source_scalars_1);

  float4 rgba;
  rgba = map_rgba(max(0.0f, c.s4), 0.0f, 1e5f);
  write_f4(pos.x, pos.y, pos.z, &rgba, b_target);
}
