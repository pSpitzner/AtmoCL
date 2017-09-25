__private float4 map_rgba(float var, float offset, float factor) {
  // background blue
  float r = 0.0f;
  float g = 90.0f;
  float b = 255.0f;
  // abusing unused alpha channel to get the actual variable in full precision into host memory for further processing
  float a = var;
  float scale = 0.0f;

  scale = (var+offset)*factor;

  r += (255.0f-r)*scale;
  g += (255.0f-g)*scale;

  if (r > 255.0f) r = 255.0f;
  if (g > 255.0f) g = 255.0f;
  if (b > 255.0f) b = 255.0f;
  if (r < 0.0f) r = 0.0f;
  if (g < 0.0f) g = 0.0f;
  if (b < 0.0f) b = 0.0f;
  float4 rgba = (float4)(r, g, b, a);
  return rgba;
}

__kernel void ke_int_iwp_kernel_main(__private parameters par,
                                     __private uint ref,
                                     __private uint dim,
                                     __read_only image3d_t b_source_scalars_0,
                                     __read_only image3d_t b_source_scalars_1,
                                     __read_only image3d_t b_source_scalars_2,
                                     __read_only image3d_t b_source_momenta,
                                     __write_only image3d_t b_target)
{
  position pos = get_pos_bc(par, get_global_id(0), get_global_id(1), get_global_id(2));


  // ice water path
  float rho_i = 0.0f;

  // XY
  if (dim == 2) {
    for (int z = 0; z < par.sz; z++) {
      // rho_i += read_f4(pos.x, pos.y, z, b_source_scalars_0).s3*par.dz; // rho_c
      // rho_i += read_f4(pos.x, pos.y, z, b_source_scalars_1).s0*par.dz; // rho_r
      rho_i += read_f4(pos.x, pos.y, z, b_source_scalars_2).s0*par.dz; // rho_i
      // rho_i += read_f4(pos.x, pos.y, z, b_source_scalars_2).s1*par.dz; // rho_s
    }
  }


  float4 rgba = map_rgba(rho_i, 0.0f, 2.0f*1e1f);
  write_f4(pos.x, pos.y, pos.z, rgba, b_target);
}
