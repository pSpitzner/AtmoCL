__private float3 add_rgba_to_rgb(float3 bg, float4 src) {
  float3 target;
  target.s0 = ((1.0f - src.s3) * bg.s0) + (src.s3 * src.s0);
  target.s1 = ((1.0f - src.s3) * bg.s1) + (src.s3 * src.s1);
  target.s2 = ((1.0f - src.s3) * bg.s2) + (src.s3 * src.s2);
  return target;
}

__kernel void ke_int_rho_i_kernel_main(__private parameters par,
                                       __private uint ref,
                                       __private uint dim,
                                       __read_only image3d_t b_source_scalars_0,
                                       __read_only image3d_t b_source_scalars_1,
                                       __read_only image3d_t b_source_scalars_2,
                                       __read_only image3d_t b_source_momenta,
                                       __write_only image3d_t b_target)
{
  position pos = get_pos_bc(par, get_global_id(0), get_global_id(1), get_global_id(2));

  float rho_i = 0.0f;
  // YZ
  if (dim == 0) {
    for (int x = 0; x < par.sx; x++) {
      rho_i += read_f4(x, pos.y, pos.z, b_source_scalars_2).s0*par.dz;
    }
    rho_i/=(float)(par.sx);
  }
  // XZ
  else if (dim == 1) {
    for (int y = 0; y < par.sy; y++) {
      rho_i += read_f4(pos.x, y, pos.z, b_source_scalars_2).s0*par.dz;
    }
    rho_i/=(float)(par.sy);
  }
  // XY
  else if (dim == 2) {
    for (int z = 0; z < par.sz; z++) {
      rho_i += read_f4(pos.x, pos.y, z, b_source_scalars_2).s0*par.dz;
    }
    rho_i/=(float)(par.sz);
  }

  float ref_i = 5e-5f;

  float a_i = max(0.0f, (rho_i < ref_i ? rho_i/ref_i*255.0f : 255.0f));

  float3 bg = (float3)(0.0f, 0.0f, 0.0f);
  float4 rgba_i = (float4)(255.0f, 255.0f, 255.0f, a_i)/255.0f;

  float3 result = bg;
  result = add_rgba_to_rgb(result, rgba_i);
  result = result*255.0f;
  // if (a_c != 0.0f) printf("%d %d | %e %e | %f %f %f\n", pos.z, pos.x, a_c, rho_c, result.r, result.g, result.b);

  float4 rgba = (float4)(result, rho_i);
  write_f4(pos.x, pos.y, pos.z, rgba, b_target);
}
