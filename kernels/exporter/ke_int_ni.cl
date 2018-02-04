__kernel void ke_int_ni_kernel_main(__private parameters par,
                                    __private uint ref,
                                    __private uint dim,
                                    __read_only image3d_t b_source_scalars_0,
                                    __read_only image3d_t b_source_scalars_1,
                                    __read_only image3d_t b_source_scalars_2,
                                    __read_only image3d_t b_source_momenta,
                                    __write_only image3d_t b_target)
{
  position pos = get_pos_bc(&par);

  float ni = 0.0f;
  float n_temp = 0.0f;

  // YZ
  if      (dim == 0) {
    for (int x = 0; x < par.sx; x++) {
      n_temp = read_f4(x, pos.y, pos.z, b_source_scalars_2).s2;
      ni = ((float)(x)*ni + n_temp)/(float)(x+1);
    }
  }
  // XZ
  else if (dim == 1) {
    for (int y = 0; y < par.sy; y++) {
      n_temp = read_f4(pos.x, y, pos.z, b_source_scalars_2).s2;
      ni = ((float)(y)*ni + n_temp)/(float)(y+1);
    }
  }

  ni /= 1000.0f; // number density from per m**3 to per litre

  // abusing unused alpha channel to get the actual variable in full precision into host memory for further processing
  float4 rgba = (float4)(0.0f, 0.0f, 0.0f, ni);
  write_f4(pos.x, pos.y, pos.z, &rgba, b_target);
}
