__kernel void ke_int_velocity_variance_kernel_main(__private parameters par,
                                                   __private uint ref,
                                                   __private uint dim,
                                                   __read_only image3d_t b_source_scalars_0,
                                                   __read_only image3d_t b_source_scalars_1,
                                                   __read_only image3d_t b_source_scalars_2,
                                                   __read_only image3d_t b_source_momenta,
                                                   __write_only image3d_t b_target)
{
  position pos = get_pos_bc(&par);

  // vertical profile for velocity variance. this is not super neat since we only calculate a mean on the host side

  float w = 0.0f;
  float w_temp = 0.0f;
  float var = 0.0f;

  // YZ host averages over Y, both should give same results, not tested, only using dim = 1
  if      (dim == 0) {
  }

  // XZ host averages over X
  else if (dim == 1) {
    for (int y = 0; y < par.sy; y++) {
      // check again: average over two vc densities to get fc?
      w_temp = (read_f4(pos.x, y, pos.z, b_source_momenta)).s2 /
               (read_f4(pos.x, y, pos.z, b_source_scalars_0)).s1;
      w = ((float)(y)*w + w_temp)/(float)(y+1);
    }
    for (int y = 0; y < par.sy; y++) {
      w_temp = (read_f4(pos.x, y, pos.z, b_source_momenta)).s2 /
               (read_f4(pos.x, y, pos.z, b_source_scalars_0)).s1;
      var = ((float)(y)*var + (w_temp-w)*(w_temp-w))/(float)(y+1);
    }
  }

  // abusing unused alpha channel to get the actual variable in full precision into host memory for further processing
  float4 rgba = (float4)(0.0f, 0.0f, 0.0f, var);
  write_f4(pos.x, pos.y, pos.z, rgba, b_target);
}
