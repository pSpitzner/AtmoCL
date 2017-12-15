__kernel void ks_prep_fast_0_kernel_main(__private parameters par,
                                         __read_only image3d_t bs_source_0,
                                         __write_only image3d_t bf_target_a)
{
  position pos = get_pos_bc(&par);

  if (par.timescheme == 0) {
    // only needed for rk3, plain copy
    float4 Temp = read_imagef(bs_source_0, (int4)(pos.x, pos.y, pos.z, 0));

    write_imagef(bf_target_a, (int4)(pos.x, pos.y, pos.z, 0), Temp);
  }
}
