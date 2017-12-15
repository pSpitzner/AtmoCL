__kernel void k_damping_kernel_main(__private parameters par,
                                    __private uint damping_strength,
                                    __read_only image3d_t b_source_momenta,
                                    __write_only image3d_t b_target_momenta)
{
  position pos = get_pos_bc(&par);

  float4 TempMomenta = read_imagef(b_source_momenta, (int4)(pos.x, pos.y, pos.z, 0));

  // reduce initial pressure waves due to deviations of the
  TempMomenta *= (1.0f - 0.3f*exp(-(float)(damping_strength)));

  write_imagef(b_target_momenta,   (int4)(pos.x, pos.y, pos.z, 0), TempMomenta);
}
