__kernel void k_damping_kernel_main(__private parameters par,
                                    __private uint damping_strength,
                                    __read_only image3d_t b_source_momenta,
                                    __write_only image3d_t b_target_momenta)
{
  position pos = get_pos_bc(par, get_global_id(0), get_global_id(1), get_global_id(2));

  float4 TempMomenta = read_imagef(b_source_momenta, (int4)(pos.x, pos.y, pos.z, 0));

  // avoid waved due to discretisation when importing fresh wrf file
  TempMomenta.s0 = 0.0f;
  TempMomenta.s1 = 0.0f;
  TempMomenta.s2 *= (1.0f - 0.3f*exp(-(float)(damping_strength)/100.0f));

  write_imagef(b_target_momenta,   (int4)(pos.x, pos.y, pos.z, 0), TempMomenta);
}
