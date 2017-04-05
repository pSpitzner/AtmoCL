__kernel void k_copy_two_kernel_main(__private parameters par,
                                     __read_only image3d_t b_source_0,
                                     __read_only image3d_t b_source_1,
                                     __write_only image3d_t b_target_0,
                                     __write_only image3d_t b_target_1)
{
  position pos = get_pos_bc(par, get_global_id(0), get_global_id(1), get_global_id(2));

  float4 temp_0 = read_imagef(b_source_0, (int4)(pos.x, pos.y, pos.z, 0));
  float4 temp_1 = read_imagef(b_source_1, (int4)(pos.x, pos.y, pos.z, 0));
  write_imagef(b_target_0, (int4)(pos.x, pos.y, pos.z, 0), temp_0);
  write_imagef(b_target_1, (int4)(pos.x, pos.y, pos.z, 0), temp_1);
}
