__kernel void k_clone_kernel_main(__private parameters par,
                                  __read_only image3d_t b_source,
                                  __write_only image3d_t b_target)
{
  position pos = get_pos_bc(&par);

  float4 temp = read_imagef(b_source, (int4)(0, 0, pos.z, 0));
  write_imagef(b_target, (int4)(pos.x, pos.y, pos.z, 0), temp);
}
