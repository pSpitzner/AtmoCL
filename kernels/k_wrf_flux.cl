__kernel void k_ext_forcings_isdac_kernel_main(__private parameters par,
                                               __private uint frame_index,
                                               __read_only image3d_t b_source_scalars_0,
                                               __read_only image3d_t b_source_scalars_1,
                                               __read_only image3d_t b_source_scalars_2,
                                               __read_only image3d_t b_source_momenta,
                                               __write_only image3d_t b_target_scalars_0,
                                               __write_only image3d_t b_target_scalars_1,
                                               __write_only image3d_t b_target_scalars_2,
                                               __write_only image3d_t b_target_momenta)
{
  position pos = get_pos_bc(par, get_global_id(0), get_global_id(1), get_global_id(2));

}
