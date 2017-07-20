__kernel void k_debug_kernel_main(__private parameters par,
                                  __read_only image3d_t bf_scalars_vc_a_0,
                                  __read_only image3d_t bf_scalars_vc_a_1,
                                  __read_only image3d_t bf_scalars_vc_a_2,
                                  __read_only image3d_t bf_momenta_fc_a)
{
  position pos = get_pos_bc(par, get_global_id(0), get_global_id(1), get_global_id(2));

}
