

__kernel void k_test_kernel_main(__private parameters par,
                                 __private wrfparameters wrf)
{
  position pos = get_pos_bc(par, get_global_id(0), get_global_id(1), get_global_id(2));

  printf("kernel: %d %d %d\n", wrf.sx, wrf.sy, wrf.sz);
}
