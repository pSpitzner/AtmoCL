__kernel void k_wrf_nesting_kernel_main(__private parameters par,
                                        __private uint frame_index,
                                        __read_only image3d_t bwrf_tgt_new,
                                        __read_only image3d_t bwrf_tgt_old,
                                        __read_only image3d_t b_source,
                                        __write_only image3d_t b_target)
{
  position pos = get_pos_bc(par, get_global_id(0), get_global_id(1), get_global_id(2));

  float4 sys     = read_f4(pos.x , pos.y , pos.z , b_source);
  float4 wrf_old = read_f4(pos.x , pos.y , pos.z , bwrf_tgt_old);
  float4 wrf_new = read_f4(pos.x , pos.y , pos.z , bwrf_tgt_new);

  float wrfdt = 15.0f*60.0f;
  float arg = fmod((float)(frame_index)*par.dT,wrfdt)/wrfdt;
  float weight = 0.5f*cospi(arg)+0.5f;

  // if (pos.x < 4 || pos.y < 4 || pos.z < 4 ||
  //     pos.x > par.sx-5 || pos.y > par.sy-5 || pos.z > par.sz-5 ) {
  //   sys = (1.0f-weight)*wrf_new + (weight)*wrf_old;
  // }
  sys = (1.0f-weight)*wrf_new + (weight)*wrf_old;

  write_f4(pos.x, pos.y, pos.z, sys, b_target);
}

