__kernel void ks_prep_fast_2_kernel_main(__private parameters par,
                                         __read_only image3d_t bs_source_0,
                                         __read_only image3d_t bs_source_1,
                                         __read_only image3d_t bs_source_2,
                                         __write_only image3d_t bf_target_a)
{
  position pos = get_pos_bc(&par);

  int s = 2;
  float4 scl, y_s[3];

  y_s[0] = read_imagef(bs_source_0, (int4)(pos.x, pos.y, pos.z, 0));
  y_s[1] = read_imagef(bs_source_1, (int4)(pos.x, pos.y, pos.z, 0));
  y_s[2] = read_imagef(bs_source_2, (int4)(pos.x, pos.y, pos.z, 0));

  scl = y_s[0];
  scl += par.a[s][1]*(y_s[1]-y_s[0]);
  scl += par.a[s][2]*(y_s[2]-y_s[0]);

  write_imagef(bf_target_a,   (int4)(pos.x, pos.y, pos.z, 0), scl);
}
