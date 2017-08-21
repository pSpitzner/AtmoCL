// ----------------------------------------------------------------- //
// Add momenta back onto thermalised wrf_tgt //
// ----------------------------------------------------------------- //

__kernel void k_wrf_finish_kernel_main(__private parameters par,
                                       __read_only image3d_t bwrf_src_scalars_0,
                                       __read_only image3d_t bwrf_src_momenta,
                                       __read_only image3d_t bsys_src_scalars_0,
                                       __read_only image3d_t bsys_src_momenta,
                                       __write_only image3d_t bwrf_tgt_momenta)
{
  position pos = get_pos_bc(par, get_global_id(0), get_global_id(1), get_global_id(2));

  float wrf_rho    = read_f4(pos.x , pos.y , pos.z , bwrf_src_scalars_0).s1;
  float wrf_rho_xl = read_f4(pos.xl, pos.y , pos.z , bwrf_src_scalars_0).s1;
  float wrf_rho_yl = read_f4(pos.x , pos.yl, pos.z , bwrf_src_scalars_0).s1;

  float sys_rho    = read_f4(pos.x , pos.y , pos.z , bsys_src_scalars_0).s1;
  float sys_rho_xl = read_f4(pos.xl, pos.y , pos.z , bsys_src_scalars_0).s1;
  float sys_rho_yl = read_f4(pos.x , pos.yl, pos.z , bsys_src_scalars_0).s1;

  float sys_rho_w  = read_f4(pos.x , pos.y , pos.z , bsys_src_momenta).s2;
  float4 wrf_mom   = read_f4(pos.x , pos.y , pos.z , bwrf_src_momenta);

  float u = 2.0f*wrf_mom.s0/(wrf_rho_xl+wrf_rho);
  float v = 2.0f*wrf_mom.s1/(wrf_rho_yl+wrf_rho);


  float tgt_rho_u = 0.5f*(sys_rho_xl+sys_rho)*u;
  float tgt_rho_v = 0.5f*(sys_rho_yl+sys_rho)*v;

  // if (pos.x == par.sx/2 && pos.y == par.sy/2 && pos.z == par.sz/2) printf("%f %f %f\n", tgt_rho_u, tgt_rho_v, sys_rho_w);

  write_f4(pos.x, pos.y, pos.z, (float4)(tgt_rho_u, tgt_rho_v, sys_rho_w, 0.0f),  bwrf_tgt_momenta);
}

