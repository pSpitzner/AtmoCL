__kernel void ks_adv_scalars_kernel_main(__private parameters par,
                                         __read_only image3d_t bf_scalars_vc_a,
                                         __read_only image3d_t bf_density_vc_a,
                                         __read_only image3d_t bf_momenta_fc_a,
                                         __write_only image3d_t bs_scalars_fc_x_s,
                                         __write_only image3d_t bs_scalars_fc_y_s,
                                         __write_only image3d_t bs_scalars_fc_z_s)
{

  position pos = get_pos_bc(&par);

  // ----------------------------------------------------------------- //
  // OS X does not support writing into >8 textures. need a workaround //
  // ----------------------------------------------------------------- //


  float4 c, cxl, cxr, cyl, cyr, czl, czr;
  float4 cxll, cxrr, cyll, cyrr, czll, czrr;
  float4 Fx, Fy, Fz;                  // sit on faces, former Fxl, Fyl WITHOUT *u


  // divide by rho to get dimension right in fast step
  c    = read_f4(pos.x  , pos.y  , pos.z  , bf_scalars_vc_a) /
         read_f4(pos.x  , pos.y  , pos.z  , bf_density_vc_a).s1;
  cxl  = read_f4(pos.xl , pos.y  , pos.z  , bf_scalars_vc_a) /
         read_f4(pos.xl , pos.y  , pos.z  , bf_density_vc_a).s1;
  cxr  = read_f4(pos.xr , pos.y  , pos.z  , bf_scalars_vc_a) /
         read_f4(pos.xr , pos.y  , pos.z  , bf_density_vc_a).s1;
  cyl  = read_f4(pos.x  , pos.yl , pos.z  , bf_scalars_vc_a) /
         read_f4(pos.x  , pos.yl , pos.z  , bf_density_vc_a).s1;
  cyr  = read_f4(pos.x  , pos.yr , pos.z  , bf_scalars_vc_a) /
         read_f4(pos.x  , pos.yr , pos.z  , bf_density_vc_a).s1;
  czl  = read_f4(pos.x  , pos.y  , pos.zl , bf_scalars_vc_a) /
         read_f4(pos.x  , pos.y  , pos.zl , bf_density_vc_a).s1;
  czr  = read_f4(pos.x  , pos.y  , pos.zr , bf_scalars_vc_a) /
         read_f4(pos.x  , pos.y  , pos.zr , bf_density_vc_a).s1;

  cxll = read_f4(pos.xll, pos.y  , pos.z  , bf_scalars_vc_a) /
         read_f4(pos.xll, pos.y  , pos.z  , bf_density_vc_a).s1;
  cxrr = read_f4(pos.xrr, pos.y  , pos.z  , bf_scalars_vc_a) /
         read_f4(pos.xrr, pos.y  , pos.z  , bf_density_vc_a).s1;
  cyll = read_f4(pos.x  , pos.yll, pos.z  , bf_scalars_vc_a) /
         read_f4(pos.x  , pos.yll, pos.z  , bf_density_vc_a).s1;
  cyrr = read_f4(pos.x  , pos.yrr, pos.z  , bf_scalars_vc_a) /
         read_f4(pos.x  , pos.yrr, pos.z  , bf_density_vc_a).s1;
  czll = read_f4(pos.x  , pos.y  , pos.zll, bf_scalars_vc_a) /
         read_f4(pos.x  , pos.y  , pos.zll, bf_density_vc_a).s1;
  czrr = read_f4(pos.x  , pos.y  , pos.zrr, bf_scalars_vc_a) /
         read_f4(pos.x  , pos.y  , pos.zrr, bf_density_vc_a).s1;

  // momenta, only for direction of upwind
  pos.ulf = read_imagef(bf_momenta_fc_a, (int4)(pos.x, pos.y, pos.z, 0)).s0;
  pos.vlf = read_imagef(bf_momenta_fc_a, (int4)(pos.x, pos.y, pos.z, 0)).s1;
  pos.wlf = read_imagef(bf_momenta_fc_a, (int4)(pos.x, pos.y, pos.z, 0)).s2;

  // fixed bc
  if (pos.s_xl==-1) pos.ulf=0.0f;
  if (pos.s_yl==-1) pos.vlf=0.0f;
  if (pos.s_zl==-1) pos.wlf=0.0f;

  Fx = 0.5f*(cxl+c);
  Fy = 0.5f*(cyl+c);
  Fz = 0.5f*(czl+c);

  if (pos.ulf>0.0f) Fx = (  c/3.0f+5.0f/6.0f*cxl-cxll/6.0f);
  if (pos.ulf<0.0f) Fx = (cxl/3.0f+5.0f/6.0f*  c- cxr/6.0f);

  if (pos.vlf>0.0f) Fy = (  c/3.0f+5.0f/6.0f*cyl-cyll/6.0f);
  if (pos.vlf<0.0f) Fy = (cyl/3.0f+5.0f/6.0f*  c- cyr/6.0f);

  if (pos.wlf>0.0f) Fz = (  c/3.0f+5.0f/6.0f*czl-czll/6.0f);
  if (pos.wlf<0.0f) Fz = (czl/3.0f+5.0f/6.0f*  c- czr/6.0f);

  write_f4(pos.x, pos.y, pos.z, &Fx, bs_scalars_fc_x_s);
  write_f4(pos.x, pos.y, pos.z, &Fy, bs_scalars_fc_y_s);
  write_f4(pos.x, pos.y, pos.z, &Fz, bs_scalars_fc_z_s);
}
