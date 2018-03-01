__kernel void ks_adv_momenta_kernel_main(__private parameters par,
                                         __read_only image3d_t bf_scalars_vc_a_0,
                                         __read_only image3d_t bf_momenta_fc_a,
                                         __read_only image3d_t b_m_vc_uv,
                                         __read_only image3d_t b_m_vc_w,
                                         __write_only image3d_t bRhs_m_vc_uv,
                                         __write_only image3d_t bRhs_m_vc_w)
{
  position pos = get_pos_bc(par, get_global_id(0), get_global_id(1), get_global_id(2));

  float8 c, cxl, cxr, cyl, cyr, czl, czr;
  float8 cxll, cxrr, cyll, cyrr, czll, czrr;
  float8 Fxlp, Fxrp, Fxln, Fxrn; // upwind: Flux_x_left_positive, Flux_x_right_positive
  float8 Fylp, Fyrp, Fyln, Fyrn;
  float8 Fzlp, Fzrp, Fzln, Fzrn;
  float8 Fx, Fy, Fz;

  float rho    = read_imagef(bf_scalars_vc_a_0, (int4)(pos.x , pos.y , pos.z , 0)).s1;
  float rho_xl = read_imagef(bf_scalars_vc_a_0, (int4)(pos.xl, pos.y , pos.z , 0)).s1;
  float rho_xr = read_imagef(bf_scalars_vc_a_0, (int4)(pos.xr, pos.y , pos.z , 0)).s1;
  float rho_yl = read_imagef(bf_scalars_vc_a_0, (int4)(pos.x , pos.yl, pos.z , 0)).s1;
  float rho_yr = read_imagef(bf_scalars_vc_a_0, (int4)(pos.x , pos.yr, pos.z , 0)).s1;
  float rho_zl = read_imagef(bf_scalars_vc_a_0, (int4)(pos.x , pos.y , pos.zl, 0)).s1;
  float rho_zr = read_imagef(bf_scalars_vc_a_0, (int4)(pos.x , pos.y , pos.zr, 0)).s1;

  // momenta are saved in channels, need to obtain velocities
  pos.ulf = read_imagef(bf_momenta_fc_a, (int4)(pos.x , pos.y , pos.z , 0)).s0*2.0f/(rho_xl+rho);
  pos.urf = read_imagef(bf_momenta_fc_a, (int4)(pos.xr, pos.y , pos.z , 0)).s0*2.0f/(rho_xr+rho);
  pos.vlf = read_imagef(bf_momenta_fc_a, (int4)(pos.x , pos.y , pos.z , 0)).s1*2.0f/(rho_yl+rho);
  pos.vrf = read_imagef(bf_momenta_fc_a, (int4)(pos.x , pos.yr, pos.z , 0)).s1*2.0f/(rho_yr+rho);
  pos.wlf = read_imagef(bf_momenta_fc_a, (int4)(pos.x , pos.y , pos.z , 0)).s2*2.0f/(rho_zl+rho);
  pos.wrf = read_imagef(bf_momenta_fc_a, (int4)(pos.x , pos.y , pos.zr, 0)).s2*2.0f/(rho_zr+rho);

  // fixed bc
  if (pos.s_xr==-1) pos.urf=0.0f;
  if (pos.s_xl==-1) pos.ulf=0.0f;
  if (pos.s_yr==-1) pos.vrf=0.0f;
  if (pos.s_yl==-1) pos.vlf=0.0f;
  if (pos.s_zr==-1) pos.wrf=0.0f;
  if (pos.s_zl==-1) pos.wlf=0.0f;

  // float8 foo = (float8)(float4, float4)
  // advection of momenta with velocities
  c   =  read_f8(pos.x  , pos.y  , pos.z  , b_m_vc_uv, b_m_vc_w);
  cxl =  read_f8(pos.xl , pos.y  , pos.z  , b_m_vc_uv, b_m_vc_w);
  cxr =  read_f8(pos.xr , pos.y  , pos.z  , b_m_vc_uv, b_m_vc_w);
  cyl =  read_f8(pos.x  , pos.yl , pos.z  , b_m_vc_uv, b_m_vc_w);
  cyr =  read_f8(pos.x  , pos.yr , pos.z  , b_m_vc_uv, b_m_vc_w);
  czl =  read_f8(pos.x  , pos.y  , pos.zl , b_m_vc_uv, b_m_vc_w);
  czr =  read_f8(pos.x  , pos.y  , pos.zr , b_m_vc_uv, b_m_vc_w);

  cxll = read_f8(pos.xll, pos.y  , pos.z  , b_m_vc_uv, b_m_vc_w);
  cxrr = read_f8(pos.xrr, pos.y  , pos.z  , b_m_vc_uv, b_m_vc_w);
  cyll = read_f8(pos.x  , pos.yll, pos.z  , b_m_vc_uv, b_m_vc_w);
  cyrr = read_f8(pos.x  , pos.yrr, pos.z  , b_m_vc_uv, b_m_vc_w);
  czll = read_f8(pos.x  , pos.y  , pos.zll, b_m_vc_uv, b_m_vc_w);
  czrr = read_f8(pos.x  , pos.y  , pos.zrr, b_m_vc_uv, b_m_vc_w);

  // eg. (float2)(0.0f) gives (0.0f, 0.0f)
  // upwind scheme ansatz
  Fxlp = (float8)(0.0f);  Fylp = (float8)(0.0f);  Fzlp = (float8)(0.0f);
  Fxrp = (float8)(0.0f);  Fyrp = (float8)(0.0f);  Fzrp = (float8)(0.0f);
  Fxln = (float8)(0.0f);  Fyln = (float8)(0.0f);  Fzln = (float8)(0.0f);
  Fxrn = (float8)(0.0f);  Fyrn = (float8)(0.0f);  Fzrn = (float8)(0.0f);


  if (pos.ulf>0.0f) Fxlp = (  c/3.0f+5.0f/6.0f*cxl-cxll/6.0f)*pos.ulf;
  if (pos.urf>0.0f) Fxrp = (cxr/3.0f+5.0f/6.0f*  c- cxl/6.0f)*pos.urf;
  if (pos.ulf<0.0f) Fxln = (cxl/3.0f+5.0f/6.0f*  c- cxr/6.0f)*pos.ulf;
  if (pos.urf<0.0f) Fxrn = (  c/3.0f+5.0f/6.0f*cxr-cxrr/6.0f)*pos.urf;
  //Fx = (Fxlp - Fxrp + Fxln - Fxrn) / par.dx; /* negative flux *-1 */
  Fx = (Fxlp/par.dx - Fxrp/par.dx + Fxln/par.dx - Fxrn/par.dx) ; /* negative flux *-1 */

  if (pos.vlf>0.0f) Fylp = (  c/3.0f+5.0f/6.0f*cyl-cyll/6.0f)*pos.vlf;
  if (pos.vrf>0.0f) Fyrp = (cyr/3.0f+5.0f/6.0f*  c- cyl/6.0f)*pos.vrf;
  if (pos.vlf<0.0f) Fyln = (cyl/3.0f+5.0f/6.0f*  c- cyr/6.0f)*pos.vlf;
  if (pos.vrf<0.0f) Fyrn = (  c/3.0f+5.0f/6.0f*cyr-cyrr/6.0f)*pos.vrf;
  //Fy = (Fylp - Fyrp + Fyln - Fyrn) / par.dy;
  Fy = (Fylp/par.dz - Fyrp/par.dz + Fyln/par.dz - Fyrn/par.dz);

  if (pos.wlf>0.0f) Fzlp = (  c/3.0f+5.0f/6.0f*czl-czll/6.0f)*pos.wlf;
  if (pos.wrf>0.0f) Fzrp = (czr/3.0f+5.0f/6.0f*  c- czl/6.0f)*pos.wrf;
  if (pos.wlf<0.0f) Fzln = (czl/3.0f+5.0f/6.0f*  c- czr/6.0f)*pos.wlf;
  if (pos.wrf<0.0f) Fzrn = (  c/3.0f+5.0f/6.0f*czr-czrr/6.0f)*pos.wrf;
  //Fz = (Fzlp - Fzrp + Fzln - Fzrn) / par.dz;
  Fz = (Fzlp/par.dz - Fzrp/par.dz + Fzln/par.dz - Fzrn/par.dz);
  // end upwind

  float8 Temp = Fx + Fy + Fz;
  // should be zero anyway
  Temp.s6 = 0.0f;
  Temp.s7 = 0.0f;

  write_f8(pos.x, pos.y, pos.z, Temp, bRhs_m_vc_uv, bRhs_m_vc_w);
}
