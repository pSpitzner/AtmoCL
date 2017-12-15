__kernel void kf_pressure_kernel_main(__private parameters par,
                                      __read_only image3d_t bf_scalars_vc_a_0,
                                      __read_only image3d_t bf_scalars_vc_a_1,
                                      __read_only image3d_t bf_scalars_vc_a_2,
                                      __write_only image3d_t bRhs_p_fc)
{
  position pos = get_pos_bc(&par);

  // ----------------------------------------------------------------- //
  // ice needs treatement, fields are already linked //
  // ----------------------------------------------------------------- //

  float8 c   = read_f8(pos.x , pos.y , pos.z , bf_scalars_vc_a_0, bf_scalars_vc_a_1);
  float8 cxl = read_f8(pos.xl, pos.y , pos.z , bf_scalars_vc_a_0, bf_scalars_vc_a_1);
  float8 cyl = read_f8(pos.x , pos.yl, pos.z , bf_scalars_vc_a_0, bf_scalars_vc_a_1);
  float8 czl = read_f8(pos.x , pos.y , pos.zl, bf_scalars_vc_a_0, bf_scalars_vc_a_1);

  float4 c_ice   = read_f4(pos.x , pos.y , pos.z , bf_scalars_vc_a_2);
  float4 cxl_ice = read_f4(pos.xl, pos.y , pos.z , bf_scalars_vc_a_2);
  float4 cyl_ice = read_f4(pos.x , pos.yl, pos.z , bf_scalars_vc_a_2);
  float4 czl_ice = read_f4(pos.x , pos.y , pos.zl, bf_scalars_vc_a_2);

  float rhosig    = c  .s0;
  float rhosig_xl = cxl.s0;
  float rhosig_yl = cyl.s0;
  float rhosig_zl = czl.s0;

  float rho       = c  .s1;
  float rho_xl    = cxl.s1;
  float rho_yl    = cyl.s1;
  float rho_zl    = czl.s1;

  float rho_v     = c  .s2;
  float rho_v_xl  = cxl.s2;
  float rho_v_yl  = cyl.s2;
  float rho_v_zl  = czl.s2;

  float rho_c     = c  .s3;
  float rho_c_xl  = cxl.s3;
  float rho_c_yl  = cyl.s3;
  float rho_c_zl  = czl.s3;

  float rho_r     = c  .s4;
  float rho_r_xl  = cxl.s4;
  float rho_r_yl  = cyl.s4;
  float rho_r_zl  = czl.s4;

  float rho_i     =   c_ice.s0;
  float rho_i_xl  = cxl_ice.s0;
  float rho_i_yl  = cyl_ice.s0;
  float rho_i_zl  = czl_ice.s0;

  float rho_s     =   c_ice.s1;
  float rho_s_xl  = cxl_ice.s1;
  float rho_s_yl  = cyl_ice.s1;
  float rho_s_zl  = czl_ice.s1;

  float rho_l    = rho_c    + rho_r   ;
  float rho_l_xl = rho_c_xl + rho_r_xl;
  float rho_l_yl = rho_c_yl + rho_r_yl;
  float rho_l_zl = rho_c_zl + rho_r_zl;

  float rho_d    = rho    - rho_v    - rho_l    - rho_i    - rho_s   ;
  float rho_d_xl = rho_xl - rho_v_xl - rho_l_xl - rho_i_xl - rho_s_xl;
  float rho_d_yl = rho_yl - rho_v_yl - rho_l_yl - rho_i_yl - rho_s_yl;
  float rho_d_zl = rho_zl - rho_v_zl - rho_l_zl - rho_i_zl - rho_s_zl;

  float cpml    = rho_d   *par.cpd+rho_v   *par.cpv+rho_l   *par.cpl+rho_i   *par.cpi+rho_s   *par.cpi;
  float cpml_xl = rho_d_xl*par.cpd+rho_v_xl*par.cpv+rho_l_xl*par.cpl+rho_i_xl*par.cpi+rho_s_xl*par.cpi;
  float cpml_yl = rho_d_yl*par.cpd+rho_v_yl*par.cpv+rho_l_yl*par.cpl+rho_i_yl*par.cpi+rho_s_yl*par.cpi;
  float cpml_zl = rho_d_zl*par.cpd+rho_v_zl*par.cpv+rho_l_zl*par.cpl+rho_i_zl*par.cpi+rho_s_zl*par.cpi;

  float rml    = rho_d   *par.rd+rho_v   *par.rv;
  float rml_xl = rho_d_xl*par.rd+rho_v_xl*par.rv;
  float rml_yl = rho_d_yl*par.rd+rho_v_yl*par.rv;
  float rml_zl = rho_d_zl*par.rd+rho_v_zl*par.rv;

  // gravity

  float Gx = 0.0f;
  float Gy = 0.0f;
  float Gz = 0.0f;
  if (pos.s_zl==1) Gz = -(rho+rho_zl)*0.5f*par.gr;

  float Px = 0.0f;
  float Py = 0.0f;
  float Pz = 0.0f;
  float lnPs    = rhosig   /(cpml   -rml   );
  float lnPs_xl = rhosig_xl/(cpml_xl-rml_xl);
  float lnPs_yl = rhosig_yl/(cpml_yl-rml_yl);
  float lnPs_zl = rhosig_zl/(cpml_zl-rml_zl);
  float lnPr    = log(rml   )/(1.0f-rml   /cpml   );
  float lnPr_xl = log(rml_xl)/(1.0f-rml_xl/cpml_xl);
  float lnPr_yl = log(rml_yl)/(1.0f-rml_yl/cpml_yl);
  float lnPr_zl = log(rml_zl)/(1.0f-rml_zl/cpml_zl);
  float lnP_x   = (lnPs+lnPr+lnPs_xl+lnPr_xl)/2.0f;
  float lnP_y   = (lnPs+lnPr+lnPs_yl+lnPr_yl)/2.0f;
  float lnP_z   = (lnPs+lnPr+lnPs_zl+lnPr_zl)/2.0f;

  if (pos.s_xl==1) {
    Px = -((lnPs-lnPs_xl)+(lnPr-lnPr_xl))*exp(lnP_x)/par.dx;
  }
  if (pos.s_yl==1) {
    Py = -((lnPs-lnPs_yl)+(lnPr-lnPr_yl))*exp(lnP_y)/par.dy;
  }
  if (pos.s_zl==1) {
    Pz = -((lnPs-lnPs_zl)+(lnPr-lnPr_zl))*exp(lnP_z)/par.dz;
  }

  // final sources
  float Sx = Gx+Px;
  float Sy = Gy+Py;
  float Sz = Gz+Pz;
  // Sx = 0.0f;
  // Sy = 0.0f;
  // Sz = 0.0f;
  // if (pos.x == (int)(par.sx*0.5f) && pos.y == (int)(par.sy*0.5f)) printf("p: %2.2f %2.2f %2.7f\t\t%2.2f\n", Sx, Sy, Sz, rhosig);

  write_imagef(bRhs_p_fc, (int4)(pos.x, pos.y, pos.z, 0), (float4)(Sx, Sy, Sz, 0.0f));
}
