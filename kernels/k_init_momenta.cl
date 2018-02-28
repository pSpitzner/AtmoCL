__kernel void k_init_momenta_kernel_main(__private parameters par,
                                         __read_only image3d_t bf_scalars_vc_a_0,
                                         __write_only image3d_t bf_momenta_fc_a)
{

  //periodic bc
  position pos = get_pos_bc(par, get_global_id(0), get_global_id(1), get_global_id(2));

  // setting output
  float4 output = (float4)(0.0f,0.0f,0.0f,0.0f);
  // rho_profiled = 1.0f;
  float rho_xl, rho_yl, rho_zl, rho;

  rho    = read_imagef(bf_scalars_vc_a_0, (int4)(pos.x , pos.y , pos.z , 0)).y;
  rho_xl = read_imagef(bf_scalars_vc_a_0, (int4)(pos.xl, pos.y , pos.z , 0)).y;
  rho_yl = read_imagef(bf_scalars_vc_a_0, (int4)(pos.x , pos.yl, pos.z , 0)).y;
  rho_zl = read_imagef(bf_scalars_vc_a_0, (int4)(pos.x , pos.y , pos.zl, 0)).y;

  output.x = par.ui*0.5f*(rho_xl+rho);
  output.y = par.vi*0.5f*(rho_yl+rho);
  output.z = par.wi*0.5f*(rho_zl+rho);

  output = (float4)(0.0);

  write_imagef(bf_momenta_fc_a, (int4)(pos.x, pos.y, pos.z, 0), output);
}
