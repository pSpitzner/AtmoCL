float rhovs(float T, parameters par) {
  // density of vapour at saturation
  float sv = par.svr*pow(T/par.tr, (par.cpv-par.cpl)/par.rv)*exp(par.lre0/par.rv*(1.0f/par.tr-1.0f/T));
  return sv/par.rv/T;
}

__kernel void k_init_scalars_dycoms_kernel_main(__private parameters par,
                                                __write_only image3d_t bf_scalars_vc_a_0,
                                                __write_only image3d_t bf_scalars_vc_a_1,
                                                __write_only image3d_t bf_scalars_vc_a_2,
                                                __read_only image3d_t bf_scalars_vc_b_0,
                                                __read_only image3d_t bf_scalars_vc_b_1,
                                                __read_only image3d_t bf_scalars_vc_b_2)
{
  position pos = get_pos_bc(&par);

  float len, r, offset, T;
  float4 arg;

  float8 c = read_f8(pos.x, pos.y, pos.z, bf_scalars_vc_b_0, bf_scalars_vc_b_1);
  float rhosig = c.s0;
  float rho    = c.s1;
  float rho_v  = c.s2;
  float rho_l  = c.s3;
  float rho_d  = rho-rho_v-rho_l;

  //random noise
  unsigned long seed = get_global_id(0) + get_global_id(1)*par.sx + get_global_id(2)*par.sx*par.sy;
  // if (pos.z <= 200.0f/par.dz) rhosig-=(rhosig/1.0e7f)*(rand(seed)-0.5f);
  // if (pos.z == 1 && pos.x == 32 && pos.y == 0) rhosig-=(rhosig/1.0e7f);
  // if (pos.z == 1 && pos.x == 31 && pos.y == 0) rhosig-=(rhosig/1.0e7f);
  // printf("%2.2f\t%d\n", rand(seed), seed);

  float8 output;
  output.s0 = rhosig;     // rho*sig
  output.s1 = rho;        // rho_total
  output.s2 = rho_v;      // rho_vapour
  output.s3 = rho_l;      // (rho_c in microphys) rho_cloud_liquid, (old) bulk: rho_liquid

  output.s4 = 0.0f;       // rho_rain_liquid
  output.s5 = 60.0e6f;    // n_dirt
  output.s6 = (rho_l/8.5e-4f)*55.0e6f;       // n_cloud
  output.s7 = 0.0f;       // n_rain

  float4 output_ice;
  output_ice.s0 = 0.0f;   // rho_ice
  output_ice.s1 = 0.0f;   // rho_snow
  output_ice.s2 = 5e6f;   // n_ice
  output_ice.s3 = 5e6f;   // n_snow

  // test advection of new texture
  if (pos.z == 0 || pos.z == 1) {
    output_ice.s2 = 20e6f;
    output_ice.s3 = 30e6f;
  }


  write_f8(pos.x, pos.y, pos.z, output,     bf_scalars_vc_a_0, bf_scalars_vc_a_1);
  write_f4(pos.x, pos.y, pos.z, output_ice, bf_scalars_vc_a_2);
}

