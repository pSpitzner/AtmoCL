float rhovs(float T, parameters par) {
  // density of vapour at saturation
  float sv = par.svr*pow(T/par.tr, (par.cpv-par.cpl)/par.rv)*exp(par.lre0/par.rv*(1.0f/par.tr-1.0f/T));
  return sv/par.rv/T;
}

__kernel void k_perturb_moistbubble_kernel_main(__private parameters par,
                                                __read_only image3d_t b_source_scalars_0,
                                                __read_only image3d_t b_source_scalars_1,
                                                __read_only image3d_t b_source_scalars_2,
                                                __read_only image3d_t b_source_momenta,
                                                __write_only image3d_t b_target_momenta,
                                                __write_only image3d_t b_target_scalars_0,
                                                __write_only image3d_t b_target_scalars_1,
                                                __write_only image3d_t b_target_scalars_2)
{
  position pos = get_pos_bc(par, get_global_id(0), get_global_id(1), get_global_id(2));

  float len, r, offset, T;
  float4 arg;

  float8 c = read_f8(pos.x, pos.y, pos.z, b_source_scalars_0, b_source_scalars_1);
  float4 cice = read_f4(pos.x, pos.y, pos.z, b_source_scalars_2);
  float rhosig = c.s0;
  float rho    = c.s1;
  float rho_v  = c.s2;
  float rho_l  = c.s3;
  float rho_d  = rho-rho_v-rho_l;

  // if (pos.x == 0 && pos.y == 0) printf("foo %v8hlf\n", c);

  float cpml = rho_d*par.cpd+rho_v*par.cpv+rho_l*par.cpl;
  float rml  = rho_d*par.rd+rho_v*par.rv;
  float theta = exp(rhosig/cpml + rml/cpml*log(par.pr));
  float p = exp(rhosig/(cpml-rml)+log(rml)/(1.0f-rml/cpml));

  int i = -1;
  if (pos.x==par.sx/2 && pos.y== par.sy/2 && pos.z==2) printf("%d %2.2f %2.2f %2.2f %2.2f || %2.2f %2.2f %2.2f %2.2f\n", i, rho_d, rho_v, rho_l, rho, rml, cpml, theta, rhosig);

  // printf("%f %f %f\n", rhosig, p, theta);

  //cosine squared, r in meters
  r = 2000.0f;
  arg.s1 = ((pos.y)*par.dy/r-(par.sy*0.5f)*par.dy/r);
  arg.s0 = ((pos.x)*par.dx/r-(par.sx*0.5f)*par.dx/r);
  arg.s2 = ((pos.z)*par.dz/r-r                   /r);
  arg.s3 = 0.0f;

  float theta_vp;
  offset = 0.0f;
  len = length(arg);
  if (len < 1.0f) {
    offset = 2.0f*pow(cospi(len*0.5f), 2.0f);
    // change to dry potential temperature
    theta = theta*pow(p/par.pr,rml/cpml-par.rd/par.cpd);

    // (virtual) density potential temperature
    theta_vp = theta*(1.0f+(par.rv/par.rd)*rho_v/rho)/(1.0f+(rho_v+rho_l)/rho);

    // 320 due to background equil. initialization
    theta_vp = theta_vp*(1.0f+offset/320.0f);

    // back to 'normal' potential temp
    float theta_temp = theta_vp*(1.0f+(rho_v+rho_l)/rho)/(1.0f+(par.rv/par.rd)*rho_v/rho);
    // printf("1 (%d %d %d) %f | %f\n", pos.x, pos.y, pos.z, theta, theta_temp);
    theta = theta_temp;

    // change densities to keep pressure of bubble and surrounding constant
    for (int i = 0; i < 200; i++) {
      // if (pos.x==par.sx/2 && pos.y== par.sy/2 && pos.z==2) printf("%d %2.2f %2.2f %2.2f %2.2f || %2.2f %2.2f %2.2f\n", i, rho_d, rho_v, rho_l, rho, rml, cpml, theta);
      T     = theta*pow(p/par.pr,par.rd/par.cpd);
      rho_v = rhovs(T, par);
      rho_d = (p-rho_v*par.rv*T)/par.rd/T;
      rho_l = 0.02f*rho_d-rho_v;
      rho   = rho_d+rho_v+rho_l;
      rml   = rho_d*par.rd+rho_v*par.rv;
      cpml  = rho_d*par.cpd+rho_v*par.cpv+rho_l*par.cpl;
      theta = theta_vp*(1.0f+(rho_v+rho_l)/rho)/(1.0f+(par.rv/par.rd)*rho_v/rho);
    }

    // back to moist potential temperature
    theta = theta*pow(p/par.pr,par.rd/par.cpd-rml/cpml);
    // cpml i actually rho*cpml...
    rhosig = cpml*log(theta)-rml*log(par.pr);
  }

  float pnew = exp(rhosig/(cpml-rml)+log(rml)/(1.0f-rml/cpml));
  // if (pos.x==par.sx/2 && pos.y== par.sy/2) printf("%d %2.5f %2.5f %2.7f\t\t%2.2v4hlf\t%2.2f %2.2f\n", pos.z, theta, offset, rhosig, arg, p, pnew);

  float8 output;
  output.s0 = rhosig;     // rho*sig
  output.s1 = rho;        // rho_total
  output.s2 = rho_v;      // rho_vapour
  output.s3 = rho_l;      // rho_cloud_liquid, (old) bulk: rho_liquid

  output.s4 = c.s4;       // rho_rain_liquid
  output.s5 = c.s5;       // n_dirt
  output.s6 = c.s6;       // n_cloud
  output.s7 = c.s7;       // n_rain

  write_f8(pos.x, pos.y, pos.z, output, b_target_scalars_0, b_target_scalars_1);
  write_f4(pos.x, pos.y, pos.z, cice, b_target_scalars_2);

  float4 mom = read_f4(pos.x, pos.y, pos.z, b_source_momenta);
  write_f4(pos.x, pos.y, pos.z, mom, b_target_momenta);
}

