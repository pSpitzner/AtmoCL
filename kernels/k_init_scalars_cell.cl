__kernel void k_init_scalars_cell_kernel_main(__private parameters par,
                                              __write_only image3d_t bf_scalars_vc_a_0,
                                              __write_only image3d_t bf_scalars_vc_a_1,
                                              __write_only image3d_t bf_scalars_vc_a_2,
                                              __read_only image3d_t bf_scalars_vc_b_0,
                                              __read_only image3d_t bf_scalars_vc_b_1,
                                              __read_only image3d_t bf_scalars_vc_b_2)
{
position pos = get_pos_bc(&par);


// hydrostatic equilibrium
float z = pos.z*par.dz+0.5f*par.dy;
float p = par.pr*pow(1.0f-par.rd/par.cpd*par.gr*z/(par.rd*300.0f) , par.cpd/par.rd);
float rho = p/(pow(p/par.pr, par.rd/par.cpd)*par.rd*300.0f);

float s = par.cpd*log(300.0f) - par.rd*log(par.pr);

float8 output;
output.s0 = s*rho;      // rho*sig
output.s1 = rho;        // rho_total
output.s2 = 0.0f;      // rho_vapour
output.s3 = 0.0f;      // rho_cloud_liquid, (old) bulk: rho_liquid

output.s4 = 0.0f;       // rho_rain_liquid
output.s5 = 600.0e6f;   // n_dirt
output.s6 = 0.0f;       // n_cloud
output.s7 = 0.0f;       // n_rain

float4 output_ice = (float4)(0.0f);

write_f8(pos.x, pos.y, pos.z, &output, bf_scalars_vc_a_0, bf_scalars_vc_a_1);
write_f4(pos.x, pos.y, pos.z, &output_ice, bf_scalars_vc_a_2);
}

