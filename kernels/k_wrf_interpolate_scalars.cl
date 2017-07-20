__private float4 interpolate(float dxl, float dyl, float dzl, float dxr, float dyr, float dzr, int xwrf, int ywrf, int zwrf, __read_only image3d_t buf) {

  return ((1.0f-dxl)*(1.0f-dyl)*(1.0f-dzl)*read_f4(xwrf  ,ywrf  ,zwrf  , buf)+
          (1.0f-dxr)*(1.0f-dyl)*(1.0f-dzl)*read_f4(xwrf+1,ywrf  ,zwrf  , buf)+
          (1.0f-dxl)*(1.0f-dyr)*(1.0f-dzl)*read_f4(xwrf  ,ywrf+1,zwrf  , buf)+
          (1.0f-dxr)*(1.0f-dyr)*(1.0f-dzl)*read_f4(xwrf+1,ywrf+1,zwrf  , buf)+
          (1.0f-dxl)*(1.0f-dyl)*(1.0f-dzr)*read_f4(xwrf  ,ywrf  ,zwrf+1, buf)+
          (1.0f-dxr)*(1.0f-dyl)*(1.0f-dzr)*read_f4(xwrf+1,ywrf  ,zwrf+1, buf)+
          (1.0f-dxl)*(1.0f-dyr)*(1.0f-dzr)*read_f4(xwrf  ,ywrf+1,zwrf+1, buf)+
          (1.0f-dxr)*(1.0f-dyr)*(1.0f-dzr)*read_f4(xwrf+1,ywrf+1,zwrf+1, buf));

}

// problem with our string replace in kernel code:
// if one swaps the order of arguments to kernel main, 'b_wrf_source_scalars' will be found when searching for 'wrf'
// binding to 'wrfparameters wrf' instead solves that. Maybe implement a bind with additional 'type' argument

__kernel void k_wrf_interpolate_scalars_kernel_main(__private parameters par,
                                                    __private wrfparameters wrf,
                                                    __read_only image3d_t b_wrf_source_scalars_0,
                                                    __read_only image3d_t b_wrf_source_scalars_1,
                                                    __write_only image3d_t b_target_scalars_0,
                                                    __write_only image3d_t b_target_scalars_1,
                                                    __write_only image3d_t b_target_scalars_2)
{
  position pos = get_pos_bc(par, get_global_id(0), get_global_id(1), get_global_id(2));

  float dxl,dxr,dyl,dyr,dzl,dzr;
  float f_x, f_y, f_z, f_zl, f_zr;
  int  xwrf, ywrf, zwrf;

  zwrf = 0;
  while (wrf.hz[zwrf]<=pos.z*par.dz+0.5*par.dz+wrf.hz[wrf.zoffset]) {
    // if ( pos.x == par.sx/2 && pos.y == par.sy/2 && pos.z == 1) printf("%d %f %f\n", zwrf, wrf.hz[zwrf],pos.z*par.dz+wrf.hz[wrf.zoffset]);
    zwrf++;
  }
  zwrf--;
  // if ( pos.x == par.sx/2 && pos.y == par.sy/2 && pos.z == 1) printf("%f %f %f\n", wrf.hz[zwrf], wrf.hz[wrf.zoffset] + pos.z*par.dz+0.5*par.dz, wrf.hz[zwrf+1]);

  // printf("lat0: %.2f\n", wrf.lat0);
  // printf("lon0: %.2f\n", wrf.lon0);
  // printf("dlat: %.2f\n", wrf.dlat);
  // printf("dlon: %.2f\n", wrf.dlon);
  // printf("dx: %.2f\n", wrf.dx);
  // printf("dy: %.2f\n", wrf.dy);
  // printf("domaincenterx: %.2f\n", wrf.domaincenterx);
  // printf("domaincentery: %.2f\n", wrf.domaincentery);
  // printf("domainsizex: %.2f | %.2f\n", wrf.domainsizex, wrf.dx*wrf.sx);
  // printf("domainsizey: %.2f | %.2f\n", wrf.domainsizey, wrf.dy*wrf.sy);
  // printf("dsx: %.2f\n", wrf.dsx);
  // printf("dsy: %.2f\n", wrf.dsy);
  // printf("dx0: %.2f\n", wrf.dx0);
  // printf("dy0: %.2f\n", wrf.dy0);

  f_x = (wrf.dx0*wrf.dx + pos.x*par.dx + 0.5*par.dx)/wrf.dx;
  f_y = (wrf.dy0*wrf.dy + pos.y*par.dy + 0.5*par.dy)/wrf.dy;

  dxl = f_x-floor(f_x);
  dyl = f_y-floor(f_y);
  dxr = ceil(f_x)-f_x;
  dyr = ceil(f_y)-f_y;

  xwrf = (int)(floor(f_x));
  ywrf = (int)(floor(f_y));

  dzl = ((wrf.hz[wrf.zoffset]+pos.z*par.dz+0.5*par.dz) - wrf.hz[zwrf])
        /(wrf.hz[zwrf+1] - wrf.hz[zwrf]);
  dzr = (wrf.hz[zwrf+1] - (wrf.hz[wrf.zoffset]+pos.z*par.dz+0.5*par.dz))
        /(wrf.hz[zwrf+1] - wrf.hz[zwrf]);

  float4 ip_0, ip_1;

  ip_0 = interpolate(dxl, dyl, dzl, dxr, dyr, dzr, xwrf, ywrf, zwrf, b_wrf_source_scalars_0);
  ip_1 = interpolate(dxl, dyl, dzl, dxr, dyr, dzr, xwrf, ywrf, zwrf, b_wrf_source_scalars_1);

  // printf("pos: %.2f(+%.2f) %.2f(+%.2f) %.2f(+%.2f+%.2f)\n", pos.x*par.dx,0.5f*par.dx, pos.y*par.dy,0.5f*par.dx, pos.z*par.dz,0.5f*par.dx, wrf.hz[wrf.zoffset]);
  // printf("hz: %.2f %.2f\n", wrf.hz[zwrf+wrf.zoffset+1], wrf.hz[zwrf+wrf.zoffset]);
  // if (pos.x == par.sx/2 && pos.y == par.sy/2) printf("%.1f %.2f %.2f\n", pos.z*par.dz, dzl, dzr);
  // printf("f_x|y: %f %f\n", f_x, f_y);
  // printf("x|y|zwrf: %d %d %d\n", xwrf, ywrf,zwrf);
  // printf("s_i: %.2f\nfrom:\n", ip_0.s0);
  // //test
  if ((dzl < 0.0 || dzl > 1.0 )) printf("%4d %4d %4d | %4d %4d %4d | %.2f %.2f %.2f | %.2f\n", xwrf, ywrf, zwrf, pos.x, pos.y, pos.z, (1.0f-dxl),(1.0f-dyl),(1.0f-dzl), read_f4(xwrf  ,ywrf  ,zwrf  , b_wrf_source_scalars_0).s0);
  // if ( pos.x == par.sx/2 && pos.y == par.sy/2 && (dzl < 0.0 || dzl > 1.0 )) printf("%4d %4d %4d | %4d %4d %4d | %.2f %.2f | %.2f %.2f %.2f\n", xwrf, ywrf, zwrf, pos.x, pos.y, pos.z, dzl, dzr, (wrf.hz[wrf.zoffset]+pos.z*par.dz+0.5*par.dz), wrf.hz[zwrf+wrf.zoffset], wrf.hz[zwrf+wrf.zoffset+1]);
  // printf("%4d %4d %4d | %4d %4d %4d | %.2f %.2f %.2f | %.2f\n", xwrf, ywrf, zwrf, pos.x, pos.y, pos.z, (1.0f-dxr),(1.0f-dyl),(1.0f-dzl), read_f4(xwrf+1,ywrf  ,zwrf  , b_wrf_source_scalars_0).s0);
  // printf("%4d %4d %4d | %4d %4d %4d | %.2f %.2f %.2f | %.2f\n", xwrf, ywrf, zwrf, pos.x, pos.y, pos.z, (1.0f-dxl),(1.0f-dyr),(1.0f-dzl), read_f4(xwrf  ,ywrf+1,zwrf  , b_wrf_source_scalars_0).s0);
  // printf("%4d %4d %4d | %4d %4d %4d | %.2f %.2f %.2f | %.2f\n", xwrf, ywrf, zwrf, pos.x, pos.y, pos.z, (1.0f-dxr),(1.0f-dyr),(1.0f-dzl), read_f4(xwrf+1,ywrf+1,zwrf  , b_wrf_source_scalars_0).s0);
  // printf("%4d %4d %4d | %4d %4d %4d | %.2f %.2f %.2f | %.2f\n", xwrf, ywrf, zwrf, pos.x, pos.y, pos.z, (1.0f-dxl),(1.0f-dyl),(1.0f-dzr), read_f4(xwrf  ,ywrf  ,zwrf+1, b_wrf_source_scalars_0).s0);
  // printf("%4d %4d %4d | %4d %4d %4d | %.2f %.2f %.2f | %.2f\n", xwrf, ywrf, zwrf, pos.x, pos.y, pos.z, (1.0f-dxr),(1.0f-dyl),(1.0f-dzr), read_f4(xwrf+1,ywrf  ,zwrf+1, b_wrf_source_scalars_0).s0);
  // printf("%4d %4d %4d | %4d %4d %4d | %.2f %.2f %.2f | %.2f\n", xwrf, ywrf, zwrf, pos.x, pos.y, pos.z, (1.0f-dxl),(1.0f-dyr),(1.0f-dzr), read_f4(xwrf  ,ywrf+1,zwrf+1, b_wrf_source_scalars_0).s0);
  // printf("%4d %4d %4d | %4d %4d %4d | %.2f %.2f %.2f | %.2f\n", xwrf, ywrf, zwrf, pos.x, pos.y, pos.z, (1.0f-dxr),(1.0f-dyr),(1.0f-dzr), read_f4(xwrf+1,ywrf+1,zwrf+1, b_wrf_source_scalars_0).s0);

  float8 c;
  float4 cice;

  c = (float8)(0.0f);
  cice = (float4)(0.0f);

  c.s0 = ip_0.s0;
  c.s1 = ip_0.s1;
  c.s2 = ip_0.s2;
  c.s3 = ip_0.s3;
  c.s4 = ip_1.s0;
  c.s5 = ip_1.s1;
  c.s6 = ip_1.s2;
  c.s7 = ip_1.s3;
  cice.s0 = 0.0f;
  cice.s1 = 0.0f;
  cice.s2 = 0.0f;
  cice.s3 = 0.0f;


  write_f8(pos.x, pos.y, pos.z, c,    b_target_scalars_0, b_target_scalars_1);
  write_f4(pos.x, pos.y, pos.z, cice, b_target_scalars_2);
}

