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

__kernel void k_wrf_interpolate_momenta_kernel_main(__private parameters par,
                                                    __private wrfparameters wrf,
                                                    __read_only image3d_t b_source_velocites,
                                                    __read_only image3d_t b_scalars_vc,
                                                    __write_only image3d_t b_target_momenta)
{
  position pos = get_pos_bc(par, get_global_id(0), get_global_id(1), get_global_id(2));

  float dxl,dxr,dyl,dyr,dzl,dzr;
  float f_x, f_y, f_z, f_zl, f_zr;
  int  xwrf, ywrf, zwrf;

  zwrf = 0;
  while (pos.z*par.dz+wrf.hz[wrf.zoffset]>=wrf.hz[zwrf+wrf.zoffset]) zwrf++;
  zwrf--;

  // no 0.5dx shift to get fc not vc
  f_x = (wrf.dx0*wrf.dx + pos.x*par.dx)/wrf.dx;
  f_y = (wrf.dy0*wrf.dy + pos.y*par.dy)/wrf.dy;

  dxl = f_x-floor(f_x);
  dyl = f_y-floor(f_y);
  dxr = ceil(f_x)-f_x;
  dyr = ceil(f_y)-f_y;

  xwrf = (int)(floor(f_x));
  ywrf = (int)(floor(f_y));

  dzl = ((wrf.hz[wrf.zoffset]+pos.z*par.dz) - wrf.hz[zwrf+wrf.zoffset])
        /(wrf.hz[zwrf+wrf.zoffset+1] - wrf.hz[zwrf+wrf.zoffset]);
  dzr = (wrf.hz[zwrf+wrf.zoffset+1] - (wrf.hz[wrf.zoffset]+pos.z*par.dz))
        /(wrf.hz[zwrf+wrf.zoffset+1] - wrf.hz[zwrf+wrf.zoffset]);


  // interpolation done already
  float rho    = read_f4(pos.x , pos.y , pos.z , b_scalars_vc).s1;
  float rho_xl = read_f4(pos.xl, pos.y , pos.z , b_scalars_vc).s1;
  float rho_yl = read_f4(pos.x , pos.yl, pos.z , b_scalars_vc).s1;
  float rho_zl = read_f4(pos.x , pos.yl, pos.zl, b_scalars_vc).s1;

  float4 ip = interpolate(dxl, dyl, dzl, dxr, dyr, dzr, xwrf, ywrf, zwrf, b_source_velocites);

  float4 mom = (float4)(0.0f);

  if (pos.x>0) mom.s0 = (rho+rho_xl)*0.5f*ip.s0;
  else         mom.s0 = rho*ip.s0;
  if (pos.y>0) mom.s1 = (rho+rho_yl)*0.5f*ip.s1;
  else         mom.s1 = rho*ip.s1;
  // if (pos.z>0) mom.s2 = (rho+rho_zl)*0.5f*ip.s2;
  // else         mom.s2 = rho*ip.s2;
  mom.s2 = 0.0f;

  // if (pos.x == par.sx/2 && pos.y == par.sy/2) printf("%d | %.2f %.2f %.2f->%.2f\n", pos.z, rho, rho_xl, ip.s0,mom.s0);

  write_f4(pos.x, pos.y, pos.z, mom,  b_target_momenta);
}

