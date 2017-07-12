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

__kernel void k_wrf_interpolate_kernel_main(__private parameters par,
                                            __private wrfparameters wrf,
                                            __read_only image3d_t b_wrf_srouce_0,
                                            __read_only image3d_t b_wrf_srouce_1,
                                            __read_only image3d_t b_wrf_srouce_2,
                                            __write_only image3d_t b_target_scalars_0,
                                            __write_only image3d_t b_target_scalars_1,
                                            __write_only image3d_t b_target_scalars_2,
                                            __write_only image3d_t b_target_momenta)
{
  position pos = get_pos_bc(par, get_global_id(0), get_global_id(1), get_global_id(2));
  // pos in target domain

  float8 c;
  float4 cice, mom;

  c = (float8)(0.0f);
  cice = (float8)(0.0f);
  mom = (float8)(0.0f);

  float dxl,dxr,dyl,dyr,dzl,dzr;
  float u,v,w,mytheta,rho,rhov,rhoc,nv,nc;
  float f_x, f_y, f_z, f_zl, f_zr;
  int  xwrf, ywrf, zwrf;

  zwrf = 0;
  while (pos.z*par.dz+wrf.hz[zoffset]>=hz[zwrf+zoffset]) zwrf++;
  zwrf--;

  f_x = (wrf.dx0*wrf.dx + pos.x*par.dx + 0.5*par.dx)/wrf.dx;
  f_y = (wrf.dy0*wrf.dy + pos.y*par.dy + 0.5*par.dy)/wrf.dy;

  dxl = f_x-floor(f_x);
  dyl = f_y-floor(f_y);
  dxr = ceil(f_x)-f_x;
  dyr = ceil(f_y)-f_y;

  xwrf = (int)(floor(f_x));
  ywrf = (int)(floor(f_y));

  dzl = ((wrf.hz[zoffset]+pos.z*par.dz+0.5*par.dz) - wrf.hz[zwrf+zoffset])
        /(wrf.hz[zwrf+zoffset+1] - wrf.hz[zwrf+zoffset]);
  dzr = (wrf.hz[zwrf+zoffset+1] - (wrf.hz[zoffset]+pos.z*par.dz+0.5*par.dz))
        /(wrf.hz[zwrf+zoffset+1] - wrf.hz[zwrf+zoffset]);

  float4 c_0, c_1, c_2;

  c_0 = interpolate(dxl, dyl, dzl, dxr, dyr, dzr, xwrf, ywrf, zwrf, b_wrf_srouce_0);
  c_1 = interpolate(dxl, dyl, dzl, dxr, dyr, dzr, xwrf, ywrf, zwrf, b_wrf_srouce_1);
  c_2 = interpolate(dxl, dyl, dzl, dxr, dyr, dzr, xwrf, ywrf, zwrf, b_wrf_srouce_2);


  write_f8(pos.x, pos.y, pos.z, c,    b_target_scalars_0, b_target_scalars_1);
  write_f4(pos.x, pos.y, pos.z, cice, b_target_scalars_2);
  write_f4(pos.x, pos.y, pos.z, mom,  b_target_momenta);
}

