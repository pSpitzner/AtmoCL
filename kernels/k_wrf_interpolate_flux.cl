__private float4 interpolate(float dxl, float dyl, float dxr, float dyr, int xwrf, int ywrf, __read_only image3d_t buf) {

  return ((1.0f-dxl)*(1.0f-dyl)*read_f4(xwrf  ,ywrf  , 0, buf)+
          (1.0f-dxr)*(1.0f-dyl)*read_f4(xwrf+1,ywrf  , 0, buf)+
          (1.0f-dxl)*(1.0f-dyr)*read_f4(xwrf  ,ywrf+1, 0, buf)+
          (1.0f-dxr)*(1.0f-dyr)*read_f4(xwrf+1,ywrf+1, 0, buf));

}

__kernel void k_wrf_interpolate_flux_kernel_main(__private parameters par,
                                                 __private wrfparameters wrf,
                                                 __read_only image3d_t b_wrf_flux,
                                                 __write_only image3d_t b_sys_flux)
{
  position pos = get_pos_bc(par, get_global_id(0), get_global_id(1), get_global_id(2));

  float dxl,dxr,dyl,dyr;
  float f_x, f_y;
  int  xwrf, ywrf;

  f_x = (wrf.dx0*wrf.dx + pos.x*par.dx + 0.5*par.dx)/wrf.dx;
  f_y = (wrf.dy0*wrf.dy + pos.y*par.dy + 0.5*par.dy)/wrf.dy;

  dxl = f_x-floor(f_x);
  dyl = f_y-floor(f_y);
  dxr = ceil(f_x)-f_x;
  dyr = ceil(f_y)-f_y;

  xwrf = (int)(floor(f_x));
  ywrf = (int)(floor(f_y));

  float4 ip;

  ip = interpolate(dxl, dyl, dxr, dyr, xwrf, ywrf, b_wrf_flux);

  write_f4(pos.x, pos.y, pos.z, ip, b_sys_flux);
}

