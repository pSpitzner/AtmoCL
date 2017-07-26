__kernel void ke_uv_kernel_main(__private parameters par,
                                __private uint ref,
                                __private uint dim,
                                __read_only image3d_t b_source_scalars_0,
                                __read_only image3d_t b_source_scalars_1,
                                __read_only image3d_t b_source_scalars_2,
                                __read_only image3d_t b_source_momenta,
                                __write_only image3d_t b_target)
{
  position pos = get_pos_bc(par, get_global_id(0), get_global_id(1), get_global_id(2));

  float u = 0.0f;
  float v = 0.0f;

  if      (dim == 0) u = (read_f4(pos.x+ref, pos.y,     pos.z,     b_source_momenta)).s0/(read_f4(pos.x+ref, pos.y,     pos.z,     b_source_scalars_0)).s1;
  else if (dim == 1) u = (read_f4(pos.x    , pos.y+ref, pos.z,     b_source_momenta)).s0/(read_f4(pos.x    , pos.y+ref, pos.z,     b_source_scalars_0)).s1;
  else if (dim == 2) u = (read_f4(pos.x    , pos.y,     pos.z+ref, b_source_momenta)).s0/(read_f4(pos.x    , pos.y,     pos.z+ref, b_source_scalars_0)).s1;

  if      (dim == 0) v = (read_f4(pos.x+ref, pos.y,     pos.z,     b_source_momenta)).s1/(read_f4(pos.x+ref, pos.y,     pos.z,     b_source_scalars_0)).s1;
  else if (dim == 1) v = (read_f4(pos.x    , pos.y+ref, pos.z,     b_source_momenta)).s1/(read_f4(pos.x    , pos.y+ref, pos.z,     b_source_scalars_0)).s1;
  else if (dim == 2) v = (read_f4(pos.x    , pos.y,     pos.z+ref, b_source_momenta)).s1/(read_f4(pos.x    , pos.y,     pos.z+ref, b_source_scalars_0)).s1;

  // u = (float)(pos.x - par.sx/2)/(float)(par.sx);
  // v = (float)(pos.y - par.sy/2)/(float)(par.sy);
  // u = 0.5;
  // v = 1;

  float umax = 10.0f; // [m/s]
  float vmax = 10.0f;
  u= (u < umax ? u/umax : 1.0f);
  v= (v < umax ? v/vmax : 1.0f);

  float sat = sqrt(v*v+u*u)/1.414;
  // sat = 1.0f;
  float hue = atan2pi(v,u);
  if (hue < 0) hue += 2;
  hue *= 180;

  float value = 1.0f;
  float chroma = value * sat;
  float hue1 = hue / 60.0f;
  float x = chroma * (1.0f-fabs((fmod(hue1,2)) - 1.0f));
  float r1, g1, b1;
  if (hue1 >= 0 && hue1 <= 1) {
    r1 = chroma;
    g1 = x;
    b1 = 0;
  } else if (hue1 >= 1 && hue1 <= 2) {
    r1 = x;
    g1 = chroma;
    b1 = 0;
  } else if (hue1 >= 2 && hue1 <= 3) {
    r1 = 0;
    g1 = chroma;
    b1 = x;
  } else if (hue1 >= 3 && hue1 <= 4) {
    r1 = 0;
    g1 = x;
    b1 = chroma;
  } else if (hue1 >= 4 && hue1 <= 5) {
    r1 = x;
    g1 = 0;
    b1 = chroma;
  } else if (hue1 >= 5 && hue1 <= 6) {
    r1 = chroma;
    g1 = 0;
    b1 = x;
  }

  float m = value - chroma;
  float r, g, b;
  r = r1+m;
  g = g1+m;
  b = b1+m;

  // if (pos.x == 0 && pos.y == 0) printf("%f %f %f %f\n", sat, r, g, b);

  float4 rgba = (float4)(fabs(r)*255.0f,fabs(g)*255.0f,fabs(b)*255.0f,0.0f);


  write_f4(pos.x, pos.y, pos.z, rgba, b_target);
}
