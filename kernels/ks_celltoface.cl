__kernel void ks_celltoface_kernel_main(__private parameters par,
                                        __read_only image3d_t bRhs_m_vc_uv,
                                        __read_only image3d_t bRhs_m_vc_w,
                                        __write_only image3d_t bsRhs_momenta_fc_s)
{
  position pos = get_pos_bc(par, get_global_id(0), get_global_id(1), get_global_id(2));

  // pos_ul + posxl_ur
  float u  = (read_imagef(bRhs_m_vc_uv, (int4)(pos.x, pos.y, pos.z, 0)).s0 + read_imagef(bRhs_m_vc_uv, (int4)(pos.xl, pos.y , pos.z , 0)).s1) * 0.5f;
  // pos_vl + posyl_vr
  float v  = (read_imagef(bRhs_m_vc_uv, (int4)(pos.x, pos.y, pos.z, 0)).s2 + read_imagef(bRhs_m_vc_uv, (int4)(pos.x , pos.yl, pos.z , 0)).s3) * 0.5f;
  // pos_wl + poszl_wr
  float w  = (read_imagef(bRhs_m_vc_w , (int4)(pos.x, pos.y, pos.z, 0)).s0 + read_imagef(bRhs_m_vc_w , (int4)(pos.x , pos.y , pos.zl, 0)).s1) * 0.5f;

  // fixed bc
  if (pos.z == 0) {
    w = 0.0f;
  }

  write_imagef(bsRhs_momenta_fc_s, (int4)(pos.x, pos.y, pos.z, 0), (float4)(u, v, w, 0.0f));
}
