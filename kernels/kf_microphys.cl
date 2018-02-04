// template
// (*cRhs).s0 += dsig;
// (*cRhs).s1 += drho;
// (*cRhs).s2 += drho_v;
// (*cRhs).s3 += drho_c;
// (*cRhs).s4 += drho_r;
// (*cRhs).s5 += dn_d;
// (*cRhs).s6 += dn_c;
// (*cRhs).s7 += dn_r;
// (*cRhs_ice).s0 += drho_i;
// (*cRhs_ice).s1 += drho_s;
// (*cRhs_ice).s2 += dn_i;
// (*cRhs_ice).s3 += dn_s;

typedef struct {
  float nu;    //..Breiteparameter der Verteil.
  float mu;    //..Exp.-parameter der Verteil.
  float x_max; //..maximale Teilchenmasse
  float x_min; //..minimale Teilchenmasse
  float a_geo; //..Koeff. Geometrie
  float b_geo; //..Koeff. Geometrie = 1/3
  float a_vel; //..Koeff. Fallgesetz
  float b_vel; //..Koeff. Fallgesetz
  float a_ven; //..Koeff. Ventilationsparam.
  float b_ven; //..Koeff. Ventilationsparam.
  float cap;   //..Koeff. Kapazitaet
  float k_au;  //..Koeff. Autokonversion
  float k_sc;  //..Koeff. Selfcollection
  float k_c;   //..Long-Kernel
  float k_1;   //..Parameter fuer Phi-Fkt.
  float k_2;   //..Parameter fuer Phi-Fkt.
} particle;

__constant float cloudRelaxation = 0.1f;
__constant particle pt_cloud = {
  .nu    = 0.333333f,
  .mu    = 0.666666f,
  .x_max = 2.60e-10f,
  .x_min = 4.20e-15f,
  .a_geo = 1.24e-01f,
  .b_geo = 0.333333f,
  .a_vel = 3.75e+05f,
  .b_vel = 0.666667f,
  .a_ven = 0.780000f,
  .b_ven = 0.308000f,
  .cap   = 2.0f,
  .k_au  = 1.99047e+19f,
  .k_sc  = 2.05131e+10f,
  .k_c   = 4.44e+9f,
  .k_1   = 4.00e+2f,
  .k_2   = 0.7f
};
__constant particle pt_rain = {
  .nu    = -0.666666f,
  .mu    =  0.333333f,
  .x_max = 3.00e-06f,
  .x_min = 2.60e-10f,
  .a_geo = 1.24e-01f,
  .b_geo = 0.333333f,
  .a_vel = 1.59e+02f,
  .b_vel = 0.266667f,
  .a_ven = 0.780000f,
  .b_ven = 0.308000f,
  .cap   = 2.0f,
  .k_au  = 0.0f,
  .k_sc  = 0.0f,
  .k_c   = 0.0f,
  .k_1   = 0.0f,
  .k_2   = 2.0f
};
__constant particle pt_ice = {
  // HK87 'hex plates'
  .nu    = -0.333333f,
  .mu    =  0.333333f,
  .x_max = 1.00e-07f,
  .x_min = 1.00e-12f,
  .a_geo = 2.17e-01f,
  .b_geo = 0.302115f,
  .a_vel = 4.19e+01f,
  .b_vel = 0.260000f,
  .a_ven = 0.860000f,
  .b_ven = 0.280000f,
  .cap   = 3.141593f,
  .k_au  = 0.0f,
  .k_sc  = 0.0f,
  .k_c   = 0.0f,
  .k_1   = 0.0f,
  .k_2   = 2.0f
};
__constant particle pt_snow = {
  //  nach Locatelli und Hobbs
  .nu    = 0.500000f,
  .mu    = 0.50000f,
  .x_max = 2.00e-06f,
  .x_min = 1.73e-09f,
  .a_geo = 8.16e-00f,
  .b_geo = 0.526316f,
  .a_vel = 2.77e+01f,
  .b_vel = 0.215790f,
  .a_ven = 0.780000f,
  .b_ven = 0.308000f,
  .cap   = 2.0f,
  .k_au  = 0.0f,
  .k_sc  = 0.0f,
  .k_c   = 0.0f,
  .k_1   = 0.0f,
  .k_2   = 2.0f
};

void bulk(parameters *par, state *st, float8 *cRhs) {

  float forcing = st->rho_v-(st->sv/st->T/par->rv);
  float limiter = st->rho_l;
  float drho_l  = cloudRelaxation*(forcing-limiter+sqrt(forcing*forcing+limiter*limiter));
  float dsig    = st->lv/st->T*drho_l+(st->lnT*(par->cpl-par->cpv)+st->lnP*par->rv)*drho_l;

  //using rho_c as rho_l;
  (*cRhs).s0 += dsig;
  (*cRhs).s2 -= drho_l;
  (*cRhs).s3 += drho_l;
}

void nucleation(parameters *par, position *pos, state *st, state *st_zr, float8 *cRhs) {

  // Saturation gradient, Nucleation parameterisation according to Seifert (16), (17)
  float dSdzw = (st_zr->sat-st->sat)/par->dz*pos->wrf;

  float kccn = 0.462f; // how much is water attracted by dirt
  float nuc  = 0.0f;   // nucleation keim

  // T>233, sonst eis
  if (st->sat>0.0f && dSdzw>0.0f && st->T>233.15f && pos->wrf>0.0f) {
    // nuc=st->n_d*kccn/pow(st->sat,kccn-1.0f)*dSdzw;  // nuc = rhs for n_c but not limited yet
    nuc=st->n_d*kccn/st->sat*dSdzw;  // nuc = rhs for n_c but not limited yet

    // if (pos->x == 9)
    // printf("%d %d\tpv:%2.3e\tsv:%2.3e\tsat:%2.3e\tw:%2.2e\tdsdzw:%2.2e\tnuc:%2.2e\n", pos->z, pos->x, st->pv, st->sv, st->sat, pos->wrf, dSdzw,nuc);
  }

  float forcing = nuc;
  float limiter = st->n_d;
  float dn_c    = cloudRelaxation*(forcing+limiter-sqrt(forcing*forcing+limiter*limiter)); // rhs for n_c

  // change in densities due to nucleation, 1e-12 is smallest (nuemerically usable) drop mass, activation with smallest possible mass
  float drho_c  = dn_c*1.0e-12f;
  // rhs for rho*sig
  float dsig =st->lv/st->T*drho_c+(st->lnT*(par->cpl-par->cpv)+st->lnP*par->rv)*drho_c;

  // if (pos->x == (int)(par->sx*0.5f) && pos->y == (int)(par->sy*0.5f)) printf("%d\tS_%2.2f\tdS_%2.2f\tnd_%2.2f\n", pos->z, S, dSdzw, n_d);

  (*cRhs).s0 += dsig;
  (*cRhs).s2 -= drho_c;
  (*cRhs).s3 += drho_c;
  (*cRhs).s5 -= dn_c;
  (*cRhs).s6 += dn_c;
}

void condensation(parameters *par, position *pos, state *st, float8 *cRhs) {
  // wasserdampf condensiert auf schon existierende tropfen
  float forcing = st->rho_v-(st->sv/st->T/par->rv);
  float limiter = st->rho_c;
  float drho_c  = cloudRelaxation*(forcing-limiter+sqrt(forcing*forcing+limiter*limiter));
  float dsig    = st->lv/st->T*drho_c+(st->lnT*(par->cpl-par->cpv)+st->lnP*par->rv)*drho_c;

  float dn_c = 0.0f;
  float dn_r = 0.0f;
  // float dn_i = 0.0f;
  // float dn_s = 0.0f;

  //Seifert (94) bissle huebscher -> zu kleine und zu grosse tropfen loeschen und erstellen
  forcing = 1.0e-1f*min(0.0f, (st->rho_c/pt_cloud.x_min-st->n_c));
  limiter = st->n_c;
  dn_c   += cloudRelaxation*(forcing-limiter+sqrt(forcing*forcing+limiter*limiter));

  forcing = 1.0e-1f*max(0.0f, (st->rho_c/pt_cloud.x_max-st->n_c));
  limiter = st->n_c;
  dn_c   += cloudRelaxation*(forcing+limiter-sqrt(forcing*forcing+limiter*limiter));

  // limiter wenn keine dichte, dann keine tropfen
  dn_r =0.0f;
  dn_r+=0.01f*min(0.0f, (st->rho_r/pt_rain.x_min-st->n_r));
  dn_r+=0.01f*max(0.0f, (st->rho_r/pt_rain.x_max-st->n_r));

  (*cRhs).s0 += dsig;
  (*cRhs).s2 -= drho_c;
  (*cRhs).s3 += drho_c;
  (*cRhs).s6 += dn_c;
  (*cRhs).s7 += dn_r;
}

void autoconversion_and_selfcollection_cloud_rain(parameters *par, position *pos, state *st, float8 *cRhs) {
  // auto: wolkentropfen anwachsen bis regentropfen, creates rain
  // self: zusammenstossende tropfen, werden groesser, nc wird vernichtet

  // SB 2005
  float au    = 0.0f;;     // autoconversion
  float sc;                // self collection
  float k_rr  = 7.12f;
  float x_i = 1.0f;      // mittlere masse in SI
  float tau;               // dimensionless internal time scale
  float phi;               // similarity function
  float forcing, limiter;
  float drho_r  = 0.0f;
  float dn_r_au = 0.0f;
  float dn_r_sc = 0.0f;
  float dn_c_sc = 0.0f;

  if (st->rho_c>1.0e-9f) {
    // avoid divisions by zero, use 1.0e-33 instead!
    x_i = min(max(st->rho_c/(st->n_c+1.0e-33f), pt_cloud.x_min), pt_cloud.x_max);
    au = pt_cloud.k_au*st->rho_c*st->rho_c*x_i*x_i*par->rr/(st->rho+1.0e-33f);
  }
  if (st->rho_c>1.0e-6f) {
    tau = min(max(1.0f-st->rho_c/(st->rho_c+st->rho_r+1.0e-33f), 1.0e-33f), 0.9f);
    phi = 400.0f*pow(tau, 0.7f)*pow(1.0f-pow(tau, 0.7f), 3.0f);
    au  = au*(1.0f+phi/pow(1.0f-tau,2.0f));
  }

  forcing = au/x_i;
  limiter = st->n_c; // need ensured positve values
  dn_r_au = cloudRelaxation*(forcing+limiter-sqrt(forcing*forcing+limiter*limiter));

  forcing = au;
  limiter = st->rho_r;
  drho_r   = cloudRelaxation*(forcing+limiter-sqrt(forcing*forcing+limiter*limiter));

  // self collection cloud
  sc      = pt_cloud.k_sc*st->rho_c*st->rho_c*par->rr/(st->rho+1.0e-33f);
  forcing = -sc;
  limiter = st->n_c;
  dn_c_sc = cloudRelaxation*(forcing-limiter+sqrt(forcing*forcing+limiter*limiter));

  // self collection rain
  x_i   = min(max(st->rho_r/(st->n_r+1.0e-33f), pt_rain.x_min), pt_rain.x_max);
  sc      = k_rr*st->n_r*st->rho_r*sqrt(par->rr/(st->rho+1.0e-33f));
  forcing = -sc;
  limiter = st->n_r;
  dn_r_sc = cloudRelaxation*(forcing-limiter+sqrt(forcing*forcing+limiter*limiter));

  (*cRhs).s3 -= drho_r;
  (*cRhs).s4 += drho_r;
  (*cRhs).s6 += dn_c_sc - dn_r_au;
  (*cRhs).s7 += dn_r_sc + dn_r_au;
}

void rain_sedimentation_and_subsidence(parameters *par, position *pos, state *st, state *st_zr, state *st_zl, float8 *cRhs) {
  float a_qc   = 8420.0f;
  float b_qc   = 0.8f;
  float rho_qc = 1000.0f;
  float n_qc   = 8.0f;

  float x_i_r, dm, mue, dr, vn, vq;
  float flux, dwrainflx, subsidence;
  float h;
  float theta, theta_zr, theta_zl;
  float drho_r, drho_v, dsig, dn_r;
  drho_r = 0.0f;
  drho_v = 0.0f;
  dsig   = 0.0f;
  dn_r   = 0.0f;

  float alf,bet,gam;
  alf=9.65e+0f;
  bet=1.03e+1f;
  gam=6.00e+2f;

  //reinfallen regen
  if (pos->s_zr==1 && (st_zr->rho_r>1.0e-10f)) {
    x_i_r = st_zr->rho_r/(st_zr->n_r+1.0e-33f);                                  //..mittlere Masse in SI
    x_i_r = min(max(x_i_r,pt_rain.x_min), pt_rain.x_max);
    dm = pow(6.0f/(1000.0f*(float)(M_PI_F))*x_i_r, 1.0f/3.0f);
    mue= (pt_rain.nu+1.0f)/pt_rain.b_geo - 1.0f;
    dr = pow(pow(dm,3.0f)/((mue+3.0f)*(mue+2.0f)*(mue+1.0f)),1.0f/3.0f);
    vn = alf-bet/pow(1.0f+gam*dr,mue+1.0f);
    vq = alf-bet/pow(1.0f+gam*dr,mue+4.0f);
    vn = vn*pow(par->rr/(st_zr->rho+1.0e-33f),0.5f);
    vq = vq*pow(par->rr/(st_zr->rho+1.0e-33f),0.5f);
    vn = max(vn,1.0e-1f);
    vq = max(vq,1.0e-1f);
    vn = min(vn,2.0e+1f);
    vq = min(vq,2.0e+1f);

    flux    = min(st_zr->rho_r*vq/par->dz, st_zr->rho_r/2.0f/par->dt);
    drho_r += flux;
    dsig   += st->lnT*par->cpl*flux;

    flux  = min(st_zr->n_r*vn/par->dz, st_zr->n_r/2.0f/par->dt);
    dn_r += flux;
  }

  //rausfallen regen
  if (st->rho_r>1.0e-10f) {
    x_i_r = st->rho_r/(st->n_r+1.0e-33f);                                      //..mittlere Masse in SI
    x_i_r = min(max(x_i_r, pt_rain.x_min), pt_rain.x_max);
    dm = pow(6.0f/(1000.0f*(float)(M_PI_F))*x_i_r,1.0f/3.0f);
    mue= (pt_rain.nu+1.0f)/pt_rain.b_geo - 1.0f;
    dr = pow(pow(dm, 3.0f)/((mue+3.0f)*(mue+2.0f)*(mue+1.0f)), 1.0f/3.0f);
    vn = alf-bet/pow(1.0f+gam*dr,mue+1.0f);
    vq = alf-bet/pow(1.0f+gam*dr,mue+4.0f);
    vn = vn*pow(par->rr/(st->rho+1.0e-33f), 0.5f);
    vq = vq*pow(par->rr/(st->rho+1.0e-33f), 0.5f);
    vn = max(vn,1.0e-1f);
    vq = max(vq,1.0e-1f);
    vn = min(vn,2.0e+1f);
    vq = min(vq,2.0e+1f);

    dwrainflx = min(st->rho_r*vq/par->dz, st->rho_r/2.0f/par->dt);
    drho_r -= dwrainflx;
    dsig   -= st->lnT*par->cpl*dwrainflx;

    flux  = min(st->n_r*vn/par->dz, st->n_r/2.0f/par->dt);
    dn_r -= flux;
  }

  // subsidenz - anwachsen der grenzschicht verhindern durch neuen absinkende-kaltluft term
  if (pos->s_zr==1 && pos->s_zl==1) {
    h=((float)(pos->z)+0.5f)*par->dz;
    // grenzschicht given at 825m dycoms
    // subsidence = 0.0f;
    // if (h>=825.0f) subsidence = -3.75e-6f*h;

    // grenzschicht isdac, test the sign
    subsidence = -0.412510e-2f;
    if (h<=825.0f) subsidence = -5e-6f*h;

    theta    = exp((   st->sig+   st->rml*log(par->pr))/   st->cpml);
    theta_zr = exp((st_zr->sig+st_zr->rml*log(par->pr))/st_zr->cpml);
    theta_zl = exp((st_zl->sig+st_zl->rml*log(par->pr))/st_zl->cpml);

    flux  = -subsidence*(theta_zr-theta_zl)/2.0f/par->dz;
    dsig += st->cpml/theta*flux;

    flux   = -subsidence*(st_zr->rho_v/st_zr->rho-st_zl->rho_v/st_zl->rho)/2.0f/par->dz;
    drho_v = flux*st->rho;
    dsig  += (st->lnT*(par->cpv-par->cpd)-st->lnP*(par->rv-par->rd))*flux*st->rho;
  }

  (*cRhs).s0 += dsig;
  (*cRhs).s2 += drho_v;
  (*cRhs).s4 += drho_r;
  (*cRhs).s7 += dn_r;
}

void rain_evaporation(parameters *par, position *pos, state *st, float8 *cRhs) {
  float g_d,x_r,d_r,v_r,N_re,f_v,f_q,eva,eva_q,eva_n;
  float a_q,b_q,c_r,D_vtp;
  float K_T  = 2.500e-2f;
  float nu_l = 1.460e-5f;    //..Kinem. Visc. von Luft
  float N_sc = 0.710f;        //..Schmidt-Zahl (PK, S.541)
  float n_f  = 0.333f;        //..Exponent von N_sc im Vent-koeff. (PK, S.541)
  float m_f  = 0.500f;        //..Exponent von N_re im Vent-koeff. (PK, S.541)
  float q_krit = 1.000e-9f;
  float drho_r, drho_v, dsig, dn_r;
  float precevap, precevapcooling;
  float forcing, limiter;

  a_q = 0.429250753972163723f;
  b_q = 0.180944125457962318f;
  c_r = 1.0f/pt_rain.cap;
  drho_r = 0.0f;
  drho_v = 0.0f;
  dsig   = 0.0f;
  dn_r   = 0.0f;
  eva_q  = 0.0f;
  eva_n  = 0.0f;
  precevap        = 0.0f;
  precevapcooling = 0.0f;

  if (st->sat<0.0f && st->rho_r>0.0f) {
    D_vtp = pow(8.7602e-5f*st->T,1.81f)/st->P;
    g_d   = 2.0f*(float)(M_PI_F)/(st->lv*st->lv/(K_T*par->rd*st->T*st->T)+par->rd*st->T/(D_vtp*exp(st->sv)));
    x_r   = min(max(st->rho_r/(st->n_r+1e-10f),pt_rain.x_min),pt_rain.x_max);
    d_r   = pt_rain.a_geo * pow(x_r,pt_rain.b_geo);
    v_r   = pt_rain.a_vel * pow(x_r,pt_rain.b_vel) * pow(par->rr/(st->rho+1e-33f),0.5f);

    N_re  = v_r*d_r/nu_l;
    f_v   = pow(N_sc,n_f)*pow(N_re,m_f);

    f_q   = a_q+b_q*f_v;

    eva   = g_d*max(st->n_r,1e2f)*f_q*c_r*d_r*st->sat;
    eva_q = min(-5e-12f,eva);
    precevap=eva;
    eva_n = eva_q/x_r;
  }
  forcing = eva_q;
  limiter = st->rho_r;
  drho_r  = cloudRelaxation*(forcing-limiter+sqrt(forcing*forcing+limiter*limiter));
  drho_v  = -drho_r;
  forcing = eva_n;
  limiter = st->n_r;
  dn_r    = cloudRelaxation*(forcing-limiter+sqrt(forcing*forcing+limiter*limiter));
  dsig    = -st->lv/st->T*drho_v-(st->lnT*(par->cpl-par->cpv)+st->lnP*par->rv)*drho_v;

  precevapcooling=dsig;

  (*cRhs).s0 += dsig;
  (*cRhs).s2 += drho_v;
  (*cRhs).s4 -= drho_v;
  (*cRhs).s5 -= dn_r;
  (*cRhs).s7 += dn_r;
}

void ice_freeze_cloud(parameters *par, position *pos, state *st, float8 *cRhs, float4 *cRhs_ice) {
  float fr_q = 0.0f;
  float fr_n = 0.0f;
  float drho_i = 0.0f;
  float dn_i = 0.0f;
  float dsig = 0.0f;
  float forcing, limiter;
  float x_i = 1.0f; // mittlere Teilchenmasse
  float j_hom = 0.0f; // temperature function for freezing
  float j_het = 0.0f;
  float j_tot = 0.0f;

  float a_het = 0.2f;
  float b_het = 0.65f;
  float facg = 2.173f; // factor from gamma distribution

  float T_c = st->T - 273.15f; // temperature in celsius
  if (T_c < 273.15f) {
    if (T_c < -50.0f) {
      // homogeneous freezing
      // below -50 C everything freezes homogene
      fr_q = st->rho_c;
      fr_n = st->n_c;
    } else {
      // homgeneous freezing of cloud droplets after Cotton and Field (2002) equ. 12 in SI [1/m3 s]
      x_i = min(max(st->rho_c/(st->n_c), pt_cloud.x_min), pt_cloud.x_max);
      if (T_c <= -30.0f) {
        j_hom = 1.0e6f*exp(-243.4f - 14.75f*T_c - 0.307f*pow(T_c, 2.0f) - 0.00287f*pow(T_c, 3.0f) - 0.0000102f*pow(T_c,4.0f));
      } else {
        j_hom = 1.0e6f*exp(-7.63f - 2.996f*(T_c + 30.0f));
      }
      // heterogeneous
      // Seifert 2005 (44)
      j_het = max(a_het*(exp(b_het*(273.15f-st->T) - 1.0f)), 0.0f);
      j_tot = (j_hom/1000.0f)+j_het; // different units for the two j
      // j_tot = 0.0f;
      fr_q = j_tot*st->rho_c*x_i*facg; // rate of change of mass densities (liquid water content)
      fr_n = j_tot*st->rho_c;            // rate of change of number distribution
      // if(fabs(fr_q) > 1e-12f)printf("%g %g %e %g %g\n", j_tot, st->rho_c, st->rho_i, fr_q, fr_n);
    }
    // fr_q = min(fr_q, st->rho_c);
    // fr_n = min(fr_n, st->n_c);

    forcing = fr_q;
    limiter = max(st->rho_c, 0.0f);
    drho_i  = cloudRelaxation*(forcing+limiter-sqrt(forcing*forcing+limiter*limiter));
    // if (fabs(drho_i) >= 0) printf("%g %g\n", drho_i, st->rho_c);
    // if (pos->x == 1 && pos->y == 1 && pos->z == 35) printf("%g %g\n", drho_i, st->rho_c);

    forcing = fr_n;
    limiter = st->n_c;
    dn_i    = cloudRelaxation*(forcing+limiter-sqrt(forcing*forcing+limiter*limiter));


    dsig = par->lvf/st->T*drho_i+(st->lnT*(par->cpi-par->cpl))*drho_i;

    (*cRhs).s0 += dsig;
    (*cRhs).s3 -= drho_i;
    (*cRhs).s6 -= dn_i;
    (*cRhs_ice).s0 += drho_i;
    (*cRhs_ice).s2 += dn_i;
  }
}

void ice_melt(parameters *par, position *pos, state *st, float8 *cRhs, float4 *cRhs_ice) {
  float x_i;
  float d_i, v_i;
  float N_Re, D_T;
  float fv_q, fh_q;
  float melt, melt_h, melt_v, melt_q, melt_n;
  float forcing, limiter;
  float dsig = 0.0f;
  float drho_c = 0.0f;
  float dn_c = 0.0f;

  float nu_l = 1.460e-5f;    //..Kinem. Visc. von Luft
  float K_T  = 2.500e-2f;
  float a_melt_q = 0.67625f;
  float b_melt_q = 0.293159f;
  float N_sc = 0.710f;        //..Schmidt-Zahl (PK, S.541)
  float n_f  = 0.333f;        //..Exponent von N_sc im Vent-koeff. (PK, S.541)
  float m_f  = 0.500f;        //..Exponent von N_re im Vent-koeff. (PK, S.541)
  float D_V = 3e-5f; // diffusivity of water vapor

  if (st->T > 273.15f && st->rho_i > 0.0f && st->n_i > 0.0f) {
    x_i = min(max(st->rho_i/(st->n_i), pt_ice.x_min), pt_ice.x_max); // mass of ice [kg]
    d_i = pt_ice.a_geo*pow(x_i, pt_ice.b_geo);                     // diameter-mass relation [m]
    v_i = pt_ice.a_vel*pow(x_i, pt_ice.b_vel)*sqrt(par->rr/st->rho); // velocity-mass relation [m s^-1]
    N_Re = (v_i*d_i)/nu_l;                                           // Reynoldsnumber
    D_T = K_T/(par->cpd*par->rr);                                      // diffusivity of heat
    // D_T  = 2.5e-2/(1005.7*rho); // Stefan
    // Seifert (73)
    fv_q = a_melt_q + b_melt_q*pow(N_sc, n_f)*pow(N_Re, m_f);          // averaged ventilation coeff vapour
    fh_q = D_T/D_V*fv_q;                                            // averaged vent coeff heat

    melt = 2.0f*M_PI_F/par->lvf*d_i*st->n_i;
    melt_h = melt * K_T*(st->T-273.15f);
    melt_v = melt * D_V*st->lv/par->rv*(st->sv/st->T-par->svr/par->tr);
    melt_q = (melt_h*fh_q + melt_v*fv_q);
    melt_n = min(max((melt_q - st->rho_i)/x_i + st->n_i, 0.0f), st->n_i);

    if (st->T>283.15f) {
      melt_q = st->rho_i;
      melt_n = st->n_i;
    }

    forcing = melt_q;
    limiter = st->rho_i;
    drho_c  = cloudRelaxation*(forcing+limiter-sqrt(forcing*forcing+limiter*limiter));

    forcing = melt_n;
    limiter = st->n_i;
    dn_c    = cloudRelaxation*(forcing+limiter-sqrt(forcing*forcing+limiter*limiter));

    dsig = -par->lvf/st->T*drho_c-(st->lnT*(par->cpi-par->cpl))*drho_c;

    (*cRhs).s0 += dsig;
    (*cRhs).s3 += drho_c;
    (*cRhs).s6 += dn_c;
    (*cRhs_ice).s0 -= drho_c;
    (*cRhs_ice).s2 -= dn_c;
  }
}


void selfcollection_ice_to_snow(parameters *par, position *pos, state *st, float8 *cRhs, float4 *cRhs_ice) {
  float x_i,D_i,e_coll,v_i,x_s,D_s,v_s;
  float ice_s_vel=0.25f; // dispersion der fallges.
  float snow_s_vel=0.25f; // dispersion der fallges.

  float delta_n=3.56212f;
  float theta_n=0.100138f;
  float drho_s=0.0f;
  float drho_i=0.0f;
  float dn_s=0.0f;
  float dn_i=0.0f;
  float sn=0.0f;
  float sq=0.0f;

  float forcing, limiter;

  if (st->rho_i>1e-9f) {
    //..Selfcollection Ice
    x_i = min(max(st->rho_i/(st->n_i+1e-33f),pt_ice.x_min),pt_ice.x_max);
    D_i = pow(pt_ice.a_geo*x_i,pt_ice.b_geo);                     //..mittlerer Durchmesser

    if (st->T>273.15f)
      e_coll = 1.0f;
    else {
      //.. Temperaturabhaengige Efficiency nach Cotton et al. (1986)
      //   (siehe auch Straka, 1989; S. 53)
      e_coll = min(pow(10.0f,0.035f*(st->T-273.15f)-0.7f),0.2e0f);
    }

    v_i = pow(pt_ice.a_vel*x_i,pt_ice.b_vel)*par->rr/(st->rho+1e-33f); //..mittlere Sedimentationsgeschw.

    sn = M_PI_F*0.25e0f*e_coll*st->n_i*  st->n_i*D_i*D_i*sqrt(1.0f*v_i*v_i+2.0f*ice_s_vel*ice_s_vel);
    sq = M_PI_F*0.25e0f*e_coll*st->n_i*st->rho_i*D_i*D_i*sqrt(1.0f*v_i*v_i+2.0f*ice_s_vel*ice_s_vel);

  }
  forcing  = -sq;
  limiter  = max(st->rho_i, 0.0f);
  drho_i   = cloudRelaxation*(forcing-limiter+sqrt(forcing*forcing+limiter*limiter));
  // if (pos->x == 1 && pos->y == 1 && pos->z == 35 && drho_i != 0.0f) printf("%g\t%g\t %g\n", sq, drho_i, st->T);
  // if (drho_i != 0.0f) printf("%g\t%g\t %g\n", sq, drho_i, st->T);

  forcing  = -sn;
  limiter  = st->n_i;
  dn_i     = cloudRelaxation*(forcing-limiter+sqrt(forcing*forcing+limiter*limiter));

  // selfcollection snow
  sn=0.0f;
  if (st->rho_s>0.0f) {
    if (st->T>273.15f) e_coll = 1.0f;
    else e_coll = max(0.1e0f,min(exp(0.09f*(st->T-273.15f)),1.0e0f));

    x_s = min(max(st->rho_s/(st->n_s+1e-33f),pt_snow.x_min),pt_snow.x_max);    //..mittlere Masse in SI
    D_s = pow(pt_snow.a_geo*x_s,pt_snow.b_geo);                     //..mittlerer Durchmesser
    v_s = pow(pt_snow.a_vel*x_s,pt_snow.b_vel)*par->rr/(st->rho+1e-33f);    //..mittlere Sedimentationsgeschw.

    sn = M_PI_F*0.125e0f*e_coll*st->n_s*st->n_s*delta_n*D_s*D_s*sqrt(theta_n*v_s*v_s+2.0f*snow_s_vel*snow_s_vel);
  }

  forcing  = -sn;

  limiter  = st->n_s;
  dn_s     = cloudRelaxation*(forcing-limiter+sqrt(forcing*forcing+limiter*limiter));

  (*cRhs_ice).s0 += drho_i;
  (*cRhs_ice).s1 -= drho_i;
  (*cRhs_ice).s2 += dn_i;
  (*cRhs_ice).s3 += dn_s-dn_i;
}

void snow_sedimentation(parameters *par, position *pos, state *st, state *st_zr, float8 *cRhs, float4 *cRhs_ice) {

  float xr,xs,n0,lam,Dm,Dr,mue,vn,vq,qr;

  float alf=9.65e+0f;
  float bet=1.03e+1f;
  float gam=6.00e+2f;
  float drho_s=0.0f;
  float dsig=0.0f;
  float dn_s=0.0f;
  float Flux=0.0f;

  // reinfallen
  if (pos->s_zr==1 && (st_zr->rho_s>1.0e-10f)) {
    xs = st_zr->rho_s/st_zr->n_s;
    xs = min(max(xs,pt_snow.x_min),pt_snow.x_max);
    lam= pow(0.083333f*xs,pt_snow.b_vel)*sqrt(par->rr/(st->rho+1e-33f));
    vn = 42.7154f*lam;
    vq = 54.1323f*lam;
    vn = max(vn,0.1e+0f);
    vq = max(vq,0.1e+0f);
    vn = min(vn,3.0e+0f);
    vq = min(vq,3.0e+0f);

    Flux   = st_zr->rho_s*vq/par->dz;
    drho_s += Flux;
    dsig   += st->lnT*par->cpi*Flux;

    Flux   = st_zr->n_s*vn/par->dz;
    dn_s   += Flux;
  }

  if (st->rho_s>1.0e-10f) {
    xs = st->rho_s/st->n_s;
    xs = min(max(xs,pt_snow.x_min),pt_snow.x_max);
    lam= pow(0.083333f*xs,pt_snow.b_vel)*sqrt(par->rr/(st->rho+1e-33f));
    vn = 42.7154f*lam;
    vq = 54.1323f*lam;
    vn = max(vn,0.1e+0f);
    vq = max(vq,0.1e+0f);
    vn = min(vn,3.0e+0f);
    vq = min(vq,3.0e+0f);

    Flux   = st->rho_s*vq/par->dz;
    drho_s -= Flux;
    dsig   -= st->lnT*par->cpi*Flux;

    Flux   = st->n_s*vn/par->dz;
    dn_s   -= Flux;
  }

  (*cRhs).s0 += dsig;
  (*cRhs_ice).s1 += drho_s;
  (*cRhs_ice).s3 += dn_s;

}

void ice_deposition(parameters *par, position *pos, state *st, float8 *cRhs, float4 *cRhs_ice) {
  // growth of ice particles by water vapour deposition
  float K_T  = 2.500e-2f;
  float nu_l = 1.460e-5f;     //..Kinem. Visc. von Luft
  float N_sc = 0.710f;        //..Schmidt-Zahl (PK, S.541)
  float n_f  = 0.333f;        //..Exponent von N_sc im Vent-koeff. (PK, S.541)
  float m_f  = 0.500f;        //..Exponent von N_re im Vent-koeff. (PK, S.541)
  float D_V = 3e-5f;          // diffusivity of water vapor

  float fr_q = 0.0f;
  float D_i, v_i, G_iv, N_Re, F_vi, x_i;
  float limiter, forcing;
  float drho_i = 0.0f;
  float dsig = 0.0f;

  float S_i=((st->rho_v*par->rv*st->T/st->svi)-1.0f);                //supersaturation wrt. ice
  if (st->rho_i>0.0f && st->n_i > 0.0f && st->T < par->tr) {
    x_i=min(max(st->rho_i/st->n_i, pt_ice.x_min),pt_ice.x_max);  //mass of ice [kg]
    D_i=pow(pt_ice.a_geo*x_i, pt_ice.b_geo);                   //diameter-mass relation [m]
    v_i=pow(pt_ice.a_vel*x_i, pt_ice.b_vel)*sqrt(par->rr/st->rho);     //velocity-mass relation [m s^-1]
    G_iv=1.0f/(((par->rv*st->T)/(st->svi*D_V))+(par->lrs0/(K_T*st->T))*((par->lrs0/(par->rv*st->T))-1.0f));
    N_Re=(v_i*D_i)/nu_l;                                       //Reynoldsnumber
    F_vi=pt_ice.a_ven+pt_ice.b_ven*pow(N_sc, n_f)*pow(N_Re, m_f);
    fr_q=4.0f*G_iv*D_i*F_vi*S_i*st->n_i;                        //mass density of ice
  }


  forcing = fr_q;
  limiter = max(st->rho_v, 0.0f);
  drho_i  = cloudRelaxation*(forcing+limiter-sqrt(forcing*forcing+limiter*limiter));

  dsig   = (par->lrs0+(par->cpv-par->cpi)*st->T)/st->T*drho_i+(st->lnT*(par->cpi-par->cpv)+st->lnP*par->rv)*drho_i;

  (*cRhs).s0 += dsig;
  (*cRhs).s2 -= drho_i;
  (*cRhs_ice).s0 += drho_i;

}



__kernel void kf_microphys_kernel_main(__private parameters par,
                                       __private uint phys,
                                       __read_only image3d_t bf_scalars_vc_a_0,
                                       __read_only image3d_t bf_scalars_vc_a_1,
                                       __read_only image3d_t bf_scalars_vc_a_2,
                                       __read_only image3d_t bf_momenta_fc_b,
                                       __write_only image3d_t bRhs_mp_vc_0,
                                       __write_only image3d_t bRhs_mp_vc_1,
                                       __write_only image3d_t bRhs_mp_vc_2)
{
  position pos = get_pos_bc(&par);

  float8 c, czr, czl;
  float4 c_ice, czr_ice, czl_ice;
  c   = read_f8(pos.x , pos.y , pos.z , bf_scalars_vc_a_0, bf_scalars_vc_a_1);
  czr = read_f8(pos.x , pos.y , pos.zr, bf_scalars_vc_a_0, bf_scalars_vc_a_1);
  czl = read_f8(pos.x , pos.y , pos.zl, bf_scalars_vc_a_0, bf_scalars_vc_a_1);

  //velocity at right/upper surface
  pos.wrf = read_f4(pos.x, pos.y, pos.zr, bf_momenta_fc_b).s2*2.0f/(czr.s1+c.s1);
  //bc
  if (pos.s_zr==-1) pos.wrf=0.0f;

  state st, st_zr, st_zl;
  if (include_ice) {
    c_ice   = read_f4(pos.x , pos.y , pos.z , bf_scalars_vc_a_2);
    czr_ice = read_f4(pos.x , pos.y , pos.zr, bf_scalars_vc_a_2);
    czl_ice = read_f4(pos.x , pos.y , pos.zl, bf_scalars_vc_a_2);
    st    = init_state_with_ice(&par, &c,   &c_ice);
    st_zr = init_state_with_ice(&par, &czr, &czr_ice);
    st_zl = init_state_with_ice(&par, &czl, &czl_ice);
  } else {
    st    = init_state(&par, &c);
    st_zr = init_state(&par, &czr);
    st_zl = init_state(&par, &czl);
  }

  float8 cRhs = (float8)(0.0f);
  float4 cRhs_ice = (float4)(0.0f);

  // using bitwise and to enable/disable functions. sum up the integers of the ones to be enabled
  if (phys &   1) nucleation(&par, &pos, &st, &st_zr, &cRhs);
  if (phys &   2) condensation(&par, &pos, &st, &cRhs);
  if (phys &   4) autoconversion_and_selfcollection_cloud_rain(&par, &pos, &st, &cRhs);
  if (phys &   8) rain_sedimentation_and_subsidence(&par, &pos, &st, &st_zr, &st_zl, &cRhs);
  if (phys &  16) rain_evaporation(&par, &pos, &st, &cRhs);
  if (phys &  32) ice_freeze_cloud(&par, &pos, &st, &cRhs, &cRhs_ice);
  if (phys &  64) selfcollection_ice_to_snow(&par, &pos, &st, &cRhs, &cRhs_ice);
  if (phys & 128) ice_melt(&par, &pos, &st, &cRhs, &cRhs_ice);
  if (phys & 256) snow_sedimentation(&par, &pos, &st, &st_zr, &cRhs, &cRhs_ice);
  if (phys & 512) ice_deposition(&par, &pos, &st, &cRhs, &cRhs_ice);

  // bulk(par, c, &cRhs);
  // printf("%f %f %f %f\n", cRhs_ice.s0, cRhs_ice.s1, cRhs_ice.s2, cRhs_ice.s3);


  write_f8(pos.x, pos.y, pos.z, &cRhs, bRhs_mp_vc_0, bRhs_mp_vc_1);
  write_f4(pos.x, pos.y, pos.z, &cRhs_ice, bRhs_mp_vc_2);
}
