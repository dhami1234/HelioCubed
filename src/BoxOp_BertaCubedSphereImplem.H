// ====================================================================================
// Berta-style 4th-order FV pointwise-WENOZ implementation.
// Corrected for face indexing convention:
//   face index i = lower face of cell i = x_{i-1/2}
// so the face lies between cells i-1 and i.
//
// Fixes included here:
//   1. thetaCell is defined on the full cell box and defaults to 0 at boundaries.
//   2. thetaReconFace and thetaFluxFace are separated:
//        thetaReconFace: full WENOZ stencil flag, used to choose WENOZ vs fallback.
//        thetaFluxFace : adjacent-cell flag, used for flux point-to-average correction.
//   3. Default fallback is WENO3, not PLM.
// ====================================================================================

// ====================================================================================
// 1. Primitive conversion helpers.
// ====================================================================================

template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::sphPointToPrim(
                                          State& a_WPoint,
                                          const State& a_UPoint,
                                          uint a_block) const
{
  PR_TIMERS("BoxOp_EulerCubedSphere::sphPointToPrim()");

  // Pure pointwise conserved -> pointwise primitive conversion.
  // No deconvolution here.
  // No primitive-average correction here.
  a_WPoint = forall<T,NUMCOMPS,MEM>(f_sphToPrim,a_UPoint);
  const bool magneticTest =
      isConstantMagneticFieldTest(ParseInputs::get_init_condition_type());
  BoxData<T,1,MEM> deltaMe(a_WPoint.box());
  deltaMe.setVal(0.0);
  if (magneticTest)
  {
    testProblemDeltaMagneticEnergy(deltaMe,a_WPoint,a_block);
    addTestProblemBackgroundB(a_WPoint,a_block);
  }

  BoxData<T,1,MEM> kePoint = kineticEnergyW(a_WPoint,a_block,-1);
  BoxData<T,1,MEM> mePoint = magneticEnergyW(a_WPoint,a_block,-1);

  forallInPlace(
      [ ] PROTO_LAMBDA
      (
        Var<T,NUMCOMPS>& a_w,
        const Var<T,1,MEM>& a_ke,
        const Var<T,1,MEM>& a_me,
        const Var<T,1,MEM>& a_deltaMe,
        bool a_magneticTest,
        T a_gamma)
      {
        if (a_magneticTest)
        {
          a_w(iP) = c_TEST_P0 + (a_gamma - 1.0)*
              (a_w(iP) - a_ke(0) - a_deltaMe(0));
        }
        else
        {
          a_w(iP) = (a_w(iP) - a_ke(0) - a_me(0))*(a_gamma - 1.0);
        }
        if (a_w(iP) < 0.0) a_w(iP) = 0.0;
      },
      a_WPoint,
      kePoint,
      mePoint,
      deltaMe,
      magneticTest,
      m_gamma);
}


template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::sphAvgToPrim(
                                          State& a_WAvg,
                                          const State& a_UAvg,
                                          uint a_block) const
{
  PR_TIMERS("BoxOp_EulerCubedSphere::sphAvgToPrim()");

  // Direct V(<U>) conversion.
  // Used only for Berta detector, not for reconstruction.
  a_WAvg = forall<T,NUMCOMPS,MEM>(f_sphToPrim,a_UAvg);
  const bool magneticTest =
      isConstantMagneticFieldTest(ParseInputs::get_init_condition_type());
  BoxData<T,1,MEM> deltaMe(a_WAvg.box());
  deltaMe.setVal(0.0);
  if (magneticTest)
  {
    testProblemDeltaMagneticEnergy(deltaMe,a_WAvg,a_block);
    addTestProblemBackgroundB(a_WAvg,a_block);
  }

  BoxData<T,1,MEM> keAvg = kineticEnergyW(a_WAvg,a_block,-1);
  BoxData<T,1,MEM> meAvg = magneticEnergyW(a_WAvg,a_block,-1);

  forallInPlace(
      [ ] PROTO_LAMBDA
      (
        Var<T,NUMCOMPS>& a_w,
        const Var<T,1,MEM>& a_ke,
        const Var<T,1,MEM>& a_me,
        const Var<T,1,MEM>& a_deltaMe,
        bool a_magneticTest,
        T a_gamma)
      {
        if (a_magneticTest)
        {
          a_w(iP) = c_TEST_P0 + (a_gamma - 1.0)*
              (a_w(iP) - a_ke(0) - a_deltaMe(0));
        }
        else
        {
          a_w(iP) = (a_w(iP) - a_ke(0) - a_me(0))*(a_gamma - 1.0);
        }
        if (a_w(iP) < 0.0) a_w(iP) = 0.0;
      },
      a_WAvg,
      keAvg,
      meAvg,
      deltaMe,
      magneticTest,
      m_gamma);
}


// ====================================================================================
// 2. Berta Jameson-type discontinuity detector.
// ====================================================================================

template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_berta_jameson_theta_rawF(
    Var<T,1,MEM>&              a_theta,
    const Var<T,NUMCOMPS,MEM>& a_Wm_x,
    const Var<T,NUMCOMPS,MEM>& a_W0,
    const Var<T,NUMCOMPS,MEM>& a_Wp_x,
    const Var<T,NUMCOMPS,MEM>& a_Wm_y,
    const Var<T,NUMCOMPS,MEM>& a_Wp_y,
    const Var<T,NUMCOMPS,MEM>& a_Wm_z,
    const Var<T,NUMCOMPS,MEM>& a_Wp_z,
    const T                    a_eta_d,
    const T                    a_eps)
{
  auto eta_one_dir = [&] PROTO_LAMBDA
    (
      const T qm,
      const T q0,
      const T qp)
    {
      T num = abs(qp - 2.0*q0 + qm);
      T den = abs(qp) + 2.0*abs(q0) + abs(qm) + a_eps;
      return num/den;
    };

  auto magp = [&] PROTO_LAMBDA
    (
      const Var<T,NUMCOMPS,MEM>& W)
    {
      T bx = W(iBX);
      T by = W(iBY);
      T bz = W(iBZ);
      return 0.5*(bx*bx + by*by + bz*bz);
    };

  T rho0 = a_W0(iRHO);
  T p0   = a_W0(iP);
  T mp0  = magp(a_W0);

  T etax_rho = eta_one_dir(a_Wm_x(iRHO),rho0,a_Wp_x(iRHO));
  T etax_p   = eta_one_dir(a_Wm_x(iP),  p0,  a_Wp_x(iP));
  T etax_mp  = eta_one_dir(magp(a_Wm_x),mp0, magp(a_Wp_x));
  T etax     = max(etax_rho,max(etax_p,etax_mp));

  T etay_rho = eta_one_dir(a_Wm_y(iRHO),rho0,a_Wp_y(iRHO));
  T etay_p   = eta_one_dir(a_Wm_y(iP),  p0,  a_Wp_y(iP));
  T etay_mp  = eta_one_dir(magp(a_Wm_y),mp0, magp(a_Wp_y));
  T etay     = max(etay_rho,max(etay_p,etay_mp));

  T etaz_rho = eta_one_dir(a_Wm_z(iRHO),rho0,a_Wp_z(iRHO));
  T etaz_p   = eta_one_dir(a_Wm_z(iP),  p0,  a_Wp_z(iP));
  T etaz_mp  = eta_one_dir(magp(a_Wm_z),mp0, magp(a_Wp_z));
  T etaz     = max(etaz_rho,max(etaz_p,etaz_mp));

  T eta = sqrt(etax*etax + etay*etay + etaz*etaz);

  // theta = 1: smooth, use 4th-order correction.
  // theta = 0: troubled, reduce order.
  a_theta(0) = (eta < a_eta_d) ? 1.0 : 0.0;
}
PROTO_KERNEL_END(f_berta_jameson_theta_rawF, f_berta_jameson_theta_raw)


template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_berta_dilate_thetaF(
    Var<T,1,MEM>&       a_theta,
    const Var<T,1,MEM>& a_t0,
    const Var<T,1,MEM>& a_tm_x,
    const Var<T,1,MEM>& a_tp_x,
    const Var<T,1,MEM>& a_tm_y,
    const Var<T,1,MEM>& a_tp_y,
    const Var<T,1,MEM>& a_tm_z,
    const Var<T,1,MEM>& a_tp_z)
{
  T tmin = a_t0(0);

  tmin = min(tmin,a_tm_x(0));
  tmin = min(tmin,a_tp_x(0));
  tmin = min(tmin,a_tm_y(0));
  tmin = min(tmin,a_tp_y(0));
  tmin = min(tmin,a_tm_z(0));
  tmin = min(tmin,a_tp_z(0));

  // If this cell or any first neighbor is troubled, theta=0.
  a_theta(0) = tmin;
}
PROTO_KERNEL_END(f_berta_dilate_thetaF, f_berta_dilate_theta)


template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::BertaComputeThetaJameson(
    BoxData<T,1,MEM>& a_theta,
    const State&      a_WAvgForDetector,
    T                 a_eta_d) const
{
  PR_TIMERS("BoxOp_EulerCubedSphere::BertaComputeThetaJameson()");

  Point ex = Point::Basis(0);
  Point ey = Point::Basis(1);
  Point ez = Point::Basis(2);

  // Raw detector needs first neighbors.
  Box bxRaw = a_WAvgForDetector.box().grow(-1);
  BoxData<T,1,MEM> thetaRaw(bxRaw);

  const T eps = 1.0e-40;

  State Wm_x = (1.0*Shift(-ex))(a_WAvgForDetector);
  State W0   = (1.0*Shift(Point::Zeros()))(a_WAvgForDetector);
  State Wp_x = (1.0*Shift( ex))(a_WAvgForDetector);
  State Wm_y = (1.0*Shift(-ey))(a_WAvgForDetector);
  State Wp_y = (1.0*Shift( ey))(a_WAvgForDetector);
  State Wm_z = (1.0*Shift(-ez))(a_WAvgForDetector);
  State Wp_z = (1.0*Shift( ez))(a_WAvgForDetector);

  forallInPlace(
      f_berta_jameson_theta_raw,
      thetaRaw,
      Wm_x,
      W0,
      Wp_x,
      Wm_y,
      Wp_y,
      Wm_z,
      Wp_z,
      a_eta_d,
      eps);

  // Dilation needs first neighbors of thetaRaw.
  Box bxInterior = thetaRaw.box().grow(-1);
  BoxData<T,1,MEM> thetaInterior(bxInterior);

  BoxData<T,1,MEM> theta0   = (1.0*Shift(Point::Zeros()))(thetaRaw);
  BoxData<T,1,MEM> thetam_x = (1.0*Shift(-ex))(thetaRaw);
  BoxData<T,1,MEM> thetap_x = (1.0*Shift( ex))(thetaRaw);
  BoxData<T,1,MEM> thetam_y = (1.0*Shift(-ey))(thetaRaw);
  BoxData<T,1,MEM> thetap_y = (1.0*Shift( ey))(thetaRaw);
  BoxData<T,1,MEM> thetam_z = (1.0*Shift(-ez))(thetaRaw);
  BoxData<T,1,MEM> thetap_z = (1.0*Shift( ez))(thetaRaw);

  forallInPlace(
      f_berta_dilate_theta,
      thetaInterior,
      theta0,
      thetam_x,
      thetap_x,
      thetam_y,
      thetap_y,
      thetam_z,
      thetap_z);

  // Full-box theta. Cells without enough detector stencil default to fallback.
  a_theta.define(a_WAvgForDetector.box());
  a_theta.setVal(0.0);
  thetaInterior.copyTo(a_theta);
}


// ====================================================================================
// 3. Berta derivative-ratio discontinuity detector.
// ====================================================================================

template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_berta_derivative_ratio_theta_rawF(
    Var<T,1,MEM>&              a_theta,
    const Var<T,NUMCOMPS,MEM>& a_Wmm_x,
    const Var<T,NUMCOMPS,MEM>& a_Wm_x,
    const Var<T,NUMCOMPS,MEM>& a_W0,
    const Var<T,NUMCOMPS,MEM>& a_Wp_x,
    const Var<T,NUMCOMPS,MEM>& a_Wpp_x,
    const Var<T,NUMCOMPS,MEM>& a_Wmm_y,
    const Var<T,NUMCOMPS,MEM>& a_Wm_y,
    const Var<T,NUMCOMPS,MEM>& a_Wp_y,
    const Var<T,NUMCOMPS,MEM>& a_Wpp_y,
    const Var<T,NUMCOMPS,MEM>& a_Wmm_z,
    const Var<T,NUMCOMPS,MEM>& a_Wm_z,
    const Var<T,NUMCOMPS,MEM>& a_Wp_z,
    const Var<T,NUMCOMPS,MEM>& a_Wpp_z,
    const T                    a_eta_d,
    const T                    a_gamma,
    const T                    a_eps)
{
  auto mag = [&] PROTO_LAMBDA
    (
      const Var<T,NUMCOMPS,MEM>& W)
    {
      T bx = W(iBX);
      T by = W(iBY);
      T bz = W(iBZ);
      return sqrt(bx*bx + by*by + bz*bz);
    };

  auto ref_value = [&] PROTO_LAMBDA
    (
      const int comp)
    {
      if (comp == iRHO || comp == iP)
      {
        return abs(a_W0(comp));
      }
      if (comp == iVX || comp == iVY || comp == iVZ)
      {
        T rho = max(abs(a_W0(iRHO)),a_eps);
        T p   = max(abs(a_W0(iP)),a_eps);
        return sqrt(abs(a_gamma)*p/rho);
      }
      return mag(a_W0);
    };

  auto eta_comp = [&] PROTO_LAMBDA
    (
      const T qmm,
      const T qm,
      const T q0,
      const T qp,
      const T qpp,
      const T qref)
    {
      T d1 = 0.5*(qp - qm);
      T d2 = qp - 2.0*q0 + qm;
      T d3 = 0.5*(qpp - 2.0*qp + 2.0*qm - qmm);
      T d4 = qpp - 4.0*qp + 6.0*q0 - 4.0*qm + qmm;

      T eta_odd  = abs(d3)/(abs(qref) + abs(d1) + abs(d3) + a_eps);
      T eta_even = abs(d4)/(abs(qref) + abs(d2) + abs(d4) + a_eps);
      return max(eta_odd,eta_even);
    };

  auto eta_dir = [&] PROTO_LAMBDA
    (
      const Var<T,NUMCOMPS,MEM>& Wmm,
      const Var<T,NUMCOMPS,MEM>& Wm,
      const Var<T,NUMCOMPS,MEM>& Wp,
      const Var<T,NUMCOMPS,MEM>& Wpp)
    {
      T eta = 0.0;
      for (int comp = 0; comp < NUMCOMPS; comp++)
      {
        eta = max(eta,eta_comp(Wmm(comp),Wm(comp),a_W0(comp),Wp(comp),Wpp(comp),ref_value(comp)));
      }
      return eta;
    };

  T etax = eta_dir(a_Wmm_x,a_Wm_x,a_Wp_x,a_Wpp_x);
  T etay = eta_dir(a_Wmm_y,a_Wm_y,a_Wp_y,a_Wpp_y);
  T etaz = eta_dir(a_Wmm_z,a_Wm_z,a_Wp_z,a_Wpp_z);
  T eta  = sqrt(etax*etax + etay*etay + etaz*etaz);

  a_theta(0) = (eta < a_eta_d) ? 1.0 : 0.0;
}
PROTO_KERNEL_END(f_berta_derivative_ratio_theta_rawF, f_berta_derivative_ratio_theta_raw)


template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::BertaComputeThetaDerivativeRatio(
    BoxData<T,1,MEM>& a_theta,
    const State&      a_WAvgForDetector,
    T                 a_eta_d) const
{
  PR_TIMERS("BoxOp_EulerCubedSphere::BertaComputeThetaDerivativeRatio()");

  Point ex = Point::Basis(0);
  Point ey = Point::Basis(1);
  Point ez = Point::Basis(2);

  Box bxRaw = a_WAvgForDetector.box().grow(-2);
  BoxData<T,1,MEM> thetaRaw(bxRaw);

  const T eps = 1.0e-40;

  State Wmm_x = (1.0*Shift(-2*ex))(a_WAvgForDetector);
  State Wm_x  = (1.0*Shift(-ex))(a_WAvgForDetector);
  State W0    = (1.0*Shift(Point::Zeros()))(a_WAvgForDetector);
  State Wp_x  = (1.0*Shift( ex))(a_WAvgForDetector);
  State Wpp_x = (1.0*Shift( 2*ex))(a_WAvgForDetector);

  State Wmm_y = (1.0*Shift(-2*ey))(a_WAvgForDetector);
  State Wm_y  = (1.0*Shift(-ey))(a_WAvgForDetector);
  State Wp_y  = (1.0*Shift( ey))(a_WAvgForDetector);
  State Wpp_y = (1.0*Shift( 2*ey))(a_WAvgForDetector);

  State Wmm_z = (1.0*Shift(-2*ez))(a_WAvgForDetector);
  State Wm_z  = (1.0*Shift(-ez))(a_WAvgForDetector);
  State Wp_z  = (1.0*Shift( ez))(a_WAvgForDetector);
  State Wpp_z = (1.0*Shift( 2*ez))(a_WAvgForDetector);

  forallInPlace(
      f_berta_derivative_ratio_theta_raw,
      thetaRaw,
      Wmm_x,
      Wm_x,
      W0,
      Wp_x,
      Wpp_x,
      Wmm_y,
      Wm_y,
      Wp_y,
      Wpp_y,
      Wmm_z,
      Wm_z,
      Wp_z,
      Wpp_z,
      a_eta_d,
      m_gamma,
      eps);

  // Flag the first neighbors too, as in the paper's troubled-cell propagation.
  Box bxInterior = thetaRaw.box().grow(-1);
  BoxData<T,1,MEM> thetaInterior(bxInterior);

  BoxData<T,1,MEM> theta0   = (1.0*Shift(Point::Zeros()))(thetaRaw);
  BoxData<T,1,MEM> thetam_x = (1.0*Shift(-ex))(thetaRaw);
  BoxData<T,1,MEM> thetap_x = (1.0*Shift( ex))(thetaRaw);
  BoxData<T,1,MEM> thetam_y = (1.0*Shift(-ey))(thetaRaw);
  BoxData<T,1,MEM> thetap_y = (1.0*Shift( ey))(thetaRaw);
  BoxData<T,1,MEM> thetam_z = (1.0*Shift(-ez))(thetaRaw);
  BoxData<T,1,MEM> thetap_z = (1.0*Shift( ez))(thetaRaw);

  forallInPlace(
      f_berta_dilate_theta,
      thetaInterior,
      theta0,
      thetam_x,
      thetap_x,
      thetam_y,
      thetap_y,
      thetam_z,
      thetap_z);

  a_theta.define(a_WAvgForDetector.box());
  a_theta.setVal(0.0);
  thetaInterior.copyTo(a_theta);
}


// ====================================================================================
// 4. Theta-aware average-to-point recovery.
// ====================================================================================

template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_berta_avg_to_point_thetaF(
    Var<T,NUMCOMPS,MEM>&       a_UPoint,
    const Var<T,NUMCOMPS,MEM>& a_UAvg,
    const Var<T,NUMCOMPS,MEM>& a_LapUAvg,
    const Var<T,1,MEM>&        a_theta)
{
  for (int comp = 0; comp < NUMCOMPS; comp++)
  {
    if (a_theta(0) > 0.5)
    {
      a_UPoint(comp) = a_UAvg(comp) - a_LapUAvg(comp)/24.0;
    }
    else
    {
      a_UPoint(comp) = a_UAvg(comp);
    }
  }
}
PROTO_KERNEL_END(f_berta_avg_to_point_thetaF, f_berta_avg_to_point_theta)


template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::BertaAverageToPoint(
    State&                  a_UPoint,
    const State&            a_UAvg,
    const BoxData<T,1,MEM>& a_theta) const
{
  PR_TIMERS("BoxOp_EulerCubedSphere::BertaAverageToPoint(theta)");

  State lapUAvg = Stencil<T>::Laplacian()(a_UAvg);

  a_UPoint.define(a_UAvg.box());

  forallInPlace(
      f_berta_avg_to_point_theta,
      a_UPoint,
      a_UAvg,
      lapUAvg,
      a_theta);
}


// ====================================================================================
// 4. WENOZ, WENO3, and PLM face reconstruction.
// Corrected convention:
//   face index i = lower face of cell i = x_{i-1/2}
//   adjacent cells are i-1 and i.
// ====================================================================================

template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_wenoz_pointwise_faceF(
    Var<T,NUMCOMPS,MEM>&       a_WL,
    Var<T,NUMCOMPS,MEM>&       a_WR,
    const Var<T,NUMCOMPS,MEM>& a_Wm2,
    const Var<T,NUMCOMPS,MEM>& a_Wm1,
    const Var<T,NUMCOMPS,MEM>& a_W0,
    const Var<T,NUMCOMPS,MEM>& a_Wp1,
    const Var<T,NUMCOMPS,MEM>& a_Wp2,
    const Var<T,NUMCOMPS,MEM>& a_Wp3)
{
  const T eps = 1.0e-40;
  const T d0  = 1.0/16.0;
  const T d1  = 5.0/8.0;
  const T d2  = 5.0/16.0;

  for (int comp = 0; comp < NUMCOMPS; comp++)
  {
    // Left state at target face, using a_W0 as the left cell.
    T vm2 = a_Wm2(comp);
    T vm1 = a_Wm1(comp);
    T v0  = a_W0(comp);
    T vp1 = a_Wp1(comp);
    T vp2 = a_Wp2(comp);

    T q0 = ( 3.0*vm2 - 10.0*vm1 + 15.0*v0  ) / 8.0;
    T q1 = (    -vm1 +  6.0*v0  +  3.0*vp1 ) / 8.0;
    T q2 = ( 3.0*v0  +  6.0*vp1 -      vp2 ) / 8.0;

    T b0a = vm2 - 2.0*vm1 + v0;
    T b0b = vm2 - 4.0*vm1 + 3.0*v0;

    T b1a = vm1 - 2.0*v0 + vp1;
    T b1b = vm1 - vp1;

    T b2a = v0 - 2.0*vp1 + vp2;
    T b2b = 3.0*v0 - 4.0*vp1 + vp2;

    T beta0 = (13.0/12.0)*b0a*b0a + 0.25*b0b*b0b;
    T beta1 = (13.0/12.0)*b1a*b1a + 0.25*b1b*b1b;
    T beta2 = (13.0/12.0)*b2a*b2a + 0.25*b2b*b2b;

    T tau5 = abs(beta0 - beta2);

    T alpha0 = d0*(1.0 + tau5/(beta0 + eps));
    T alpha1 = d1*(1.0 + tau5/(beta1 + eps));
    T alpha2 = d2*(1.0 + tau5/(beta2 + eps));

    T asum = alpha0 + alpha1 + alpha2;

    a_WL(comp) = (alpha0*q0 + alpha1*q1 + alpha2*q2)/asum;

    // Right state at target face: reversed stencil, using a_Wp1 as the right cell.
    T rm2 = a_Wp3(comp);
    T rm1 = a_Wp2(comp);
    T r0  = a_Wp1(comp);
    T rp1 = a_W0(comp);
    T rp2 = a_Wm1(comp);

    T rq0 = ( 3.0*rm2 - 10.0*rm1 + 15.0*r0  ) / 8.0;
    T rq1 = (    -rm1 +  6.0*r0  +  3.0*rp1 ) / 8.0;
    T rq2 = ( 3.0*r0  +  6.0*rp1 -      rp2 ) / 8.0;

    T rb0a = rm2 - 2.0*rm1 + r0;
    T rb0b = rm2 - 4.0*rm1 + 3.0*r0;

    T rb1a = rm1 - 2.0*r0 + rp1;
    T rb1b = rm1 - rp1;

    T rb2a = r0 - 2.0*rp1 + rp2;
    T rb2b = 3.0*r0 - 4.0*rp1 + rp2;

    T rbeta0 = (13.0/12.0)*rb0a*rb0a + 0.25*rb0b*rb0b;
    T rbeta1 = (13.0/12.0)*rb1a*rb1a + 0.25*rb1b*rb1b;
    T rbeta2 = (13.0/12.0)*rb2a*rb2a + 0.25*rb2b*rb2b;

    T rtau5 = abs(rbeta0 - rbeta2);

    T ralpha0 = d0*(1.0 + rtau5/(rbeta0 + eps));
    T ralpha1 = d1*(1.0 + rtau5/(rbeta1 + eps));
    T ralpha2 = d2*(1.0 + rtau5/(rbeta2 + eps));

    T rasum = ralpha0 + ralpha1 + ralpha2;

    a_WR(comp) = (ralpha0*rq0 + ralpha1*rq1 + ralpha2*rq2)/rasum;
  }
}
PROTO_KERNEL_END(f_wenoz_pointwise_faceF, f_wenoz_pointwise_face)


template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_weno3_pointwise_faceF(
    Var<T,NUMCOMPS,MEM>&       a_WL,
    Var<T,NUMCOMPS,MEM>&       a_WR,
    const Var<T,NUMCOMPS,MEM>& a_Wm1,
    const Var<T,NUMCOMPS,MEM>& a_W0,
    const Var<T,NUMCOMPS,MEM>& a_Wp1,
    const Var<T,NUMCOMPS,MEM>& a_Wp2,
    const T                    a_dx2)
{
  const T d0  = 3.0/4.0;
  const T d1  = 1.0/4.0;
  const T eps = a_dx2;

  for (int comp = 0; comp < NUMCOMPS; comp++)
  {
    // Left state using a_W0 as the left cell of the face.
    T vm1 = a_Wm1(comp);
    T v0  = a_W0(comp);
    T vp1 = a_Wp1(comp);

    T Dm = v0  - vm1;
    T Dp = vp1 - v0;

    T beta0 = Dm*Dm;
    T beta1 = Dp*Dp;
    T tau   = (Dp - Dm)*(Dp - Dm);

    T a0 = d0*(1.0 + tau/(beta0 + eps));
    T a1 = d1*(1.0 + tau/(beta1 + eps));

    T q0 = 0.5*(vp1 + v0);
    T q1 = 0.5*(-vp1 + 3.0*v0);

    a_WL(comp) = (a0*q0 + a1*q1)/(a0 + a1);

    // Right state using a_Wp1 as the right cell of the face.
    T rm1 = a_Wp2(comp);
    T r0  = a_Wp1(comp);
    T rp1 = a_W0(comp);

    T rDm = r0  - rm1;
    T rDp = rp1 - r0;

    T rbeta0 = rDm*rDm;
    T rbeta1 = rDp*rDp;
    T rtau   = (rDp - rDm)*(rDp - rDm);

    T ra0 = d0*(1.0 + rtau/(rbeta0 + eps));
    T ra1 = d1*(1.0 + rtau/(rbeta1 + eps));

    T rq0 = 0.5*(rp1 + r0);
    T rq1 = 0.5*(-rp1 + 3.0*r0);

    a_WR(comp) = (ra0*rq0 + ra1*rq1)/(ra0 + ra1);
  }
}
PROTO_KERNEL_END(f_weno3_pointwise_faceF, f_weno3_pointwise_face)


template<typename T>
inline T berta_minmod3_scalar(T a, T b, T c)
{
  if (a > 0.0 && b > 0.0 && c > 0.0) return min(a,min(b,c));
  if (a < 0.0 && b < 0.0 && c < 0.0) return max(a,max(b,c));
  return 0.0;
}


template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_plm_pointwise_faceF(
    Var<T,NUMCOMPS,MEM>&       a_WL,
    Var<T,NUMCOMPS,MEM>&       a_WR,
    const Var<T,NUMCOMPS,MEM>& a_Wm1,
    const Var<T,NUMCOMPS,MEM>& a_W0,
    const Var<T,NUMCOMPS,MEM>& a_Wp1,
    const Var<T,NUMCOMPS,MEM>& a_Wp2)
{
  const T theta = 1.5;

  for (int comp = 0; comp < NUMCOMPS; comp++)
  {
    // Left state from cell represented by a_W0.
    T dm = a_W0(comp)  - a_Wm1(comp);
    T dp = a_Wp1(comp) - a_W0(comp);
    T dc = 0.5*(a_Wp1(comp) - a_Wm1(comp));

    T slopeL = berta_minmod3_scalar(theta*dm,dc,theta*dp);

    a_WL(comp) = a_W0(comp) + 0.5*slopeL;

    // Right state from cell represented by a_Wp1.
    T rdm = a_Wp1(comp) - a_W0(comp);
    T rdp = a_Wp2(comp) - a_Wp1(comp);
    T rdc = 0.5*(a_Wp2(comp) - a_W0(comp));

    T slopeR = berta_minmod3_scalar(theta*rdm,rdc,theta*rdp);

    a_WR(comp) = a_Wp1(comp) - 0.5*slopeR;
  }
}
PROTO_KERNEL_END(f_plm_pointwise_faceF, f_plm_pointwise_face)


template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_face_theta_adjacentF(
    Var<T,1,MEM>&       a_thetaFace,
    const Var<T,1,MEM>& a_thetaLeft,
    const Var<T,1,MEM>& a_thetaRight)
{
  // Face i lies between cells i-1 and i.
  a_thetaFace(0) = min(a_thetaLeft(0),a_thetaRight(0));
}
PROTO_KERNEL_END(f_face_theta_adjacentF, f_face_theta_adjacent)


template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_face_theta_recon_stencilF(
    Var<T,1,MEM>&       a_thetaFace,
    const Var<T,1,MEM>& a_thetam3,
    const Var<T,1,MEM>& a_thetam2,
    const Var<T,1,MEM>& a_thetam1,
    const Var<T,1,MEM>& a_theta0,
    const Var<T,1,MEM>& a_thetap1,
    const Var<T,1,MEM>& a_thetap2)
{
  // WENOZ at lower face of cell i uses cells i-3 ... i+2.
  // If any cell in that stencil is troubled, use fallback reconstruction.
  T tmin = a_thetam3(0);
  tmin = min(tmin,a_thetam2(0));
  tmin = min(tmin,a_thetam1(0));
  tmin = min(tmin,a_theta0(0));
  tmin = min(tmin,a_thetap1(0));
  tmin = min(tmin,a_thetap2(0));
  a_thetaFace(0) = tmin;
}
PROTO_KERNEL_END(f_face_theta_recon_stencilF, f_face_theta_recon_stencil)


template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_select_reconstructionF(
    Var<T,NUMCOMPS,MEM>&       a_WL,
    Var<T,NUMCOMPS,MEM>&       a_WR,
    const Var<T,NUMCOMPS,MEM>& a_WL4,
    const Var<T,NUMCOMPS,MEM>& a_WR4,
    const Var<T,NUMCOMPS,MEM>& a_WLF,
    const Var<T,NUMCOMPS,MEM>& a_WRF,
    const Var<T,1,MEM>&        a_thetaReconFace)
{
  T th = a_thetaReconFace(0);

  if (th > 0.5)
  {
    for (int comp = 0; comp < NUMCOMPS; comp++)
    {
      a_WL(comp) = a_WL4(comp);
      a_WR(comp) = a_WR4(comp);
    }
  }
  else
  {
    for (int comp = 0; comp < NUMCOMPS; comp++)
    {
      a_WL(comp) = a_WLF(comp);
      a_WR(comp) = a_WRF(comp);
    }
  }

  // Positivity safety.
  a_WL(iRHO) = max(a_WL(iRHO),RHO00);
  a_WR(iRHO) = max(a_WR(iRHO),RHO00);

  a_WL(iP) = max(a_WL(iP),1.0e-12);
  a_WR(iP) = max(a_WR(iP),1.0e-12);
}
PROTO_KERNEL_END(f_select_reconstructionF, f_select_reconstruction)


template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::BertaReconstructWithFallback(
    State&                  a_WFace_L,
    State&                  a_WFace_R,
    BoxData<T,1,MEM>&       a_thetaFluxFace,
    BoxData<T,1,MEM>&       a_thetaReconFace,
    const State&            a_WPoint,
    const BoxData<T,1,MEM>& a_thetaCell,
    uint                    a_dir,
    int                     a_fallback_type) const
{
  PR_TIMERS("BoxOp_EulerCubedSphere::BertaReconstructWithFallback()");

  Point e = Point::Basis(a_dir);

  // Face index i = lower face of cell i = x_{i-1/2}.
  // WENOZ needs cells i-3 ... i+2.
  // grow(-3) is conservative. You may later optimize asymmetrically.
  Box faceBox = a_WPoint.box().extrude(e).grow(-3);

  State WL_wenoz(faceBox);
  State WR_wenoz(faceBox);
  State WL_fallback(faceBox);
  State WR_fallback(faceBox);

  a_WFace_L.define(faceBox);
  a_WFace_R.define(faceBox);
  a_thetaFluxFace.define(faceBox);
  a_thetaReconFace.define(faceBox);

  // Adjacent-cell theta for flux point-to-average correction.
  // Face i lies between cells i-1 and i.
  BoxData<T,1,MEM> thetam1 = (1.0*Shift(-1*e))(a_thetaCell);
  BoxData<T,1,MEM> theta0  = (1.0*Shift(Point::Zeros()))(a_thetaCell);

  forallInPlace(
      f_face_theta_adjacent,
      a_thetaFluxFace,
      thetam1,
      theta0);

  // Full-stencil theta for reconstruction fallback.
  BoxData<T,1,MEM> thetam3 = (1.0*Shift(-3*e))(a_thetaCell);
  BoxData<T,1,MEM> thetam2 = (1.0*Shift(-2*e))(a_thetaCell);
  BoxData<T,1,MEM> thetap1 = (1.0*Shift( 1*e))(a_thetaCell);
  BoxData<T,1,MEM> thetap2 = (1.0*Shift( 2*e))(a_thetaCell);

  forallInPlace(
      f_face_theta_recon_stencil,
      a_thetaReconFace,
      thetam3,
      thetam2,
      thetam1,
      theta0,
      thetap1,
      thetap2);

  State Wm3 = (1.0*Shift(-3*e))(a_WPoint);
  State Wm2 = (1.0*Shift(-2*e))(a_WPoint);
  State Wm1 = (1.0*Shift(-1*e))(a_WPoint);
  State W0  = (1.0*Shift(Point::Zeros()))(a_WPoint);
  State Wp1 = (1.0*Shift( 1*e))(a_WPoint);
  State Wp2 = (1.0*Shift( 2*e))(a_WPoint);

  // WENOZ candidate at lower face of cell i.
  forallInPlace(
      f_wenoz_pointwise_face,
      WL_wenoz,
      WR_wenoz,
      Wm3,
      Wm2,
      Wm1,
      W0,
      Wp1,
      Wp2);

  if (a_fallback_type == 0)
  {
    // WENO3 fallback at lower face of cell i.
    T dx2 = this->dx()[a_dir]*this->dx()[a_dir];

    forallInPlace(
        f_weno3_pointwise_face,
        WL_fallback,
        WR_fallback,
        Wm2,
        Wm1,
        W0,
        Wp1,
        dx2);
  }
  else
  {
    // PLM fallback at lower face of cell i.
    forallInPlace(
        f_plm_pointwise_face,
        WL_fallback,
        WR_fallback,
        Wm2,
        Wm1,
        W0,
        Wp1);
  }

  forallInPlace(
      f_select_reconstruction,
      a_WFace_L,
      a_WFace_R,
      WL_wenoz,
      WR_wenoz,
      WL_fallback,
      WR_fallback,
      a_thetaReconFace);
}


// ====================================================================================
// 5. Theta-aware point-to-average flux correction.
// ====================================================================================

template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_flux_point_to_avg_thetaF(
    Var<T,NUMCOMPS,MEM>&       a_fluxAvg,
    const Var<T,NUMCOMPS,MEM>& a_fluxPoint,
    const Var<T,NUMCOMPS,MEM>& a_lapFluxPoint,
    const Var<T,1,MEM>&        a_thetaFluxFace)
{
  for (int comp = 0; comp < NUMCOMPS; comp++)
  {
    if (a_thetaFluxFace(0) > 0.5)
    {
      a_fluxAvg(comp) = a_fluxPoint(comp) + a_lapFluxPoint(comp)/24.0;
    }
    else
    {
      a_fluxAvg(comp) = a_fluxPoint(comp);
    }
  }
}
PROTO_KERNEL_END(f_flux_point_to_avg_thetaF, f_flux_point_to_avg_theta)


template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::BertaPrimPointToFluxAverage(
           State&                  a_physFlux,
           const State&            a_WFacePoint,
           const BoxData<T,1,MEM>& a_thetaFluxFace,
           const BoxData<T,1,MEM>& a_adjDdir,
           uint                    a_block,
           uint                    a_dir) const
{
  PR_TIMERS("BoxOp_EulerCubedSphere::BertaPrimPointToFluxAverage()");

  Stencil<T> deconvolveStencilRadial =
      (1.0*Shift(Point::Zeros())
       - (1.0/24.0)*Stencil<T>::Derivative(2,0,2));

  BoxData<T,1,MEM> kePoint(a_WFacePoint.box());
  BoxData<T,1,MEM> mePoint(a_WFacePoint.box());
  BoxData<T,1,MEM> vbPoint(a_WFacePoint.box());

  EnergyW(kePoint,mePoint,vbPoint,a_WFacePoint,a_block,a_dir);

  BoxData<T,1,MEM> adjDdirPoint;
  if (a_dir != 0)
  {
    adjDdirPoint = deconvolveStencilRadial(a_adjDdir);
  }
  else
  {
    adjDdirPoint.define(a_adjDdir.box());
    a_adjDdir.copyTo(adjDdirPoint);
  }

  State fluxPoint =
      forall<T,NUMCOMPS>(
          f_primflux,
          a_WFacePoint,
          kePoint,
          mePoint,
          vbPoint,
          adjDdirPoint,
          m_gamma,
          a_block,
          a_dir);

  State lapFluxPoint = Stencil<T>::LaplacianFace(a_dir)(fluxPoint);

  a_physFlux.define(fluxPoint.box());

  forallInPlace(
      f_flux_point_to_avg_theta,
      a_physFlux,
      fluxPoint,
      lapFluxPoint,
      a_thetaFluxFace);
}


// ====================================================================================
// 6. Mapped flux from pointwise face primitive state.
// ====================================================================================

template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::BertaPrimPointToMappedFlux(
    State&                         a_fluxdir,
    const State&                   a_physFluxAvg,
    const State&                   a_WFacePoint,
    const BoxData<T,1,MEM>&        a_thetaFluxFace,
    const BoxData<T,1,MEM>&        a_adjDrDir,
    const BoxData<T,DIM,MEM,DIM>&  a_deta_a,
    const BoxData<T,1,MEM>&        a_deta,
    const BoxData<T,DIM,MEM>&      a_ntrow,
    uint                           a_block,
    uint                           a_dir) const
{
  PR_TIMERS("BoxOp_EulerCubedSphere::BertaPrimPointToMappedFlux()");

  BoxData<T,DIM,MEM> fluxdirMom =
      slice<T,NUMCOMPS,DIM,MEM>(a_fluxdir,iMOMX);

  BoxData<T,DIM,MEM> physFluxMom =
      slice<T,NUMCOMPS,DIM,MEM>(a_physFluxAvg,iMOMX);

  BoxData<T,DIM,MEM> tmp =
      Operator::_faceMatrixProductAB(
          a_deta_a,physFluxMom,
          a_deta_a,physFluxMom,
          a_dir);

  tmp.copyTo(fluxdirMom);

  BoxData<T,DIM,MEM> fluxdirB =
      slice<T,NUMCOMPS,DIM,MEM>(a_fluxdir,iBX);

  BoxData<T,DIM,MEM> physFluxB =
      slice<T,NUMCOMPS,DIM,MEM>(a_physFluxAvg,iBX);

  BoxData<T,DIM,MEM> tmp2 =
      Operator::_faceMatrixProductAB(
          a_deta_a,physFluxB,
          a_deta_a,physFluxB,
          a_dir);

  tmp2.copyTo(fluxdirB);

  forallInPlace(
      [ ] PROTO_LAMBDA
      (
        Var<T,NUMCOMPS,MEM>&       a_advflux,
        const Var<T,NUMCOMPS,MEM>& a_fluxin,
        const Var<T,1,MEM>&        a_deta)
      {
        a_advflux(iRHO) = a_deta(0)*a_fluxin(iRHO);
        a_advflux(iE  ) = a_deta(0)*a_fluxin(iE);
      },
      a_fluxdir,
      a_physFluxAvg,
      a_deta);

  for (int comp = 0; comp < NUMCOMPS; comp++)
  {
    BoxData<T,1,MEM> advfluxcomp = slice(a_fluxdir,comp);
    BoxData<T,1,MEM> advfluxprod =
        faceProduct(advfluxcomp,a_adjDrDir,a_dir);

    advfluxprod.copyTo(advfluxcomp);
  }

  // Pressure force from point total pressure, then theta-aware face averaging.
  BoxData<T,1,MEM> pPoint  = slice(a_WFacePoint,iP);
  BoxData<T,1,MEM> mePoint = magneticEnergyW(a_WFacePoint,a_block,a_dir);
  BoxData<T,1,MEM> pTotPoint = pPoint + mePoint;

  BoxData<T,1,MEM> lapPTotPoint =
      Stencil<T>::LaplacianFace(a_dir)(pTotPoint);

  BoxData<T,1,MEM> pTotAvg(pTotPoint.box());

  forallInPlace(
      [ ] PROTO_LAMBDA
      (
        Var<T,1,MEM>&       a_pAvg,
        const Var<T,1,MEM>& a_pPoint,
        const Var<T,1,MEM>& a_lapP,
        const Var<T,1,MEM>& a_theta)
      {
        if (a_theta(0) > 0.5)
        {
          a_pAvg(0) = a_pPoint(0) + a_lapP(0)/24.0;
        }
        else
        {
          a_pAvg(0) = a_pPoint(0);
        }
      },
      pTotAvg,
      pTotPoint,
      lapPTotPoint,
      a_thetaFluxFace);

  BoxData<T,1,MEM> tmp3 =
      faceProduct(pTotAvg,a_adjDrDir,a_dir);

  BoxData<T,DIM,MEM> pforce(a_ntrow.box());
  pforce.setVal(0.0);

  for (int dirvec = 0; dirvec < DIM; dirvec++)
  {
    BoxData<T,1,MEM> ntSlice = slice(a_ntrow,dirvec);
    BoxData<T,1,MEM> pforceSlice = slice(pforce,dirvec);

    BoxData<T,1,MEM> tmpPF =
        faceProduct(tmp3,ntSlice,a_dir);

    tmpPF.copyTo(pforceSlice);
  }

  fluxdirMom += pforce;
}


// ====================================================================================
// 7. Rusanov dissipation from pointwise face states.
// ====================================================================================

template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::BertaPrimPointToMappedRusanovFlux(
    State&                         a_fluxdir,
    const State&                   a_WFacePoint_L,
    const State&                   a_WFacePoint_R,
    const BoxData<T,1,MEM>&        a_thetaFluxFace,
    const BoxData<T,1,MEM>&        a_adjDdir,
    const BoxData<T,DIM,MEM,DIM>&  a_deta_a,
    const BoxData<T,1,MEM>&        a_deta,
    const BoxData<T,DIM,MEM,DIM>&  a_A,
    const BoxData<T,DIM,MEM>&      a_ntrow,
    uint                           a_block,
    uint                           a_dir) const
{
  PR_TIMERS("BoxOp_EulerCubedSphere::BertaPrimPointToMappedRusanovFlux()");

  Stencil<T> deconvolveStencilRadial =
      (1.0*Shift(Point::Zeros())
       - (1.0/24.0)*Stencil<T>::Derivative(2,0,2));

  BoxData<T,1,MEM> adjDdirPoint;
  if (a_dir != 0)
  {
    adjDdirPoint = deconvolveStencilRadial(a_adjDdir);
  }
  else
  {
    adjDdirPoint.define(a_adjDdir.box());
    a_adjDdir.copyTo(adjDdirPoint);
  }

  BoxData<T,1,MEM> kePoint_L(a_WFacePoint_L.box());
  BoxData<T,1,MEM> mePoint_L(a_WFacePoint_L.box());
  EnergyW(kePoint_L,mePoint_L,a_WFacePoint_L,a_block,a_dir);

  State qPoint_L =
      forall<T,NUMCOMPS>
      ([ ] PROTO_LAMBDA
       (
        Var<T,NUMCOMPS>& a_u,
        const Var<T,NUMCOMPS>& a_wface,
        const Var<T,1>& a_ke,
        const Var<T,1>& a_me,
        const Var<T>& a_adjddir,
        const T a_gamma)
       {
          a_u(iRHO ) = a_adjddir(0)*a_wface(iRHO);
          a_u(iMOMX) = a_adjddir(0)*a_wface(iRHO)*a_wface(iVX);
          a_u(iMOMY) = a_adjddir(0)*a_wface(iRHO)*a_wface(iVY);
          a_u(iMOMZ) = a_adjddir(0)*a_wface(iRHO)*a_wface(iVZ);
          a_u(iE   ) = a_adjddir(0)*(a_ke(0) + a_me(0)
                          + a_wface(iP)/(a_gamma - 1.0));
          a_u(iBX  ) = a_adjddir(0)*a_wface(iBX);
          a_u(iBY  ) = a_adjddir(0)*a_wface(iBY);
          a_u(iBZ  ) = a_adjddir(0)*a_wface(iBZ);
       },
       a_WFacePoint_L,
       kePoint_L,
       mePoint_L,
       adjDdirPoint,
       m_gamma);

  BoxData<T,1,MEM> kePoint_R(a_WFacePoint_R.box());
  BoxData<T,1,MEM> mePoint_R(a_WFacePoint_R.box());
  EnergyW(kePoint_R,mePoint_R,a_WFacePoint_R,a_block,a_dir);

  State qPoint_R =
      forall<T,NUMCOMPS>
      ([ ] PROTO_LAMBDA
       (
        Var<T,NUMCOMPS>& a_u,
        const Var<T,NUMCOMPS>& a_wface,
        const Var<T,1>& a_ke,
        const Var<T,1>& a_me,
        const Var<T>& a_adjddir,
        const T a_gamma)
       {
          a_u(iRHO ) = a_adjddir(0)*a_wface(iRHO);
          a_u(iMOMX) = a_adjddir(0)*a_wface(iRHO)*a_wface(iVX);
          a_u(iMOMY) = a_adjddir(0)*a_wface(iRHO)*a_wface(iVY);
          a_u(iMOMZ) = a_adjddir(0)*a_wface(iRHO)*a_wface(iVZ);
          a_u(iE   ) = a_adjddir(0)*(a_ke(0) + a_me(0)
                          + a_wface(iP)/(a_gamma - 1.0));
          a_u(iBX  ) = a_adjddir(0)*a_wface(iBX);
          a_u(iBY  ) = a_adjddir(0)*a_wface(iBY);
          a_u(iBZ  ) = a_adjddir(0)*a_wface(iBZ);
       },
       a_WFacePoint_R,
       kePoint_R,
       mePoint_R,
       adjDdirPoint,
       m_gamma);

  State lapQ_L = Stencil<T>::LaplacianFace(a_dir)(qPoint_L);
  State lapQ_R = Stencil<T>::LaplacianFace(a_dir)(qPoint_R);

  State qAvg_L(qPoint_L.box());
  State qAvg_R(qPoint_R.box());

  forallInPlace(
      f_flux_point_to_avg_theta,
      qAvg_L,
      qPoint_L,
      lapQ_L,
      a_thetaFluxFace);

  forallInPlace(
      f_flux_point_to_avg_theta,
      qAvg_R,
      qPoint_R,
      lapQ_R,
      a_thetaFluxFace);

  qAvg_R -= qAvg_L;

  T one = 1.0;

  BoxData<T,DIM,MEM> fluxdirMom =
      slice<T,NUMCOMPS,DIM,MEM>(a_fluxdir,iMOMX);

  BoxData<T,DIM,MEM> Con_diff_Mom =
      slice<T,NUMCOMPS,DIM,MEM>(qAvg_R,iMOMX);

  forallInPlace(f_matVecProd,fluxdirMom,a_deta_a,Con_diff_Mom,one);

  BoxData<T,DIM,MEM> fluxdirB =
      slice<T,NUMCOMPS,DIM,MEM>(a_fluxdir,iBX);

  BoxData<T,DIM,MEM> Con_diff_B =
      slice<T,NUMCOMPS,DIM,MEM>(qAvg_R,iBX);

  forallInPlace(f_matVecProd,fluxdirB,a_deta_a,Con_diff_B,one);

  forallInPlace(
      [ ] PROTO_LAMBDA
      (
        Var<T,NUMCOMPS,MEM>&       a_Rusflux,
        const Var<T,NUMCOMPS,MEM>& a_fluxin,
        const Var<T,1,MEM>&        a_deta)
      {
        a_Rusflux(iRHO) = a_deta(0)*a_fluxin(iRHO);
        a_Rusflux(iE  ) = a_deta(0)*a_fluxin(iE);
      },
      a_fluxdir,
      qAvg_R,
      a_deta);

  forallInPlace(
      [ ] PROTO_LAMBDA
      (
        Var<T,NUMCOMPS,MEM>& a_Rusflux,
        const Var<T,NUMCOMPS>& a_wface_L,
        const Var<T,NUMCOMPS>& a_wface_R,
        const Var<T,DIM,MEM,DIM>& a_Amatrix,
        const Var<T,1,MEM>& a_detA,
        const Var<T,DIM,MEM>& a_adjArow,
        const T a_gamma,
        const int a_d)
      {
        double rho = 0.5*(a_wface_L(iRHO) + a_wface_R(iRHO));
        double w_d = 0.5*(a_wface_L(iVX+a_d) + a_wface_R(iVX+a_d));

        double B_x = 0.5*(a_wface_L(iBX) + a_wface_R(iBX));
        double B_y = 0.5*(a_wface_L(iBY) + a_wface_R(iBY));
        double B_z = 0.5*(a_wface_L(iBZ) + a_wface_R(iBZ));

        double p = 0.5*(a_wface_L(iP) + a_wface_R(iP));
        if (p < 0.0) p = 0.0;

        rho = max(rho,RHO00);

        vecType bSph;
        vecType bCart;

        bSph[0] = B_x;
        bSph[1] = B_y;
        bSph[2] = B_z;

        bCart[0] = 0.0;
        bCart[1] = 0.0;
        bCart[2] = 0.0;

        double B_n   = 0.0;
        double B_mag = 0.0;

        for (int dir0 = 0; dir0 < DIM; dir0++)
        {
          for (int dir1 = 0; dir1 < DIM; dir1++)
          {
            bCart[dir0] += a_Amatrix(dir0,dir1)*bSph[dir1];
          }

          B_n   += a_Amatrix(dir0,0)*bCart[dir0];
          B_mag += bCart[dir0]*bCart[dir0];
        }

        B_mag = sqrt(B_mag);

        double ce = sqrt(a_gamma*p/rho);
        double aa = ce*ce + B_mag*B_mag/(4.0*M_PI*rho);
        double bb = abs(B_n)*ce/sqrt(M_PI*rho);

        double af =
            0.5*(sqrt(max(0.0,aa + bb)) + sqrt(max(0.0,aa - bb)));

        double A_d_mag =
            sqrt(a_Amatrix(0,a_d)*a_Amatrix(0,a_d)
               + a_Amatrix(1,a_d)*a_Amatrix(1,a_d)
               + a_Amatrix(2,a_d)*a_Amatrix(2,a_d));

        double speed = abs(w_d) + af*A_d_mag;

        for (int comp = 0; comp < NUMCOMPS; comp++)
        {
          a_Rusflux(comp) *= speed;
        }
      },
      a_fluxdir,
      a_WFacePoint_L,
      a_WFacePoint_R,
      a_A,
      a_deta,
      a_ntrow,
      m_gamma,
      a_dir);
}


// ====================================================================================
// 8. Main Berta-WENOZ operator.
// ====================================================================================

template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::BertaOperator
  (
   State&                          a_Rhs,
   Array<State, DIM>&              a_fluxes,
   const State&                    a_USph,
   T                               a_dt,
   T                               a_scale,
   int                             a_counter) const
{
  PR_TIMERS("BoxOp_EulerCubedSphere::operator()::BertaWENOZ");

  int a_block = this->block();

  Box bxRad = a_USph.box().grow(Point::Basis(0,2));

  BoxData<T,DIM,MEM> Dr(bxRad);
  BoxData<T,1,MEM>   radius(bxRad);
  BoxData<T,DIM,MEM> adjDr(bxRad);
  BoxData<T,1,MEM>   dVolr(bxRad);

  radialMetrics(radius,Dr,adjDr,dVolr,bxRad);

  Box bx = a_USph.box();

  // Detector/fallback settings.
  // Increase eta_d for less fallback. Decrease eta_d for more fallback.
  const T eta_d = ParseInputs::get_Berta_eta_d();
  const int detector_type = ParseInputs::get_Berta_detector_type(); // 0 = Jameson, 1 = derivative ratio.
  const int fallback_type = ParseInputs::get_Berta_fallback_type(); // 0 = WENO3, 1 = PLM.

  // Detector uses V(<U>), not V(UPoint).
  State WAvgForDetector;
  sphAvgToPrim(WAvgForDetector,a_USph,a_block);

  BoxData<T,1,MEM> thetaCell;
  if (detector_type == 1)
  {
    BertaComputeThetaDerivativeRatio(thetaCell,WAvgForDetector,eta_d);
  }
  else
  {
    BertaComputeThetaJameson(thetaCell,WAvgForDetector,eta_d);
  }

  // Berta average-to-point recovery.
  State UPoint;
  BertaAverageToPoint(UPoint,a_USph,thetaCell);

  // Pointwise primitive variables.
  State WPoint;
  sphPointToPrim(WPoint,UPoint,a_block);

  T dxiPerp = this->dx()[1];
  T dxiR    = this->dx()[0];

  T half = 0.5;

  BoxData<T,1,MEM>       deta0(bx);
  BoxData<T,DIM,MEM>     ntrow0(bx);
  BoxData<T,DIM,MEM,DIM> A_matrix0(bx);

  BoxData<T,1,MEM> minusJDivB(a_Rhs.box());
  minusJDivB.setVal(0.0);

  for (int dir = 0; dir < DIM; dir++)
  {
    PR_TIME("BoxOp_EulerCubedSphere::operator()::Berta directional flux");

    Point e = Point::Basis(dir);

    // Pointwise WENOZ reconstruction with WENO3/PLM fallback.
    State WFace_L;
    State WFace_R;
    BoxData<T,1,MEM> thetaFluxFace;
    BoxData<T,1,MEM> thetaReconFace;

    BertaReconstructWithFallback(
        WFace_L,
        WFace_R,
        thetaFluxFace,
        thetaReconFace,
        WPoint,
        thetaCell,
        dir,
        fallback_type);

    Box bxface = WFace_L.box();

    // Mapping variables.
    BoxData<T,1,MEM> adjDrDir = slice(adjDr,dir);

    BoxData<T,DIM,MEM,DIM> deta_a =
        CubedSphereShell::detA_A<T,MEM>(bxface,dxiPerp,a_block,dir);

    BoxData<T,1,MEM> deta =
        CubedSphereShell::detA<T,MEM>(bxface,dxiPerp,dir);

    if (dir == 0) deta.copyTo(deta0);

    Box bxCell = WPoint.box();
    Box bxdir  = bxCell.extrude(Point::Basis(dir));

    BoxData<T,DIM,MEM> ntrow =
        CubedSphereShell::adjA<T,MEM>(bxdir,dxiPerp,a_block,dir);

    if (dir == 0) ntrow.copyTo(ntrow0);

    Array<Array<uint,DIM>,6> permute =
        {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};

    Array<Array<int,DIM>,6> sign =
        {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}};

    Point high = bxCell.high();
    Point low  = bxCell.low();
    high[0] = low[0];

    Box bx0(low,high);

    BoxData<T,DIM,MEM,DIM> A_matrix(bxCell);

    T offseta = half;
    T offsetb = half;

    if (dir == 1) offseta = 0.0;
    if (dir == 2) offsetb = 0.0;

    forallInPlace_p(
        f_Amatrix,
        bx0,
        A_matrix,
        permute[a_block],
        sign[a_block],
        dxiPerp,
        offseta,
        offsetb);

    spreadSlice(A_matrix);

    if (dir == 0) A_matrix.copyTo(A_matrix0);

    // Powell term. Use centered reconstructed face state.
    State WFace_C(WFace_L.box());
    WFace_L.copyTo(WFace_C);
    WFace_C += WFace_R;
    WFace_C *= 0.5;

    DivBPowell1(minusJDivB,WFace_C,deta,adjDrDir,dir,bxface);

    // Physical fluxes:
    // pointwise primitive -> pointwise flux -> theta-aware face-average flux.
    // Use thetaFluxFace, not thetaReconFace, for flux averaging.
    BoxData<T,NUMCOMPS> physFlux_L;
    BoxData<T,NUMCOMPS> physFlux_R;

    BertaPrimPointToFluxAverage(
        physFlux_L,
        WFace_L,
        thetaFluxFace,
        adjDrDir,
        a_block,
        dir);

    BertaPrimPointToFluxAverage(
        physFlux_R,
        WFace_R,
        thetaFluxFace,
        adjDrDir,
        a_block,
        dir);

    // Mapped fluxes and Rusanov dissipation.
    Box bxface0 = a_fluxes[dir].box();

    State fluxdir_L(bxface);
    State fluxdir_R(bxface);
    State fluxdir_Rus(bxface);
    State fluxdir(bxface0);
    fluxdir_L.setVal(0.0);
    fluxdir_R.setVal(0.0);
    fluxdir_Rus.setVal(0.0);
    fluxdir.setVal(0.0);

    BertaPrimPointToMappedFlux(
        fluxdir_L,
        physFlux_L,
        WFace_L,
        thetaFluxFace,
        adjDrDir,
        deta_a,
        deta,
        ntrow,
        a_block,
        dir);

    BertaPrimPointToMappedFlux(
        fluxdir_R,
        physFlux_R,
        WFace_R,
        thetaFluxFace,
        adjDrDir,
        deta_a,
        deta,
        ntrow,
        a_block,
        dir);

    BertaPrimPointToMappedRusanovFlux(
        fluxdir_Rus,
        WFace_L,
        WFace_R,
        thetaFluxFace,
        adjDrDir,
        deta_a,
        deta,
        A_matrix,
        ntrow,
        a_block,
        dir);

    forallInPlace(
      [ ] PROTO_LAMBDA
      (
        Var<T,NUMCOMPS,MEM>&       a_fluxout,
        const Var<T,NUMCOMPS,MEM>& a_fluxin_L,
        const Var<T,NUMCOMPS,MEM>& a_fluxin_R,
        const Var<T,NUMCOMPS,MEM>& a_fluxin_Rus)
      {
        for (int comp = 0; comp < NUMCOMPS; comp++)
        {
          a_fluxout(comp) =
              0.5*(a_fluxin_L(comp) + a_fluxin_R(comp) - a_fluxin_Rus(comp));
        }
      },
      fluxdir,
      fluxdir_L,
      fluxdir_R,
      fluxdir_Rus);

    fluxdir.copyTo(a_fluxes[dir]);

    T dxInv = -a_scale/(this->dx()[dir]);

    BoxData<T,NUMCOMPS,MEM> rhsIncr;
    rhsIncr = Stencil<T>::FluxDivergence(dir)(fluxdir,dxInv);

    forallInPlace(
      [ ] PROTO_LAMBDA
      (
        Var<T,NUMCOMPS,MEM>&       a_rhs,
        const Var<T,NUMCOMPS,MEM>& a_rhsIncr)
      {
        for (int comp = 0; comp < NUMCOMPS; comp++)
        {
          a_rhs(comp) += a_rhsIncr(comp);
        }
      },
      a_Rhs,
      rhsIncr);
  }

  if (ParseInputs::get_divergence_cleaning() == 1)
  {
    DivBPowell2(a_Rhs,minusJDivB,WPoint,A_matrix0);
  }
}

// ====================================================================================
// Berta-style 4th-order FV pointwise-WENOZ implementation ends here.
// ====================================================================================
