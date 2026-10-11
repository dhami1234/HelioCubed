#include "MHD_Constants.H"

template <typename T, MemType MEM>
PROTO_KERNEL_START
void f_radialoutflowInit_F(
    Point a_pt,
    Var<T, NUMCOMPS, MEM> &a_W,
    Var<T> &a_radius,
    Var<T, DIM> &a_X_cart,
    T a_gamma,
    int a_nradius,
    T a_strength,
    double r_BC)
{
  T r = sqrt(a_X_cart(0) * a_X_cart(0) + a_X_cart(1) * a_X_cart(1) + a_X_cart(2) * a_X_cart(2));
  T theta = atan2(sqrt(a_X_cart(0) * a_X_cart(0) + a_X_cart(1) * a_X_cart(1)), a_X_cart(2)); // polar angle
  T phi = atan2(a_X_cart(1), a_X_cart(0)); // azimuthal angle
  // Compute spherical initial data.
  T p0 = 1.0e-7; // p near 21.5 c_SR is about 1e-7 dyne/cm2
  T rho0 = 700*c_MP; // rho at 21.5 c_SR is about 700/cm3
  T u0 = 500.0*1e5; // v at 21.5 c_SR is about 500 km/s
  T rho = rho0*pow(r_BC/r,2.0);
  T p = p0 * pow(r_BC/r,2.0*a_gamma);
  // T ur = u0;
  // T ur = u0 - u0/10 * sin(5.0*phi)*sin(theta);
  const T perturbationScale = a_strength/100.0;
  const T sinTheta = sin(theta);
  // The original longitude-only profile is multi-valued at the poles.
  // Taper each azimuthal mode with sin(theta)^m, which leaves the
  // equatorial profile unchanged and makes it smooth on the sphere.
  const T angularPattern =
      pow(sinTheta,5)*sin(5*phi) + pow(sinTheta,10)*sin(10*phi);
  T ur = u0 - perturbationScale*u0/3
      * angularPattern;
  a_W(iRHO) = rho;
  a_W(iVX) = ur;
  a_W(iVY) = 0.0;
  a_W(iVZ) = 0.0;
  a_W(iP) = p;
  a_W(iBX) = 0.0;
  a_W(iBY) = 0.0;
  a_W(iBZ) = 0.0;
}

PROTO_KERNEL_END(f_radialoutflowInit_F, f_radialoutflowInit)

template<
    typename T,
    typename MAP,
    MemType MEM>
inline void BoxOp_EulerCubedSphere<T,MAP,MEM>::
initialize(State& a_W,
          const State& a_BC,
          const Scalar& a_radius,
          const Vector& a_Xcart,
          T a_gamma,
          int a_nradius,
          T a_h,
          uint a_block) const
{
  PR_TIME("BoxOp_EulerCubedSphere::initialize()");
  int init_condition_type = ParseInputs::get_init_condition_type();
  double innerRadius = ParseInputs::get_inner_radius();
  double outerRadius = ParseInputs::get_outer_radius();

  // Calculate A_matrix_inv
  BoxData<T,DIM,MEM,DIM> amat = CubedSphereShell::detA_A<T,MEM>(a_W.box(),a_h,a_block,0);
  BoxData<T,1,MEM> deta = CubedSphereShell::detA<T,MEM>(a_W.box(),a_h,0);
  forallInPlace([ ] PROTO_LAMBDA
    (Var<T, DIM, MEM, DIM> &a_amat,
    Var<T, 1, MEM> &a_deta)
    {
      double inv_deta = 1.0 / a_deta(0);
      for (int i = 0; i < DIM; i++){
        for (int j = 0; j < DIM; j++){
          a_amat(i,j) = a_amat(i,j) * inv_deta;
        }
      }
    },amat,deta);
  BoxData<T,DIM,MEM,DIM> A_matrix_inv(a_W.box());
  A_matrix_inv.setToZero();
  forallInPlace(f_matinv3by3,A_matrix_inv,amat);


  T strength = ParseInputs::get_strength();
  T testBx = testProblemBackgroundBx(init_condition_type);
  T testBy = 0.0;
  T testBz = 0.0;
  if (init_condition_type == 2 || init_condition_type == 3)
  {
    forallInPlace_p(f_smoothSphericalWaveInit, a_W, a_gamma, a_nradius);
    if (isConstantMagneticFieldTest(init_condition_type))
    {
      forallInPlace_p(
          f_setConstantCartesianBPrimitive,
          a_W,
          A_matrix_inv,
          testBx,
          testBy,
          testBz);
    }
  }
  if (init_condition_type == 4 || init_condition_type == 5)
  {
    forallInPlace_p(
        f_nonRadialSphericalWaveInit,
        a_W,
        a_Xcart,
        a_gamma,
        innerRadius,
        outerRadius);
    if (isConstantMagneticFieldTest(init_condition_type))
    {
      forallInPlace_p(
          f_setConstantCartesianBPrimitive,
          a_W,
          A_matrix_inv,
          testBx,
          testBy,
          testBz);
    }
  }
  if (init_condition_type == 8 || init_condition_type == 9)
  {
    forallInPlace_p(f_nonRadialBlastInit, a_W, a_Xcart,
                   innerRadius, outerRadius);
    if (isConstantMagneticFieldTest(init_condition_type))
    {
      forallInPlace_p(
          f_setConstantCartesianBPrimitive,
          a_W,
          A_matrix_inv,
          testBx,
          testBy,
          testBz);
    }
  }
  if (init_condition_type == 1)
  {
    forallInPlace_p(
        f_radialoutflowInit,
        a_W,
        a_radius,
        a_Xcart,
        a_gamma,
        a_nradius,
        strength,
        innerRadius);
  }

  if (init_condition_type == 0)
  {
    solarWindPrimitive(a_W,a_BC,a_radius,a_gamma,a_block);
  }

  if (init_condition_type == 6 || init_condition_type == 7)
  {
    forallInPlace_p(
        f_symmetricStrongShockInit, a_W, a_gamma, strength, a_nradius);
    if (isConstantMagneticFieldTest(init_condition_type))
    {
      forallInPlace_p(
          f_setConstantCartesianBPrimitive,
          a_W,
          A_matrix_inv,
          testBx,
          testBy,
          testBz);
    }
  }
  
}

// One pointwise conversion for initialization and the boundary ghost fill.
// The latter averages the resulting conserved state before it enters the
// finite-volume stencil. Never invert <det(A) A>/<det(A)> at a point.
template<typename T, typename MAP, MemType MEM>
inline void BoxOp_EulerCubedSphere<T,MAP,MEM>::solarWindPrimitive(
    State& out, const State& boundary, const Scalar& radius, T gamma, uint block) const
{
  const T innerRadius = ParseInputs::get_inner_radius();
  Point low = boundary.box().low(), high = boundary.box().high();
  high[0] = low[0];
  Point requiredLow = out.box().low(), requiredHigh = out.box().high();
  requiredLow[0] = requiredHigh[0] = low[0];
  PROTO_ASSERT(boundary.box().containsBox(Box(requiredLow,requiredHigh))
      && radius.box().containsBox(out.box()),
      "solarWindPrimitive | Insufficient boundary or radius coverage.");
  for (int r = out.box().low()[0]; r <= out.box().high()[0]; ++r)
    boundary.copyTo(out,Box(low,high),Point::Basis(0,r-low[0]));
  auto basis = pointwiseSphericalBasis(out.box(),block,-1);
  BoxData<T,DIM,MEM,DIM> inverse(out.box());
  forallInPlace(f_matinv3by3,inverse,basis);
  forallInPlace([] PROTO_LAMBDA(
      Var<T,NUMCOMPS,MEM>& w, const Var<T,DIM,MEM,DIM>& A,
      const Var<T,DIM,MEM,DIM>& invA, const Var<T,1,MEM>& rad, T r0, T gamma)
    {
      const double x = A(0,0), y = A(1,0), z = A(2,0);
      const double s = hypot(x,y), cp = s > 0 ? x/s : 1.0, sp = s > 0 ? y/s : 0.0;
      const double er[3] = {x,y,z}, et[3] = {z*cp,z*sp,-s}, ep[3] = {-sp,cp,0};
      const double scale = r0/rad(0);
      // File ordering is r, phi, theta. Cartesian rotation and the cubed
      // basis are evaluated at exactly the same cell-center location.
      for (int vector = 0; vector < 2; ++vector)
      {
        const int first = vector == 0 ? iVX : iBX;
        double cart[3] = {}, local[3] = {};
        for (int i = 0; i < 3; ++i)
          cart[i] = er[i]*w(first)+ep[i]*w(first+1)+et[i]*w(first+2);
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j) local[i] += invA(i,j)*cart[j];
        for (int i = 0; i < 3; ++i)
          w(first+i) = local[i]*(first == iBX ? 1.e-6*scale*(i == 0 ? scale : 1.0) : 1.0);
      }
      w(iRHO) *= c_MP*scale*scale;
      w(iP) *= 1.e-12*pow(scale,2.0*gamma);
    },out,basis,inverse,radius,innerRadius,gamma);
}
