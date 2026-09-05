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
    T a_strength)
{
  T r = sqrt(a_X_cart(0) * a_X_cart(0) + a_X_cart(1) * a_X_cart(1) + a_X_cart(2) * a_X_cart(2));
  T theta = atan2(sqrt(a_X_cart(0) * a_X_cart(0) + a_X_cart(1) * a_X_cart(1)), a_X_cart(2)); // polar angle
  T phi = atan2(a_X_cart(1), a_X_cart(0)); // azimuthal angle
  // Compute spherical initial data.
  double r_BC = CUBED_SPHERE_SHELL_R0; //in cm
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
    template<MemType> typename MAP,
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
  T testBx = c_TEST_BX0;
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
        a_gamma);
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
        strength);
  }

  if (init_condition_type == 0) {
    Point source_lo = Point(a_BC.box().low()[0],a_BC.box().low()[1],a_BC.box().low()[2]);
    Point source_hi = Point(a_BC.box().low()[0],a_BC.box().high()[1],a_BC.box().high()[2]);
    Box sourceBox(source_lo,source_hi);
    int low_index = a_W.box().low()[0];
    int high_index = a_W.box().high()[0];
    for (int i = -NGHOST; i <= (high_index-low_index)-NGHOST; i++){
      a_BC.copyTo(a_W,sourceBox,Point::Basis(0)*(i+1));
    }

    // Convert from classical spherical to Cartesian (confirmed correct)
    forallInPlace([ ] PROTO_LAMBDA
    (Var<double, NUMCOMPS, HOST> &a_W,
    Var<double, DIM> &a_X_cart)
    {  
      double x = a_X_cart(0);
      double y = a_X_cart(1);
      double z = a_X_cart(2);
      double r = sqrt(x*x + y*y + z*z);
      double theta = acos(z/r);
      double phi = atan2(y,x);
      // Convert spherical to Cartesian velocities
      double u_r = a_W(iVX);
      double u_theta = a_W(iVZ); //BC files are written as r,phi,theta
      double u_phi = a_W(iVY);
      a_W(iVX) = sin(theta)*cos(phi)*u_r + cos(theta)*cos(phi)*u_theta - sin(phi)*u_phi;
      a_W(iVY) = sin(theta)*sin(phi)*u_r + cos(theta)*sin(phi)*u_theta + cos(phi)*u_phi;
      a_W(iVZ) = cos(theta)*u_r - sin(theta)*u_theta;
      // Convert spherical to Cartesian magnetic fields
      double B_r = a_W(iBX);
      double B_theta = a_W(iBZ); //BC files are written as r,phi,theta
      double B_phi = a_W(iBY);
      a_W(iBX) = sin(theta)*cos(phi)*B_r + cos(theta)*cos(phi)*B_theta - sin(phi)*B_phi;
      a_W(iBY) = sin(theta)*sin(phi)*B_r + cos(theta)*sin(phi)*B_theta + cos(phi)*B_phi;
      a_W(iBZ) = cos(theta)*B_r - sin(theta)*B_theta;
    },a_W,a_Xcart);

    // Convert from cartesian to cubed-spherical Primitive variables w=A_inv*u
    BoxData<T,DIM,MEM> V_cart = slice<T,NUMCOMPS,DIM,MEM>(a_W,iVX);
    BoxData<T,DIM,MEM> tmpV(a_W.box());
    T one = 1.0;  
    forallInPlace(f_matVecProd,tmpV,A_matrix_inv,V_cart,one);
    tmpV.copyTo(V_cart);

    BoxData<T,DIM,MEM> B_cart = slice<T,NUMCOMPS,DIM,MEM>(a_W,iBX);
    BoxData<T,DIM,MEM> tmpB(a_W.box());
    forallInPlace(f_matVecProd,tmpB,A_matrix_inv,B_cart,one);
    tmpB.copyTo(B_cart);

    // In primitive spherical. To correct units. scale by radius.
    forallInPlace([ ] PROTO_LAMBDA
    (Var<double, NUMCOMPS, HOST> &a_WSph,
    Var<double> &a_radius,
    double a_gamma)
    {  
      double r_BC = CUBED_SPHERE_SHELL_R0; //in cm
      double rad = a_radius(0);
      double rho = a_WSph(0)*c_MP*pow(r_BC/rad,2.0);
      double u = a_WSph(1);
      double v = a_WSph(2);
      double w = a_WSph(3);
      double p = a_WSph(4)*1.0e-12*pow(r_BC/rad,2.0*a_gamma); // From picodyne to dyne 
      double Bx = a_WSph(5)*1.0e-6*pow(r_BC/rad,2.0); // MicroGauss to Gauss
      double By = a_WSph(6)*1.0e-6*pow(r_BC/rad,1.0);
      double Bz = a_WSph(7)*1.0e-6*pow(r_BC/rad,1.0);
      a_WSph(0) = rho; //rho
      a_WSph(1) = u; //V-x
      a_WSph(2) = v; //V-y
      a_WSph(3) = w; //V-z
      a_WSph(4) = p; //Pressure
      a_WSph(5) = Bx; //Bx
      a_WSph(6) = By; //By
      a_WSph(7) = Bz; //Bz
    },a_W,a_radius,a_gamma);
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
