#include "MHDReader.H"
#include "MHD_Constants.H"
#include <memory>  
extern MHDReader BC_global;
template <typename T, MemType MEM>
PROTO_KERNEL_START
void f_radialoutflowIBCs_F(
    Point a_pt,
    Var<T, NUMCOMPS, MEM> &a_USph,
    Var<T> &a_radius,
    Var<T, DIM> &a_X_cart,
    T a_gamma,
    int a_nradius,
    T a_time,
    T a_strength,
    double r_BC)
{
  T r = sqrt(a_X_cart(0) * a_X_cart(0) + a_X_cart(1) * a_X_cart(1) + a_X_cart(2) * a_X_cart(2));
  T theta = atan2(sqrt(a_X_cart(0) * a_X_cart(0) + a_X_cart(1) * a_X_cart(1)), a_X_cart(2)); // polar angle
  T phi = atan2(a_X_cart(1), a_X_cart(0)); // azimuthal angle
  // Compute spherical BCs.
  T p0 = 1.0e-7; // p near 21.5 c_SR is about 1e-7 dyne/cm2
  T rho0 = 700*c_MP; // rho at 21.5 c_SR is about 700/cm3
  T u0 = 500.0*1e5; // v at 21.5 c_SR is about 500 km/s
  T rho = rho0*pow(r_BC/r,2.0);
  T p = p0*pow(r_BC/r,2.0*a_gamma);
  double carr_rot_time = 25.38 * 24 * 60 * 60; // Seconds
  T omega = 2.0 * M_PI / carr_rot_time; 
  T phit  = phi - omega * a_time; // or + for opposite rotation
  const T perturbationScale = a_strength/100.0;
  const T sinTheta = sin(theta);
  // Match the smooth-on-the-sphere azimuthal profile used at initialization.
  const T angularPattern =
      pow(sinTheta,5)*sin(5*phit) + pow(sinTheta,10)*sin(10*phit);
  T ur = u0 - perturbationScale*u0/3
      * angularPattern;
  a_USph(iRHO) = rho;
  a_USph(iMOMX) = ur * rho;
  a_USph(iMOMY) = 0.0;
  a_USph(iMOMZ) = 0.0;
  a_USph(iBX) = 0.0;
  a_USph(iBY) = 0.0;
  a_USph(iBZ) = 0.0;
  T umag = sqrt(a_USph(iMOMX) * a_USph(iMOMX) + a_USph(iMOMY) * a_USph(iMOMY) + a_USph(iMOMZ) * a_USph(iMOMZ)) / a_USph(iRHO);
  T Bmag = sqrt(a_USph(iBX) * a_USph(iBX) + a_USph(iBY) * a_USph(iBY) + a_USph(iBZ) * a_USph(iBZ));
  a_USph(iE) = p / (a_gamma - 1.0) + 0.5 * a_USph(iRHO) * umag * umag + Bmag * Bmag/8.0/M_PI;
}
PROTO_KERNEL_END(f_radialoutflowIBCs_F, f_radialoutflowIBCs)

// Radial boundary conditions here.
template<
  typename T,
  typename MAP,
  MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::
PreStagePatch(BoxData<T,NUMCOMPS,MEM>& a_USphStage,
              BoxData<T,NUMCOMPS,MEM>& a_BC,
              BoxData<T,1,MEM>& a_dVolr,
              Box a_blockBox,
              T a_time)
{
  PR_TIMERS("BoxOp_EulerCubedSphere::PreStagePatch()");
  HDF5Handler h5;

  // if (procID() == 0) h5.writePatch(1,a_USphStage,"a_USphStage_noBC");

  int a_block = this->block();
  Box bxGhosted = a_USphStage.box();
  BoxData<double> radius(a_dVolr.box());
  BoxData<double, DIM, HOST> Dr(a_dVolr.box());
  BoxData<double, DIM, HOST> adjDr(a_dVolr.box());
  radialMetrics(radius, Dr, adjDr, a_dVolr, Dr.box());
  Array<double, DIM> dx;
  dx = this->dx();
  double dxradius = dx[0];
  double a_h = dx[2];
  int thickness = std::round(1.0 / dxradius);
  // Set radial boundary condtions in radial ghost cells.
  int init_condition_type = ParseInputs::get_init_condition_type();
  double innerRadius = ParseInputs::get_inner_radius();
  Point low = bxGhosted.low();
  low[0] = a_blockBox.high()[0] + 1;
  Box bdryBoxHigh(low, bxGhosted.high());
  Point high = bxGhosted.high();
  high[0] = a_blockBox.low()[0] - 1;
  Box bdryBoxLow(bxGhosted.low(), high);

  // Outer BCs
  if (a_blockBox.high()[0] < bxGhosted.high()[0])
  {
    if (init_condition_type == 2
        || init_condition_type == 4
        || init_condition_type == 6
        || init_condition_type == 8){
      forallInPlace_p(
          f_testProblemHydroBoundary,
          bdryBoxHigh,
          a_USphStage,
          m_gamma);
    } 
    if (init_condition_type == 0 || init_condition_type == 1)
      {
        // Copy last layer to ghost cells and scale by radius.
        // Extrapolation can make -ve P in ghost cells.
        Point source_low = Point(a_blockBox.high()[0],
                                 bxGhosted.low()[1],bxGhosted.low()[2]);
        Point source_high = Point(a_blockBox.high()[0],bxGhosted.high()[1],
                                  bxGhosted.high()[2]);
        Box sourceBox(source_low,source_high);
        State U_ghost(bdryBoxHigh);
        Scalar radius_ghost(bdryBoxHigh);
        Stencil<double> m_exterp_f_2nd;
        for (int i = 1; i <= NGHOST; i++ ) {
          a_USphStage.copyTo(U_ghost,sourceBox,Point::Basis(0)*(i));
          // Using shifting option of copyTo
          radius.copyTo(radius_ghost,sourceBox,Point::Basis(0)*(i));
        }
        U_ghost.copyTo(a_USphStage,bdryBoxHigh);
        //scale by radius.
        forallInPlace_p([ ] PROTO_LAMBDA
                        (Point a_pt,
                         Var<T, NUMCOMPS, MEM> &a_Usph,
                         Var<T> &a_radius,
                         Var<T> &a_radius_ghost)
                        {  
                          double rad = a_radius(0);
                          double rad_last = a_radius_ghost(0);
                          a_Usph(0) *= pow(rad_last/rad,2.0);
                          a_Usph(1) *= pow(rad_last/rad,2.0);
                          a_Usph(2) *= pow(rad_last/rad,2.0);
                          a_Usph(3) *= pow(rad_last/rad,2.0);
                          a_Usph(4) *= pow(rad_last/rad,2.0);
                          a_Usph(5) *= pow(rad_last/rad,2.0);
                          a_Usph(6) *= pow(rad_last/rad,1.0); 
                          a_Usph(7) *= pow(rad_last/rad,1.0);
                        },bdryBoxHigh,a_USphStage,radius,radius_ghost);
    }

    if (isConstantMagneticFieldTest(init_condition_type))
    {
      // Preserve the cell-averaged multiblock state at radial/angular
      // corners. Reconstructing an analytic point value here is inconsistent
      // with the fourth-order averages in a_USphStage and seeds visible block
      // seams in low-beta magnetic tests.
      Point sourceLow = bxGhosted.low();
      Point sourceHigh = bxGhosted.high();
      sourceLow[0] = a_blockBox.high()[0];
      sourceHigh[0] = a_blockBox.high()[0];
      const Box sourceBox(sourceLow, sourceHigh);
      for (int ii = 1; ii <= NGHOST; ++ii)
      {
        a_USphStage.copyTo(
            a_USphStage, sourceBox, Point::Basis(0)*ii);
      }
    }
  }
  // Inner BCs
  if (a_blockBox.low()[0] > bxGhosted.low()[0])
  {

    // Find X_cart
    double half = 0.5;
    BoxData<double, DIM, HOST> XCart;
    if (init_condition_type == 1)
    {
      // Only the inner radial ghost cells consume these coordinates.
      XCart = forall_p<double,DIM,HOST>
        (f_cubedSphereMap3,bdryBoxLow,radius,dx,half,half,a_block);
    }

    if (init_condition_type == 2
        || init_condition_type == 4
        || init_condition_type == 6
        || init_condition_type == 8){
      forallInPlace_p(
          f_testProblemHydroBoundary,
          bdryBoxLow,
          a_USphStage,
          m_gamma);
    }
    if (init_condition_type == 1){
      T strength = ParseInputs::get_strength();
      forallInPlace_p(
          f_radialoutflowIBCs,
          bdryBoxLow,
          a_USphStage,
          radius,
          XCart,
          m_gamma,
          thickness,
          a_time,
          strength,
          innerRadius);
    }
    if (init_condition_type == 0)
    {
      // The reader supplies native conserved cell averages, including the
      // radial/angular corner ghosts required by reconstruction.
      PROTO_ASSERT(a_BC.box().containsBox(bdryBoxLow),
          "Inner boundary averages do not cover the radial ghost region");
      a_BC.copyTo(a_USphStage,bdryBoxLow);
    }

    if (isConstantMagneticFieldTest(init_condition_type))
    {
      Point sourceLow = bxGhosted.low();
      Point sourceHigh = bxGhosted.high();
      sourceLow[0] = a_blockBox.low()[0];
      sourceHigh[0] = a_blockBox.low()[0];
      const Box sourceBox(sourceLow, sourceHigh);
      for (int ii = 1; ii <= NGHOST; ++ii)
      {
        a_USphStage.copyTo(
            a_USphStage, sourceBox, Point::Basis(0)*(-ii));
      }
    }
  }
  // if (procID() == 0) h5.writePatch(1,a_USphStage,"a_USphStage_withBC");
}


/// Exchange, block boundary conditions, transform from JU to USph. Set state variables in boundary registers.
template<
    typename T,
    typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::
preStageLevel(MBLevelBoxData<T,NUMCOMPS,MEM>& a_UStage,
              MBLevelBoxData<T,NUMCOMPS,MEM>& a_k,
              MBLevelBoxData<T,1,MEM>& a_dVolrLev,
              MBInterpOp& a_iop,
              MBBoundaryRegister<T,NUMREG,MEM>& a_blockreg,
              Array<Reduction<T,Operation::Sum,MEM>,NUMCOMPS>& a_consSum,
              T a_time,
              T a_dt,
              T a_weight,
              int a_kstage,
              bool a_flag_viscosity)
{
  
  PR_TIMERS("BoxOp_EulerCubedSphere::PreStageLevel()");
  int init_condition_type = ParseInputs::get_init_condition_type();
  double gamma = ParseInputs::get_gamma();
  auto layout = a_UStage.layout();

  auto& geometry = levelGeometry(layout);
  auto& map = *geometry.map;
  auto& eulerOp = *geometry.op;

  // Transform JU -> USph.
  HelioCubed::consToSphInterpEuler(a_UStage,a_iop,a_dVolrLev, 4, 0);

  // initialize a_k to zero. 
  a_k.setVal(0.);
  MBLevelBoxData<double, 8, HOST> dstData;
  BoxData<T,NUMCOMPS,MEM> unusedBoundary;
  if (init_condition_type == 0)
  {
    // The BC interpolation writes patch storage directly; this temporary is
    // never exchanged and needs no interblock boundary buffers.
    Array<Point,DIM+1> bcGhost;
    bcGhost.fill(Point::Zeros());
    bcGhost[0] = NGHOST*Point::Ones();
    dstData.define(layout, bcGhost);
    BC_global.BoxData_to_BC<T,MEM>(dstData, map, a_time, true);
  }

  for (auto dit : a_UStage.layout())
  {
    unsigned int block = layout.block(dit);
    Box blockBox = layout.getBlock(block).domain().box();
    auto &rhs_i = a_k[dit];
    auto &USph_i =a_UStage[dit];
    
    // Conversion requires the fixed volume metrics initialized before stepping.
    // Only halos crossing a physical radial boundary need boundary work.
    if (USph_i.box().low()[0] < blockBox.low()[0]
        || USph_i.box().high()[0] > blockBox.high()[0])
    {
      eulerOp[dit].PreStagePatch(USph_i,init_condition_type == 0 ? dstData[dit] : unusedBoundary,a_dVolrLev[dit],blockBox, a_time);
    }
    auto dx = eulerOp[dit].dx();    
    // Compute viscous operators here, and increment a_rhsTotal.
    const bool viscosityEnabled =
      (ParseInputs::get_linear_viscosity() == 1 ||
       ParseInputs::get_nonlinear_viscosity() == 1);
    if (a_flag_viscosity)
      {
        PR_TIMERS("artificial viscosity");
        unsigned int block = layout.block(dit);
        Box blockBox = layout.getBlock(block).domain().box();
        // auto mbdata = a_blockreg.bounds(dit);
        Array<BoxData<double, NUMCOMPS>, DIM> fluxes;
        fluxes[0].define(rhs_i.box().extrude(0), 0.0);
        fluxes[1].define(rhs_i.box().extrude(1), 0.0);
        fluxes[2].define(rhs_i.box().extrude(2), 0.0);
        if (viscosityEnabled) eulerOp[dit].viscFlux(rhs_i, fluxes, USph_i,a_dt);
        if (ParseInputs::get_Sun_gravity() == 1) eulerOp[dit].Gravity_Euler(rhs_i, USph_i);
        // viscFlux updates the state with +div(Fvisc). The register and
        // radial balance use the conservative -div(F) convention.
        a_blockreg.increment(fluxes, dit, -a_dt);
        CubedSphereShell::radialCons<T,NUMCOMPS,MEM>(fluxes[0],a_consSum,dx[0],blockBox,-a_dt);
      }
  }
}
