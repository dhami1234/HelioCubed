#include "MHD_Constants.H"

template<typename T, typename MAP, MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::Gravity(
    State& a_Rhs,
    const State& a_USph,
    const BoxData<T,1,MEM>& a_drDxi,
    const BoxData<T,DIM,MEM,DIM>& a_detA_A,
    const BoxData<T,1,MEM>& a_detA) const
{
  PR_TIME("BoxOp_EulerCubedSphere::Gravity()");
  if (ParseInputs::get_spatial_order() == 2)
  {
    forallInPlace([] PROTO_LAMBDA(
        Var<T,NUMCOMPS,MEM>& rhs,const Var<T,NUMCOMPS,MEM>& u,
        const Var<T,1,MEM>& dr,const Var<T,DIM,MEM,DIM>& mat,
        const Var<T,1,MEM>& det)
    {
      const T factor=-c_G*c_M0*dr(0);
      for(int d=0;d<DIM;++d) rhs(iMOMX+d)+=factor*mat(d,0)*u(iRHO);
      rhs(iE)+=factor*det(0)*u(iMOMX);
    },a_Rhs,a_USph,a_drDxi,a_detA_A,a_detA);
    return;
  }
  // J = det(A) r^2 dr/dxi and g = -GM/r^2 e_r. Cancel r^2
  // analytically before averaging. This also preserves the angular metric
  // factorization of a radial state. USph stores <rho> and <rho V_r>.
  const Box inputBox = a_Rhs.box().grow(2);
  PROTO_ASSERT(a_USph.box().containsBox(inputBox)
      && a_drDxi.box().containsBox(inputBox),
      "Gravity radial product lacks stencil coverage");
  const Box angularBox = a_Rhs.box().grow(Point(0,1,1));
  PROTO_ASSERT(a_detA_A.box().containsBox(angularBox)
      && a_detA.box().containsBox(angularBox),
      "Gravity angular product lacks stencil coverage");

  auto rho = slice(a_USph,iRHO);
  auto radialMomentum = slice(a_USph,iMOMX);
  auto radialDensity = cellProduct(a_drDxi,rho);
  auto radialEnergyWork = cellProduct(a_drDxi,radialMomentum);
  const T minusGM = -c_G*c_M0;

  for (int dir = 0; dir < DIM; ++dir)
  {
    // Column zero is the averaged det(A) times the radial unit vector.
    auto radialDirection = slice(a_detA_A,dir,0);
    auto source = faceProduct(radialDirection,radialDensity,0);
    PROTO_ASSERT(source.box().containsBox(a_Rhs.box()),
        "Gravity momentum source does not cover the requested RHS");
    source *= minusGM;
    auto rhs = slice(a_Rhs,iMOMX+dir);
    rhs += source;
  }
  // The angular metric is independent of radius, so only angular
  // covariance is needed for this final product (face normal = radial).
  auto energySource = faceProduct(a_detA,radialEnergyWork,0);
  PROTO_ASSERT(energySource.box().containsBox(a_Rhs.box()),
      "Gravity energy source does not cover the requested RHS");
  energySource *= minusGM;
  auto energyRhs = slice(a_Rhs,iE);
  energyRhs += energySource;
}

template<typename T, typename MAP, MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::Gravity_Euler(
    State& a_Rhs, const State& a_USph) const
{
  PR_TIME("BoxOp_EulerCubedSphere::Gravity_Euler()");
  const Box radialBox = a_Rhs.box().grow(2);
  BoxData<T,DIM,MEM> Dr(radialBox), adjDr(radialBox);
  BoxData<T,1,MEM> radius(radialBox), dVolr(radialBox);
  radialMetrics(radius,Dr,adjDr,dVolr,radialBox);
  auto drDxi = slice(Dr,0);

  const Box angularBox = a_Rhs.box().grow(Point(0,1,1));
  const T h = this->dx()[1];
  auto detA = CubedSphereShell::detA<T,MEM>(angularBox,h,0);
  auto detA_A = CubedSphereShell::detA_A<T,MEM>(
      angularBox,h,this->block(),0);
  Gravity(a_Rhs,a_USph,drDxi,detA_A,detA);
}
