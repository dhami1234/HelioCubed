#include "MHD_Constants.H"

template<typename T, typename MAP, MemType MEM>
inline void BoxOp_EulerCubedSphere<T,MAP,MEM>::gravitySecondOrderPrimitives(
    State& rhs,const State& w) const
{
  PR_TIMERS("HelioCubed::gravitySecondOrderPrimitives");
  const Box box=rhs.box();
  if constexpr (MEM == HOST)
  {
    if (!m_gravityGeometry || m_gravityGeometry->box != box)
    {
      PR_TIMERS("HelioCubed::gravityGeometrySetup");
      auto cache=std::make_shared<SecondOrderGravityGeometry>();
      cache->box=box;
      Point high=box.high();high[0]=box.low()[0];
      const Box plane(box.low(),high);
      auto basis=pointwiseSphericalBasis(plane,this->block(),-1);
      auto det=CubedSphereShell::detA<T,MEM>(plane,this->dx()[1],0);
      auto mat=CubedSphereShell::detA_A<T,MEM>(plane,this->dx()[1],this->block(),0);
      for (auto p : plane)
      {
        std::array<T,7> a;
        for (int d=0; d<DIM; ++d) {a[d]=basis(p,d,0);a[4+d]=mat(p,d,0);}
        a[3]=det(p,0);cache->angular.push_back(a);
      }
      high=box.low();high[0]=box.high()[0];
      const Box line(box.low(),high);
      Scalar radius(line),volume(line);
      Vector dr(line),adj(line);
      radialMetrics(radius,dr,adj,volume,line);
      for (auto p : line) cache->radial.push_back(dr(p,0));
      m_gravityGeometry=std::move(cache);
    }
    const auto& cache=*m_gravityGeometry;
    const auto low=box.low(),high=box.high();
    const auto sourceBox=w.box();
    const auto sizes=sourceBox.sizes();
    const size_t sourceStride=sourceBox.size(),destStride=box.size();
    const T* source=w.data();T* dest=rhs.data();
    size_t angular=0,offset=0;
    for (int k=low[2]; k<=high[2]; ++k) for (int j=low[1]; j<=high[1]; ++j)
    {
      const auto& a=cache.angular[angular++];
      const size_t row=(low[0]-sourceBox.low()[0])+size_t(sizes[0])*
          ((j-sourceBox.low()[1])+size_t(sizes[1])*(k-sourceBox.low()[2]));
      for (int i=low[0]; i<=high[0]; ++i,++offset)
      {
        const size_t s=row+i-low[0];
        const T factor=-c_G*c_M0*cache.radial[i-low[0]]*source[s];
        T radialVelocity=0;
        for (int d=0; d<DIM; ++d)
        {
          dest[offset+(iMOMX+d)*destStride]+=factor*a[4+d];
          radialVelocity+=a[d]*source[s+(iVX+d)*sourceStride];
        }
        dest[offset+iE*destStride]+=factor*a[3]*radialVelocity;
      }
    }
    return;
  }
  Scalar radius(box),volume(box);
  Vector dr(box),adj(box);
  radialMetrics(radius,dr,adj,volume,box);
  auto basis=pointwiseSphericalBasis(box,this->block(),-1);
  auto det=CubedSphereShell::detA<T,MEM>(box,this->dx()[1],0);
  auto mat=CubedSphereShell::detA_A<T,MEM>(box,this->dx()[1],this->block(),0);
  forallInPlace([] PROTO_LAMBDA(Var<T,NUMCOMPS,MEM>& rhs,
      const Var<T,NUMCOMPS,MEM>& w,const Var<T,DIM,MEM>& dr,
      const Var<T,DIM,MEM,DIM>& A,const Var<T,1,MEM>& det,
      const Var<T,DIM,MEM,DIM>& mat)
  {
      const T factor=-c_G*c_M0*dr(0)*w(iRHO);
      T radialVelocity=0;
      for (int d=0; d<DIM; ++d)
      {
          rhs(iMOMX+d)+=factor*mat(d,0);
          radialVelocity+=A(d,0)*w(iVX+d);
      }
      rhs(iE)+=factor*det(0)*radialVelocity;
  },rhs,w,dr,basis,det,mat);
}

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
