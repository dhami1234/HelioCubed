template<typename T, typename MAP, MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::BertaPowellMagneticAverages(
    BoxData<T,DIM,MEM>& a_BAvg,
    const State& a_USph,
    uint a_block) const
{
  // Magnetic test states contain perturbations; reconstructed primitives
  // contain the full field. Include the background in the smooth source too,
  // so switching to the troubled-face flux cannot introduce a background jump.
  const Box pointBox = a_USph.box().grow(1);
  State backgroundState(pointBox);
  backgroundState.setVal(0.0);
  addTestProblemBackgroundB(backgroundState,a_block);
  auto backgroundPoint = forall<T,DIM,MEM>([] PROTO_LAMBDA(
      Var<T,DIM,MEM>& a_B,
      const Var<T,NUMCOMPS,MEM>& a_W)
    {
      for (int d = 0; d < DIM; ++d) a_B(d) = a_W(iBX+d);
    },backgroundState);
  a_BAvg = (1.0*Shift(Point::Zeros())
      + (1.0/24.0)*Stencil<T>::Laplacian())(backgroundPoint,a_USph.box());
  forallInPlace([] PROTO_LAMBDA(
      Var<T,DIM,MEM>& a_B,
      const Var<T,NUMCOMPS,MEM>& a_U)
    {
      for (int d = 0; d < DIM; ++d) a_B(d) += a_U(iBX+d);
    },a_BAvg,a_USph);
}

template<
    typename T,
    typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::DivBPowell1(
                        BoxData<T,1,MEM>& a_minusJDivB,
                        const Scalar& a_Bnormal,
                        const BoxData<T,1,MEM>& a_deta,
                        const BoxData<T,1,MEM>& a_adjDrDir,
                        int a_dir,
                        const Scalar* a_BnormalPoint,
                        const Scalar* a_thetaFluxFace) const
{
  // The transverse covariance is required for fourth-order div(B) on the
  // curved grid. A face point multiplied by an averaged metric is only second order.
  auto Bfluxdir = faceProduct(a_deta,a_Bnormal,a_dir);
  if (a_thetaFluxFace && a_BnormalPoint)
    forallInPlace([ ] PROTO_LAMBDA
      (Var<T,1,MEM>& a_flux,
       const Var<T,1,MEM>& a_theta,
       const Var<T,1,MEM>& a_det,
       const Var<T,1,MEM>& a_Bnormal)
      {
        // Retain the existing local flux near troubled faces; use the same
        // point-to-average switch as Berta's physical flux.
        if (a_theta(0) <= T(0.5)) a_flux(0) = a_det(0)*a_Bnormal(0);
      },Bfluxdir,*a_thetaFluxFace,a_deta,*a_BnormalPoint);

  BoxData<T,1,MEM> Bfluxprod = faceProduct(Bfluxdir,a_adjDrDir,a_dir);
  PROTO_ASSERT(Bfluxprod.box().containsBox(a_minusJDivB.box().extrude(a_dir)),
      "Powell magnetic face product does not cover the requested divergence");

  BoxData<T,1,MEM> BrhsIncr;
  T BdxInv = -1.0/(this->dx()[a_dir]);
  BrhsIncr = Stencil<T>::FluxDivergence(a_dir)(Bfluxprod,BdxInv);

  forallInPlace([ ] PROTO_LAMBDA
    (Var<T,1,MEM>& a_rhs,
      const Var<T,1,MEM>& a_rhsIncr)
    {  
      a_rhs(0) += a_rhsIncr(0);
    },a_minusJDivB,BrhsIncr);
}

template<
    typename T,
    typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::DivBPowell2(
                        State& a_Rhs,
                        const BoxData<T,1,MEM>& a_minusJDivB,
                        const State& a_W,
                        const BoxData<T,DIM,MEM,DIM>& a_A_matrix0,
                        bool a_cellAverages) const
{
  PR_TIME("BoxOp_EulerCubedSphere::DivBPowell2()");
  if (ParseInputs::get_spatial_order() == 2)
  {
    forallInPlace([] PROTO_LAMBDA(
        Var<T,NUMCOMPS,MEM>& rhs,const Var<T,1,MEM>& minusDiv,
        const Var<T,NUMCOMPS,MEM>& w,const Var<T,DIM,MEM,DIM>& A)
    {
      T v[DIM]={},b[DIM]={},vb=0;
      for (int i=0;i<DIM;++i)
      {
        for (int j=0;j<DIM;++j) {v[i]+=A(i,j)*w(iVX+j);b[i]+=A(i,j)*w(iBX+j);}
        rhs(iMOMX+i)+=minusDiv(0)*b[i]/(4*M_PI);
        rhs(iBX+i)+=minusDiv(0)*v[i];
        vb+=v[i]*b[i];
      }
      rhs(iE)+=minusDiv(0)*vb/(4*M_PI);
    },a_Rhs,a_minusJDivB,a_W,a_A_matrix0);
    return;
  }
  const Box sourceBox = a_Rhs.box().grow(1);
  const Box pointBox = a_Rhs.box().grow(2);
  PROTO_ASSERT(a_minusJDivB.box().containsBox(sourceBox),
      "Powell source product requires one divergence ghost layer");
  PROTO_ASSERT(a_A_matrix0.box().containsBox(pointBox),
      "Powell source coefficient lacks pointwise basis coverage");

  // The standard operator supplies <W>; Berta already supplies point values.
  // Convert the complete coefficient at points, including V dot B, then
  // average it. Multiplying a point basis by <W> would be only second order.
  State deconvolved;
  const State* WPoint = &a_W;
  if (a_cellAverages)
  {
    PROTO_ASSERT(a_W.box().containsBox(pointBox.grow(1)),
        "Powell primitive deconvolution lacks stencil coverage");
    deconvolved = (1.0*Shift(Point::Zeros())
        - (1.0/24.0)*Stencil<T>::Laplacian())(a_W,pointBox);
    WPoint = &deconvolved;
  }
  else
  {
    PROTO_ASSERT(a_W.box().containsBox(pointBox),
        "Powell point primitives lack stencil coverage");
  }

  State coefficientPoint = forallOverwrite<T,NUMCOMPS,MEM>([] PROTO_LAMBDA(
      Var<T,NUMCOMPS,MEM>& a_coefficient,
      const Var<T,NUMCOMPS,MEM>& a_w,
      const Var<T,DIM,MEM,DIM>& a_A)
    {
      T velocity[DIM] = {}, magnetic[DIM] = {};
      for (int i = 0; i < DIM; ++i)
      for (int j = 0; j < DIM; ++j)
      {
        velocity[i] += a_A(i,j)*a_w(iVX+j);
        magnetic[i] += a_A(i,j)*a_w(iBX+j);
      }
      a_coefficient(iRHO) = 0;
      a_coefficient(iE) = 0;
      for (int i = 0; i < DIM; ++i)
      {
        a_coefficient(iMOMX+i) = magnetic[i]/(4.0*M_PI);
        a_coefficient(iBX+i) = velocity[i];
        a_coefficient(iE) += velocity[i]*magnetic[i]/(4.0*M_PI);
      }
    },pointBox,*WPoint,a_A_matrix0);
  State coefficientAvg = (1.0*Shift(Point::Zeros())
      + (1.0/24.0)*Stencil<T>::Laplacian())(coefficientPoint,sourceBox);

  // <S K> = <S><K> + sum_d delta_d(S) delta_d(K)/12 + O(h^4),
  // with S = -J div(B). Only the source needs the additional one-cell halo.
  for (int comp = 1; comp < NUMCOMPS; ++comp)
  {
    auto coefficient = slice(coefficientAvg,comp);
    auto source = cellProduct(a_minusJDivB,coefficient);
    PROTO_ASSERT(source.box().containsBox(a_Rhs.box()),
        "Powell source product does not cover the requested RHS");
    auto rhs = slice(a_Rhs,comp);
    rhs += source;
  }
}

template<
    typename T,
    typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::computeDivB(
                        BoxData<T,1,MEM>& a_divB,
                        const State& a_USph) const
{
  PR_TIME("BoxOp_EulerCubedSphere::computeDivB()");

  const int block = this->block();
  const T dxiPerp = this->dx()[1];
  Box bxRad = a_USph.box().grow(Point::Basis(0,2));
  BoxData<T,DIM,MEM> Dr(bxRad);
  BoxData<T,1,MEM> radius(bxRad);
  BoxData<T,DIM,MEM> adjDr(bxRad);
  BoxData<T,1,MEM> dVolr(bxRad);
  radialMetrics(radius,Dr,adjDr,dVolr,bxRad);

  BoxData<T,1,MEM> minusJDivB(a_divB.box());
  minusJDivB.setVal(0.0);

  if (ParseInputs::get_spatial_order() == 2)
  {
    State native,w;
    const bool magneticTest=isConstantMagneticFieldTest(ParseInputs::get_init_condition_type());
    sphToPrim_2O(native,a_USph,block,!magneticTest);
    auto basis=pointwiseSphericalBasis(native.box(),block,-1);
    HelioCubed::nativePrimitivesToCartesian(w,native,basis);
    for (int dir=0;dir<DIM;++dir)
    {
      const Box faces=a_divB.box().extrude(dir);
      State left,right;
      reconstructPLM(left,right,w,faces,dir);
      if (magneticTest)
      {
        const T b0=testProblemBackgroundBx(ParseInputs::get_init_condition_type());
        forallInPlace([] PROTO_LAMBDA(Var<T,NUMCOMPS,MEM>& l,
            Var<T,NUMCOMPS,MEM>& r,T b0)
        {l(iBX)+=b0;r(iBX)+=b0;},left,right,b0);
      }
      auto normal=CubedSphereShell::adjA<T,MEM>(faces,dxiPerp,block,dir);
      auto area=slice(adjDr,dir);
      auto flux=forallOverwrite<T,1,MEM>([] PROTO_LAMBDA(
          Var<T,1,MEM>& f,const Var<T,NUMCOMPS,MEM>& l,
          const Var<T,NUMCOMPS,MEM>& r,const Var<T,DIM,MEM>& normal,
          const Var<T,1,MEM>& area)
      {f(0)=0;for(int j=0;j<DIM;++j)f(0)+=T(0.5)*(l(iBX+j)+r(iBX+j))*normal(j)*area(0);},left,right,normal,area);
      minusJDivB+=Stencil<T>::FluxDivergence(dir)(flux,-T(1)/this->dx()[dir]);
    }
  }
  else if (ParseInputs::get_spatial_operator() == 2)
  {
    const T eta_d = ParseInputs::get_Berta_eta_d();
    const int detector_type = ParseInputs::get_Berta_detector_type();
    const int fallback_type = ParseInputs::get_Berta_fallback_type();

    State WAvgForDetector;
    sphAvgToPrim(WAvgForDetector,a_USph,block);

    BoxData<T,1,MEM> thetaCell;
    if (detector_type == 1)
    {
      BertaComputeThetaDerivativeRatio(thetaCell,WAvgForDetector,eta_d);
    }
    else
    {
      BertaComputeThetaJameson(thetaCell,WAvgForDetector,eta_d);
    }

    State UPoint;
    BertaAverageToPoint(UPoint,a_USph,thetaCell);
    State WPoint;
    sphPointToPrim(WPoint,UPoint,block);

    const bool magneticTest =
        isConstantMagneticFieldTest(ParseInputs::get_init_condition_type());
    BoxData<T,DIM,MEM> magneticAvg;
    if (magneticTest) BertaPowellMagneticAverages(magneticAvg,a_USph,block);

    for (int dir = 0; dir < DIM; dir++)
    {
      State WFace_L;
      State WFace_R;
      BoxData<T,1,MEM> thetaFluxFace;
      BoxData<T,1,MEM> thetaReconFace;
      BertaReconstructWithFallback(
          WFace_L,WFace_R,thetaFluxFace,thetaReconFace,
          WPoint,thetaCell,dir,fallback_type);

      State WFace_C(WFace_L.box());
      WFace_L.copyTo(WFace_C);
      WFace_C += WFace_R;
      WFace_C *= 0.5;

      const Box bxface = minusJDivB.box().extrude(dir).grow(
          2*(Point::Ones()-Point::Basis(dir)));
      BoxData<T,1,MEM> deta =
          CubedSphereShell::detA<T,MEM>(bxface,dxiPerp,dir);
      BoxData<T,1,MEM> adjDrDir = slice(adjDr,dir);
      auto Bnormal = magneticTest ? slice(magneticAvg,dir) : slice(a_USph,iBX + dir);
      Scalar Bface = Stencil<T>::CellToFace(dir)(Bnormal);
      auto Bpoint = slice(WFace_C,iBX + dir);
      DivBPowell1(minusJDivB,Bface,deta,adjDrDir,dir,&Bpoint,&thetaFluxFace);
    }
  }
  else
  {
    State W;
    State WBar;
    sphToPrim(W,WBar,a_USph,block);

    int limiter = ParseInputs::get_limiter();
    if (limiter != 0 && limiter != 1 && limiter != 2) limiter = 1;

    for (int dir = 0; dir < DIM; dir++)
    {
      State WFace = Stencil<T>::CellToFace(dir)(W);
      Box bxface = WFace.box().grow(-2);
      if (limiter == 2) bxface = WFace.box();

      BoxData<T,1,MEM> deta =
          CubedSphereShell::detA<T,MEM>(bxface,dxiPerp,dir);
      BoxData<T,1,MEM> adjDrDir = slice(adjDr,dir);
      auto Bnormal = slice(WFace,iBX + dir);
      DivBPowell1(minusJDivB,Bnormal,deta,adjDrDir,dir);
    }
  }

  BoxData<T,1,MEM> detaCell =
      CubedSphereShell::detA<T,MEM>(a_divB.box(),dxiPerp,0);
  forallInPlace([ ] PROTO_LAMBDA
      (Var<T,1,MEM>& a_divB,
       const Var<T,1,MEM>& a_minusJDivB,
       const Var<T,1,MEM>& a_deta,
       const Var<T,1,MEM>& a_dVolr)
      {
        a_divB(0) = -a_minusJDivB(0)/(a_deta(0)*a_dVolr(0));
      },a_divB,minusJDivB,detaCell,dVolr);
}
