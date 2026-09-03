template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::DivBPowell1(
                        BoxData<T,1,MEM>& a_minusJDivB,
                        const State& a_WFace,
                        const BoxData<T,1,MEM>& a_deta,
                        const BoxData<T,1,MEM>& a_adjDrDir,
                        int a_dir,
                        const Box& a_bxface) const  
{
  BoxData<T,1,MEM> Bfluxdir(a_bxface);
  forallInPlace([ ] PROTO_LAMBDA
    (Var<T,1,MEM>& a_fluxout,
      const Var<T,NUMCOMPS,MEM>& a_Wsph,
      const Var<T,1,MEM>& a_deta,
      int a_dir)
    {  
      a_fluxout(0) = a_deta(0)*a_Wsph(iBX + a_dir);
    },Bfluxdir,a_WFace,a_deta,a_dir);

  BoxData<T,1,MEM> Bfluxprod = faceProduct(Bfluxdir,a_adjDrDir,a_dir);
  Bfluxprod.copyTo(Bfluxdir);

  BoxData<T,1,MEM> BrhsIncr;
  T BdxInv = -1.0/(this->dx()[a_dir]);
  BrhsIncr = Stencil<T>::FluxDivergence(a_dir)(Bfluxdir,BdxInv);

  forallInPlace([ ] PROTO_LAMBDA
    (Var<T,1,MEM>& a_rhs,
      const Var<T,1,MEM>& a_rhsIncr)
    {  
      a_rhs(0) += a_rhsIncr(0);
    },a_minusJDivB,BrhsIncr);
}

template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::DivBPowell2(
                        State& a_Rhs,
                        const BoxData<T,1,MEM>& a_minusJDivB,
                        const State& a_W,
                        const BoxData<T,DIM,MEM,DIM>& a_A_matrix0) const
{
  PR_TIME("BoxOp_EulerCubedSphere::DivBPowell2()");
  BoxData<T,DIM,MEM> Vsph = slice<T,NUMCOMPS,DIM,MEM>(a_W,iVX);
  BoxData<T,DIM,MEM> Vcart = Operator::_cellMatrixProductAB(a_A_matrix0,Vsph,a_A_matrix0,Vsph);

  BoxData<T,DIM,MEM> Bsph = slice<T,NUMCOMPS,DIM,MEM>(a_W,iBX);
  BoxData<T,DIM,MEM> Bcart = Operator::_cellMatrixProductAB(a_A_matrix0,Bsph,a_A_matrix0,Bsph);

  forallInPlace([ ] PROTO_LAMBDA
      (Var<T,NUMCOMPS,MEM>& a_rhs,
        const Var<T,1,MEM>& a_minusJDivB,
        const Var<T,DIM,MEM>& a_Vcart,
        const Var<T,DIM,MEM>& a_Bcart)
      {  
        T VdotB = a_Vcart(0)*a_Bcart(0) + a_Vcart(1)*a_Bcart(1) + a_Vcart(2)*a_Bcart(2);
        T source = a_minusJDivB(0);
        a_rhs(iMOMX) += source*a_Bcart(0)/4.0/M_PI;
        a_rhs(iMOMY) += source*a_Bcart(1)/4.0/M_PI;
        a_rhs(iMOMZ) += source*a_Bcart(2)/4.0/M_PI;
        a_rhs(iE) += source*VdotB/4.0/M_PI;
        a_rhs(iBX) += source*a_Vcart(0);
        a_rhs(iBY) += source*a_Vcart(1);
        a_rhs(iBZ) += source*a_Vcart(2);
      },a_Rhs,a_minusJDivB,Vcart,Bcart);
}

template<
    typename T,
    template<MemType> typename MAP,
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

  if (ParseInputs::get_spatial_operator() == 2)
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

      Box bxface = WFace_C.box();
      BoxData<T,1,MEM> deta =
          CubedSphereShell::detA<T,MEM>(bxface,dxiPerp,dir);
      BoxData<T,1,MEM> adjDrDir = slice(adjDr,dir);
      DivBPowell1(minusJDivB,WFace_C,deta,adjDrDir,dir,bxface);
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
      DivBPowell1(minusJDivB,WFace,deta,adjDrDir,dir,bxface);
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
