template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::DivBPowell1(
                        State& a_Rhs,
                        const State& a_WFace,
                        const BoxData<T,1,MEM>& a_deta,
                        const BoxData<T,1,MEM>& a_adjDrDir,
                        int a_dir,
                        const Box& a_bxface) const  
{
  State Bfluxdir(a_bxface);
  forallInPlace([ ] PROTO_LAMBDA
    (Var<T,NUMCOMPS,MEM>& a_fluxout,
      const Var<T,NUMCOMPS,MEM>& a_Wsph,
      const Var<T,1,MEM>& a_deta,
      int a_dir)
    {  
      for (int comp = 0; comp < NUMCOMPS;comp++)
      {                  
        a_fluxout(comp) = a_deta(0)*a_Wsph(iBX + a_dir);
      }
    },Bfluxdir,a_WFace,a_deta,a_dir);

  for (int comp = 0; comp < NUMCOMPS;comp++)
  {
    BoxData<T,1,MEM> advfluxcomp = slice(Bfluxdir,comp);
    BoxData<T,1,MEM> advfluxprod = faceProduct(advfluxcomp,a_adjDrDir,a_dir);
    advfluxprod.copyTo(advfluxcomp);
  }
  BoxData<T,NUMCOMPS,MEM> BrhsIncr;
  T BdxInv = -1.0/(this->dx()[a_dir]);
  BrhsIncr = Stencil<T>::FluxDivergence(a_dir)(Bfluxdir,BdxInv);

  forallInPlace([ ] PROTO_LAMBDA
    (Var<T,NUMCOMPS,MEM>& a_rhs,
      const Var<T,NUMCOMPS,MEM>& a_rhsIncr)
    {  
      for (int comp = 0; comp < NUMCOMPS;comp++)
      {                  
        a_rhs(comp) += a_rhsIncr(comp);
      }
    },a_Rhs,BrhsIncr);
}

template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::DivBPowell2(
                        State& a_Rhs,
                        State& a_Brhs,
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
        const Var<T,DIM,MEM>& a_Vcart,
        const Var<T,DIM,MEM>& a_Bcart)
      {  
        T VdotB = a_Vcart(0)*a_Bcart(0) + a_Vcart(1)*a_Bcart(1) + a_Vcart(2)*a_Bcart(2);                
        a_rhs(iRHO) *= 0;
        a_rhs(iMOMX) *= a_Bcart(0)/4.0/M_PI;
        a_rhs(iMOMY) *= a_Bcart(1)/4.0/M_PI;
        a_rhs(iMOMZ) *= a_Bcart(2)/4.0/M_PI;
        a_rhs(iE) *= VdotB/4.0/M_PI;
        a_rhs(iBX) *= a_Vcart(0);
        a_rhs(iBY) *= a_Vcart(1);
        a_rhs(iBZ) *= a_Vcart(2);
        
      },a_Brhs,Vcart,Bcart);

  a_Rhs += a_Brhs;
}