#include "MHD_Constants.H"

template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::Gravity(
                        State& a_Rhs,
                        const State& a_W,
                        const BoxData<T,1,MEM>& a_radius,
                        const BoxData<T,DIM,MEM,DIM>& a_A_matrix0,
                        const BoxData<T,1,MEM>& a_deta0,
                        const BoxData<T,1,MEM>& a_dVolr) const
{
    PR_TIME("BoxOp_EulerCubedSphere::Gravity()");

    BoxData<T,DIM,MEM> Gsph(a_W.box());
    forallInPlace([ ] PROTO_LAMBDA
    (Var<T,DIM,MEM>& a_Gsph,
    const Var<T,1,MEM>& a_radius)
    {  
        a_Gsph(0) = -c_G*c_M0/pow(a_radius(0),2);
        a_Gsph(1) = 0;
        a_Gsph(2) = 0;
    },Gsph,a_radius);

    BoxData<T,1,MEM> rho = slice<T,NUMCOMPS,1,MEM>(a_W,iRHO);
    BoxData<T,DIM,MEM> Vsph = slice<T,NUMCOMPS,DIM,MEM>(a_W,iVX);

    BoxData<T,DIM,MEM> Gcart = Operator::_cellMatrixProductAB(a_A_matrix0,Gsph,a_A_matrix0,Gsph);
  
    forallInPlace([ ] PROTO_LAMBDA
        (Var<T,NUMCOMPS,MEM>& a_rhs,
        const Var<T,DIM,MEM>& a_Vsph,
        const Var<T,DIM,MEM>& a_Gsph,
        const Var<T,DIM,MEM>& a_Gcart,
        const Var<T,1,MEM>& a_rho,
        const Var<T,1,MEM>& l_deta0,
        const Var<T,1,MEM>& l_dVolr)
        {               
          a_rhs(iMOMX) += l_deta0(0)*l_dVolr(0)*a_rho(0)*a_Gcart(0);
          a_rhs(iMOMY) += l_deta0(0)*l_dVolr(0)*a_rho(0)*a_Gcart(1);
          a_rhs(iMOMZ) += l_deta0(0)*l_dVolr(0)*a_rho(0)*a_Gcart(2);
          a_rhs(iE) += l_deta0(0)*l_dVolr(0)*a_rho(0)*a_Vsph(0)*a_Gsph(0);
        },a_Rhs,Vsph,Gsph,Gcart,rho,a_deta0,a_dVolr);
}

template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::Gravity_Euler(
                        State& a_Rhs,
                        const State& a_USph) const
{
    PR_TIME("BoxOp_EulerCubedSphere::Gravity_Euler()");

    int block = this->block();
    Box bxRad = a_USph.box().grow(Point::Basis(0,2));
    BoxData<T,DIM,MEM> Dr(bxRad);
    BoxData<T,1,MEM> radius(bxRad);
    BoxData<T,DIM,MEM> adjDr(bxRad);
    BoxData<T,1,MEM> dVolr(bxRad);
    radialMetrics(radius,Dr,adjDr,dVolr,bxRad);

    State W;
    State WBar;
    sphToPrim(W,WBar,a_USph,block);

    Box bx = a_USph.box();
    T dxiPerp = this->dx()[1];
    BoxData<T,1,MEM> deta0(bx);
    BoxData<T,DIM,MEM,DIM> A_matrix0(bx);
    T half = 0.5;

    BoxData<T,1,MEM> deta = CubedSphereShell::detA<T,MEM>(bx,dxiPerp,0);
    deta.copyTo(deta0);

    Array<Array<uint,DIM>,6> permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
    Array<Array<int,DIM>,6> sign = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}};
    Point high = bx.high();
    Point low = bx.low();
    high[0] = low[0];
    Box bx0(low,high);
    BoxData<T,DIM,MEM,DIM> A_matrix(bx);
    forallInPlace_p(f_Amatrix,bx0,A_matrix,permute[block],sign[block],
                    dxiPerp,half,half);
    spreadSlice(A_matrix);
    A_matrix.copyTo(A_matrix0);

    Gravity(a_Rhs,W,radius,A_matrix0,deta0,dVolr);
}
