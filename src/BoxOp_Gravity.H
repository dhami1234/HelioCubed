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
                        const BoxData<T,DIM,MEM,DIM>& a_A_matrix0) const
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
        const Var<T,1,MEM>& a_rho)
        {               
        a_rhs(iRHO) += 0;
        a_rhs(iMOMX) += a_rho(0)*a_Gcart(0);
        a_rhs(iMOMY) += a_rho(0)*a_Gcart(1);
        a_rhs(iMOMZ) += a_rho(0)*a_Gcart(2);
        a_rhs(iE) += a_rho(0)*a_Vsph(0)*a_Gsph(0);
        a_rhs(iBX) += 0;
        a_rhs(iBY) += 0;
        a_rhs(iBZ) += 0;

        },a_Rhs,Vsph,Gsph,Gcart,rho);
}