#pragma once

#include "MHD_Roe8Wave.H"

template<typename T, typename MAP, MemType MEM>
inline void BoxOp_EulerCubedSphere<T,MAP,MEM>::repairRoeFaceStates(
    State& a_WFace_L, State& a_WFace_R,
    const State& a_WFaceBar, const State& a_WBar, uint a_dir) const
{
  PR_TIMERS("BoxOp_EulerCubedSphere::repairRoeFaceStates()");
  State pointLeft = _deconvolveFace(a_WFace_L,a_WFaceBar,a_dir);
  State pointRight = _deconvolveFace(a_WFace_R,a_WFaceBar,a_dir);
  State correction = Stencil<T>::LaplacianFace(a_dir)(a_WFaceBar,T(1.0/24.0));
  State donorLeft = (1.0*Shift(Point::Basis(a_dir,-1)))(a_WBar);
  forallInPlace([] PROTO_LAMBDA(
      Var<T,NUMCOMPS,MEM>& a_left, Var<T,NUMCOMPS,MEM>& a_right,
      const Var<T,NUMCOMPS,MEM>& a_pointLeft,
      const Var<T,NUMCOMPS,MEM>& a_pointRight,
      const Var<T,NUMCOMPS,MEM>& a_cellLeft,
      const Var<T,NUMCOMPS,MEM>& a_cellRight,
      const Var<T,NUMCOMPS,MEM>& a_correction)
    {
      MHD_Reconstruction::repairRiemannFaceStates<NUMCOMPS>(
          a_left,a_right,a_pointLeft,a_pointRight,a_cellLeft,a_cellRight,
          a_correction,RHO00);
    },a_WFace_L,a_WFace_R,pointLeft,pointRight,donorLeft,a_WBar,correction);
}

template<typename T, typename MAP, MemType MEM>
inline void BoxOp_EulerCubedSphere<T,MAP,MEM>::primToMappedRoeFlux(
    State& a_fluxdir,
    const State& a_WFace_L,
    const State& a_WFace_R,
    const State& a_WFaceBar,
    const BoxData<T,1,MEM>& a_adjDdir,
    const BoxData<T,DIM,MEM,DIM>& a_detA_A,
    const BoxData<T,1,MEM>& a_detA,
    const BoxData<T,DIM,MEM,DIM>& a_A,
    uint a_dir,
    const State* a_leftPoint,
    const State* a_rightPoint,
    const BoxData<T,1,MEM>* a_radialPoint) const
{
  PR_TIMERS("BoxOp_EulerCubedSphere::primToMappedRoeFlux()");
  // Standalone callers can still supply face averages only.
  State ownedLeft,ownedRight;
  if (!a_leftPoint) { ownedLeft = _deconvolveFace(a_WFace_L,a_WFaceBar,a_dir); a_leftPoint = &ownedLeft; }
  if (!a_rightPoint) { ownedRight = _deconvolveFace(a_WFace_R,a_WFaceBar,a_dir); a_rightPoint = &ownedRight; }
  const State& left = *a_leftPoint;
  const State& right = *a_rightPoint;
  BoxData<T,1,MEM> ownedRadial;
  if (!a_radialPoint)
  {
    if (a_dir == 0)
    {
      ownedRadial.define(a_adjDdir.box());
      a_adjDdir.copyTo(ownedRadial);
    }
    else
    {
      ownedRadial = (1.0*Shift(Point::Zeros())
          - (1.0/24.0)*Stencil<T>::Derivative(2,0,2))(a_adjDdir);
    }
    a_radialPoint = &ownedRadial;
  }
  const auto& radialAreaPoint = *a_radialPoint;

  const int direction = a_dir;
  // Convolution and the face products need one transverse layer around the
  // requested flux. Avoid solving the Riemann problem in unused outer layers.
  const Box available = left.box() & right.box() & a_A.box() & radialAreaPoint.box();
  const Box needed = a_fluxdir.box().grow(Point::Ones()-Point::Basis(a_dir));
  const Box pointBox = a_fluxdir.box().empty() ? available : (needed & available);
  const double entropyFix = ParseInputs::get_entropy_fix_coeff();
  State dissPoint = forall<T,NUMCOMPS,MEM>(
      [] PROTO_LAMBDA(
          Var<T,NUMCOMPS,MEM>& out,
          const Var<T,NUMCOMPS,MEM>& lo,
          const Var<T,NUMCOMPS,MEM>& hi,
          const Var<T,DIM,MEM,DIM>& A,
          const Var<T,1,MEM>& radialArea,
          int dir, T gamma, double entropyFix)
      {
        double wL[8], wR[8], matrix[3][3], diss[8];
        for (int c = 0; c < 8; ++c)
        {
          wL[c] = lo(c);
          wR[c] = hi(c);
        }
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j)
            matrix[i][j] = A(i,j);
        MHD_Riemann_Solvers::factoredRoe8WaveDissipation(
            diss,wL,wR,matrix,radialArea(0),dir,gamma,entropyFix,RHO00);
        for (int c = 0; c < 8; ++c) out(c) = diss[c];
      },pointBox,left,right,a_A,radialAreaPoint,direction,m_gamma,entropyFix);

  // Average the coefficients before applying the angular metric averages.
  // For radial hydro flow these coefficients are constant across each
  // radial face, so the same <det(A) A> used by the rest of the operator
  // preserves radial momentum. Face products retain the transverse gradient
  // terms needed for fourth-order accuracy when the coefficients vary.
  State dissAvg = (1.0*Shift(Point::Zeros())
      + (1.0/24.0)*Stencil<T>::LaplacianFace(a_dir))(dissPoint);
  for (int c : {iRHO,iE})
  {
    auto avg = slice(dissAvg,c);
    auto point = slice(dissPoint,c);
    auto mapped = Operator::_faceProduct(a_detA,avg,a_detA,point,a_dir);
    PROTO_ASSERT(mapped.box().containsBox(a_fluxdir.box()),
        "primToMappedRoeFlux | Insufficient scalar face stencil coverage.");
    auto out = slice(a_fluxdir,c);
    mapped.copyTo(out);
  }
  for (int c : {iMOMX,iBX})
  {
    auto avg = slice<T,NUMCOMPS,DIM,MEM>(dissAvg,c);
    auto point = slice<T,NUMCOMPS,DIM,MEM>(dissPoint,c);
    auto mapped = Operator::_faceMatrixProductAB(
        a_detA_A,avg,a_detA_A,point,a_dir);
    PROTO_ASSERT(mapped.box().containsBox(a_fluxdir.box()),
        "primToMappedRoeFlux | Insufficient vector face stencil coverage.");
    auto out = slice<T,NUMCOMPS,DIM,MEM>(a_fluxdir,c);
    mapped.copyTo(out);
  }
}
