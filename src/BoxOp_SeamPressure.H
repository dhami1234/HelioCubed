#pragma once

#include <limits>

template<typename T, typename MAP, MemType MEM>
inline void BoxOp_EulerCubedSphere<T,MAP,MEM>::fillRadialBoundaries(
    MBLevelBoxData<T,NUMCOMPS,MEM>& a_state,
    MBLevelBoxData<T,1,MEM>& a_dVolr, T a_time)
{
  PR_TIMERS("BoxOp_EulerCubedSphere::fillRadialBoundaries()");
  const auto& layout = a_state.layout();
  auto& geometry = levelGeometry(layout);
  const bool fileBoundary = ParseInputs::get_init_condition_type() == 0;
  MBLevelBoxData<T,NUMCOMPS,MEM> boundary;
  BoxData<T,NUMCOMPS,MEM> unused;
  if (fileBoundary)
  {
    // The reader writes patch storage, including radial/angular corners.
    Array<Point,DIM+1> boundaryGhost;
    boundaryGhost.fill(Point::Zeros());
    boundaryGhost[0] = a_state.ghost()[0];
    boundary.define(layout,boundaryGhost);
    BC_global.BoxData_to_BC<T,MEM>(boundary,*geometry.map,a_time,true);
  }
  for (auto index : layout)
  {
    const Box block = layout.getBlock(layout.block(index)).domain().box();
    const Box storage = a_state[index].box();
    if (storage.low()[0] < block.low()[0] || storage.high()[0] > block.high()[0])
      (*geometry.op)[index].PreStagePatch(a_state[index],
          fileBoundary ? boundary[index] : unused,a_dVolr[index],block,a_time);
  }
}

template<typename T, typename MAP, MemType MEM>
inline void BoxOp_EulerCubedSphere<T,MAP,MEM>::correctSeamPressure(
    MBLevelBoxData<T,NUMCOMPS,MEM>& a_state, MBInterpOp& a_iop)
{
  PR_TIMERS("BoxOp_EulerCubedSphere::correctSeamPressure()");
  const auto& layout = a_state.layout();
  auto& ops = *levelGeometry(layout).op;
  // Match the source halo used by native state conversion. The scalar uses
  // exactly the same weights and donor stencil as the conserved interpolation.
  const int padding = ParseInputs::get_state_conversion() == 1 ? 3 : 4;
  const auto pressureGhost = a_state.ghost() + Point(0,padding,padding);
  MBLevelBoxData<T,1,MEM> pressure(layout,pressureGhost);
  MBLevelBoxData<T,1,MEM> target(layout,pressureGhost);
  pressure.setVal(std::numeric_limits<T>::quiet_NaN());
  target.setVal(std::numeric_limits<T>::quiet_NaN());
  for (auto index : layout)
  {
    State w,bar;
    ops[index].sphToPrim(w,bar,a_state[index],layout.block(index));
    auto p = slice(w,iP);
    p.copyTo(pressure[index],layout[index]);
  }
  pressure.exchange();
  pressure.copyTo(target);
  a_iop.apply(target,pressure);
  target.exchange0();

  const T inverseGammaMinusOne = T(1)/(ParseInputs::get_gamma()-T(1));
  for (auto index : layout)
  {
    const Box block = layout.getBlock(layout.block(index)).domain().box();
    const Box storage = a_state[index].box();
    if (storage.low()[1] >= block.low()[1] && storage.high()[1] <= block.high()[1]
        && storage.low()[2] >= block.low()[2] && storage.high()[2] <= block.high()[2])
      continue;
    // The first pass suffices in smooth states. Repeat against the fixed
    // target when a correction changes the existing positivity-fallback branch.
    for (int pass = 0; pass < 3; ++pass)
    {
      State w,bar;
      ops[index].sphToPrim(w,bar,a_state[index],layout.block(index));
      forallInPlace_p([] PROTO_LAMBDA(
          Point p, Var<T,NUMCOMPS,MEM>& u,
          const Var<T,NUMCOMPS,MEM>& primitive, const Var<T,1,MEM>& desired,
          Point low, Point high, T inverseGammaMinusOne)
      {
        const bool angular = p[1] < low[1] || p[1] > high[1]
                          || p[2] < low[2] || p[2] > high[2];
        // Leave physical radial BCs, including corners, and all owned cells
        // untouched. A nonpositive interpolated pressure cannot be a target.
        if (angular && p[0] >= low[0] && p[0] <= high[0]
            && std::isfinite(desired(0)) && desired(0) > T(0)
            && std::isfinite(primitive(iP)))
          u(iE) += (desired(0)-primitive(iP))*inverseGammaMinusOne;
      },w.box(),a_state[index],w,target[index],block.low(),block.high(),
          inverseGammaMinusOne);
    }
  }
}

template<typename T, typename MAP, MemType MEM>
inline void BoxOp_EulerCubedSphere<T,MAP,MEM>::prepareState(
    MBLevelBoxData<T,NUMCOMPS,MEM>& a_state, MBInterpOp& a_iop,
    MBLevelBoxData<T,1,MEM>& a_dVolr, T a_time)
{
  PR_TIMERS("BoxOp_EulerCubedSphere::prepareState()");
  if (ParseInputs::get_seam_pressure_correction() == 0)
  {
    HelioCubed::consToSphInterpEuler(a_state,a_iop,a_dVolr,4,0);
    fillRadialBoundaries(a_state,a_dVolr,a_time);
    return;
  }

  const auto& layout = a_state.layout();
  // sphToPrim consumes one halo. Prepare one extra layer in every direction
  // so energy correction covers the entire ordinary halo, including locations
  // where radial or angular patch cuts meet a block seam.
  const Point nativeGhost = a_state.ghost()[0] + Point::Ones();
  const int padding = ParseInputs::get_state_conversion() == 1 ? 3 : 4;
  const Point sourceGhost = nativeGhost + Point(0,padding,padding);
  struct InterpolationCache
  {
    LayoutT layout;
    Point ghost;
    MBInterpOp op;
    InterpolationCache(const LayoutT& a_layout, Point a_ghost)
      : layout(a_layout), ghost(a_ghost),
        op(CubedSphereShell::MHDInterpOp<MEM>(a_layout,a_ghost,4)) {}
  };
  // Retain the expensive interpolation plan, not full-level state buffers.
  // The saved layout keeps its partition identity alive across grid changes.
  static std::unique_ptr<InterpolationCache> interpolation;
  if (!interpolation || &interpolation->layout.partition() != &layout.partition()
      || interpolation->ghost != sourceGhost)
    interpolation = std::make_unique<InterpolationCache>(layout,sourceGhost);

  MBLevelBoxData<T,NUMCOMPS,MEM> native(layout,nativeGhost);
  MBLevelBoxData<T,1,MEM> volume(layout,nativeGhost+Point::Basis(0,2));
  auto& ops = *levelGeometry(layout).op;
  for (auto index : layout)
  {
    Scalar radius(volume[index].box());
    Vector dr(volume[index].box()),adj(volume[index].box());
    ops[index].radialMetrics(radius,dr,adj,volume[index],dr.box());
  }
  native.setVal(T(0));
  a_state.copyTo(native);
  HelioCubed::consToSphInterpEuler(native,interpolation->op,volume,4,0);
  fillRadialBoundaries(native,volume,a_time);
  correctSeamPressure(native,interpolation->op);
  // Copy patch storage explicitly: the corrected seam ghosts must survive.
  for (auto index : layout) native[index].copyTo(a_state[index]);
}
