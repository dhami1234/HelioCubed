#pragma once

#include "Proto_accel.H"
#include <cmath>

namespace MHD_Reconstruction
{
// This also works on mapped conserved states: the cell-volume factor is
// positive, while total energy may be stored relative to a background.
template<unsigned int C, typename ConservedState>
inline ACCEL_DECORATION bool finiteConservedStateWithPositiveDensity(
    const ConservedState& a_state)
{
  for (unsigned int c = 0; c < C; ++c)
    if (!std::isfinite(a_state(c))) return false;
  return a_state(0) > 0.0;
}

// Deconvolution can give a negative point density from positive cell averages.
// Revert the entire reconstructed conserved state, rather than flooring only
// density and dividing the unrepaired momentum by that floor. The cell average
// itself is never modified. A bad cell average requires a separate remedy.
template<unsigned int C, typename PointState, typename CellState>
inline ACCEL_DECORATION bool useCellAverageForInvalidDensity(
    PointState& a_point, const CellState& a_cell, double a_densityFloor)
{
  if (std::isfinite(a_point(0)) && a_point(0) >= a_densityFloor) return false;
  if (!std::isfinite(a_cell(0)) || a_cell(0) < a_densityFloor) return false;
  for (unsigned int c = 0; c < C; ++c) a_point(c) = a_cell(c);
  return true;
}

template<unsigned int C, typename PrimitiveState>
inline ACCEL_DECORATION bool usablePrimitiveState(
    const PrimitiveState& a_state, double a_densityFloor)
{
  static_assert(C > 4, "Primitive states require density and gas pressure.");
  for (unsigned int c = 0; c < C; ++c)
    if (!std::isfinite(a_state(c))) return false;
  return a_state(0) >= a_densityFloor && a_state(4) >= 0.0;
}

// Apply the same complete-state recovery when deconvolution gives an invalid
// gas pressure. Clipping pressure alone leaves the inconsistent velocity.
template<unsigned int C, typename PointState, typename CellState>
inline ACCEL_DECORATION bool useCellAverageForInvalidPrimitive(
    PointState& a_point, const CellState& a_cell, double a_densityFloor)
{
  if (usablePrimitiveState<C>(a_point,a_densityFloor)) return false;
  if (!usablePrimitiveState<C>(a_cell,a_densityFloor)) return false;
  for (unsigned int c = 0; c < C; ++c) a_point(c) = a_cell(c);
  return true;
}

// Repair the common inputs to the physical flux and Riemann dissipation.
// The consumers subtract the transverse correction to obtain point states,
// so add it here to recover the donor states in both consumers. Repair both
// sides together when either point state is invalid.
template<unsigned int C, typename FaceState, typename PointState,
         typename CellState, typename Correction>
inline ACCEL_DECORATION bool repairRiemannFaceStates(
    FaceState& a_left, FaceState& a_right,
    const PointState& a_pointLeft, const PointState& a_pointRight,
    const CellState& a_cellLeft, const CellState& a_cellRight,
    const Correction& a_correction, double a_densityFloor)
{
  if (usablePrimitiveState<C>(a_pointLeft,a_densityFloor) &&
      usablePrimitiveState<C>(a_pointRight,a_densityFloor)) return false;
  if (!usablePrimitiveState<C>(a_cellLeft,a_densityFloor) ||
      !usablePrimitiveState<C>(a_cellRight,a_densityFloor)) return false;
  for (unsigned int c = 0; c < C; ++c)
  {
    a_left(c) = a_cellLeft(c) + a_correction(c);
    a_right(c) = a_cellRight(c) + a_correction(c);
  }
  return true;
}

// A viscosity coefficient must not acquire an enormous wave speed from an
// interpolated near-zero/negative density between two positive cells. Bound
// thermodynamic face values by the neighboring cell values. This limits the
// coefficient of artificial viscosity, not the reconstructed Riemann states.
inline ACCEL_DECORATION double boundedFaceValue(
    double a_value, double a_left, double a_right, double a_floor)
{
  const double low = fmax(a_floor,fmin(a_left,a_right));
  const double high = fmax(low,fmax(a_left,a_right));
  return fmax(low,fmin(high,a_value));
}
}
