#pragma once

#include "Proto.H"
#include <algorithm>
#include <cmath>

namespace HelioCubed
{
// The file contains spherical primitive point values. Sample at angular CELL
// centers: Proto's BCBilinearOp samples mesh nodes, which shifts the forcing
// by half a cell in a different physical direction on each cube panel.
// Source data are replicated on every rank, as in MHDReader's existing cache.
template<typename T, Proto::MemType MEM>
Proto::MBInterpOp solarWindBoundaryOp(
    const Proto::MBLevelMap<Proto::MBMap_CubedSphereShell<MEM>,MEM>& map,
    Proto::MBDisjointBoxLayout& sourceLayout,
    const std::vector<T>& theta, const std::vector<T>& phi, int ghost)
{
  using namespace Proto;
  const auto& layout = map.layout();
  MBInterpOp op;
  PROTO_ASSERT(theta.size() >= 2 && phi.size() >= 2,
      "Solar-wind boundary interpolation needs at least two samples per angle.");
  for (auto index : layout)
  {
    const int block = layout.block(index);
    const Box domain = layout.domain().getBlock(block).box();
    Box target = layout[index].face(0,Side::Lo,1).shift(0,-1)
        .grow(Point(0,ghost,ghost));
    for (auto p : target)
    {
      Array<double,DIM> xi, x;
      xi[0] = 0.0;
      for (int d = 1; d < DIM; ++d) xi[d] = (p[d]+0.5)/domain.size(d);
      double jac;
      f_cubedSphereMap(x,jac,xi,block,cubedSphereShellGeometry());
      const double polar = atan2(hypot(x[0],x[1]),x[2]);
      double azimuth = atan2(x[1],x[0]);
      if (azimuth < phi.front()) azimuth += 2.0*M_PI;
      // Clamp at the polar endpoints; wrap the longitude interpolation.
      int it = std::upper_bound(theta.begin(),theta.end(),polar)-theta.begin()-1;
      it = std::max(0,std::min(it,int(theta.size())-2));
      const double wt = std::max(0.0,std::min(1.0,
          (polar-theta[it])/(theta[it+1]-theta[it])));
      int ip = std::upper_bound(phi.begin(),phi.end(),azimuth)-phi.begin()-1;
      ip = std::max(0,std::min(ip,int(phi.size())-1));
      const int jp = (ip+1)%phi.size();
      const double nextPhi = jp == 0 ? phi[0]+2.0*M_PI : phi[jp];
      const double wp = (azimuth-phi[ip])/(nextPhi-phi[ip]);
      std::vector<T> weights = {(1-wt)*(1-wp),(1-wt)*wp,wt*(1-wp),wt*wp};
      std::vector<MBDataPoint> sources;
      for (int k = 0; k < 4; ++k)
      {
        Point q(0,it+k/2,k%2 == 0 ? ip : jp);
        sources.emplace_back(*sourceLayout.begin(),q,sourceLayout);
      }
      op.insert(std::make_shared<MBPointInterpOp>(
          MBDataPoint(index,p,layout),sources,weights));
    }
  }
  return op;
}
}
