#include "BoxOp_EulerCubedSphere.H"
#include <iostream>
#include <limits>

MHDReader BC_global;
using OP = BoxOp_EulerCubedSphere<double,MBMap_CubedSphereShell<HOST>,HOST>;
using Level = MBLevelBoxData<double,NUMCOMPS,HOST>;

// A sharp radial velocity/Br front with exactly constant thermal pressure.
// No SW file, flux evaluation or time advance is required.
int main(int argc,char** argv)
{
#ifdef PR_MPI
  MPI_Init(&argc,&argv);
#endif
  int result = 0;
  {
    ParseInputs::getInstance().parsenow(argc,argv);
    const int n = ParseInputs::get_domainSize(), nr = ParseInputs::get_thickness();
    if (ParseInputs::get_init_condition_type() != 1
        || ParseInputs::get_seam_pressure_correction() != 1 || n < 48 || nr < 20)
      result = 1;
    auto domain = CubedSphereShell::Domain(n,nr,0);
    MBDisjointBoxLayout layout(domain,Point(ParseInputs::get_boxSize_rad(),
        ParseInputs::get_boxSize_nonrad(),ParseInputs::get_boxSize_nonrad()));
    auto& ops = *OP::levelGeometry(layout).op;
    Level initial(layout,OP::ghost()), before(layout,OP::ghost()), after(layout,OP::ghost());
    MBLevelBoxData<double,1,HOST> volume(layout,OP::ghost()+Point::Basis(0,2));
    constexpr double rho = 700*c_MP, pressure = 8e-8;
    initial.setVal(0.);
    for (auto index : layout)
    {
      const int block = layout.block(index);
      Scalar radius(volume[index].box());
      Vector dr(volume[index].box()),adj(volume[index].box());
      ops[index].radialMetrics(radius,dr,adj,volume[index],dr.box());
      State w(initial[index].box()),ju;
      auto basis = ops[index].pointwiseSphericalBasis(w.box(),block,-1);
      w.setVal(0.);
      for (auto p : w.box())
      {
        const double transition = std::tanh((basis(p,1,0)-basis(p,0,0))
            /std::sqrt(2.)/(3*M_PI/180));
        w(p,iRHO) = rho;
        w(p,iVX) = 5e7+2e7*transition;
        w(p,iP) = pressure;
        w(p,iBX) = 6e-4*transition;
      }
      ops[index].primToCons(ju,w,volume[index],ParseInputs::get_gamma(),1./n,block);
      ju.copyTo(initial[index],layout[index]);
    }
    const int padding = ParseInputs::get_state_conversion() == 1 ? 3 : 4;
    auto interpolation = CubedSphereShell::MHDInterpOp<HOST>(
        layout,OP::ghost()+Point(0,padding,padding),4);
    initial.copyTo(before);
    HelioCubed::consToSphInterpEuler(before,interpolation,volume,4,0);
    OP::fillRadialBoundaries(before,volume,0.);
    initial.copyTo(after);
    OP::prepareState(after,interpolation,volume,0.);

    double unchangedError = 0., sums[3] = {0.,0.,0.};
    long changed = 0, radialChanges[2] = {0,0};
    const double scale[NUMCOMPS] = {rho,rho*7e7,rho*7e7,rho*7e7,1e-6,6e-4,6e-4,6e-4};
    for (auto index : layout)
    {
      State wb,barb,wa,bara;
      ops[index].sphToPrim(wb,barb,before[index],layout.block(index));
      ops[index].sphToPrim(wa,bara,after[index],layout.block(index));
      auto basis = ops[index].pointwiseSphericalBasis(wa.box(),layout.block(index),-1);
      for (auto p : after[index].box())
      {
        const bool angular = p[1] < 0 || p[1] >= n || p[2] < 0 || p[2] >= n;
        const bool corrected = angular && p[0] >= 0 && p[0] < nr;
        for (int c = 0; c < NUMCOMPS; ++c)
        {
          const double delta = std::abs(after[index](p,c)-before[index](p,c))/scale[c];
          if (!std::isfinite(delta)) unchangedError = std::numeric_limits<double>::infinity();
          else if (c != iE || !corrected) unchangedError = std::max(unchangedError,delta);
          else if (delta > 1e-12)
          {
            ++changed;
            if (p[0] == 0) ++radialChanges[0];
            if (p[0] == nr-1) ++radialChanges[1];
          }
        }
      }
      // Compare with known physical pressure, independently of the transferred
      // target. Avoid radial boundaries, whose analytic BC is a different flow.
      for (auto p : wa.box())
      {
        const int distance = std::min(std::min(p[1],p[2]),std::min(n-1-p[1],n-1-p[2]));
        if (distance != -1 || p[0] != nr/2) continue;
        // Count each seam location once even where patch halos overlap.
        bool duplicate=false;
        for (int d=0;d<DIM;++d)
        {
          const int extent=d==0 ? nr : n;
          if (p[d]>=0 && p[d]<extent
              && (p[d]<layout[index].low()[d] || p[d]>layout[index].high()[d]))
            duplicate=true;
        }
        if (duplicate) continue;
        const double q = (basis(p,0,0)+basis(p,1,0))/std::sqrt(2.)*std::cos(12*M_PI/180)
                       + basis(p,2,0)*std::sin(12*M_PI/180);
        if (q < std::cos(18*M_PI/180)) continue;
        const double weight = std::exp(2*(q-1)/std::pow(6*M_PI/180,2));
        sums[0] += weight*std::pow(wb(p,iP)/pressure-1,2);
        sums[1] += weight*std::pow(wa(p,iP)/pressure-1,2);
        sums[2] += weight;
      }
    }
    // Reuse the cached interpolation plan with fresh state buffers, as the
    // RK stages and output path do. No previous ghost contents may leak in.
    Level repeat(layout,OP::ghost());
    initial.copyTo(repeat);
    OP::prepareState(repeat,interpolation,volume,0.);
    for (auto index : layout) for (auto p : repeat[index].box())
      for (int c=0;c<NUMCOMPS;++c)
      {
        const double delta=std::abs(repeat[index](p,c)-after[index](p,c))/scale[c];
        if (!std::isfinite(delta)) unchangedError=std::numeric_limits<double>::infinity();
        else unchangedError=std::max(unchangedError,delta);
      }
#ifdef PR_MPI
    double globalError,globalSums[3]; long globalChanged,globalRadial[2];
    MPI_Allreduce(&unchangedError,&globalError,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(sums,globalSums,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&changed,&globalChanged,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(radialChanges,globalRadial,2,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
    unchangedError=globalError; changed=globalChanged;
    for (int i=0;i<3;++i) sums[i]=globalSums[i];
    for (int i=0;i<2;++i) radialChanges[i]=globalRadial[i];
#endif
    const double errorBefore=std::sqrt(sums[0]/sums[2]),errorAfter=std::sqrt(sums[1]/sums[2]);
    if (!std::isfinite(unchangedError) || unchangedError > 5e-13 || !changed
        || !radialChanges[0] || !radialChanges[1] || !(sums[2]>0)
        || !std::isfinite(errorBefore) || !std::isfinite(errorAfter)
        || errorAfter > .005 || errorAfter > .1*errorBefore)
      result=1;
    if (!procID()) std::cout << "Seam pressure RMS: " << errorBefore << " -> " << errorAfter
      << "; protected state maximum scaled change: " << unchangedError
      << "; corrected cells: " << changed << "; radial end layers: "
      << radialChanges[0] << ", " << radialChanges[1] << std::endl;
  }
#ifdef PR_MPI
  MPI_Finalize();
#endif
  return result;
}
