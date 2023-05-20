
#ifndef _PROTO_AMR_Euler__
#define _PROTO_AMR_Euler__
#define RKORDER 1
#define RKSTAGES 1
#include "Proto.H"

extern int TIME_STEP;
extern int PLOT_NUM;

template <
    template <typename, MemType> class OpType,
    typename T,
    template <typename, unsigned int, MemType, Centering> class BCType = PeriodicBC,
    MemType MEM = MEMTYPE_DEFAULT>
class AMREuler
{
public:
    typedef OpType<T, MEM> OP;
    typedef BCType<T, OP::numState(), MEM, PR_CELL> BC;
    typedef LevelOp<OpType, double, BCType, MEM> LOP;
    typedef AMRData<T, OP::numState(), MEM> AMRStateData;

    inline AMREuler() {}
    inline AMREuler(AMRStateData &a_data, const Array<T, DIM> &a_dx,
                    int a_regridInterval = 1, int a_bufferSize = 2);
    inline void define(AMRStateData &a_data, const Array<T, DIM> &a_dx,
                       int a_regridInterval = 1, int a_bufferSize = 2);
    inline void advance(T a_dt);

private:
    inline void define(AMRStateData &a_data, const Array<T, DIM> &a_dx, int a_level,
                       int a_regridInterval = 1, int a_bufferSize = 2);
    inline void advance(int a_level, T a_dt, T a_chi);
    inline void interpBoundaries(int a_level, int a_stage, T a_chi);
    inline void stageUpdate(int a_level, int a_stage, T a_dt);
    inline void regrid(int a_baseLevel);

    Array<double, RKSTAGES> m_stageWeights = {1.};
    Array<double, RKSTAGES + 1> m_stageTimes = {0., 1.};

    AMRStateData *m_data;

    std::vector<int> m_timeRefRatios;
    std::vector<Point> m_spaceRefRatios;
    std::vector<LOP> m_ops;
    std::vector<int> m_numSteps;
    std::vector<int> m_regridInterval;
    std::vector<int> m_bufferSize;
    std::vector<InterpStencil<T>> m_cfInterp;

    std::vector<std::shared_ptr<LevelBoxData<T, OP::numState(), MEM>>> m_UStage;
    std::vector<std::shared_ptr<LevelBoxData<T, OP::numState(), MEM>>> m_RHSTotal;
    std::vector<std::shared_ptr<LevelBoxData<T, OP::numState() * RKORDER, MEM>>> m_UTaylor;
    std::vector<std::shared_ptr<LevelBoxData<T, OP::numState() * RKORDER, MEM>>> m_UTaylorC;
    std::vector<std::shared_ptr<LevelFluxRegister<T, OP::numState(), MEM>>> m_fluxRegisters;
};

#include "MHD_AMR_EulerImplem.H"
#endif // end include guard
