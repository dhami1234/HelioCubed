#pragma once
#ifndef _BOX_OP_divB_
#define _BOX_OP_divB_

#include "Proto.H"
extern Parsefrominputs inputs;
// constexpr MemType MEM = MEMTYPE_DEFAULT;


#include "MHD_Limiters.H"
#include "MHD_Mapping.H"
#include "MHD_Riemann_Solvers.H"
#include "MHD_Output_Writer.H"
#include "MHD_Input_Parsing.H"
#include "MHD_Constants.H"
#include "MHD_CFL.H"
#include "MHDLevelDataRK4.H"
#include "MHD_Operator.H"


using namespace Proto;



/**
 * @brief Function to calculate the Powell term, without the divB multiplied to it. 
 * @param a_P the output terms.
 * @param a_W the input primitive variables.
 */
PROTO_KERNEL_START
void PowellF(State&         a_P,
            const State&   a_W)
{

#if DIM==2
    a_P(0) = 0.;
    a_P(1) = a_W(4)/4.0/c_PI;
    a_P(2) = a_W(5)/4.0/c_PI;
    a_P(3) = a_W(1)*a_W(4)/4.0/c_PI + a_W(2)*a_W(5)/4.0/c_PI;
    a_P(4) = a_W(1);
    a_P(5) = a_W(2);
#endif

#if DIM==3
    a_P(0) = 0.;
    a_P(1) = a_W(5)/4.0/c_PI;
    a_P(2) = a_W(6)/4.0/c_PI;
    a_P(3) = a_W(7)/4.0/c_PI;
    a_P(4) = a_W(1)*a_W(5)/4.0/c_PI + a_W(2)*a_W(6)/4.0/c_PI + a_W(3)*a_W(7)/4.0/c_PI;
    a_P(5) = a_W(1);
    a_P(6) = a_W(2);
    a_P(7) = a_W(3);

#endif
}
PROTO_KERNEL_END(PowellF, Powell)


/**
 * @brief Function to average cell centered B on both sides of a face. 
 * @param a_Bavg the output BoxData.
 * @param a_W_ave the input cell averaged BoxData of primitive variables.
 * @param a_d direction to average.
 */ 
PROTO_KERNEL_START
void BavgcalcF(State& a_Bavg,
            const State& a_W_ave,
            int a_d)
{
    for (int i=0; i< NUMCOMPS; i++) {
        a_Bavg(i) = a_W_ave(2+DIM+a_d);
    }
}
PROTO_KERNEL_END(BavgcalcF, Bavgcalc)

template<typename T, MemType MEM = MEMTYPE_DEFAULT>
class BoxOp_divB : public BoxOp<T, NUMCOMPS, 1, MEM>
{
    public:

    T gamma = inputs.gamma;
    mutable T umax;

    // How many ghost cells does the operator need from the state variables
    inline static Point ghost() { return Point::Ones(NGHOST);}
    
    // How many ghost cells does the operator need from the auxiliary variables
    inline static Point auxGhost() { return Point::Zeros();}
    
    // What is the intended order of accuracy of the operator
    inline static constexpr int order() { return 4; }
    
    inline BoxOp_divB() : BoxOp<T, NUMCOMPS, 1, MEM>() {}
    inline BoxOp_divB(Array<double, DIM> a_dx) : BoxOp<T, NUMCOMPS, 1, MEM>(a_dx) {}
    inline BoxOp_divB(double a_dx) : BoxOp<T, NUMCOMPS, 1, MEM>(a_dx) {}
    
    // Initialization
    inline void init()
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            m_interp_H[dir] = Stencil<double>::CellToFaceH(dir);
            m_interp_L[dir] = Stencil<double>::CellToFaceL(dir);
            m_interp_edge[dir] = Stencil<double>::CellToFace(dir);
            m_divergence[dir] = Stencil<double>::FluxDivergence(dir);
            m_laplacian_f[dir] = Stencil<double>::LaplacianFace(dir);
        }
    }

    // Helper Function
    inline void computeFlux(
            BoxData<T, NUMCOMPS>& a_flux,
            BoxData<T, NUMCOMPS>& a_W_ave,
            int a_dir) const
    { 
        Vector W_ave_edge = m_interp_edge[a_dir](a_W_ave);
        a_flux = forall<double,NUMCOMPS>(Bavgcalc, W_ave_edge, a_dir);
    }
   

    // Apply BCs by filling ghost cells in stage values. For Euler, this is done by calling
    // exchange. For the MHD code, it will be more complicated.
    // The interface is very provisional. We expect it to evolve as we d more real problems.
    inline void bcStage(
                        LevelBoxData<T,NUMCOMPS>& a_UStage,
                        const LevelBoxData<T,NUMCOMPS>& a_U0,
                        int a_stage)
    {
      a_UStage.exchange();
    }                 
    
    // Apply Operator
    inline void operator()(
            BoxData<T, NUMCOMPS>&                   a_Rhs,
            Array<BoxData<T, NUMCOMPS>, DIM>&  a_fluxes,
            const BoxData<T, NUMCOMPS>&             a_U,
            T                                       a_scale = 1.0) const
    {
        T dx = this->dx()[0];
        PR_TIME("BoxOp_divB::operator()");        
        // COMPUTE W_AVE
        a_Rhs.setVal(0.0);    
        Vector W_bar = forall<double, NUMCOMPS>(consToPrim, a_U, gamma);
        Vector U = Operator::deconvolve(a_U);
        Vector W = forall<double, NUMCOMPS>(consToPrim, U, gamma);
        Vector W_ave = Operator::_convolve(W, W_bar);
        

        // COMPUTE DIV FLUXES
        for (int dir = 0; dir < DIM; dir++)
        {
            computeFlux(a_fluxes[dir], W_ave, dir);
            a_Rhs += m_divergence[dir](a_fluxes[dir]);
        }
        a_Rhs *= -(a_scale / dx); //Assuming isotropic grid spacing
        Vector Powell_term = forall<double,NUMCOMPS>(Powell,W_ave);
		a_Rhs *= Powell_term;
    }

    static inline void generateTags(
                                  TagData& a_tags,
                                  BoxData<T, NUMCOMPS>& a_state)
  {
    forallInPlace(f_threshold, a_tags, a_state);
  }

private:
  Stencil<T> m_interp_H[DIM];
  Stencil<T> m_interp_L[DIM];
  Stencil<T> m_interp_edge[DIM];
  Stencil<T> m_divergence[DIM];
  Stencil<T> m_laplacian_f[DIM];
};

#endif //end include guard
