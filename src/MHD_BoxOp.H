#pragma once
#ifndef _BOX_OP_MHD_
#define _BOX_OP_MHD_

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

#define NUMCOMPS 2 * DIM + 2

using namespace Proto;

typedef BoxData<double> Scalar;
typedef BoxData<double, NUMCOMPS> Vector;
/**
 * @brief Function to covert conserevd variables to primitive variables.
 * @param a_W the output primitive variables.
 * @param a_U the input conserved variables.
 * @param a_gamma gamma.
 */
PROTO_KERNEL_START
void consToPrimF(State &a_W,
				 const State &a_U,
				 double a_gamma)
{
	double rho = a_U(0);
	double v2 = 0.0;
	double B2 = 0.0;
	double gamma = a_gamma;
	a_W(0) = rho;

	for (int i = 1; i <= DIM; i++)
	{
		double v, B;
		v = a_U(i) / rho;
		B = a_U(DIM + 1 + i);
		a_W(i) = v;
		a_W(DIM + 1 + i) = a_U(DIM + 1 + i);
		v2 += v * v;
		B2 += B * B;
	}

	a_W(NUMCOMPS - 1 - DIM) = (a_U(NUMCOMPS - 1 - DIM) - .5 * rho * v2 - B2 / 8.0 / c_PI) * (gamma - 1.0);
}
PROTO_KERNEL_END(consToPrimF, consToPrim)

PROTO_KERNEL_START
void f_thresholdF(
	Var<short> &a_tags,
	Var<double, NUMCOMPS> &a_U)
{
	double thresh = 1.001;
	if (a_U(0) > thresh)
	{
		a_tags(0) = 1;
	}
	else
	{
		a_tags(0) = 0;
	};
};
PROTO_KERNEL_END(f_thresholdF, f_threshold);

/**
 * @brief Function to transfer data from lower dimensional BoxData to a higher dimensional BoxData.
 * @param a_F the output BoxData.
 * @param a_F_temp the input BoxData.
 * @param a_s where to transfer.
 */
PROTO_KERNEL_START
void Fill_flux_calcF(const Point &a_pt,
					 Var<double, DIM, MEM, NUMCOMPS> &a_F,
					 const Var<double, NUMCOMPS> &a_F_temp,
					 const int a_s)
{
	for (int i = 0; i < NUMCOMPS; i++)
	{
		a_F(a_s, i) = a_F_temp(i);
	}
}
PROTO_KERNEL_END(Fill_flux_calcF, Fill_flux_calc)

/**
 * @brief Function to transpose BoxData.
 * @param a_F the output BoxData.
 * @param a_F_temp the input BoxData.
 */
PROTO_KERNEL_START
void Transpose_calcF(const Point &a_pt,
					 Var<double, NUMCOMPS> &a_F,
					 const Var<double, 1, MEM, NUMCOMPS> &a_F_temp)
{
	for (int i = 0; i < NUMCOMPS; i++)
	{
		a_F(i) = a_F_temp(0, i);
	}
}
PROTO_KERNEL_END(Transpose_calcF, Transpose_calc)

template <typename T, MemType MEM = MEMTYPE_DEFAULT>
class BoxOp_MHD : public BoxOp<T, NUMCOMPS, 1, MEM>
{
public:
	using BoxOp<T, NUMCOMPS, 1, MEM>::BoxOp;
	T gamma = inputs.gamma;
	mutable T umax;
	
	// inline static T min_dt() { return 1e-10; }
	// How many ghost cells does the operator need from the state variables
	inline static Point ghost() { return Point::Ones(NGHOST); }

	// How many ghost cells does the operator need from the auxiliary variables
	inline static Point auxGhost() { return Point::Zeros(); }

	// What is the intended order of accuracy of the operator
	inline static constexpr int order() { return 4; }

	MHD_Mapping::mapping_variables m_map_vars;

	// Initialization
	inline void
	init()
	{
		m_laplacian = Stencil<double>::Laplacian();
		m_deconvolve = (-1.0 / 24.0) * m_laplacian + (1.0) * Shift(Point::Zeros());
		m_copy = 1.0 * Shift(Point::Zeros());
		for (int dir = 0; dir < DIM; dir++)
		{
			m_laplacian_f[dir] = Stencil<double>::LaplacianFace(dir);
			m_deconvolve_f[dir] = (-1.0 / 24.0) * m_laplacian_f[dir] + 1.0 * Shift(Point::Zeros());
			m_convolve_f[dir] = (1.0 / 24.0) * m_laplacian_f[dir] + 1.0 * Shift(Point::Zeros());
			m_interp_H[dir] = Stencil<double>::CellToFaceH(dir);
			m_interp_L[dir] = Stencil<double>::CellToFaceL(dir);
			m_interp_edge[dir] = Stencil<double>::CellToFace(dir);
			m_divergence[dir] = Stencil<double>::FluxDivergence(dir);
			m_derivative[dir] = Stencil<double>::Derivative(1, dir, 2);
		}

		HDF5Handler h5;
		using namespace std;
		double a_dx, a_dy, a_dz;
		a_dx = this->dx()[0];
#if DIM == 2
		a_dy = this->dx()[1];
#endif
#if DIM == 3
		a_dy = this->dx()[1];
		a_dz = this->dx()[2];
#endif
		if (DIM == 2)
			a_dz = 1.0;

		Box bxmap = this->box().grow(NGHOST);
		for (int dir = 0; dir < DIM; dir++)
		{
			bxmap = bxmap.extrude(dir);
		}

		MHD_Mapping::Regular_map_filling_func(m_map_vars, bxmap, a_dx, a_dy, a_dz);

	}

	// Apply BCs by filling ghost cells in stage values. For Euler, this is done by calling
	// exchange. For the MHD code, it will be more complicated.
	// The interface is very provisional. We expect it to evolve as we d more real problems.
	inline void bcStage(
		LevelBoxData<T, NUMCOMPS> &a_UStage,
		const LevelBoxData<T, NUMCOMPS> &a_U0,
		int a_stage)
	{
		a_UStage.exchange();
	}

	// Apply Operator
	inline void operator()(
		BoxData<T, NUMCOMPS> &a_Rhs,
		Array<BoxData<T, NUMCOMPS>, DIM> &a_fluxes,
		const BoxData<T, NUMCOMPS> &a_JU_ave,
		T a_scale = 1.0) const
	{
		HDF5Handler h5;
		using namespace std;
		double a_dx, a_dy, a_dz;
		a_dx = this->dx()[0];
#if DIM == 2
		a_dy = this->dx()[1];
#endif
#if DIM == 3
		a_dy = this->dx()[1];
		a_dz = this->dx()[2];
#endif
		if (DIM == 2)
			a_dz = 1.0;
		double volume = a_dx * a_dy * a_dz;
		double dxd[3] = {a_dx, a_dy, a_dz};

		double dt_new;

		Box dbx0 = a_JU_ave.box();
		Box dbx1 = dbx0;
		Box dbx2 = dbx0.grow(0 - NGHOST);
		a_Rhs.setVal(0.0);
		auto a_U_ave = Operator::_cellTensorQuotient(a_JU_ave, m_map_vars.m_J, a_JU_ave, m_map_vars.m_J);
		Vector W_bar = forall<double, NUMCOMPS>(consToPrim, a_U_ave, gamma);
		Vector U = m_deconvolve(a_U_ave);
		Vector W = forall<double, NUMCOMPS>(consToPrim, U, gamma);
		Vector W_ave = m_laplacian(W_bar, 1.0 / 24.0);
		W_ave += W;
		// h5.writePatch(1,W_ave,"W_ave");
		for (int d = 0; d < DIM; d++)
		{
			Vector W_ave_low_temp(dbx0), W_ave_high_temp(dbx0);
			Vector W_ave_low(dbx0), W_ave_high(dbx0);
			W_ave_low_temp = m_interp_L[d](W_ave);
			W_ave_high_temp = m_interp_H[d](W_ave);
			Vector W_ave_edge = m_interp_edge[d](W_ave);
			MHD_Limiters::MHD_Limiters_4O(W_ave_low, W_ave_high, W_ave_low_temp, W_ave_high_temp, W_ave, W_bar, d, a_dx, a_dy, a_dz);
			Vector W_low = m_deconvolve_f[d](W_ave_low);
			Vector W_high = m_deconvolve_f[d](W_ave_high);
			BoxData<double, DIM, MEM, NUMCOMPS> F_f(dbx1);
			BoxData<double, DIM, MEM, NUMCOMPS> F_ave_f(dbx1);
			BoxData<double, NUMCOMPS> F_f_temp(dbx1);
			Vector F_f_mapped(dbx1);
			F_f_mapped.setVal(0.0);
			BoxData<double, DIM, MEM, NUMCOMPS> BF_f(dbx1), BF_ave_f(dbx1);
			Vector BF_f_mapped(dbx1);
			double dx_d = dxd[d];
			for (int s = 0; s < DIM; s++)
			{
				if (inputs.Riemann_solver_type == 1)
				{
					MHD_Riemann_Solvers::Rusanov_Solver(F_f_temp, W_low, W_high, s, gamma);
				}
				if (inputs.Riemann_solver_type == 2)
				{
					MHD_Riemann_Solvers::Roe8Wave_Solver(F_f_temp, W_low, W_high, s, gamma);
				}
				forallInPlace_p(Fill_flux_calc, F_f, F_f_temp, s);
			}
			F_ave_f = m_convolve_f[d](F_f);
			BoxData<double, 1, MEM, NUMCOMPS> fluxdir = Operator::_faceMatrixProductATB(m_map_vars.m_NT[d], F_ave_f, m_map_vars.m_NT[d], F_ave_f, d);
			forallInPlace_p(Transpose_calc, F_f_mapped, fluxdir);
			// h5.writePatch(1,F_f_mapped,"F_f_mapped"+to_string(d));
			Vector Rhs_d = m_divergence[d](F_f_mapped);
			Rhs_d *= -1. / volume;
			a_Rhs += Rhs_d;
		}
	}

	static inline void generateTags(
		TagData &a_tags,
		BoxData<T, NUMCOMPS> &a_state)
	{
		forallInPlace(f_threshold, a_tags, a_state);
	}

private:
	Stencil<T> m_interp_H[DIM];
	Stencil<T> m_interp_L[DIM];
	Stencil<T> m_divergence[DIM];
	Stencil<T> m_laplacian_f[DIM];
	Stencil<T> m_laplacian;
	Stencil<T> m_deconvolve;
	Stencil<T> m_copy;
	Stencil<T> m_deconvolve_f[DIM];
	Stencil<T> m_convolve_f[DIM];
	Stencil<T> m_interp_edge[DIM];
	Stencil<T> m_derivative[DIM];
};

#endif // end include guard
