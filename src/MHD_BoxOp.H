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
#include "PolarExchangeCopier.H"

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

/**
 * @brief Function to covert conserevd variables to primitive variables, used in spherical mapping.
 * @param a_W_sph the output primitive variables.
 * @param a_U_sph the input conserved variables, but the v and B vectors are scaled with row magnitudes of A matrix (see overleaf document).
 * @param a_U_sph_actual the input conserved variables.
 * @param a_gamma gamma.
 */
PROTO_KERNEL_START
void consToPrimSphF(State &a_W_sph,
					const State &a_U_sph,
					const State &a_U_sph_actual,
					double a_gamma)
{
	double rho = a_U_sph(0);
	double v2 = 0.0;
	double B2 = 0.0;
	double gamma = a_gamma;
	a_W_sph(0) = rho;

	for (int i = 1; i <= DIM; i++)
	{
		double v, v_actual, B, B_actual;
		v = a_U_sph(i) / rho;
		v_actual = a_U_sph_actual(i) / rho;
		B = a_U_sph(DIM + 1 + i);
		B_actual = a_U_sph_actual(DIM + 1 + i);
		a_W_sph(i) = v;
		a_W_sph(DIM + 1 + i) = a_U_sph(DIM + 1 + i);
		v2 += v_actual * v_actual;
		B2 += B_actual * B_actual;
	}

	a_W_sph(NUMCOMPS - 1 - DIM) = (a_U_sph(NUMCOMPS - 1 - DIM) - .5 * rho * v2 - B2 / 8.0 / c_PI) * (gamma - 1.0);
	// a_W_sph(NUMCOMPS-1-DIM) = a_U_sph(NUMCOMPS-1-DIM);
}
PROTO_KERNEL_END(consToPrimSphF, consToPrimSph)

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

/**
 * @brief Function to scale flux with face area.
 * @param a_F_scaled the output scaled flux.
 * @param a_F the input flux.
 * @param a_face_area the input face area.
 * @param a_d direction index.
 */
PROTO_KERNEL_START
void Scale_with_A_Ff_calcF(const Point &a_pt,
						   Var<double, NUMCOMPS> &a_F_scaled,
						   Var<double, NUMCOMPS> &a_F,
						   Var<double, DIM> &a_face_area,
						   int a_d)
{
	double area = a_face_area(a_d);

	for (int i = 0; i < NUMCOMPS; i++)
	{
		a_F_scaled(i) = -a_F(i) * area;
	}
}
PROTO_KERNEL_END(Scale_with_A_Ff_calcF, Scale_with_A_Ff_calc)

/**
 * @brief Function to scale flux at a face with cell volume.
 * @param a_F_scaled the output scaled flux.
 * @param a_F the input flux.
 * @param a_cell_volume the input cell volume.
 */
PROTO_KERNEL_START
void Scale_with_V_calcF(const Point &a_pt,
						Var<double, NUMCOMPS> &a_F_scaled,
						Var<double, NUMCOMPS> &a_F,
						Var<double, 1> &a_cell_volume)
{
	double volume = a_cell_volume(0);

	for (int i = 0; i < NUMCOMPS; i++)
	{
		a_F_scaled(i) = a_F(i) / volume;
	}
}
PROTO_KERNEL_END(Scale_with_V_calcF, Scale_with_V_calc)

PROTO_KERNEL_START
void QuotientF(State &a_RHSbyJ,
			   const State &a_RHS,
			   const Var<double, 1> &a_J)
{
	for (int i = 0; i < NUMCOMPS; i++)
	{
		a_RHSbyJ(i) = a_RHS(i) / a_J(0);
	}
}
PROTO_KERNEL_END(QuotientF, Quotient)

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

	mutable MHD_Mapping::mapping_variables m_map_vars;

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

		// MHD_Mapping::Regular_map_filling_func(m_map_vars, bxmap, a_dx, a_dy, a_dz);

		if (inputs.grid_type_global == 2)
		{
			MHD_Mapping::Spherical_2O_map_filling_func(m_map_vars, bxmap, a_dx, a_dy, a_dz);
			MHD_Mapping::Spherical_map_filling_func(m_map_vars, bxmap, a_dx, a_dy, a_dz);
		}
		else
		{
			MHD_Mapping::Regular_map_filling_func(m_map_vars, bxmap, a_dx, a_dy, a_dz);
		}
	}

	// Apply BCs by filling ghost cells in stage values. For Euler, this is done by calling
	// exchange. For the MHD code, it will be more complicated.
	// The interface is very provisional. We expect it to evolve as we d more real problems.
	inline void bcStage(
		LevelBoxData<T, NUMCOMPS> &a_UStage,
		const LevelBoxData<T, NUMCOMPS> &a_U0,
		int a_stage)
	{


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


		if (inputs.grid_type_global == 2)
		{
			a_UStage.template defineExchange<PolarExchangeCopier>(2, 1);
		}
		a_UStage.exchange();

		// DisjointBoxLayout layout = a_UStage.layout();

		// if (inputs.LowBoundType != 0 || inputs.HighBoundType != 0)
		// {
		// 	if (inputs.Spherical_2nd_order == 0)
		// 	{
		// 		MHD_Set_Boundary_Values::Set_Boundary_Values(a_UStage, a_dx, a_dy, a_dz);
		// 	}
		// 	if (inputs.Spherical_2nd_order == 1)
		// 	{
		// 		MHD_Set_Boundary_Values::Set_Boundary_Values_Spherical_2O(a_UStage, a_dx, a_dy, a_dz);
		// 	}
		// }
	}

	// Apply Operator
	inline void operator()(
		BoxData<T, NUMCOMPS> &a_Rhs,
		Array<BoxData<T, NUMCOMPS>, DIM> &a_fluxes,
		BoxData<T, NUMCOMPS> &a_JU_ave,
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

		if (inputs.grid_type_global == 2)
		{
			if (inputs.Spherical_2nd_order == 0)
			{

				HDF5Handler h5;
				Box dbx0 = a_JU_ave.box();
				Box dbx1 = dbx0.grow(0 - NGHOST);
				a_Rhs.setVal(0.0);
				Vector a_U_Sph_ave(dbx0), a_U_Sph_actual_ave(dbx0), a_U_cart_ave(dbx0), a_U_cart_ave_temp(dbx0), Powell_term(dbx0), F_ave_f(dbx0);

				MHD_Mapping::JU_to_U_Sph_ave_calc_func(a_U_Sph_ave, a_JU_ave, m_map_vars.m_detAA_inv_avg, m_map_vars.m_r2rdot_avg, m_map_vars.m_detA_avg, m_map_vars.m_A_row_mag_avg, false, 4);
				MHD_Mapping::JU_to_U_Sph_ave_calc_func(a_U_Sph_actual_ave, a_JU_ave, m_map_vars.m_detAA_inv_avg, m_map_vars.m_r2rdot_avg, m_map_vars.m_detA_avg, m_map_vars.m_A_row_mag_avg, true, 4);
				MHD_Mapping::JU_to_U_ave_calc_func(a_U_cart_ave, a_JU_ave, m_map_vars.m_r2rdot_avg, m_map_vars.m_detA_avg);
				MHD_Mapping::Correct_V_theta_phi_at_poles(a_U_Sph_ave, a_dx, a_dy, a_dz);
				MHD_Mapping::Correct_V_theta_phi_at_poles(a_U_Sph_actual_ave, a_dx, a_dy, a_dz);
				Vector W_bar = forall<double, NUMCOMPS>(consToPrimSph, a_U_Sph_ave, a_U_Sph_actual_ave, gamma);
				Vector W_bar_actual = forall<double, NUMCOMPS>(consToPrimSph, a_U_Sph_actual_ave, a_U_Sph_actual_ave, gamma);

				Vector U = m_deconvolve(a_U_Sph_ave);
				Vector U_actual = m_deconvolve(a_U_Sph_actual_ave);
				Vector W = forall<double, NUMCOMPS>(consToPrimSph, U, U_actual, gamma);
				Vector W_ave = m_laplacian(W_bar, 1.0 / 24.0);
				W_ave += W;
				Vector W_actual = forall<double, NUMCOMPS>(consToPrimSph, a_U_Sph_actual_ave, a_U_Sph_actual_ave, gamma);

				for (int d = 0; d < DIM; d++)
				{
					Vector W_ave_low_temp(dbx0), W_ave_high_temp(dbx0);
					Vector W_ave_low(dbx0), W_ave_high(dbx0);
					Vector W_ave_low_actual(dbx0), W_ave_high_actual(dbx0);

					// W_ave_low_temp = m_interp_L[d](W_ave);
					// W_ave_high_temp = m_interp_H[d](W_ave);
					// MHD_Limiters::MHD_Limiters_4O(W_ave_low,W_ave_high,W_ave_low_temp,W_ave_high_temp,W_ave,W_bar,d,a_dx, a_dy, a_dz);

					MHD_Limiters::MHD_Limiters_minmod(W_ave_low, W_ave_high, W_bar, m_map_vars.m_x_sph_cc, m_map_vars.m_dx_sph, d);

					MHD_Mapping::W_Sph_to_W_normalized_sph(W_ave_low_actual, W_ave_low, m_map_vars.m_A_row_mag_face_avg[d], d);
					MHD_Mapping::W_Sph_to_W_normalized_sph(W_ave_high_actual, W_ave_high, m_map_vars.m_A_row_mag_face_avg[d], d);
					// Vector W_low_boxed(dbx1);
					// W_ave_low_actual.copyTo(W_low_boxed);
					// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, W_low_boxed, "W_low_4O_"+to_string(d));
					// F_ave_f.setVal(0.0);

					double dx_d = dxd[d];
					MHD_Riemann_Solvers::Spherical_Riemann_Solver(F_ave_f, W_ave_low, W_ave_high, W_ave_low_actual, W_ave_high_actual, m_map_vars.m_Dr_detA_avg[d], m_map_vars.m_Dr_detA_A_avg[d], m_map_vars.m_Dr_AdjA_avg[d], m_map_vars.m_A_row_mag_face_avg[d], d, gamma, a_dx, a_dy, a_dz, 4);
					// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, F_ave_f, "F_ave_f_"+to_string(d));
					Vector Rhs_d = m_divergence[d](F_ave_f);
					Rhs_d *= -1. / dx_d;
					a_Rhs += Rhs_d;
				}

				// Vector RHSbyJ = forall<double,NUMCOMPS>(Quotient,a_Rhs,m_map_vars.m_Jacobian_ave);
				// Vector RHS(dbx1);
				// RHSbyJ.copyTo(RHS);
				// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, RHS, "RHS_4");
			}

			if (inputs.Spherical_2nd_order == 1)
			{
				Box dbx0 = a_JU_ave.box();
				Box dbx1 = dbx0.grow(0 - NGHOST);
				Vector F_f_sph(dbx0), F_f(dbx0), F_f_scaled(dbx0), Rhs_d(dbx0), RhsV(dbx0);
				RhsV.setVal(0.0);

				MHDOp::Fix_negative_P(a_JU_ave, inputs.gamma);
				Vector W_low_temp(dbx0), W_high_temp(dbx0), W_low(dbx0), W_high(dbx0);
				BoxData<double, NUMCOMPS> W_sph(dbx0);
				Vector W_cart = forall<double, NUMCOMPS>(consToPrim, a_JU_ave, gamma);
				MHD_Mapping::Cartesian_to_Spherical(W_sph, W_cart, m_map_vars.m_x_sph_cc);
				MHD_Mapping::Correct_V_theta_phi_at_poles(W_sph, a_dx, a_dy, a_dz);
				HDF5Handler h5;

				// Vector W_boxed(dbx1);
				// W_sph.copyTo(W_boxed);
				// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, W_boxed, "W_2");

				for (int d = 0; d < DIM; d++)
				{
					MHD_Limiters::MHD_Limiters_minmod(W_low, W_high, W_sph, m_map_vars.m_x_sph_cc, m_map_vars.m_dx_sph, d);
					// Vector W_low_boxed(dbx1);
					// W_low.copyTo(W_low_boxed);
					// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, W_low_boxed, "W_low_2O_"+to_string(d));
					if (inputs.Riemann_solver_type == 1)
						MHD_Riemann_Solvers::Rusanov_Solver(F_f_sph, W_low, W_high, d, gamma);
					if (inputs.Riemann_solver_type == 2)
						MHD_Riemann_Solvers::Roe8Wave_Solver(F_f_sph, W_low, W_high, d, gamma);
					if (d == 0)
						MHD_Mapping::Spherical_to_Cartesian(F_f, F_f_sph, m_map_vars.m_x_sph_fc_1);
					if (d == 1)
						MHD_Mapping::Spherical_to_Cartesian(F_f, F_f_sph, m_map_vars.m_x_sph_fc_2);
					if (d == 2)
						MHD_Mapping::Spherical_to_Cartesian(F_f, F_f_sph, m_map_vars.m_x_sph_fc_3);
					// Vector F_f_boxed(dbx1);
					// F_f.copyTo(F_f_boxed);
					// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, F_f_boxed, "F_f_2O_"+to_string(d));
					forallInPlace_p(Scale_with_A_Ff_calc, F_f_scaled, F_f, m_map_vars.m_face_area, d);
					Rhs_d = m_divergence[d](F_f_scaled);
					RhsV += Rhs_d;
				}
				forallInPlace_p(Scale_with_V_calc, a_Rhs, RhsV, m_map_vars.m_cell_volume);

				// Vector RHS(dbx1);
				// a_State.m_divB.copyTo(RHS);
				// a_Rhs.copyTo(RHS);
				// if (procID() == 0) h5.writePatch({"density","Vx","Vy","Vz", "p","Bx","By","Bz"}, 1, RHS, "RHS_2");
			}
		}
		else
		{

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
