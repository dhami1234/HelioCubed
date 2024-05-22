#include "Proto.H"
#include "MHD_Limiters.H"
#include "CommonTemplates.H"
#include "Proto_Timer.H"
// #include "Proto_WriteBoxData.H"
#include "MHD_Input_Parsing.H"
#include "MHD_Constants.H"
#include <iomanip>

typedef BoxData<double,1,HOST> Scalar;
typedef BoxData<double,NUMCOMPS,HOST> Vector;

extern Parsefrominputs inputs;
/// @brief MHD_Limiters namespace
namespace MHD_Limiters {

	PROTO_KERNEL_START
	void minmod_calcF(
		const Point& a_pt,
		State& W_low_ahead_limited,
		State& W_high_limited,
		const State& a_W,
		const State& a_dWL,
		const State& a_dWR,
		Var<double,DIM>& a_dx_sph,
		Var<double,DIM>& a_dx_sph_behind,
		Var<double,DIM>& a_dx_sph_ahead,
		int a_d)
	{
		double c1     = 1.0/(a_dx_sph(a_d) + a_dx_sph_ahead(a_d));
		double c2     = 1.0/(a_dx_sph(a_d) + a_dx_sph_behind(a_d));

		for (int i = 0; i< NUMCOMPS; i++){ 
			double WP     = a_W(i);
			double dWR    = a_dWR(i);
			double dWL    = a_dWL(i);

			// SLOPE  = dx*minmod( dWR*c1, dWL*c2 )

			// double signA = /abs(dWR*c1);
			double minmod;
			double A = dWR*c1;
			double B = dWL*c2;
			if (A*B < 0){
				minmod = 0;
			} else if (abs(A) < abs(B)){
			 	minmod = A;
			} else {
				minmod = B;
			}

			// double minmod = signA*max(0.0,min(abs(dWR*c1),signA*dWL*c2));
			double SLOPE = a_dx_sph(a_d)*minmod;

			W_low_ahead_limited(i) = WP + SLOPE;
			W_high_limited(i) = WP - SLOPE;

			// Reduce order if density or p become negative at faces
			// if (i==0 || i==4){
			if (i==4){
				if (W_low_ahead_limited(i) < 0) W_low_ahead_limited(i) = 0.;//WP;
				if (W_high_limited(i) < 0) W_high_limited(i) = 0.;//WP;
			}
		}
	}
	PROTO_KERNEL_END(minmod_calcF, minmod_calc)




	void MHD_Limiters_minmod(BoxData<double,NUMCOMPS>& a_W_low_lim,
	                  BoxData<double,NUMCOMPS>& a_W_high_lim,
	                  BoxData<double,NUMCOMPS>& a_W,
	                  BoxData<double,DIM>& a_x_sph,
	                  BoxData<double,DIM>& a_dx_sph,
	                  const int a_d)
	{
		static Stencil<double> m_difference;
		m_difference = -1.0*Shift(Point::Zeros()) + 1.0*Shift(Point::Basis(a_d)); 
		static Stencil<double> m_copy;
		m_copy = 1.0*Shift(Point::Zeros());

		Vector W_low_ahead = alias(a_W,Point::Basis(a_d)*(-1));
		Vector W_low_ahead_limited = m_copy(W_low_ahead);
		Vector W_high_limited = m_copy(a_W);

		if (inputs.limiter_apply == 1) {
			Vector dWR = m_difference(a_W);
			Vector dWL = alias(dWR,Point::Basis(a_d)*(1));
			BoxData<double,DIM> a_dx_sph_behind = alias(a_dx_sph,Point::Basis(a_d)*(1));
			BoxData<double,DIM> a_dx_sph_ahead = alias(a_dx_sph,Point::Basis(a_d)*(-1));
			forallInPlace_p(minmod_calc,W_low_ahead_limited, W_high_limited, a_W, dWL,
			                 dWR, a_dx_sph, a_dx_sph_behind, a_dx_sph_ahead, a_d);
		}
		a_W_high_lim = alias(W_high_limited);
		a_W_low_lim = alias(W_low_ahead_limited,Point::Basis(a_d)*(1));		
	}
}
