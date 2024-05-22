#include "Proto.H"
#include "MHD_Mapping.H"
#include "MHDOp.H"
#include "MHD_Initialize.H"
#include "MHD_Output_Writer.H"
// #include "Proto_WriteBoxData.H"
#include "MHD_Constants.H"
#include <algorithm>    // std::min

typedef BoxData<double,1> Scalar;
typedef BoxData<double,NUMCOMPS> Vector;
/// @brief MHD_CFL namespace
namespace MHD_CFL {

    PROTO_KERNEL_START
	void lambdacalcF(Var<double,1>& a_lambda,
	                 const State& a_W_ave,
	                 int a_d,
	                 double a_gamma)
	{
		double gamma = a_gamma;
		double rho=0., u=0., v=0., w=0., p=0., Bx=0., By=0., Bz=0., ce, af, B_mag, u_mag, Bdir, udir;

		rho = a_W_ave(0);
		u   = a_W_ave(1);
		v   = a_W_ave(2);
		w   = a_W_ave(3);
		p   = a_W_ave(4);
		Bx  = a_W_ave(5);
		By  = a_W_ave(6);
		Bz  = a_W_ave(7);

		if (a_d == 0) {
			Bdir = Bx;
			udir = u;
		};
		if (a_d == 1) {
			Bdir = By;
			udir = v;
		};
		if (a_d == 2) {
			Bdir = Bz;
			udir = w;
		};

		if (p < 0.0) p = 0.0;
		ce = sqrt(gamma*p/rho);
		B_mag = sqrt(Bx*Bx+By*By+Bz*Bz);
		af = sqrt(ce*ce + B_mag*B_mag/4.0/c_PI/rho);
		u_mag = sqrt(u*u+v*v+w*w);
		// a_lambda(0) = af + u_mag;
		a_lambda(0) = af + abs(udir);

	}
	PROTO_KERNEL_END(lambdacalcF, lambdacalc)


    PROTO_KERNEL_START
	void dt_dcalcF(Var<double,1>& a_dt_d,
	               Var<double,1>& a_Lambda_f,
                   V& a_x,
                   V& a_x_ahead)
	{
        double dx_d;
        #if DIM == 2
            dx_d = sqrt((a_x_ahead(0)-a_x(0))*(a_x_ahead(0)-a_x(0)) + (a_x_ahead(1)-a_x(1))*(a_x_ahead(1)-a_x(1)));
        #endif
        #if DIM == 3
            dx_d = sqrt((a_x_ahead(0)-a_x(0))*(a_x_ahead(0)-a_x(0)) + (a_x_ahead(1)-a_x(1))*(a_x_ahead(1)-a_x(1)) + (a_x_ahead(2)-a_x(2))*(a_x_ahead(2)-a_x(2)));
        #endif
		a_dt_d(0) = dx_d/a_Lambda_f(0);
	}
	PROTO_KERNEL_END(dt_dcalcF, dt_dcalc)


	void Min_dt_calc_func(double& a_dt,
	                    const BoxData<double,NUMCOMPS>& a_W_ave,
						const Box a_dbx1,
                        const double a_dx,
                        const double a_dy,
                        const double a_dz,
	                    const double a_gamma)
	{   
        double dt[DIM];
        for (int dir = 0; dir < DIM; dir++)
		{
		    Scalar Lambda_f = forall<double>(lambdacalc, a_W_ave, dir, a_gamma);
            Box dbx0 = a_dbx1;
            BoxData<double,DIM> x(dbx0);		
			MHD_Mapping::get_cart_coords_fc(x, dbx0, a_dx, a_dy, a_dz, dir);
            BoxData<double,DIM> x_ahead = alias(x,Point::Basis(dir)*(-1));
            Scalar dt_d = forall<double>(dt_dcalc, Lambda_f, x, x_ahead);
            dt[dir] = dt_d.min();
        }
		double a_dt_temp = min(dt[0], dt[1]);
		a_dt = min(a_dt_temp, dt[2]);
	}

}
