#if TURB == 1
#include "Proto.H"
#include "MHDOp.H"
#include "MHD_Mapping.H"
#include "MHD_Turbulence.H"
#include "MHD_Input_Parsing.H"
#include "MHD_Output_Writer.H"
#include "MHD_Constants.H"
extern Parsefrominputs inputs;
// constexpr MemType MEM = MEMTYPE_DEFAULT;
typedef BoxData<double,1,HOST> Scalar;
typedef BoxData<double,NUMCOMPS,HOST> Vector;
/// @brief MHD_Turbulence namespace
namespace MHD_Turbulence {


    PROTO_KERNEL_START
	void Turb_Flux_calcF(const Point &a_pt,
                         Var<double, NUMCOMPS> &a_F,
                         Var<double, NUMCOMPS> &a_W_lo,
                         Var<double, NUMCOMPS> &a_W_hi,
                         const int a_dir)
	{
		int inorm  = iVX +      a_dir;
        int itan1  = iVX + (( a_dir + 1) % DIM);
        int itan2  = iVX + (( a_dir + 2) % DIM);

        int inormB = iBX   +      a_dir;
        int itanB1 = iBX   + (( a_dir + 1) % DIM);
        int itanB2 = iBX   + (( a_dir + 2) % DIM);

        double sigmaDTM = inputs.SigmaD;  // Check dimensions
        double C1       = (3.0 + sigmaDTM)/12.0;
        double C2       = sigmaDTM/3.0;

        double RAV    = 0.5*(a_W_lo(iRHO ) + a_W_hi(iRHO ));
        double UAV    = 0.5*(a_W_lo(inorm) + a_W_hi(inorm));
        double VAV    = 0.5*(a_W_lo(itan1) + a_W_hi(itan1));
        double WAV    = 0.5*(a_W_lo(itan2) + a_W_hi(itan2));

        double BXAV   = 0.5*(a_W_lo(inormB) + a_W_hi(inormB));
        double BYAV   = 0.5*(a_W_lo(itanB1) + a_W_hi(itanB1));
        double BZAV   = 0.5*(a_W_lo(itanB2) + a_W_hi(itanB2));

        double Z2AV   = 0.5*(a_W_lo(iZ2) + a_W_hi(iZ2));
        double SCAV   = 0.5*(a_W_lo(iSIGMA) + a_W_hi(iSIGMA));

        double RZ2    = RAV*Z2AV;
        double FL     = C1*RZ2;

        double BB     = BXAV*BXAV + BYAV*BYAV + BZAV*BZAV;

        if( BB > 1.0e-10 ){
            double UB     = BXAV*UAV + BYAV*VAV + BZAV*WAV;

            double RZ2B     = RZ2*BXAV;
            double CRZ2B_B  = C2*RZ2B/BB;

            a_F(inorm) = a_F(inorm) - CRZ2B_B*BXAV;
            a_F(itan1) = a_F(itan1) - CRZ2B_B*BYAV;
            a_F(itan2) = a_F(itan2) - CRZ2B_B*BZAV;
            a_F(iP) = a_F(iP) - (C2*UB/BB + SCAV/(4.0*sqrt(M_PI*RAV)))*RZ2B;
        }

        double RhoU   = a_F(iRHO);

        double Rho, Z2, SC, LM;
        if( RhoU > 0.0 ){
            Rho    = a_W_lo(iRHO);
            Z2     = a_W_lo(iZ2);
            SC     = a_W_lo(iSIGMA);
            LM     = a_W_lo(iLAMBDA);
        } else {
            Rho    = a_W_hi(iRHO);
            Z2     = a_W_hi(iZ2);
            SC     = a_W_hi(iSIGMA);
            LM     = a_W_hi(iLAMBDA);
        }

        a_F(inorm) = a_F(inorm) + C1*Rho *Z2;
        a_F(iP) = a_F(iP) + C1*RhoU*Z2;

        a_F(iRHOZ2)  = RhoU*Z2;
        a_F(iRHOZ2SIGMA)  = RhoU*Z2*SC;
        a_F(iRHOLAMBDA)  = RhoU*LM;
	}
	PROTO_KERNEL_END(Turb_Flux_calcF, Turb_Flux_calc)



    void Turb_Flux(BoxData<double,NUMCOMPS>& a_F,
			  BoxData<double,NUMCOMPS>& a_W_lo,
              BoxData<double,NUMCOMPS>& a_W_hi,
			  const int a_dir)
    {
        forallInPlace_p(Turb_Flux_calc, a_F, a_W_lo, a_W_hi, a_dir);
    }


    PROTO_KERNEL_START
    void Turb_Source_calcF(const Point &a_pt,
                         Var<double, NUMCOMPS> &a_S,
                         Var<double, NUMCOMPS> &a_W,
                         Var<double, 1> &a_divV)
    {
        double smallTM = 1.0e-8;
        double alphaTM = 0.8;
        double betaTM  = 0.4;
        double sigmaDTM    = inputs.SigmaD;
        double C3       = 0.5*(1.0 + sigmaDTM);
    
        double Rho    =               a_W(iRHO);
        double Z2     = max( smallTM, a_W(iZ2) );
        double SC     =               a_W(iSIGMA);
        double Lm     = max( smallTM, a_W(iLAMBDA) );

        double S2     = SC*SC;
        double H1     = sqrt( (1.0 - S2)*(1.0 + SC) );
        double H2     = sqrt( (1.0 - S2)*(1.0 - SC) );
        double FP     = 0.5*(H1 + H2);
        double FM     = 0.5*(H1 - H2);
        
        double RhoZ   = Rho*sqrt( Z2 );
        H1     = alphaTM*RhoZ*Z2/Lm;

        a_S(iRHO) = 0.0;
        a_S(iMOMX) = 0.0;
        a_S(iMOMY) = 0.0;
        a_S(iMOMZ) = 0.0;
        a_S(iE) = 0.5*FP*H1;
        a_S(iBX) = 0.0;
        a_S(iBY) = 0.0;
        a_S(iBZ) = 0.0;
        a_S(iRHOZ2)  = (-FP*H1 - C3  *Rho*Z2*a_divV(0));
        a_S(iRHOZ2SIGMA)  =-( FM*H1 + 0.5*Rho*Z2*a_divV(0)*SC);
        a_S(iRHOLAMBDA)  = betaTM*FP*RhoZ;        
    }
    PROTO_KERNEL_END(Turb_Source_calcF, Turb_Source_calc)


    void Turb_Source(BoxData<double,NUMCOMPS>& a_S,
			  BoxData<double,NUMCOMPS>& a_W,
              BoxData<double,1>& a_divV){
        forallInPlace_p(Turb_Source_calc, a_S, a_W, a_divV);
    }

}
#endif