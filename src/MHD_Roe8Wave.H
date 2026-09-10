#pragma once

#include "Proto_accel.H"
#include <cmath>

namespace MHD_Riemann_Solvers
{
// Ported from FullSphere b6e94ff, src/MHD_Riemann_Solvers.cpp.
// The original solver uses arithmetic primitive averages (iAveraging = 0)
// and an eight-wave entropy fix. Keep that eigensystem and return
// D = |A| (U_R - U_L), without the factor 1/2 or the centered flux.
// States are [rho, v_normal, v_tangent1, v_tangent2, p, B_n, B_t1, B_t2]
// in CGS units. False requests local Rusanov dissipation for an inadmissible
// input, wave fan, or eigensystem. Disable the fan check only for direct
// regression against the unguarded FullSphere kernel.
inline ACCEL_DECORATION bool roe8WaveDissipation(
    double a_out[8], const double a_W_lo[8], const double a_W_hi[8],
    double a_gamma, double a_entropyFix, bool a_checkIntermediate = true)
{
  double gamma = a_gamma;

  double RhoL, UL, VL, WL, PGasL, BXL, BYL, BZL;
  double RhoR, UR, VR, WR, PGasR, BXR, BYR, BZR;

  RhoL = a_W_lo[0]; RhoR = a_W_hi[0];
  UL = a_W_lo[1]; UR = a_W_hi[1];
  VL = a_W_lo[2]; VR = a_W_hi[2];
  WL = a_W_lo[3]; WR = a_W_hi[3];
  PGasL = a_W_lo[4]; PGasR = a_W_hi[4];
  BXL = a_W_lo[5]; BXR = a_W_hi[5];
  BYL = a_W_lo[6]; BYR = a_W_hi[6];
  BZL = a_W_lo[7]; BZR = a_W_hi[7];
  if (!(RhoL > 0.0 && RhoR > 0.0 && PGasL > 0.0 && PGasR > 0.0))
    return false;

  double RUL, RVL, RWL, RUR, RVR, RWR;
  double PGas, Rho,  U,  V,  W, BX,  BY,  BZ;
  double BxBL, RUxUL, PL, EL, AL2, HL, coefRL;
  double BxBR, RUxUR, PR, ER, AR2, HR, coefRR;
  double BxBH, P, HAV, BXAV, DBX, DBY, DBZ;

  double aas, as, aaf, af, alf, als, c, cc, ccmi, ccma, dfs, X;
  double hyy, hyz, hzz, szb, h, kk, sih, z, del, dell2;
  double b, bb, BX2, BY2, BZ2, BYZ2, sqPiRho, sq2PiRho;
  double Aeigen1, Aeigen2, Aeigen3, Aeigen4, Aeigen5, Aeigen6;
  double Aeigen7, Aeigen8;
  double eigen1,  eigen2,  eigen3,  eigen4,  eigen5,  eigen6;
  double eigen7,  eigen8,  eigmax;
  double du1, du2, du3, du4, du5, du6, du7, du8;
  double lu1, lu2, lu3, lu4, lu5, lu6, lu7, lu8;
  double ra1, ra2, ra3, ra4, ra5, ra6, ra7, ra8;

  double sp11, sp12, sp13, sp14, sp15, sp16, sp17, sp18;
  double sp21, sp22, sp23, sp24, sp25, sp26, sp27, sp28;
  double sp31, sp32, sp33, sp34, sp35, sp36, sp37, sp38;
  double sp41, sp42, sp43, sp44, sp45, sp46, sp47, sp48;
  double sp51, sp52, sp53, sp54, sp55, sp56, sp57, sp58;
  double sp61, sp62, sp63, sp64, sp65, sp66, sp67, sp68;
  double sp71, sp72, sp73, sp74, sp75, sp76, sp77, sp78;
  double sp81, sp82, sp83, sp84, sp85, sp86, sp87, sp88;

  double so11, so12, so13, so14, so15, so16, so17, so18;
  double so21, so22, so23, so24, so25, so26, so27, so28;
  double so31, so32, so33, so34, so35, so36, so37, so38;
  double so41, so42, so43, so44, so45, so46, so47, so48;
  double so51, so52, so53, so54, so55, so56, so57, so58;
  double so61, so62, so63, so64, so65, so66, so67, so68;
  double so71, so72, so73, so74, so75, so76, so77, so78;
  double so81, so82, so83, so84, so85, so86, so87, so88;

  RUL    = RhoL*UL;
  RVL    = RhoL*VL;
  RWL    = RhoL*WL;

  double d_1_8PI = 0.039788735772973833942220940843128591;
  double half = 0.5;
  double one = 1.0;
  double zero = 0.0;
  double two = 2.0;
  double four = 4.0;
  double d_SQRT_2 = 0.707106781186547524400844362104849039;
  double d_1_4PI  = 0.079577471545947667884441881686257181;
  double smallB = 1.0e-9;
  double d_PI     = 3.14159265358979323846264338327950288;
  double d_2PI    = 6.28318530717958647692528676655900576;
  int iAveraging = 0;
  int iLaxFriedrix = 0;
  double del2 = a_entropyFix;
  double hgamma = gamma-1.0;
  BxBL   = d_1_8PI*(BXL*BXL + BYL*BYL + BZL*BZL);
  RUxUL  = half   *(RUL*UL  + RVL*VL  + RWL*WL );

  PL     = PGasL + BxBL;
  EL     = PGasL/hgamma + RUxUL + BxBL;

  AL2    = gamma*PGasL/RhoL;

  RUR    = RhoR*UR;
  RVR    = RhoR*VR;
  RWR    = RhoR*WR;

  BxBR   = d_1_8PI*(BXR*BXR + BYR*BYR + BZR*BZR);
  RUxUR  = half   *(RUR*UR  + RVR*VR  + RWR*WR );

  PR     = PGasR + BxBR;
  ER     = PGasR/hgamma + RUxUR + BxBR;

  AR2    = gamma*PGasR/RhoR;

  if( AL2 > AR2 ) {
      ccmi   = AR2;
      ccma   = AL2;
  } else {
      ccmi   = AL2;
      ccma   = AR2;
  }

  if( iAveraging == 1 ) {
      HL     = (EL + PL)/RhoL;
      HR     = (ER + PR)/RhoR;

      coefRL = one/(one + sqrt( RhoR/RhoL ));
      coefRR = one - coefRL;

      Rho    = RhoR*coefRL + RhoL*coefRR;
      U      = UL  *coefRL + UR  *coefRR;
      V      = VL  *coefRL + VR  *coefRR;
      W      = WL  *coefRL + WR  *coefRR;
      HAV    = HL  *coefRL + HR  *coefRR;

      BX     = BXR *coefRL + BXL *coefRR;
      BY     = BYR *coefRL + BYL *coefRR;
      BZ     = BZR *coefRL + BZL *coefRR;

      BXAV   = half*(BXL  + BXR );

      BX2    = BXAV*BXAV;
      BY2    = BY  *BY;
      BZ2    = BZ  *BZ;

      BYZ2   = BY2 + BZ2;
      BxBH   = d_1_8PI*(BX2 + BYZ2);

      kk     = half*(U*U + V*V + W*W);

      DBX    = BXR - BXL;
      DBY    = BYR - BYL;
      DBZ    = BZR - BZL;

      X      = d_1_8PI*(DBX*DBX + DBY*DBY + DBZ*DBZ)
               * coefRL*coefRL/RhoL;

      cc     = hgamma*(HAV - kk - d_1_4PI*(BX*BXAV + BYZ2)/Rho)
               + (two - gamma)*X;
  } else {
      Rho    = half*(RhoL + RhoR);
      U      = half*(UL   + UR  );
      V      = half*(VL   + VR  );
      W      = half*(WL   + WR  );
      BX     = half*(BXL  + BXR );
      BY     = half*(BYL  + BYR );
      BZ     = half*(BZL  + BZR );
      P      = half*(PL   + PR  );

      BX2    = BX*BX;
      BY2    = BY*BY;
      BZ2    = BZ*BZ;

      BYZ2   = BY2 + BZ2;
      BxBH   = d_1_8PI*(BX2 + BYZ2);

      kk     = half*(U*U + V*V + W*W);

      PGas   = P - BxBH;
      cc     = gamma*PGas/Rho;
  }

  bb     = d_1_4PI*BX2/Rho;
  b      = sqrt( bb );

  if( cc < ccmi ) cc    = ccmi;
  if( cc > ccma ) cc    = ccma;
  if (!(cc > 0.0) || !std::isfinite(cc)) return false;
  c      = sqrt( cc );

  if( BYZ2 < smallB*smallB ) {
      hyy    = d_SQRT_2;
      hzz    = d_SQRT_2;

      if( bb > cc ) {
          aaf    = bb;
          aas    = cc;
          af     = b;
          as     = c;
      } else {
          aaf    = cc;
          aas    = bb;
          af     = c;
          as     = b;
      }
  } else {
      hyz    = sqrt( BYZ2 );
      hyy    = BY/hyz;
      hzz    = BZ/hyz;

      szb    = cc + (BxBH + BxBH)/Rho;
      h      = sqrt( fmax(0.0, szb*szb - four*cc*bb) );
      aas    = half*(szb - h);

      if( aas <= zero ) {
          aas    = zero;
          as     = zero;
      } else {
          if( aas >= cc ) {
              aas    = cc;
              as     = c;
          } else {
              as     = sqrt( aas );
          }
      }

      aaf    = half*(szb + h);
      if( aaf <= cc ) {
          aaf    = cc;
          af     = c;
      } else {
          af     = sqrt( aaf );
      }
  }

  if( iLaxFriedrix == 1 ) {
      eigmax   = fabs( U ) + af;

      Aeigen1  = eigmax;
      Aeigen2  = eigmax;
      Aeigen3  = eigmax;
      Aeigen4  = eigmax;
      Aeigen5  = eigmax;
      Aeigen6  = eigmax;
      Aeigen7  = eigmax;
      Aeigen8  = eigmax;
  } else {
      eigen1 = U;
      eigen2 = U;
      eigen3 = U + b;
      eigen4 = U - b;
      eigen5 = U + af;
      eigen6 = U - af;
      eigen7 = U + as;
      eigen8 = U - as;

      eigmax = fabs( U ) + fabs( V ) + fabs( W ) + af;

      del    = del2*eigmax;
      dell2  = del*del;

      z      = fabs( eigen1 );
      if( z >= del ) {
          Aeigen1  = z;
      } else {
          Aeigen1  = half*(z*z + dell2)/del;
      }

      z      = fabs( eigen2 );
      if( z >= del ) {
          Aeigen2  = z;
      } else {
          Aeigen2  = half*(z*z + dell2)/del;
      }

      z      = fabs( eigen3 );
      if( z >= del ) {
          Aeigen3  = z;
      } else {
          Aeigen3  = half*(z*z + dell2)/del;
      }

      z      = fabs( eigen4 );
      if( z >= del ) {
          Aeigen4  = z;
      } else {
          Aeigen4  = half*(z*z + dell2)/del;
      }

      z      = fabs( eigen5 );
      if( z >= del ) {
          Aeigen5  = z;
      } else {
          Aeigen5  = half*(z*z + dell2)/del;
      }

      z      = fabs( eigen6 );
      if( z >= del ) {
          Aeigen6  = z;
      } else {
          Aeigen6  = half*(z*z + dell2)/del;
      }

      z      = fabs( eigen7 );
      if( z >= del ) {
          Aeigen7  = z;
      } else {
          Aeigen7  = half*(z*z + dell2)/del;
      }

      z      = fabs( eigen8 );
      if( z >= del ) {
          Aeigen8  = z;
      } else {
          Aeigen8  = half*(z*z + dell2)/del;
      }
  }

  dfs    = aaf - aas;

  if( dfs < 1.0e-8 ) {
      alf    = one;
      als    = zero;
  } else {
      if( (cc - aas) <= zero ) {
          alf    = zero;
          als    = one;
      } else {
          if( (aaf - cc) <= zero ) {
              alf    = one;
              als    = zero;
          } else {
              alf    = sqrt( (cc  - aas)/dfs );
              als    = sqrt( (aaf - cc )/dfs );
          }
      }
  }

  sih    =-one;
  if( BX >= zero ) sih   = one;

  du1    = RhoL - RhoR;
  du2    = RUL  - RUR;
  du3    = RVL  - RVR;
  du4    = RWL  - RWR;
  du5    = EL   - ER;
  du6    = BXL  - BXR;
  du7    = BYL  - BYR;
  du8    = BZL  - BZR;

  du5    =-hgamma*(-kk*du1 +   U*du2 +  V*du3 +  W*du4 - du5
                   + (BX*du6 + BY*du7 + BZ*du8)*d_1_4PI);

  du2    = (-U*du1 + du2)/Rho;
  du3    = (-V*du1 + du3)/Rho;
  du4    = (-W*du1 + du4)/Rho;

  sqPiRho  = sqrt( d_PI *Rho );
  sq2PiRho = sqrt( d_2PI*Rho );

  sp11   = one;
  sp12   = zero;
  sp13   = zero;
  sp14   = zero;
  sp15   = Rho*alf;
  sp16   = Rho*alf;
  sp17   = Rho*als;
  sp18   = Rho*als;

  sp21   = zero;
  sp22   = zero;
  sp23   = zero;
  sp24   = zero;
  sp25   = alf*af;
  sp26   =-alf*af;
  sp27   = als*as;
  sp28   =-als*as;

  sp31   = zero;
  sp32   = zero;
  sp33   =-hzz*d_SQRT_2;
  sp34   =-hzz*d_SQRT_2;
  sp35   =-als*as*hyy*sih;
  sp36   = als*as*hyy*sih;
  sp37   = alf*af*hyy*sih;
  sp38   =-alf*af*hyy*sih;

  sp41   = zero;
  sp42   = zero;
  sp43   = hyy*d_SQRT_2;
  sp44   = hyy*d_SQRT_2;
  sp45   =-als*as*hzz*sih;
  sp46   = als*as*hzz*sih;
  sp47   = alf*af*hzz*sih;
  sp48   =-alf*af*hzz*sih;

  sp51   = zero;
  sp52   = zero;
  sp53   = zero;
  sp54   = zero;
  sp55   = alf*Rho*cc;
  sp56   = alf*Rho*cc;
  sp57   = als*Rho*cc;
  sp58   = als*Rho*cc;

  sp61   = zero;
  sp62   = one;
  sp63   = zero;
  sp64   = zero;
  sp65   = zero;
  sp66   = zero;
  sp67   = zero;
  sp68   = zero;

  sp71   = zero;
  sp72   = zero;
  sp73   = hzz*sq2PiRho*sih;
  sp74   =-hzz*sq2PiRho*sih;
  sp75   = two*als*sqPiRho*c*hyy;
  sp76   = two*als*sqPiRho*c*hyy;
  sp77   =-two*alf*sqPiRho*c*hyy;
  sp78   =-two*alf*sqPiRho*c*hyy;

  sp81   = zero;
  sp82   = zero;
  sp83   =-hyy*sq2PiRho*sih;
  sp84   = hyy*sq2PiRho*sih;
  sp85   = two*als*sqPiRho*c*hzz;
  sp86   = two*als*sqPiRho*c*hzz;
  sp87   =-two*alf*sqPiRho*c*hzz;
  sp88   =-two*alf*sqPiRho*c*hzz;

  so11   = one;
  so12   = zero;
  so13   = zero;
  so14   = zero;
  so15   =-one/cc;
  so16   = zero;
  so17   = zero;
  so18   = zero;

  so21   = zero;
  so22   = zero;
  so23   = zero;
  so24   = zero;
  so25   = zero;
  so26   = one;
  so27   = zero;
  so28   = zero;

  so31   = zero;
  so32   = zero;
  so33   =-hzz*d_SQRT_2;
  so34   = hyy*d_SQRT_2;
  so35   = zero;
  so36   = zero;
  so37   = half*hzz*sih/sq2PiRho;
  so38   =-half*hyy*sih/sq2PiRho;

  so41   = zero;
  so42   = zero;
  so43   =-hzz*d_SQRT_2;
  so44   = hyy*d_SQRT_2;
  so45   = zero;
  so46   = zero;
  so47   =-half*hzz*sih/sq2PiRho;
  so48   = half*hyy*sih/sq2PiRho;

  so51   = zero;
  so52   = half*alf*af/cc;
  so53   =-half*als*as*hyy*sih/cc;
  so54   =-half*als*as*hzz*sih/cc;
  so55   = half*alf/(cc*Rho);
  so56   = zero;
  so57   = 0.25*als*hyy/(c*sqPiRho);
  so58   = 0.25*als*hzz/(c*sqPiRho);

  so61   = zero;
  so62   =-half*alf*af/cc;
  so63   = half*als*as*hyy*sih/cc;
  so64   = half*als*as*hzz*sih/cc;
  so65   = half*alf/(cc*Rho);
  so66   = zero;
  so67   = 0.25*als*hyy/(c*sqPiRho);
  so68   = 0.25*als*hzz/(c*sqPiRho);

  so71   = zero;
  so72   = half*als*as/cc;
  so73   = half*alf*af*hyy*sih/cc;
  so74   = half*alf*af*hzz*sih/cc;
  so75   = half*als/(Rho*cc);
  so76   = zero;
  so77   =-0.25*alf*hyy/(c*sqPiRho);
  so78   =-0.25*alf*hzz/(c*sqPiRho);

  so81   = zero;
  so82   =-half*als*as/cc;
  so83   =-half*alf*af*hyy*sih/cc;
  so84   =-half*alf*af*hzz*sih/cc;
  so85   = half*als/(Rho*cc);
  so86   = zero;
  so87   =-0.25*alf*hyy/(c*sqPiRho);
  so88   =-0.25*alf*hzz/(c*sqPiRho);

  // The FullSphere arithmetic eigensystem does not guarantee admissible
  // intermediate states. Walk its conservative wave fan in speed order;
  // reject a fan with nonpositive density or negative thermal energy.
  if (a_checkIntermediate)
  {
    const double left[8][8] = {{so11,so12,so13,so14,so15,so16,so17,so18},
        {so21,so22,so23,so24,so25,so26,so27,so28},
        {so31,so32,so33,so34,so35,so36,so37,so38},
        {so41,so42,so43,so44,so45,so46,so47,so48},
        {so51,so52,so53,so54,so55,so56,so57,so58},
        {so61,so62,so63,so64,so65,so66,so67,so68},
        {so71,so72,so73,so74,so75,so76,so77,so78},
        {so81,so82,so83,so84,so85,so86,so87,so88}};
    const double right[8][8] = {{sp11,sp12,sp13,sp14,sp15,sp16,sp17,sp18},
        {sp21,sp22,sp23,sp24,sp25,sp26,sp27,sp28},
        {sp31,sp32,sp33,sp34,sp35,sp36,sp37,sp38},
        {sp41,sp42,sp43,sp44,sp45,sp46,sp47,sp48},
        {sp51,sp52,sp53,sp54,sp55,sp56,sp57,sp58},
        {sp61,sp62,sp63,sp64,sp65,sp66,sp67,sp68},
        {sp71,sp72,sp73,sp74,sp75,sp76,sp77,sp78},
        {sp81,sp82,sp83,sp84,sp85,sp86,sp87,sp88}};
    const double jump[8] = {du1,du2,du3,du4,du5,du6,du7,du8};
    const double speeds[8] = {U,U,U+b,U-b,U+af,U-af,U+as,U-as};
    const int order[8] = {5,3,7,0,1,6,2,4};
    double state[8] = {RhoL,RUL,RVL,RWL,EL,BXL,BYL,BZL};
    for (int wave = 0; wave < 7; ++wave)
    {
      const int k = order[wave];
      double amplitude = 0.0;
      for (int c = 0; c < 8; ++c) amplitude -= left[k][c]*jump[c];
      double delta[8];
      for (int c = 0; c < 8; ++c) delta[c] = amplitude*right[c][k];
      delta[4] = kk*delta[0] + Rho*(U*delta[1]+V*delta[2]+W*delta[3])
          + delta[4]/hgamma + d_1_4PI*(BX*delta[5]+BY*delta[6]+BZ*delta[7]);
      delta[1] = U*delta[0] + Rho*delta[1];
      delta[2] = V*delta[0] + Rho*delta[2];
      delta[3] = W*delta[0] + Rho*delta[3];
      for (int c = 0; c < 8; ++c) state[c] += delta[c];
      // Coincident waves have no state between them; test only after the group.
      if (speeds[order[wave+1]] > speeds[k])
      {
        if (!(state[0] > 0.0) || !std::isfinite(state[0])) return false;
        const double thermal = state[4]
            - 0.5*(state[1]*state[1]+state[2]*state[2]+state[3]*state[3])/state[0]
            - d_1_8PI*(state[5]*state[5]+state[6]*state[6]+state[7]*state[7]);
        if (!(thermal >= 0.0) || !std::isfinite(thermal)) return false;
      }
    }
  }

  lu1    = Aeigen1*(so11*du1 + so12*du2 + so13*du3 + so14*du4
                    + so15*du5 + so16*du6 + so17*du7 + so18*du8);
  lu2    = Aeigen2*(so21*du1 + so22*du2 + so23*du3 + so24*du4
                    + so25*du5 + so26*du6 + so27*du7 + so28*du8);
  lu3    = Aeigen3*(so31*du1 + so32*du2 + so33*du3 + so34*du4
                    + so35*du5 + so36*du6 + so37*du7 + so38*du8);
  lu4    = Aeigen4*(so41*du1 + so42*du2 + so43*du3 + so44*du4
                    + so45*du5 + so46*du6 + so47*du7 + so48*du8);
  lu5    = Aeigen5*(so51*du1 + so52*du2 + so53*du3 + so54*du4
                    + so55*du5 + so56*du6 + so57*du7 + so58*du8);
  lu6    = Aeigen6*(so61*du1 + so62*du2 + so63*du3 + so64*du4
                    + so65*du5 + so66*du6 + so67*du7 + so68*du8);
  lu7    = Aeigen7*(so71*du1 + so72*du2 + so73*du3 + so74*du4
                    + so75*du5 + so76*du6 + so77*du7 + so78*du8);
  lu8    = Aeigen8*(so81*du1 + so82*du2 + so83*du3 + so84*du4
                    + so85*du5 + so86*du6 + so87*du7 + so88*du8);

  ra1    = sp11*lu1 + sp12*lu2 + sp13*lu3 + sp14*lu4
           + sp15*lu5 + sp16*lu6 + sp17*lu7 + sp18*lu8;
  ra2    = sp21*lu1 + sp22*lu2 + sp23*lu3 + sp24*lu4
           + sp25*lu5 + sp26*lu6 + sp27*lu7 + sp28*lu8;
  ra3    = sp31*lu1 + sp32*lu2 + sp33*lu3 + sp34*lu4
           + sp35*lu5 + sp36*lu6 + sp37*lu7 + sp38*lu8;
  ra4    = sp41*lu1 + sp42*lu2 + sp43*lu3 + sp44*lu4
           + sp45*lu5 + sp46*lu6 + sp47*lu7 + sp48*lu8;
  ra5    = sp51*lu1 + sp52*lu2 + sp53*lu3 + sp54*lu4
           + sp55*lu5 + sp56*lu6 + sp57*lu7 + sp58*lu8;
  ra6    = sp61*lu1 + sp62*lu2 + sp63*lu3 + sp64*lu4
           + sp65*lu5 + sp66*lu6 + sp67*lu7 + sp68*lu8;
  ra7    = sp71*lu1 + sp72*lu2 + sp73*lu3 + sp74*lu4
           + sp75*lu5 + sp76*lu6 + sp77*lu7 + sp78*lu8;
  ra8    = sp81*lu1 + sp82*lu2 + sp83*lu3 + sp84*lu4
           + sp85*lu5 + sp86*lu6 + sp87*lu7 + sp88*lu8;

  ra5    = ra1*kk     + Rho    *(U *ra2 + V *ra3 + W *ra4)
           + ra5/hgamma + d_1_4PI*(BX*ra6 + BY*ra7 + BZ*ra8);
  ra2    = U*ra1 + Rho*ra2;
  ra3    = V*ra1 + Rho*ra3;
  ra4    = W*ra1 + Rho*ra4;
  a_out[0] = -ra1;
  a_out[1] = -ra2;
  a_out[2] = -ra3;
  a_out[3] = -ra4;
  a_out[4] = -ra5;
  a_out[5] = -ra6;
  a_out[6] = -ra7;
  a_out[7] = -ra8;
  for (int comp = 0; comp < 8; ++comp)
    if (!std::isfinite(a_out[comp])) return false;
  return true;
}

// Local Lax-Friedrichs dissipation for nonpositive pressure or an unusable
// eigensystem. The acoustic-plus-Alfven speed bounds the fast-wave speed.
inline ACCEL_DECORATION void rusanovDissipation(
    double a_out[8], const double a_left[8], const double a_right[8],
    double a_gamma)
{
  constexpr double inv4Pi = 0.079577471545947667884441881686257181;
  double conserved[2][8];
  double speed = 0.0;
  for (int side = 0; side < 2; ++side)
  {
    const double* w = side == 0 ? a_left : a_right;
    double v2 = 0.0, b2 = 0.0;
    conserved[side][0] = w[0];
    for (int d = 0; d < 3; ++d)
    {
      v2 += w[1+d]*w[1+d];
      b2 += w[5+d]*w[5+d];
      conserved[side][1+d] = w[0]*w[1+d];
      conserved[side][5+d] = w[5+d];
    }
    conserved[side][4] = w[4]/(a_gamma-1.0)
        + 0.5*w[0]*v2 + 0.5*inv4Pi*b2;
    speed = fmax(speed, fabs(w[1]) + sqrt((a_gamma*w[4] + inv4Pi*b2)/w[0]));
  }
  for (int c = 0; c < 8; ++c)
    a_out[c] = speed*(conserved[1][c] - conserved[0][c]);
}

// A maps the stored contravariant vectors to Cartesian components. Its
// columns are not orthogonal: the normal is a cofactor row, not column dir.
// Return Cartesian dissipation times the physical face-area Jacobian.
inline ACCEL_DECORATION void mappedRoe8WaveDissipation(
    double a_out[8], const double a_left[8], const double a_right[8],
    const double a_A[3][3], double a_radialArea, int a_dir,
    double a_gamma, double a_entropyFix, double a_densityFloor)
{
  const int t1 = (a_dir + 1) % 3;
  const int t2 = (a_dir + 2) % 3;
  double frame[3][3]; // rows: unit normal and two orthonormal tangents
  double normal2 = 0.0, tangent2 = 0.0;
  for (int i = 0; i < 3; ++i)
  {
    const int j = (i + 1) % 3, k = (i + 2) % 3;
    frame[0][i] = a_A[j][t1]*a_A[k][t2] - a_A[k][t1]*a_A[j][t2];
    frame[1][i] = a_A[i][t1];
    normal2 += frame[0][i]*frame[0][i];
    tangent2 += frame[1][i]*frame[1][i];
  }
  const double normalMag = sqrt(normal2);
  const double tangentMag = sqrt(tangent2);
  for (int i = 0; i < 3; ++i)
  {
    frame[0][i] /= normalMag;
    frame[1][i] /= tangentMag;
  }
  for (int i = 0; i < 3; ++i)
  {
    const int j = (i + 1) % 3, k = (i + 2) % 3;
    frame[2][i] = frame[0][j]*frame[1][k] - frame[0][k]*frame[1][j];
  }

  double local[2][8] = {};
  for (int side = 0; side < 2; ++side)
  {
    const double* w = side == 0 ? a_left : a_right;
    local[side][0] = fmax(w[0], a_densityFloor);
    local[side][4] = fmax(w[4], 0.0);
    double vCart[3] = {}, bCart[3] = {};
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
      {
        vCart[i] += a_A[i][j]*w[1+j];
        bCart[i] += a_A[i][j]*w[5+j];
      }
    for (int d = 0; d < 3; ++d)
      for (int i = 0; i < 3; ++i)
      {
        local[side][1+d] += frame[d][i]*vCart[i];
        local[side][5+d] += frame[d][i]*bCart[i];
      }
  }

  double diss[8];
  if (!roe8WaveDissipation(diss, local[0], local[1], a_gamma, a_entropyFix, true))
    rusanovDissipation(diss, local[0], local[1], a_gamma);

  const double area = a_radialArea*normalMag;
  a_out[0] = area*diss[0];
  a_out[4] = area*diss[4];
  for (int i = 0; i < 3; ++i)
  {
    a_out[1+i] = 0.0;
    a_out[5+i] = 0.0;
    for (int d = 0; d < 3; ++d)
    {
      a_out[1+i] += area*frame[d][i]*diss[1+d];
      a_out[5+i] += area*frame[d][i]*diss[5+d];
    }
  }
}
// Factor out the angular mapping before transverse face averaging. Scalars
// are subsequently multiplied by <det(A)>, and contravariant vectors by
// <det(A) A>, using the operator's metric averages and face-product rule.
// The radial face area remains in these coefficients. In a radial hydro
// flow they are angle-independent, so the discrete radial metric identity
// is retained rather than approximated by convolving Cartesian geometry.
inline ACCEL_DECORATION void factoredRoe8WaveDissipation(
    double a_out[8], const double a_left[8], const double a_right[8],
    const double a_A[3][3], double a_radialArea, int a_dir,
    double a_gamma, double a_entropyFix, double a_densityFloor)
{
  double mapped[8], adjugate[3][3];
  mappedRoe8WaveDissipation(mapped,a_left,a_right,a_A,a_radialArea,
      a_dir,a_gamma,a_entropyFix,a_densityFloor);
  for (int d = 0; d < 3; ++d)
  {
    const int t1 = (d + 1) % 3, t2 = (d + 2) % 3;
    for (int i = 0; i < 3; ++i)
    {
      const int j = (i + 1) % 3, k = (i + 2) % 3;
      adjugate[d][i] = a_A[j][t1]*a_A[k][t2] - a_A[k][t1]*a_A[j][t2];
    }
  }
  const double det = a_A[0][0]*adjugate[0][0]
      + a_A[1][0]*adjugate[0][1] + a_A[2][0]*adjugate[0][2];
  a_out[0] = mapped[0]/det;
  a_out[4] = mapped[4]/det;
  for (int d = 0; d < 3; ++d)
  {
    a_out[1+d] = 0.0;
    a_out[5+d] = 0.0;
    for (int i = 0; i < 3; ++i)
    {
      a_out[1+d] += adjugate[d][i]*mapped[1+i]/(det*det);
      a_out[5+d] += adjugate[d][i]*mapped[5+i]/(det*det);
    }
  }
}
} // namespace MHD_Riemann_Solvers
