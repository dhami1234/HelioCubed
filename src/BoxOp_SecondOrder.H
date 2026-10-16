#pragma once

// Face i is the low face of cell i. Reconstruction is linear in the uniform
// mapped coordinates; the geometry is evaluated at the corresponding face.
template<typename T, typename MAP, MemType MEM>
inline void BoxOp_EulerCubedSphere<T,MAP,MEM>::reconstructPLM(
    State& left, State& right, const State& w, const Box& faces, int dir) const
{
    PR_TIMERS("HelioCubed::PLM");
    left.define(faces); right.define(faces);
    State wm2=(1.0*Shift::Basis(dir,-2))(w,faces);
    State wm1=(1.0*Shift::Basis(dir,-1))(w,faces);
    State wp1=(1.0*Shift::Basis(dir,1))(w,faces);
    PROTO_ASSERT((wm2.box() & wm1.box() & w.box() & wp1.box()).containsBox(faces),
        "PLM requires two normal ghost cells");
    const int limiter=ParseInputs::get_limiter();
    forallInPlace([] PROTO_LAMBDA(
        Var<T,NUMCOMPS,MEM>& lo, Var<T,NUMCOMPS,MEM>& hi,
        const Var<T,NUMCOMPS,MEM>& mm, const Var<T,NUMCOMPS,MEM>& m,
        const Var<T,NUMCOMPS,MEM>& c, const Var<T,NUMCOMPS,MEM>& p, int limiter)
    {
        for (int k=0; k<NUMCOMPS; ++k)
        {
            T a=m(k)-mm(k), b=c(k)-m(k), d=p(k)-c(k);
            T sl=T(0.5)*(a+b), sr=T(0.5)*(b+d);
            if (limiter != 0)
            {
                sl = a*b <= 0 ? T(0) : (a>0 ? min(a,b) : max(a,b));
                sr = b*d <= 0 ? T(0) : (b>0 ? min(b,d) : max(b,d));
                // Monotonized central: minmod(2*dL, (dL+dR)/2, 2*dR).
                // Keep zero slope at extrema and on either side of a plateau.
                if (limiter == 3)
                {
                    if (a*b > 0) sl=(a>0 ? T(1) : T(-1))*
                        min(T(0.5)*fabs(a+b),T(2)*min(fabs(a),fabs(b)));
                    if (b*d > 0) sr=(b>0 ? T(1) : T(-1))*
                        min(T(0.5)*fabs(b+d),T(2)*min(fabs(b),fabs(d)));
                }
            }
            lo(k)=m(k)+T(0.5)*sl;
            hi(k)=c(k)-T(0.5)*sr;
        }
        // Fall back as a complete state if unlimited extrapolation fails.
        if (!(lo(iRHO)>=RHO00 && lo(iP)>=0))
            for (int k=0; k<NUMCOMPS; ++k) lo(k)=m(k);
        if (!(hi(iRHO)>=RHO00 && hi(iP)>=0))
            for (int k=0; k<NUMCOMPS; ++k) hi(k)=c(k);
    },left,right,wm2,wm1,w,wp1,limiter);
}

// Reconstruct physical Cartesian primitives. One integrated face normal is
// used for both the physical flux and Riemann dissipation, so uniform vectors
// do not acquire a jump from the panel's changing coordinate basis.
template<typename T, typename MAP, MemType MEM>
inline void BoxOp_EulerCubedSphere<T,MAP,MEM>::secondOrderOperator(
    State& rhs, Array<State,DIM>& fluxes, const State& state, T scale) const
{
    PR_TIMERS("HelioCubed::secondOrderOperator");
    const int block=this->block();
    const T h=this->dx()[1];
    const bool magneticTest=isConstantMagneticFieldTest(ParseInputs::get_init_condition_type());
    const bool powell=ParseInputs::get_divergence_cleaning()==1;
    State native,w;
    sphToPrim_2O(native,state,block,!magneticTest);
    auto cellBasis=pointwiseSphericalBasis(native.box(),block,-1);
    HelioCubed::nativePrimitivesToCartesian(w,native,cellBasis);
    const Box metricBox=rhs.box().grow(1);
    Scalar radius(metricBox),volume(metricBox),minusDiv(rhs.box(),T(0));
    Vector dr(metricBox),adjDr(metricBox);
    radialMetrics(radius,dr,adjDr,volume,metricBox);
    for (int dir=0; dir<DIM; ++dir)
    {
        ScopedHostStack<MEM> scratch;
        const Box faces=fluxes[dir].box();
        State left,right;
        reconstructPLM(left,right,w,faces,dir);
        if (magneticTest)
        {
            const T b0=testProblemBackgroundBx(ParseInputs::get_init_condition_type());
            forallInPlace([] PROTO_LAMBDA(Var<T,NUMCOMPS,MEM>& l,
                Var<T,NUMCOMPS,MEM>& r,T b0)
            {l(iBX)+=b0;r(iBX)+=b0;},left,right,b0);
        }
        auto normal=CubedSphereShell::adjA<T,MEM>(faces,h,block,dir);
        auto radialArea=slice(adjDr,dir);
        const bool roe=ParseInputs::get_Riemann_solver_type()==2;
        const T entropy=ParseInputs::get_entropy_fix_coeff();
        forallInPlace([] PROTO_LAMBDA(
            Var<T,NUMCOMPS,MEM>& out,
            const Var<T,NUMCOMPS,MEM>& l,const Var<T,NUMCOMPS,MEM>& r,
            const Var<T,DIM,MEM>& normal,const Var<T,1,MEM>& radialArea,
            T gamma,T entropy,bool roe)
        {
            constexpr double inv4pi=1.0/(4.0*M_PI);
            double frame[3][3]={},norm=0;
            for(int i=0;i<3;++i)norm+=normal(i)*normal(i);
            norm=sqrt(norm);
            for(int i=0;i<3;++i)frame[0][i]=normal(i)/norm;
            // Cross the normal with the least aligned Cartesian axis.
            int axis=0;
            for(int i=1;i<3;++i)if(fabs(frame[0][i])<fabs(frame[0][axis]))axis=i;
            frame[1][(axis+1)%3]=frame[0][(axis+2)%3];
            frame[1][(axis+2)%3]=-frame[0][(axis+1)%3];
            double tangent=0;
            for(int i=0;i<3;++i)tangent+=frame[1][i]*frame[1][i];
            for(int i=0;i<3;++i)frame[1][i]/=sqrt(tangent);
            for(int i=0;i<3;++i){int j=(i+1)%3,k=(i+2)%3;frame[2][i]=frame[0][j]*frame[1][k]-frame[0][k]*frame[1][j];}
            double local[2][8]={},flux[2][8]={};
            for(int side=0;side<2;++side)
            {
                local[side][0]=side==0?l(0):r(0);
                local[side][4]=side==0?l(4):r(4);
                for(int i=0;i<3;++i)for(int j=0;j<3;++j)
                {
                    local[side][1+i]+=frame[i][j]*(side==0?l(1+j):r(1+j));
                    local[side][5+i]+=frame[i][j]*(side==0?l(5+j):r(5+j));
                }
                auto* q=local[side];double v2=0,b2=0,vb=0;
                for(int i=0;i<3;++i){v2+=q[1+i]*q[1+i];b2+=q[5+i]*q[5+i];vb+=q[1+i]*q[5+i];}
                const double pressure=q[4]+0.5*inv4pi*b2;
                const double energy=q[4]/(gamma-1)+0.5*q[0]*v2+0.5*inv4pi*b2;
                flux[side][0]=q[0]*q[1];
                flux[side][4]=(energy+pressure)*q[1]-inv4pi*q[5]*vb;
                for(int i=0;i<3;++i)
                {
                    flux[side][1+i]=q[0]*q[1]*q[1+i]-inv4pi*q[5]*q[5+i]+(i==0?pressure:0);
                    flux[side][5+i]=q[1]*q[5+i]-q[1+i]*q[5];
                }
            }
            double diss[8],numerical[8];
            if(!roe || !MHD_Riemann_Solvers::roe8WaveDissipation(diss,local[0],local[1],gamma,entropy,true))
                MHD_Riemann_Solvers::rusanovDissipation(diss,local[0],local[1],gamma);
            const double area=radialArea(0)*norm;
            for(int c=0;c<8;++c)numerical[c]=0.5*area*(flux[0][c]+flux[1][c]-diss[c]);
            out(0)=numerical[0];out(4)=numerical[4];
            for(int i=0;i<3;++i)
            {
                out(1+i)=out(5+i)=0;
                for(int j=0;j<3;++j){out(1+i)+=frame[j][i]*numerical[1+j];out(5+i)+=frame[j][i]*numerical[5+j];}
            }
        },fluxes[dir],left,right,normal,radialArea,m_gamma,entropy,roe);
        rhs += Stencil<T>::FluxDivergence(dir)(fluxes[dir],-scale/this->dx()[dir]);
        if (powell)
        {
            auto bface=forallOverwrite<T,1,MEM>([] PROTO_LAMBDA(
                Var<T,1,MEM>& f,const Var<T,NUMCOMPS,MEM>& l,
                const Var<T,NUMCOMPS,MEM>& r,const Var<T,DIM,MEM>& normal,
                const Var<T,1,MEM>& area)
            {f(0)=0;for(int j=0;j<DIM;++j)f(0)+=T(0.5)*(l(iBX+j)+r(iBX+j))*normal(j)*area(0);},left,right,normal,radialArea);
            minusDiv += Stencil<T>::FluxDivergence(dir)(bface,-scale/this->dx()[dir]);
        }
    }
    if (powell)
    {
        if (magneticTest) addTestProblemBackgroundB(native,block);
        DivBPowell2(rhs,minusDiv,native,cellBasis,false);
    }
}
