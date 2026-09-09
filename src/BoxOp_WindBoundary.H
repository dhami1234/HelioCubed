#pragma once
#include "WindBoundary.H"
#include <map>

namespace WindBoundary
{
inline double densityUnit() { return c_TEST_RHO0; }
inline double velocityUnit() { return 5.e7; }
inline double timeUnit() { return CUBED_SPHERE_SHELL_R0/velocityUnit(); }
inline double magneticUnit() { return std::sqrt(4*pi*densityUnit())*velocityUnit(); }
inline double unit(int c)
{
    if (c==0) return densityUnit();
    if (c<4) return densityUnit()*velocityUnit();
    if (c==4) return densityUnit()*velocityUnit()*velocityUnit();
    return magneticUnit();
}
inline bool enabled()
{
    return ParseInputs::use_wind_boundary();
}

inline void validateSettings()
{
    const int order=ParseInputs::get_wind_boundary_order();
    if (order!=0 && order!=4) throw std::runtime_error("wind_boundary_order must be 0 or 4");
    const int type=ParseInputs::get_init_condition_type();
    if ((type==10 || type==11) && (order!=4 || ParseInputs::get_Sun_gravity()!=0))
        throw std::runtime_error("Exact wind tests 10 and 11 require wind_boundary_order=4 and Sun_gravity=0");
    if (!enabled()) return;
    if (ParseInputs::get_domainSize()<2*NGHOST || ParseInputs::get_thickness()<5)
        throw std::runtime_error("Wind BC requires at least ten angular and five radial cells per block");
    if (ParseInputs::get_divergence_cleaning()!=1 || ParseInputs::get_spatial_operator()!=1)
        throw std::runtime_error("Wind BC currently requires spatial_operator=1 and divergence_cleaning=1 (Powell)");
    if (CUBED_SPHERE_SHELL_RADIAL_COORD!=0 || CUBED_SPHERE_RADIAL_STRETCH<=0)
        throw std::runtime_error("Wind BC requires radial coordinate 0 and positive radial stretch");
    const double ghostRadius=1+(CUBED_SPHERE_SHELL_R1/CUBED_SPHERE_SHELL_R0-1)
        *std::expm1(CUBED_SPHERE_RADIAL_STRETCH*(.5-NGHOST)/ParseInputs::get_thickness())
        /std::expm1(CUBED_SPHERE_RADIAL_STRETCH);
    if(!(ghostRadius>0)) throw std::runtime_error("Wind BC: increase thickness so radial ghost centers have positive radius");
    const int temporal=ParseInputs::get_temporal_order();
    if (temporal!=1 && temporal!=3 && temporal!=4)
        throw std::runtime_error("Wind BC supports temporal_order 1, 3 and 4");
}

// Exact, smooth MHD expansion, with a density pattern advected through the
// inner boundary. u=x/(tau+t), rho=lambda^-3 f(x/lambda), p=lambda^-3gamma p0,
// B=lambda^-2 B0 solves ideal MHD without gravity (uniform p and B in space).
inline State exact(const Vec& x,double t,double gamma)
{
    const double lambda=1+t;
    const double r=std::sqrt(dot(x,x));
    // Two initial density cycles across the shell; lengths here are in units of r0.
    const double wavelength=.5*(CUBED_SPHERE_SHELL_R1/CUBED_SPHERE_SHELL_R0-1);
    State w{};
    w[0]=std::pow(lambda,-3)*(1+.1*std::sin(2*pi*(r/lambda-1)/wavelength));
    for (int d=0;d<3;++d) w[1+d]=x[d]/lambda;
    w[4]=c_TEST_P0/unit(4)*std::pow(lambda,-3*gamma);
    // Case 10 is hydrodynamic; case 11 retains the uniform Cartesian field.
    if (ParseInputs::get_init_condition_type()==11)
    {
        w[5]=.05/(lambda*lambda); w[6]=.03/(lambda*lambda); w[7]=-.02/(lambda*lambda);
    }
    return conserved(w,gamma);
}

// Stationary uniform magnetic reference for case 11, at t=0. Average its
// local components over the angular cell; scalar energy is uniform and the
// velocity is zero. This is a force reference, not the expanding solution.
inline State magneticEquilibriumAverage(int j,int k,double h,int block,double gamma)
{
    const State reference=exact({0,0,0},0,gamma);
    State result{};
    result[0]=1.;
    result[4]=reference[4];
    const double q[3]={-std::sqrt(3./5)/2,0,std::sqrt(3./5)/2};
    const double w[3]={5./18,4./9,5./18};
    for(int a=0;a<3;++a) for(int b=0;b<3;++b)
    {
        const auto local=transform(reference,inverse(basis(j+q[a],k+q[b],h,block)));
        for(int c=5;c<8;++c) result[c]+=w[a]*w[b]*local[c];
    }
    return result;
}

inline Geometry geometry(const Array<double,DIM>& dx,int block)
{
    return {CUBED_SPHERE_SHELL_R1/CUBED_SPHERE_SHELL_R0,
        CUBED_SPHERE_RADIAL_STRETCH,dx[0],dx[1],ParseInputs::get_gamma(),
        ParseInputs::get_Sun_gravity() ? c_G*c_M0/(CUBED_SPHERE_SHELL_R0*velocityUnit()*velocityUnit()) : 0.,block};
}

// Same exchange and semi-spherical representation as Proto's Euler path,
// with the missing fourth-order product in semiSphToSph restored. Keeping
// this small adapter here avoids modifying a user's external Proto checkout.
template<typename T, MemType MEM>
void consToSph(MBLevelBoxData<T,8,MEM>& u,MBInterpOp& interp,
               MBLevelBoxData<T,1,MEM>& volume,unsigned int order=4,int counter=0)
{
    if (!enabled())
    {
        CubedSphereShell::consToSphInterpEuler(u,interp,volume,order,counter);
        return;
    }
    const auto layout=u.layout();
    const auto ghosts=u.ghost();
    Array<Point,DIM+1> localGhosts={Point::Ones(3),Point::Zeros(),Point::Zeros(),Point::Zeros()};
    MBLevelBoxData<T,8,MEM> local(layout,localGhosts);
    local.setVal(0.);u.copyTo(local);local.exchange0();
    MBLevelBoxData<T,10,MEM> input(layout,ghosts+Point(0,2,2)),semi(layout,ghosts);
    input.setVal(0.);semi.setVal(0.);
    const T h=1./layout.domain().getBlock(0).box().sizes()[1];
    for(auto dit:layout)
    {
        const int block=layout.block(dit);
        CubedSphereShell::consToSemiSphNGEuler(input[dit],local[dit],volume[dit],layout[dit],
            layout.getBlock(block).domain().box(),h,block,4);
        input[dit].copyTo(semi[dit],layout[dit]);
    }
    input.exchange();input.copyTo(semi);interp.apply(semi,input);semi.exchange0();
    for(auto dit:layout)
    {
        const int block=layout.block(dit);
        const Box box=semi[dit].box();
        auto det=CubedSphereShell::detA<T,MEM>(box,h,0);
        auto a=CubedSphereShell::detA_A<T,MEM>(box,h,block,0);
        BoxData<T,DIM,MEM,DIM> inv(box);inv.setVal(0.);
        forallInPlace(f_matinv3by3,inv,a);
        for(int begin : {1,5})
        {
            auto cart=slice<T,10,DIM,MEM>(semi[dit],begin);
            auto sph=slice<T,8,DIM,MEM>(u[dit],begin);
            BoxData<T,DIM,MEM> product(box),bar(box);
            product.setVal(0.);bar.setVal(0.);
            for(int d=0;d<3;++d)
            {
                auto dst=slice(product,d),src=slice(cart,d);
                // <det(A) V> = <det(A)><V> + sum Ddet(A) DV / 12.
                // Include the outer ghost layer: the generic product's
                // centered derivatives require a smaller destination box.
                for (auto pt:box)
                {
                    T value=det(pt,0)*src(pt,0);
                    for(int axis=1;axis<3;++axis)
                    {
                        const int start=std::max(box.low()[axis],std::min(pt[axis]-1,box.high()[axis]-2));
                        const int offset=start-pt[axis];
                        const double w[3]={-offset-1.5,2.*offset+2.,-offset-.5};
                        T dd=0,dv=0;
                        for(int i=0;i<3;++i)
                        {
                            Point p=pt;p[axis]=start+i;
                            dd+=w[i]*det(p,0);dv+=w[i]*src(p,0);
                        }
                        value+=dd*dv/12.;
                    }
                    dst(pt,0)=value;
                }
            }
            CubedSphereShell::cellMatrixQuotientMatAndInv(sph,bar,a,inv,product,false);
            auto radial=slice(u[dit],begin);
            radial+=slice(semi[dit],begin==1?8:9);
        }
        for(int c : {0,4})
        {
            auto dst=slice(u[dit],c),src=slice(semi[dit],c);
            src.copyTo(dst);
        }
    }
    u.exchange0();
}

// BC files contain classical spherical primitives, ordered r,phi,theta.
inline State fromFile(const State& values,const Matrix& a,double gamma)
{
    const Vec n={a[0][0],a[1][0],a[2][0]};
    const double phi=std::atan2(n[1],n[0]), st=std::hypot(n[0],n[1]), ct=n[2];
    const Vec et={ct*std::cos(phi),ct*std::sin(phi),-st};
    const Vec ep={-std::sin(phi),std::cos(phi),0};
    State w{};
    w[0]=values[0]*c_MP/densityUnit();
    w[4]=values[4]*1.e-12/unit(4);
    for (int d=0;d<3;++d)
    {
        // Velocity files are already in cm/s, as in the legacy reader.
        w[1+d]=(values[1]*n[d]+values[2]*ep[d]+values[3]*et[d])/velocityUnit();
        w[5+d]=(values[5]*n[d]+values[6]*ep[d]+values[7]*et[d])*1.e-6/magneticUnit();
    }
    return conserved(w,gamma);
}

using CacheKey=std::array<int,7>;
inline std::map<CacheKey,Patch>& stageCache()
{ static std::map<CacheKey,Patch> cache; return cache; }
inline long& limitedColumns() { static long count=0;return count; }
inline long& limitedGhostColumns() { static long count=0;return count; }
inline long& limitedOutflowColumns() { static long count=0;return count; }

// Boundary traces use the same local conserved components as the polynomial.
template<class Data>
State boundaryState(const Data& bc,int type,int j,int k,const Geometry& geo,double time)
{
    const Matrix a=basis(j,k,geo.ha,geo.block);
    const Vec n={a[0][0],a[1][0],a[2][0]};
    if (type==10 || type==11) return transform(exact(n,time/timeUnit(),geo.gamma),inverse(a));
    if (type==1)
    {
        const double theta=std::acos(std::max(-1.,std::min(1.,n[2])));
        const double phi=std::atan2(n[1],n[0])-2*pi*time/(25.38*24*3600);
        const double pattern=std::pow(std::sin(theta),5)*std::sin(5*phi)
                            +std::pow(std::sin(theta),10)*std::sin(10*phi);
        const double vr=1-ParseInputs::get_strength()/300.*pattern;
        State w{}; w[0]=1; w[4]=c_TEST_P0/unit(4);
        w[1]=vr;
        return conserved(w,geo.gamma);
    }
    State values{};
    const Point pt(bc.box().low()[0],j,k);
    for (int c=0;c<8;++c) values[c]=bc(pt,c);
    return transform(fromFile(values,a,geo.gamma),inverse(a));
}

template<class Data>
void fill(Data& usph,const Data& bc,const Data* bcTimeSamples,
          const Box& blockBox,const Array<double,DIM>& dx,int block,
          double time,double dt,int stage)
{
    const int type=ParseInputs::get_init_condition_type();
    const Geometry geo=geometry(dx,block);
    const bool boundaryOnly=type==0 && stage<0;
    const Box box=usph.box();
    const int jl=box.low()[1],jh=box.high()[1],kl=box.low()[2],kh=box.high()[2];
    const int vjl=std::max(jl,blockBox.low()[1]),vjh=std::min(jh,blockBox.high()[1]);
    const int vkl=std::max(kl,blockBox.low()[2]),vkh=std::min(kh,blockBox.high()[2]);
    if (blockBox.low()[0]>box.low()[0])
    {
        if (box.high()[0]<blockBox.low()[0]+2)
            throw std::runtime_error("Wind BC requires three available interior radial cells");
        Patch diagnostic;
        Patch& patch=stage<0 ? diagnostic : stageCache()[{block,box.low()[0],box.high()[0],jl,jh,kl,kh}];
        if (stage<=0)
        {
            patch.define(jl,jh,kl,kh);
            for (int k=kl;k<=kh;++k) for (int j=jl;j<=jh;++j)
            {
                const Matrix a=basis(j,k,geo.ha,block);
                const State g=boundaryState(bc,type,j,k,geo,time);
                State gt{},dj{},dk{};
                // A small fixed physical interval keeps differentiation error
                // independent of the mesh. Five forward samples also work at
                // t=0; the analytic drivers are evaluated directly.
                const double eps=type==0 ? 1. : 1.e-4*timeUnit();
                const double tw[5]={-25./12,4.,-3.,4./3,-.25};
                for (int i=0;i<5;++i)
                {
                    const State u=type==0 && i>0 && bcTimeSamples
                        ? boundaryState(bcTimeSamples[i-1],type,j,k,geo,time+i*eps)
                        : boundaryState(bc,type,j,k,geo,time+i*eps);
                    for (int c=0;c<8;++c) gt[c]+=tw[i]*(u[c]-g[c])*timeUnit()/eps;
                }
                for (int dir=0;dir<2;++dir)
                {
                    const int center=dir==0?j:k,lo=dir==0?jl:kl,hi=dir==0?jh:kh;
                    const int start=std::max(lo,std::min(center-2,hi-4));
                    const auto w=weights(start-center,false,1);
                    for (int i=0;i<5;++i)
                    {
                        const State u=boundaryState(bc,type,dir==0?start+i:j,dir==1?start+i:k,geo,time);
                        for (int c=0;c<8;++c) (dir==0?dj:dk)[c]+=w[i]*u[c];
                    }
                }
                const State slope=boundaryOnly ? State{} : normalSlope(g,gt,dj,dk,geo,j,k);
                auto& polynomial=patch.initial[patch.index(j,k)];
                polynomial[0]=g; polynomial[1]=slope;
                // At a three-block corner, the distant angular ghost cells
                // are not a complete rectangular set of valid averages.
                // Reconstruct radial curvature on the physical face first.
                if (j<vjl || j>vjh || k<vkl || k>vkh) continue;
                const int js=std::max(vjl,std::min(j-2,vjh-4));
                const int ks=std::max(vkl,std::min(k-2,vkh-4));
                const auto wj=weights(js-j,true),wk=weights(ks-k,true);
                std::array<State,4> interior{};
                for (int r=0;r<(boundaryOnly?4:3);++r)
                {
                    State u{};
                    for (int b=0;b<5;++b) for (int aa=0;aa<5;++aa)
                    {
                        const Point pt(blockBox.low()[0]+r,js+aa,ks+b);
                        for (int c=0;c<8;++c) u[c]+=wj[aa]*wk[b]*usph(pt,c)/unit(c);
                    }
                    interior[r]=u;
                }
                patch.initial[patch.index(j,k)]=boundaryOnly ? reconstructValue(g,interior)
                    : reconstruct(g,slope,{interior[0],interior[1],interior[2]});
                if(limitPolynomial(patch.initial[patch.index(j,k)],geo.gamma,a)) ++limitedColumns();
            }
            for (int k=kl;k<=kh;++k) for (int j=jl;j<=jh;++j)
            {
                if (j>=vjl && j<=vjh && k>=vkl && k<=vkh) continue;
                const int js=std::max(vjl,std::min(j-2,vjh-4));
                const int ks=std::max(vkl,std::min(k-2,vkh-4));
                const auto wj=weights(js-j,false),wk=weights(ks-k,false);
                auto& p=patch.initial[patch.index(j,k)];
                for (int b=0;b<5;++b) for (int a=0;a<5;++a)
                    for (int m=(boundaryOnly?1:2);m<5;++m) for (int c=0;c<8;++c)
                        p[m][c]+=wj[a]*wk[b]*patch.initial[patch.index(js+a,ks+b)][m][c];
                if(limitPolynomial(p,geo.gamma,basis(j,k,geo.ha,block))) ++limitedGhostColumns();
            }
            patch.current=patch.initial;
        }
        else
        {
            if (patch.current.empty()) throw std::runtime_error("Wind BC: RK stage without initial boundary state");
            patch.advanceStage(geo,dt/timeUnit(),stage,ParseInputs::get_temporal_order());
            for(int k=kl;k<=kh;++k)for(int j=jl;j<=jh;++j)
            {
                const bool outside=j<vjl || j>vjh || k<vkl || k>vkh;
                if(outside)
                {
                    // A distant corner is an extrapolated stage, not an
                    // evolved physical cell. Fall back to its initial trace
                    // if that extrapolation leaves the admissible set.
                    auto& p=patch.current[patch.index(j,k)];
                    try { primitive(transform(p[0],basis(j,k,geo.ha,block)),geo.gamma); }
                    catch(const std::runtime_error&) { p=patch.initial[patch.index(j,k)];++limitedGhostColumns(); }
                }
                if(limitPolynomial(patch.current[patch.index(j,k)],geo.gamma,basis(j,k,geo.ha,block)))
                    ++(outside ? limitedGhostColumns() : limitedColumns());
            }
        }
        // Integrate the local components radially, then convolve angularly.
        // The polynomial already uses the same vector basis as usph.
        Data radial(box);
        for (int k=kl;k<=kh;++k) for (int j=jl;j<=jh;++j)
        {
            for (int r=box.low()[0];r<blockBox.low()[0];++r)
            {
                const State u=average(patch.current[patch.index(j,k)],r,r+1);
                // Only the boundary-adjacent stencil on the physical face
                // is a physical admissibility check; unused far corner slots
                // can lie outside the local chart's useful extrapolation range.
                if (stage>=0 && r>=-3 && j>=vjl && j<=vjh && k>=vkl && k<=vkh)
                {
                    try { primitive(transform(u,basis(j,k,geo.ha,block)),geo.gamma); }
                    catch(const std::runtime_error& error)
                    {
                        std::ostringstream message;
                        message<<error.what()<<" (block="<<block<<", cell="<<r<<","<<j<<","<<k<<", stage="<<stage<<")";
                        throw std::runtime_error(message.str());
                    }
                }
                for (int c=0;c<8;++c) radial(Point(r,j,k),c)=u[c]*unit(c);
            }
        }
        for (int k=kl;k<=kh;++k) for (int j=jl;j<=jh;++j)
        {
            const int js=std::max(jl,std::min(j-2,jh-4)),ks=std::max(kl,std::min(k-2,kh-4));
            auto wj=weights(js-j,false,2),wk=weights(ks-k,false,2);
            for (int r=box.low()[0];r<blockBox.low()[0];++r) for (int c=0;c<8;++c)
            {
                double u=radial(Point(r,j,k),c);
                for (int i=0;i<5;++i)
                    u+=(wj[i]*radial(Point(r,js+i,k),c)+wk[i]*radial(Point(r,j,ks+i),c))/24.;
                usph(Point(r,j,k),c)=u;
            }
        }
    }
    // Smooth outflow: degree-four extrapolation of radial cell averages.
    // The physical state and B divergence must remain compatible here too.
    if (blockBox.high()[0]<box.high()[0])
    {
        const int last=blockBox.high()[0];
        if (last-4<box.low()[0]) throw std::runtime_error("Wind BC requires five available outflow cells");
        for (int r=last+1;r<=box.high()[0];++r)
        {
            std::array<double,5> w{};
            for (int i=0;i<5;++i)
            {
                w[i]=1;
                for (int m=0;m<5;++m) if(i!=m) w[i]*=double(r-(last-4+m))/(i-m);
            }
            for (int k=kl;k<=kh;++k) for (int j=jl;j<=jh;++j) for (int c=0;c<8;++c)
            {
                double u=0;
                for (int i=0;i<5;++i) u+=w[i]*usph(Point(last-4+i,j,k),c);
                usph(Point(r,j,k),c)=u;
            }
        }
        // Quartic extrapolation can overshoot on under-resolved profiles.
        // Apply one common blend to all conserved components and ghost layers
        // in a column; resolved, admissible extrapolations remain unchanged.
        for (int k=kl;k<=kh;++k) for (int j=jl;j<=jh;++j)
        {
            const Matrix a=basis(j,k,geo.ha,block);
            State anchor{};
            for(int c=0;c<8;++c) anchor[c]=usph(Point(last,j,k),c)/unit(c);
            if(!admissible(transform(anchor,a)) && (j<vjl || j>vjh || k<vkl || k>vkh))
            {
                // Far angular corners are auxiliary chart extensions. Use a
                // physical column as the fallback, rotating its vector basis.
                const int js=std::clamp(j,vjl,vjh),ks=std::clamp(k,vkl,vkh);
                for(int c=0;c<8;++c) anchor[c]=usph(Point(last,js,ks),c)/unit(c);
                anchor=transform(transform(anchor,basis(js,ks,geo.ha,block)),inverse(a));
            }
            if(!admissible(transform(anchor,a)))
            {
                std::ostringstream message;
                message<<"Wind outflow: non-admissible anchor (block="<<block
                    <<", cell="<<last<<","<<j<<","<<k<<", stage="<<stage<<")";
                throw std::runtime_error(message.str());
            }
            bool limited=false;
            for(int attempt=0;attempt<60;++attempt)
            {
                bool valid=true;
                for(int r=last+1;r<=box.high()[0];++r)
                {
                    State u{};
                    for(int c=0;c<8;++c) u[c]=usph(Point(r,j,k),c)/unit(c);
                    valid=valid && admissible(transform(u,a));
                }
                if(valid) break;
                if(attempt==59)
                    throw std::runtime_error("Wind outflow: could not recover admissible ghost averages");
                limited=true;
                for(int r=last+1;r<=box.high()[0];++r) for(int c=0;c<8;++c)
                    usph(Point(r,j,k),c)=.5*(usph(Point(r,j,k),c)+anchor[c]*unit(c));
            }
            if(limited && stage>=0 && j>=vjl && j<=vjh && k>=vkl && k<=vkh)
                ++limitedOutflowColumns();
        }
    }
}

// Independent tensor Gauss reference for the full mapped cell average <JU>.
// No solver reconstruction or ghost values enter this diagnostic.
template<class Level>
void reportExactError(const Level& ju,double time)
{
    double sum[9]={},maximum[8]={},inner[8]={};
    int location[8]={};
    const auto layout=ju.layout();
    const auto sizes=layout.domain().getBlock(0).box().sizes();
    const Array<double,DIM> dx={1./sizes[0],1./sizes[1],1./sizes[2]};
    const double q[3]={-std::sqrt(3./5)/2,0,std::sqrt(3./5)/2};
    const double w[3]={5./18,4./9,5./18};
    for (auto dit : layout)
    {
        const auto geo=geometry(dx,layout.block(dit));
        for (auto pt : layout[dit])
        {
            State ref{}; double volume=0;
            for (int a=0;a<3;++a) for (int b=0;b<3;++b) for (int c=0;c<3;++c)
            {
                const double s=pt[0]+.5+q[a],r=geo.radius(s);
                const Matrix basisPoint=basis(pt[1]+q[b],pt[2]+q[c],geo.ha,geo.block);
                const auto& m=basisPoint;
                const double det=m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
                                -m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
                                +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
                const double jac=det*r*r*geo.drds(s)/geo.hr*w[a]*w[b]*w[c];
                const State u=exact({r*m[0][0],r*m[1][0],r*m[2][0]},time/timeUnit(),geo.gamma);
                volume+=jac;
                for (int comp=0;comp<8;++comp) ref[comp]+=jac*u[comp];
            }
            sum[8]+=volume;
            for (int c=0;c<8;++c)
            {
                const double err=std::abs(ju[dit](pt,c)/(unit(c)*std::pow(CUBED_SPHERE_SHELL_R0,3))-ref[c]);
                sum[c]+=err;
                if (err/volume>maximum[c])
                {
                    maximum[c]=err/volume;
                    location[c]=((geo.block*sizes[0]+pt[0])*sizes[1]+pt[1])*sizes[2]+pt[2];
                }
                if (pt[0]<3) inner[c]=std::max(inner[c],err/volume);
            }
        }
    }
#ifdef PR_MPI
    double localMaximum[8];std::copy(maximum,maximum+8,localMaximum);
    MPI_Allreduce(MPI_IN_PLACE,sum,9,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,maximum,8,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,inner,8,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    for(int c=0;c<8;++c) if(localMaximum[c]!=maximum[c]) location[c]=-1;
    MPI_Allreduce(MPI_IN_PLACE,location,8,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif
    for (int c=0;c<8;++c)
    {
        int id=location[c];const int k=id%sizes[2];id/=sizes[2];
        const int j=id%sizes[1];id/=sizes[1];const int r=id%sizes[0],block=id/sizes[0];
        pout(0)<<std::setprecision(12)<<"Wind exact error time="<<time<<" component="<<c
               <<" L1="<<sum[c]/sum[8]<<" Linf="<<maximum[c]<<" inner_Linf="<<inner[c]
               <<" max_cell="<<block<<":"<<r<<","<<j<<","<<k<<std::endl;
    }
}
} // namespace WindBoundary
