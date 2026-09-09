#pragma once

// Smooth super-fast inflow, in dimensionless cubed-sphere conserved variables.
// Momentum and B use the interior solver's radial/tangential basis; E is total energy.
// A quartic matches U, its ILW normal derivative, and three interior averages.
// The same RK method advances the polynomial at intermediate stages; imposing
// the exact physical solution at stage times is not stage-consistent.
// References: Tan & Shu, JCP 229 (2010), 8144;
// Zhao, Huang & Ruuth, JCP 421 (2020), 109697.
// This header has no Proto dependency so its numerical kernels can be tested
// independently of the mapped-grid solver.
#include <array>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace WindBoundary
{
using Vec = std::array<double,3>;
using State = std::array<double,8>;
using Poly = std::array<State,5>;
using Matrix = std::array<Vec,3>;
constexpr double pi = 3.14159265358979323846;

template<size_t N>
std::array<double,N> solve(std::array<std::array<double,N>,N> a,
                           std::array<double,N> b)
{
    for (size_t k=0; k<N; ++k)
    {
        size_t pivot=k;
        for (size_t i=k+1; i<N; ++i)
            if (std::abs(a[i][k])>std::abs(a[pivot][k])) pivot=i;
        if (std::abs(a[pivot][k])<1.e-13)
            throw std::runtime_error("Wind BC: singular reconstruction/Jacobian");
        std::swap(a[k],a[pivot]); std::swap(b[k],b[pivot]);
        const double diag=a[k][k];
        for (size_t j=k; j<N; ++j) a[k][j]/=diag;
        b[k]/=diag;
        for (size_t i=0; i<N; ++i) if (i!=k)
        {
            const double factor=a[i][k];
            for (size_t j=k; j<N; ++j) a[i][j]-=factor*a[k][j];
            b[i]-=factor*b[k];
        }
    }
    return b;
}

inline double dot(const Vec& a,const Vec& b)
{ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

inline Matrix orientBasis(const Matrix& m,int block)
{
    const int perm[6][3]={{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
    const int sign[6][3]={{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}};
    Matrix out{};
    for (int i=0;i<3;++i) for (int d=0;d<3;++d)
        out[i][d]=sign[block][i]*m[perm[block][i]][d];
    return out;
}

// A maps cubed-sphere vector components to Cartesian physical components.
inline Matrix basis(double j,double k,double h,int block)
{
    const double a=std::tan(-pi/4+pi/2*(j+.5)*h);
    const double b=std::tan(-pi/4+pi/2*(k+.5)*h);
    const double n=1/std::sqrt(1+a*a+b*b), n3=n*n*n;
    const Matrix m={{{n,-a*(1+a*a)*n3*pi/2,-b*(1+b*b)*n3*pi/2},
                     {a*n,(1+a*a)*(1+b*b)*n3*pi/2,-a*b*(1+b*b)*n3*pi/2},
                     {b*n,-a*b*(1+a*a)*n3*pi/2,(1+a*a)*(1+b*b)*n3*pi/2}}};
    return orientBasis(m,block);
}

// Derivative with respect to a computational angle (not the cell index).
// A = [n, d_eta n, d_zeta n], where n = (1,tan(alpha),tan(beta))/|...|.
inline Matrix basisDerivative(double j,double k,double h,int block,int dir)
{
    const double c=pi/2, a=std::tan(-pi/4+c*(j+.5)*h), b=std::tan(-pi/4+c*(k+.5)*h);
    const Vec q={1,a,b};
    const std::array<Vec,2> dq={Vec{0,c*(1+a*a),0},Vec{0,0,c*(1+b*b)}};
    const double n=1/std::sqrt(dot(q,q)), n3=n*n*n, n5=n3*n*n;
    const double d=dot(q,dq[dir]);
    Matrix m{};
    for(int x=0;x<3;++x) m[x][0]=n*dq[dir][x]-n3*q[x]*d;
    for(int axis=0;axis<2;++axis)
    {
        Vec ddq{};
        if(axis==dir) ddq[axis+1]=2*c*(axis==0?a:b)*dq[axis][axis+1];
        const double e=dot(q,dq[axis]), de=dot(dq[dir],dq[axis])+dot(q,ddq);
        for(int x=0;x<3;++x)
            m[x][axis+1]=n*ddq[x]-n3*(dq[dir][x]*e+dq[axis][x]*d+q[x]*de)
                +3*n5*q[x]*d*e;
    }
    return orientBasis(m,block);
}

inline Matrix inverse(const Matrix& a)
{
    Matrix out{};
    for (int j=0;j<3;++j)
    {
        Vec unit{}; unit[j]=1;
        const Vec col=solve<3>(a,unit);
        for (int i=0;i<3;++i) out[i][j]=col[i];
    }
    return out;
}

inline State transform(State u,const Matrix& a)
{
    const State old=u;
    for (int start : {1,5}) for (int i=0;i<3;++i)
    {
        u[start+i]=0;
        for (int j=0;j<3;++j) u[start+i]+=a[i][j]*old[start+j];
    }
    return u;
}

inline State conserved(const State& w,double gamma)
{
    State u=w;
    u[4]=w[4]/(gamma-1);
    for (int i=0;i<3;++i)
    {
        u[1+i]=w[0]*w[1+i];
        u[4]+=.5*(w[0]*w[1+i]*w[1+i]+w[5+i]*w[5+i]);
    }
    return u;
}

inline State primitive(const State& u,double gamma)
{
    if (!(u[0]>0)) throw std::runtime_error("Wind BC: nonpositive density");
    State w=u;
    for (int i=0;i<3;++i)
    {
        w[1+i]=u[1+i]/u[0];
        w[4]-=.5*(u[1+i]*w[1+i]+u[5+i]*u[5+i]);
    }
    w[4]*=gamma-1;
    if (!(w[4]>0) || !std::isfinite(w[4]))
        throw std::runtime_error("Wind BC: nonpositive/nonfinite pressure in smooth reconstruction");
    return w;
}

inline double fastMargin(const State& u,const Vec& n,double gamma)
{
    const State w=primitive(u,gamma);
    Vec v={w[1],w[2],w[3]}, b={w[5],w[6],w[7]};
    const double a2=gamma*w[4]/w[0], b2=dot(b,b)/w[0];
    const double bn2=std::pow(dot(b,n),2)/w[0];
    const double cf2=.5*(a2+b2+std::sqrt(std::max(0.,std::pow(a2+b2,2)-4*a2*bn2)));
    return dot(v,n)-std::sqrt(cf2);
}

// Cartesian ideal-MHD/Powell RHS. For divergence-free B this is conservative
// ideal MHD. This is also the eight-wave system used by HelioCubed's cleaning.
inline State rhs(const State& u,const std::array<State,3>& grad,
                 double gamma,const Vec& gravity={0,0,0})
{
    const State w=primitive(u,gamma);
    std::array<State,3> dw=grad;
    for (int d=0;d<3;++d)
    {
        for (int i=0;i<3;++i)
        {
            dw[d][1+i]=(grad[d][1+i]-w[1+i]*grad[d][0])/w[0];
            dw[d][4]-=w[1+i]*grad[d][1+i]-.5*w[1+i]*w[1+i]*grad[d][0]
                       +w[5+i]*grad[d][5+i];
        }
        dw[d][4]*=gamma-1;
    }
    double divv=0;
    for (int d=0;d<3;++d) divv+=dw[d][1+d];
    State wt{};
    for (int c=0;c<8;++c) for (int d=0;d<3;++d) wt[c]-=w[1+d]*dw[d][c];
    wt[0]-=w[0]*divv;
    wt[4]-=gamma*w[4]*divv;
    for (int i=0;i<3;++i)
    {
        double force=dw[i][4];
        for (int d=0;d<3;++d)
        {
            force+=w[5+d]*(dw[i][5+d]-dw[d][5+i]);
            wt[5+i]+=w[5+d]*dw[d][1+i];
        }
        wt[1+i]+=-force/w[0]+gravity[i];
        wt[5+i]-=w[5+i]*divv;
    }
    State ut=wt;
    ut[4]=wt[4]/(gamma-1);
    for (int i=0;i<3;++i)
    {
        ut[1+i]=w[0]*wt[1+i]+w[1+i]*wt[0];
        ut[4]+=.5*w[1+i]*w[1+i]*wt[0]+w[0]*w[1+i]*wt[1+i]+w[5+i]*wt[5+i];
    }
    return ut;
}

inline State evaluate(const Poly& p,double s,int derivative=0)
{
    State out{};
    for (int k=derivative;k<5;++k)
    {
        double factor=std::pow(s,k-derivative);
        for (int d=0;d<derivative;++d) factor*=k-d;
        for (int c=0;c<8;++c) out[c]+=factor*p[k][c];
    }
    return out;
}

inline State average(const Poly& p,double lo,double hi)
{
    State out{};
    for (int k=0;k<5;++k)
    {
        const double factor=(std::pow(hi,k+1)-std::pow(lo,k+1))/((k+1)*(hi-lo));
        for (int c=0;c<8;++c) out[c]+=factor*p[k][c];
    }
    return out;
}

// Coefficients of a quartic in s=(xi-xi_boundary)/dxi. All moments below
// are computational-cell averages, not physical-volume averages.
inline Poly reconstruct(const State& g,const State& slope,
                        const std::array<State,3>& interior)
{
    Poly p{}; p[0]=g; p[1]=slope;
    std::array<std::array<double,3>,3> m{};
    for (int i=0;i<3;++i) for (int k=2;k<5;++k)
        m[i][k-2]=(std::pow(i+1,k+1)-std::pow(i,k+1))/(k+1);
    for (int c=0;c<8;++c)
    {
        std::array<double,3> b{};
        for (int i=0;i<3;++i) b[i]=interior[i][c]-g[c]-(i+.5)*slope[c];
        const auto x=solve<3>(m,b);
        for (int k=2;k<5;++k) p[k][c]=x[k-2];
    }
    return p;
}

// Diagnostic ghost fill when only an instantaneous boundary image is
// available: match its value and four interior averages. Evolution uses ILW.
inline Poly reconstructValue(const State& g,const std::array<State,4>& interior)
{
    Poly p{};p[0]=g;
    std::array<std::array<double,4>,4> m{};
    for(int i=0;i<4;++i)for(int k=1;k<5;++k)
        m[i][k-1]=(std::pow(i+1,k+1)-std::pow(i,k+1))/(k+1);
    for(int c=0;c<8;++c)
    {
        std::array<double,4>b{};
        for(int i=0;i<4;++i)b[i]=interior[i][c]-g[c];
        const auto x=solve<4>(m,b);
        for(int k=1;k<5;++k)p[k][c]=x[k-1];
    }
    return p;
}

inline bool admissible(const State& u)
{
    if (!(u[0]>0) || !std::isfinite(u[0])) return false;
    double energy=u[4];
    for(int d=0;d<3;++d) energy-=.5*(u[1+d]*u[1+d]/u[0]+u[5+d]*u[5+d]);
    return energy>0 && std::isfinite(energy);
}

// A smooth, resolved polynomial is unchanged. For under-resolved input,
// damp its variation toward the prescribed trace until the states used by
// the PDE and ghost averages have positive density and internal energy.
inline bool limitPolynomial(Poly& p,double gamma,const Matrix& a)
{
    primitive(transform(p[0],a),gamma);
    bool limited=false;
    for(int attempt=0;attempt<60;++attempt)
    {
        bool valid=true;
        for(int s=-5;s<=3;++s) valid=valid && admissible(transform(evaluate(p,s),a));
        for(int s=-5;s<0;++s) valid=valid && admissible(transform(average(p,s,s+1),a));
        if(valid) return limited;
        limited=true;
        for(int m=1;m<5;++m)for(int c=0;c<8;++c)p[m][c]*=.5;
    }
    throw std::runtime_error("Wind BC: could not recover an admissible ghost polynomial");
}

// Five cell averages -> point value, or five point values -> derivative.
// Offset is the first stencil cell/point relative to the target center.
inline std::array<double,5> weights(int offset,bool cellAverages,int derivative=0)
{
    std::array<std::array<double,5>,5> m{};
    std::array<double,5> b{}; b[derivative]=1;
    for (int k=2;k<=derivative;++k) b[derivative]*=k;
    for (int k=0;k<5;++k) for (int i=0;i<5;++i)
    {
        const double x=offset+i;
        m[k][i]=cellAverages ? (std::pow(x+.5,k+1)-std::pow(x-.5,k+1))/(k+1)
                                             : std::pow(x,k);
    }
    return solve<5>(m,b);
}

struct Geometry
{
    double ratio, stretch, hr, ha, gamma, gravity; // r1/r0, GM/(r0 V0^2)
    int block;
    double radius(double s) const
    { return 1+(ratio-1)*std::expm1(stretch*hr*s)/std::expm1(stretch); }
    double drds(double s) const
    { return (ratio-1)*stretch*hr*std::exp(stretch*hr*s)/std::expm1(stretch); }
};

// Geometry is fixed in time and radius. Numerical angular derivatives act on
// local components only; basis derivatives supply the geometric terms exactly.
struct LocalFrame
{
    Matrix a, inv, da, db;
    LocalFrame(const Geometry& geo,int j,int k)
        : a(basis(j,k,geo.ha,geo.block)), inv(inverse(a)),
          da(basisDerivative(j,k,geo.ha,geo.block,0)),
          db(basisDerivative(j,k,geo.ha,geo.block,1)) {}

    State rate(const State& u,const State& ds,const State& dj,const State& dk,
               const Geometry& geo,double s) const
    {
        const double r=geo.radius(s);
        const State cart=transform(u,a), radial=transform(ds,a);
        State eta=transform(dj,a), zeta=transform(dk,a);
        for(int c=0;c<8;++c) { eta[c]/=geo.ha; zeta[c]/=geo.ha; }
        for(int begin : {1,5}) for(int d=0;d<3;++d) for(int q=0;q<3;++q)
        {
            eta[begin+d]+=da[d][q]*u[begin+q];
            zeta[begin+d]+=db[d][q]*u[begin+q];
        }
        std::array<State,3> grad{};
        Vec gravity{};
        for(int d=0;d<3;++d)
        {
            gravity[d]=-geo.gravity*a[d][0]/(r*r);
            for(int c=0;c<8;++c)
                grad[d][c]=inv[0][d]*radial[c]/geo.drds(s)
                    +(inv[1][d]*eta[c]+inv[2][d]*zeta[c])/r;
        }
        return transform(rhs(cart,grad,geo.gamma,gravity),inv);
    }
};

struct Patch
{
    int jlo=0,jhi=-1,klo=0,khi=-1;
    std::vector<Poly> initial, current, firstRhs;
    size_t index(int j,int k) const { return (k-klo)*(jhi-jlo+1)+(j-jlo); }
    void define(int jl,int jh,int kl,int kh)
    {
        jlo=jl;jhi=jh;klo=kl;khi=kh;
        initial.assign((jh-jl+1)*(kh-kl+1),Poly{});
        current=initial;firstRhs=initial;
    }

    Poly angularDerivative(int j,int k,int dir,int faceSize) const
    {
        const int center=dir==0?j:k, lo=std::max(0,dir==0?jlo:klo), hi=std::min(faceSize-1,dir==0?jhi:khi);
        const int start=std::max(lo,std::min(center-2,hi-4));
        const auto w=weights(start-center,false,1);
        Poly out{};
        for (int i=0;i<5;++i)
        {
            const Poly& p=current[index(dir==0?start+i:j,dir==1?start+i:k)];
            for (int m=0;m<5;++m) for (int c=0;c<8;++c) out[m][c]+=w[i]*p[m][c];
        }
        return out;
    }

    // Reconstruct the PDE RHS as a quartic as well. Nodal interpolation is
    // used only for this local polynomial; its derivative is analytic.
    std::vector<Poly> rate(const Geometry& geo) const
    {
        std::vector<Poly> out(current.size());
        std::array<std::array<double,5>,5> vand{};
        for (int i=0;i<5;++i) for (int m=0;m<5;++m) vand[i][m]=std::pow(i-2,m);
        const int faceSize=std::lround(1/geo.ha);
        const int jl=std::max(0,jlo),jh=std::min(faceSize-1,jhi);
        const int kl=std::max(0,klo),kh=std::min(faceSize-1,khi);
        for (int k=kl;k<=kh;++k) for (int j=jl;j<=jh;++j)
        {
            const LocalFrame frame(geo,j,k);
            const Poly dj=angularDerivative(j,k,0,faceSize), dk=angularDerivative(j,k,1,faceSize);
            std::array<State,5> sample{};
            for (int i=0;i<5;++i)
            {
                const double s=i-2;
                const Poly& p=current[index(j,k)];
                const State u=evaluate(p,s), ds=evaluate(p,s,1);
                const State da=evaluate(dj,s), db=evaluate(dk,s);
                sample[i]=frame.rate(u,ds,da,db,geo,s);
            }
            for (int c=0;c<8;++c)
            {
                std::array<double,5> b{};
                for (int i=0;i<5;++i) b[i]=sample[i][c];
                const auto coeff=solve<5>(vand,b);
                for (int m=0;m<5;++m) out[index(j,k)][m][c]=coeff[m];
            }
        }
        // Do not evaluate the PDE in distant three-block corner slots.
        // Continue the rate from the physical face using a quartic instead.
        for(int k=klo;k<=khi;++k) for(int j=jlo;j<=jhi;++j)
        {
            if(j>=jl && j<=jh && k>=kl && k<=kh) continue;
            const int js=std::max(jl,std::min(j-2,jh-4)),ks=std::max(kl,std::min(k-2,kh-4));
            const auto wj=weights(js-j,false),wk=weights(ks-k,false);
            for(int b=0;b<5;++b)for(int a=0;a<5;++a)
                for(int m=0;m<5;++m)for(int c=0;c<8;++c)
                    out[index(j,k)][m][c]+=wj[a]*wk[b]*out[index(js+a,ks+b)][m][c];
        }
        return out;
    }

    void advanceStage(const Geometry& geo,double dt,int stage,int order)
    {
        if (stage==0) { current=initial; return; }
        const auto f=rate(geo);
        if (stage==1) firstRhs=f;
        for (size_t i=0;i<current.size();++i) for (int m=0;m<5;++m) for (int c=0;c<8;++c)
        {
            double increment=(stage==3 ? 1. : .5)*f[i][m][c];
            if (order==3 && stage==2) increment=2*f[i][m][c]-firstRhs[i][m][c];
            current[i][m][c]=initial[i][m][c]+dt*increment;
        }
    }
};

// The unknown slope is dU/ds. Tangential derivatives and time derivatives
// are evaluated at the physical boundary, then the normal Jacobian is solved.
inline State normalSlope(const State& g,const State& gt,const State& dj,const State& dk,
                         const Geometry& geo,int j,int k)
{
    const LocalFrame frame(geo,j,k);
    const Vec n={frame.a[0][0],frame.a[1][0],frame.a[2][0]};
    if (!(fastMargin(transform(g,frame.a),n,geo.gamma)>0))
        throw std::runtime_error("Wind BC requires u_r > c_fast,r everywhere at the inner boundary");
    const State base=frame.rate(g,State{},dj,dk,geo,0);
    std::array<std::array<double,8>,8> m{};
    State b{};
    for (int col=0;col<8;++col)
    {
        State ds{};ds[col]=1.;
        const State f=frame.rate(g,ds,dj,dk,geo,0);
        for (int row=0;row<8;++row) m[row][col]=f[row]-base[row];
        b[col]=gt[col]-base[col];
    }
    return solve<8>(m,b);
}
} // namespace WindBoundary
