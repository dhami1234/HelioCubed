#pragma once
#include "Proto.H"
#include "Inputs_Parsing.H"
#include <array>
#include <algorithm>
#include <cmath>
#include <vector>
#include <map>
#include <tuple>

// The evolution operator consumes unweighted computational-cell averages of
// native conserved variables at order four. Integrate the pointwise basis
// change and nonlinear primitive-to-conserved conversion for those averages;
// at order two, evaluate both at the cell center.
// Source components are (n, Vr, Vphi, Vtheta, P, Br, Bphi, Btheta), in the
// boundary file's units: cm^-3, cm/s, picodyne/cm^2, and microgauss.
class MHD_BoundaryAverages
{
    struct Sample
    {
        std::array<Proto::Point,4> source;
        std::array<double,4> coefficient;
        double tangent[2][2]; // classical (phi,theta) -> native angular components
        double weight;
    };
    static constexpr int angularOrder = 8;
    using Samples = std::array<Sample,angularOrder*angularOrder>;
    // Only geometry and interpolation weights are cached; the rotating source
    // values are evaluated on every call. Clear these caches on file reload.
    std::map<std::tuple<int,int,int,int>,Samples> m_samples;
    std::map<std::tuple<int,int,int,int>,Sample> m_points;
    static constexpr double radialNodes[4] = {
        0.06943184420297371,0.33000947820757187,
        0.6699905217924281,0.9305681557970262};
    static constexpr double radialWeights[4] = {
        0.17392742256872693,0.32607257743127307,
        0.32607257743127307,0.17392742256872693};

    // Resolve the piecewise bilinear source map inside each angular cell.
    static constexpr double nodes[angularOrder] = {
        0.019855071751231912,0.10166676129318664,
        0.2372337950418355,0.40828267875217511,
        0.59171732124782495,0.7627662049581645,
        0.89833323870681336,0.98014492824876809};
    static constexpr double weights[angularOrder] = {
        0.050614268145188185,0.11119051722668709,
        0.15685332293894372,0.18134189168918102,
        0.18134189168918102,0.15685332293894372,
        0.11119051722668709,0.050614268145188185};

    Samples makeSamples(int block,int n,int j,int k,const Proto::Box& sourceBox,
        const std::vector<double>& theta,const std::vector<double>& phi,bool point=false) const
    {
        using namespace Proto;
        Samples samples;
        const int perm[6][3]={{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
        const int sign[6][3]={{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}};
        for(int q=0;q<(point?1:angularOrder);++q) for(int s=0;s<(point?1:angularOrder);++s)
        {
            auto& sample=samples[angularOrder*q+s];
            sample.weight=point?1:weights[q]*weights[s];
            const double x=std::tan(-M_PI/4+(j+(point?.5:nodes[q]))*M_PI/(2*n));
            const double y=std::tan(-M_PI/4+(k+(point?.5:nodes[s]))*M_PI/(2*n));
            const double z=1/std::sqrt(1+x*x+y*y),c=z*z*z*M_PI/2;
            const double local[3][3]={{z,-x*(1+x*x)*c,-y*(1+y*y)*c},
                {x*z,(1+x*x)*(1+y*y)*c,-x*y*(1+y*y)*c},
                {y*z,-x*y*(1+x*x)*c,(1+x*x)*(1+y*y)*c}};
            double A[3][3];
            for(int a=0;a<3;++a) for(int b=0;b<3;++b)
                A[a][b]=sign[block][a]*local[perm[block][a]][b];
            const double t=std::acos(std::max(-1.,std::min(1.,A[2][0])));
            double p=std::atan2(A[1][0],A[0][0]);
            p=phi.front()+std::fmod(p-phi.front()+4*M_PI,2*M_PI);
            const double ep[3]={-std::sin(p),std::cos(p),0};
            const double et[3]={std::cos(t)*std::cos(p),std::cos(t)*std::sin(p),-std::sin(t)};
            // Contravariant tangential components, with exactly zero coupling
            // from a purely radial vector to either angular component.
            double g11=0,g12=0,g22=0;
            for(int a=0;a<3;++a)
            {
                g11+=A[a][1]*A[a][1];
                g12+=A[a][1]*A[a][2];
                g22+=A[a][2]*A[a][2];
            }
            const double determinant=g11*g22-g12*g12;
            for(int b=0;b<2;++b)
            {
                double d1=0,d2=0;
                const double* e=b==0?ep:et;
                for(int a=0;a<3;++a)
                {
                    d1+=A[a][1]*e[a];
                    d2+=A[a][2]*e[a];
                }
                sample.tangent[0][b]=(g22*d1-g12*d2)/determinant;
                sample.tangent[1][b]=(g11*d2-g12*d1)/determinant;
            }
            int it=std::upper_bound(theta.begin(),theta.end(),t)-theta.begin()-1;
            int ip=std::upper_bound(phi.begin(),phi.end(),p)-phi.begin()-1;
            ip=std::max(0,ip);
            const int jp=(ip+1)%phi.size();
            double fp=(p-phi[ip])/((jp==0?phi.front()+2*M_PI:phi[jp])-phi[ip]);
            double ft=0;
            if(it<0) it=0;
            else if(it>=int(theta.size())-1) it=theta.size()-1;
            else ft=(t-theta[it])/(theta[it+1]-theta[it]);
            const int jt=std::min(it+1,int(theta.size())-1);
            for(int b=0;b<4;++b)
            {
                sample.source[b]=sourceBox.low();
                sample.source[b][1]+=(b>=2?jt:it);
                sample.source[b][2]+=(b%2?jp:ip);
                sample.coefficient[b]=(b>=2?ft:1-ft)*(b%2?fp:1-fp);
            }
        }
        return samples;
    }
public:
    void clear()
    {
        m_samples.clear();
        m_points.clear();
    }

    // Initialization consumes classical spherical primitives at cell centers.
    // Keep its source layer convention while sharing the boundary coordinates.
    template<typename T,Proto::MemType MEM>
    void fillPoints(Proto::MBLevelBoxData<T,8,MEM>& destination,
        const Proto::BoxData<double,8,Proto::HOST>& source,
        const std::vector<double>& theta,const std::vector<double>& phi,
        std::vector<Proto::BoxData<T,8,MEM>>* strips=nullptr)
    {
        using namespace Proto;
        static_assert(MEM==HOST,"Boundary input maps are host resident");
        PROTO_ASSERT(!theta.empty() && !phi.empty(),"Empty boundary coordinate arrays");
        unsigned local=0;
        for(auto index:destination.layout())
        {
            auto& out=strips ? (*strips)[local++] : destination[index];
            if (out.box().size()==0) continue;
            const int block=destination.layout().block(index);
            const int n=destination.layout().domain().getBlock(block).sizes()[1];
            Point low=out.box().low(),high=out.box().high();high[0]=low[0];
            for(auto p:Box(low,high))
            {
                auto key=std::make_tuple(block,n,p[1],p[2]);
                auto it=m_points.find(key);
                if(it==m_points.end())
                    it=m_points.emplace(key,makeSamples(block,n,p[1],p[2],
                        source.box(),theta,phi,true)[0]).first;
                const auto& sample=it->second;
                for(int c=0;c<8;++c)
                {
                    double value=0;
                    for(int q=0;q<4;++q)value+=sample.coefficient[q]*source(sample.source[q],c);
                    out(p,c)=value;
                }
            }
        }
    }

    template<typename T,Proto::MemType MEM>
    void fill(Proto::MBLevelBoxData<T,8,MEM>& destination,
        const Proto::BoxData<double,8,Proto::HOST>& source,
        const std::vector<double>& theta,const std::vector<double>& phi,
        std::vector<Proto::BoxData<T,8,MEM>>* strips=nullptr)
    {
        using namespace Proto;
        static_assert(MEM==HOST,"Boundary quadrature currently uses host-resident input maps");
        static_assert(CUBED_SPHERE_SHELL_RADIAL_COORD==0,"Native MHD boundary coordinates require radial axis 0");
        PR_TIMERS("MHDReader::conservedBoundaryAverages");
        PROTO_ASSERT(!theta.empty() && !phi.empty(),"Empty boundary coordinate arrays");
        const double r0=ParseInputs::get_inner_radius(),r1=ParseInputs::get_outer_radius();
        const double stretch=ParseInputs::get_radial_stretch(),gamma=ParseInputs::get_gamma();
        const bool midpoint=ParseInputs::get_spatial_order()==2;
        const double span=stretch==0?r1-r0:(r1-r0)/std::expm1(stretch);
        const auto& layout=destination.layout();
        std::map<int,std::array<double,4>> radial;
        unsigned local=0;
        for(auto index:layout)
        {
            auto& out=strips ? (*strips)[local++] : destination[index];
            if (out.box().size()==0) continue;
            if(out.box().low()[0]>=0) continue;
            const int block=layout.block(index),n=layout.domain().getBlock(block).sizes()[1];
            const int nr=layout.domain().getBlock(block).sizes()[0];
            PROTO_ASSERT(layout.domain().getBlock(block).sizes()[2]==n,"Boundary requires equal angular resolutions");
            Point low=out.box().low(),high=out.box().high();
            high[0]=-1;
            Box ghosts(low,high);
            for(int ir=low[0];ir<0;++ir) if(!radial.count(ir))
            {
                // <s>, <s^2>, <s^4>, <s^(2 gamma)>, where s = r_inner/r.
                std::array<double,4> mean{};
                for(int q=0;q<(midpoint?1:4);++q)
                {
                    const double xi=(ir+(midpoint?.5:radialNodes[q]))/nr;
                    const double radius=r0+span*(stretch==0?xi:std::expm1(stretch*xi));
                    PROTO_ASSERT(radius>0,"Inner radial ghost cell reaches nonpositive radius");
                    const double ratio=r0/radius,r2=ratio*ratio;
                    const double weight=midpoint?1:radialWeights[q];
                    mean[0]+=weight*ratio;
                    mean[1]+=weight*r2;
                    mean[2]+=weight*r2*r2;
                    mean[3]+=weight*std::pow(ratio,2*gamma);
                }
                radial.emplace(ir,mean);
            }
            low[0]=high[0]=-1;
            for(auto p:Box(low,high))
            {
                auto key=std::make_tuple(block,n,p[1],p[2]);
                const Sample* samples;
                if(midpoint)
                {
                    auto it=m_points.find(key);
                    if(it==m_points.end())
                        it=m_points.emplace(key,makeSamples(block,n,p[1],p[2],
                            source.box(),theta,phi,true)[0]).first;
                    samples=&it->second;
                }
                else
                {
                    auto it=m_samples.find(key);
                    if(it==m_samples.end())
                        it=m_samples.emplace(key,makeSamples(block,n,p[1],p[2],
                            source.box(),theta,phi)).first;
                    samples=it->second.data();
                }
                // Angular integrals: density, three native momenta, thermal
                // energy, kinetic + tangential magnetic energy, radial magnetic
                // energy, then the three native magnetic components. Splitting
                // the energy terms preserves their different radial scalings.
                std::array<double,10> mean{};
                for(int s=0;s<(midpoint?1:angularOrder*angularOrder);++s)
                {
                    const auto& sample=samples[s];
                    double w[8]={};
                    for(int q=0;q<4;++q) for(int c=0;c<8;++c)
                        w[c]+=sample.coefficient[q]*source(sample.source[q],c);
                    const double rho=w[0]*1.67262192e-24,weight=sample.weight;
                    mean[0]+=weight*rho;
                    mean[1]+=weight*rho*w[1];
                    mean[4]+=weight*w[4]*1e-12/(gamma-1);
                    mean[5]+=weight*(.5*rho*(w[1]*w[1]+w[2]*w[2]+w[3]*w[3])
                        +1e-12*(w[6]*w[6]+w[7]*w[7])/(8*M_PI));
                    mean[6]+=weight*1e-12*w[5]*w[5]/(8*M_PI);
                    mean[7]+=weight*1e-6*w[5];
                    for(int d=0;d<2;++d)
                    {
                        mean[2+d]+=weight*rho*(sample.tangent[d][0]*w[2]+sample.tangent[d][1]*w[3]);
                        mean[8+d]+=weight*1e-6*(sample.tangent[d][0]*w[6]+sample.tangent[d][1]*w[7]);
                    }
                }
                for(int ir=ghosts.low()[0];ir<0;++ir)
                {
                    p[0]=ir;
                    const auto& scale=radial.at(ir);
                    for(int c=0;c<4;++c) out(p,c)=mean[c]*scale[1];
                    out(p,4)=mean[4]*scale[3]+mean[5]*scale[1]+mean[6]*scale[2];
                    out(p,5)=mean[7]*scale[1];
                    out(p,6)=mean[8]*scale[0];
                    out(p,7)=mean[9]*scale[0];
                }
            }
        }
    }
};
