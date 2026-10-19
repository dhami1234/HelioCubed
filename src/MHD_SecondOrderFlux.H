#pragma once
#include "MHD_Roe8Wave.H"

namespace HelioCubed
{
struct SecondOrderFrame
{
    double frame[3][3],norm,normal[3];
};

inline ACCEL_DECORATION SecondOrderFrame secondOrderFrame(const double normal[3])
{
    SecondOrderFrame result{};
    auto& frame=result.frame;
    double& norm=result.norm;
    for(int i=0;i<3;++i)norm+=normal[i]*normal[i];
    norm=sqrt(norm);
    for(int i=0;i<3;++i)frame[0][i]=normal[i]/norm;
    // Cross the normal with the least aligned Cartesian axis.
    int axis=0;
    for(int i=1;i<3;++i)if(fabs(frame[0][i])<fabs(frame[0][axis]))axis=i;
    frame[1][(axis+1)%3]=frame[0][(axis+2)%3];
    frame[1][(axis+2)%3]=-frame[0][(axis+1)%3];
    double tangent=0;
    for(int i=0;i<3;++i)tangent+=frame[1][i]*frame[1][i];
    for(int i=0;i<3;++i)frame[1][i]/=sqrt(tangent);
    for(int i=0;i<3;++i){int j=(i+1)%3,k=(i+2)%3;frame[2][i]=frame[0][j]*frame[1][k]-frame[0][k]*frame[1][j];}
    for (int d=0; d<3; ++d) result.normal[d]=normal[d];
    return result;
}

// This is the original O2 flux arithmetic, shared by the cached HOST loop
// and the general forall path. Reconstruction and Roe/Rusanov are unchanged.
template<class Out,class Left,class Right>
inline ACCEL_DECORATION void secondOrderFaceFlux(Out& out,
    const Left& l,const Right& r,const SecondOrderFrame& geometry,
    double radialArea,double gamma,double entropy,bool roe)
{
    constexpr double inv4pi=1.0/(4.0*M_PI);
    const auto& frame=geometry.frame;
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
    const double area=radialArea*geometry.norm;
    for(int c=0;c<8;++c)numerical[c]=0.5*area*(flux[0][c]+flux[1][c]-diss[c]);
    out(0)=numerical[0];out(4)=numerical[4];
    for(int i=0;i<3;++i)
    {
        out(1+i)=out(5+i)=0;
        for(int j=0;j<3;++j){out(1+i)+=frame[j][i]*numerical[1+j];out(5+i)+=frame[j][i]*numerical[5+j];}
    }
}

template<typename T>
struct LocalFaceState
{
    T data[8];
    T& operator()(int c) {return data[c];}
    const T& operator()(int c) const {return data[c];}
};

template<typename T>
struct StridedFaceState
{
    T* data;
    size_t stride;
    T& operator()(int c) const {return data[c*stride];}
};
}
