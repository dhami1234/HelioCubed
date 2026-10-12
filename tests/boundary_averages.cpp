#include "Proto.H"
#include "Inputs_Parsing.H"
#include "MHD_BoundaryAverages.H"
#include <iostream>
#include <limits>
using namespace Proto;

// Run with inputs_boundary_averages. No external boundary file is needed.
int main(int argc,char** argv)
{
#ifdef PR_MPI
    MPI_Init(&argc,&argv);
#endif
    int result=0;
    {
        ParseInputs::getInstance().parsenow(argc,argv);
        const int n=ParseInputs::get_domainSize(),nr=ParseInputs::get_thickness();
        auto domain=CubedSphereShell::Domain(n,nr,0);
        MBDisjointBoxLayout layout(domain,Point(ParseInputs::get_boxSize_rad(),
            ParseInputs::get_boxSize_nonrad(),ParseInputs::get_boxSize_nonrad()));
        MBLevelBoxData<double,8,HOST> averages(layout,Point::Ones(8));
        MBLevelBoxData<double,8,HOST> points(layout,Point(1,8,8));
        BoxData<double,8,HOST> source(Box(Point(1,180,360)));
        std::vector<double> theta(180),phi(360);
        for(int j=0;j<180;++j)theta[j]=(j+.5)*M_PI/180;
        for(int k=0;k<360;++k)phi[k]=k*M_PI/180;
        const double values[8]={700,5e7,0,0,8e4,600,0,0};
        for(auto p:source.box())for(int c=0;c<8;++c)source(p,c)=values[c];
        MHD_BoundaryAverages bc;
        bc.fill(averages,source,theta,phi);
        bc.fillPoints(points,source,theta,phi);
        double error=0;long checks=0;
        auto check=[&](double value,double expected,double scale=1.)
        {
            const double difference=std::abs(value-expected)/scale;
            if(!std::isfinite(difference))error=std::numeric_limits<double>::infinity();
            else error=std::max(error,difference);
        };
        const double r0=ParseInputs::get_inner_radius(),r1=ParseInputs::get_outer_radius();
        const double dr=(r1-r0)/nr,rho=values[0]*1.67262192e-24;
        const double pressure=values[4]*1e-12,br=values[5]*1e-6;
        if(ParseInputs::get_radial_stretch()!=0||ParseInputs::get_gamma()!=1.5)result=1;
        const double scale[8]={rho,rho*values[1],rho*values[1],rho*values[1],
            2*pressure+.5*rho*values[1]*values[1]+br*br/(8*M_PI),br,br,br};
        for(auto index:layout)
        {
            for(auto p:averages[index].box())
            {
                if(p[0]>=0)continue;
                const double lo=r0+p[0]*dr,hi=lo+dr;
                const double s2=r0*r0/(lo*hi),s3=s2*r0*(lo+hi)/(2*lo*hi);
                const double s4=s2*s2*(lo*lo+lo*hi+hi*hi)/(3*lo*hi);
                const double expected[8]={rho*s2,rho*values[1]*s2,0,0,
                    2*pressure*s3+.5*rho*values[1]*values[1]*s2+br*br*s4/(8*M_PI),br*s2,0,0};
                for(int c=0;c<8;++c)check(averages[index](p,c),expected[c],scale[c]);
                ++checks;
            }
        }
        // Verify centered sampling independently, including periodic longitude,
        // polar cells, every patch's source strip, and angular corner ghosts.
        for(double offset:{0.,.5})
        {
          for(int k=0;k<360;++k)phi[k]=(k+offset)*M_PI/180;
          bc.clear(); // A reloaded file can have a different source-grid origin.
          for(auto p:source.box())for(int c=0;c<8;++c)
              source(p,c)=c+1+.2*std::sin(theta[p[1]])*std::cos(phi[p[2]]);
          bc.fillPoints(points,source,theta,phi);
          for(auto index:layout)
          {
            int block=layout.block(index);Point low=points[index].box().low(),high=points[index].box().high();high[0]=low[0];
            for(auto p:Box(low,high))
            {
                Array<double,DIM> xi,xyz;xi[0]=0;xi[1]=(p[1]+.5)/n;xi[2]=(p[2]+.5)/n;double jac;
                f_cubedSphereMap(xyz,jac,xi,block,cubedSphereShellGeometry());
                const double radius=xyz.norm(),t=std::acos(std::max(-1.,std::min(1.,xyz[2]/radius)));
                double lon=std::atan2(xyz[1],xyz[0]);if(lon<phi[0])lon+=2*M_PI;
                int j=std::max(0,std::min(179,int(std::floor(t*180/M_PI-.5))));
                int j1=std::min(j+1,179),k=int(std::floor((lon-phi[0])*180/M_PI))%360,k1=(k+1)%360;
                double a=j==j1?0:std::max(0.,std::min(1.,(t-theta[j])/(theta[j1]-theta[j])));
                double b=(lon-phi[k])/(M_PI/180);
                const double value=.2*((1-a)*std::sin(theta[j])+a*std::sin(theta[j1]))*((1-b)*std::cos(phi[k])+b*std::cos(phi[k1]));
                for(int c=0;c<8;++c)check(points[index](p,c),c+1+value);
                ++checks;
            }
          }
        }
#ifdef PR_MPI
        double globalError;MPI_Allreduce(&error,&globalError,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);error=globalError;
        long total;MPI_Allreduce(&checks,&total,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);checks=total;
#endif
        if(!std::isfinite(error)||error>5e-12)result=1;
        if(!procID())std::cout<<"Boundary average checks: "<<checks<<", maximum scaled error: "<<error<<std::endl;
    }
#ifdef PR_MPI
    MPI_Finalize();
#endif
    return result;
}
