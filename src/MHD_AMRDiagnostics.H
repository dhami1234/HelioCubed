#pragma once

namespace HelioCubed
{
template<bool DivB>
inline void CubedSphereAMR::writePlot(HDF5Handler& h5,const std::string& prefix,double time)
{
    prepareStates(time);
    using PlotLevel=MBLevelBoxData<double,NUMCOMPS+2+int(DivB),HOST>;
    std::vector<std::unique_ptr<PlotLevel>> storage;
    std::vector<const PlotLevel*> output;
    std::vector<const OP::MapT*> maps;
    for(auto* grid : levels())
    {
        storage.push_back(std::make_unique<PlotLevel>(grid->layout,Point::Zeros()));
        auto& plot=*storage.back();
        OP::GeometryScope geometry(grid->geometry);
        for(auto index : grid->layout)
        {
            State fourth;
            if(fourthOrder()) fourth=fourthOrderPrimitives(*grid,index);
            Scalar divergence;
            if constexpr(DivB)
            {
                divergence.define(grid->layout[index]);
                grid->op[index].computeDivB(divergence,grid->primitive[index],!fourthOrder());
            }
            for(auto p : grid->layout[index])
            {
                for(int c=0;c<NUMCOMPS;++c)
                    plot[index](p,c)=fourthOrder() ? fourth(p,c) : grid->primitive[index](p,c);
                plot[index](p,iBX)+=testProblemBackgroundBx(ParseInputs::get_init_condition_type());
                plot[index](p,NUMCOMPS)=covered(grid->number,p,grid->layout.block(index)) ? 0 : 1;
                plot[index](p,NUMCOMPS+1)=grid->number;
                if constexpr(DivB) plot[index](p,NUMCOMPS+2)=divergence(p);
            }
        }
        output.push_back(&plot);maps.push_back(&grid->map);
    }
    std::vector<std::string> names={"density","Vx","Vy","Vz","P","Bx","By","Bz","active","amr_level"};
    if constexpr(DivB) names.push_back("divB"); // Physical G/cm; AMR primitives retain their Cartesian CGS format.
    const std::vector<Point> ratios(m_ratios.begin(),m_ratios.begin()+numLevels()-1);
    h5.writeMBAMRData(names,maps,output,ratios,prefix);
}

inline std::pair<int,std::array<double,NUMCOMPS>> CubedSphereAMR::samplePreparedProbe(
    const Array<double,DIM>& position)
{
    PR_TIMERS("HelioCubed::AMR::sampleProbe");
    std::array<double,NUMCOMPS> values{};
    const double radius=std::hypot(position[0],position[1],position[2]);
    if(!std::isfinite(radius) || radius<ParseInputs::get_inner_radius()
        || radius>ParseInputs::get_outer_radius()) return {-1,values};
    // Layouts are global on every rank. Search all panels, finest first, so a
    // seam point also uses the finest covering panel and has exactly one owner.
    for(int l=int(numLevels())-1;l>=0;--l)
    {
        auto& grid=level(l);
        for(unsigned block=0;block<6;++block)
        {
            Array<double,DIM> xi;
            f_cubedSphereInverseMap(xi,position,block,cubedSphereShellGeometry());
            const auto sizes=grid.layout.domain().getBlock(block).sizes();
            Point cell,lower;Array<double,DIM> fraction;
            bool inside=true;
            for(int d=0;d<DIM;++d)
            {
                if(!std::isfinite(xi[d]) || xi[d]<-1e-12 || xi[d]>1+1e-12) {inside=false;break;}
                const double x=std::clamp(xi[d],0.0,1.0)*sizes[d];
                cell[d]=std::min(int(floor(x)),sizes[d]-1);
                lower[d]=int(floor(x-0.5));fraction[d]=x-0.5-lower[d];
            }
            if(!inside) continue;
            const auto index=grid.layout.find(cell/m_boxSizes[l],block);
            if(index==*grid.layout.end()) continue;
            const unsigned owner=grid.layout.procID(index);
            if(owner==procID())
            {
                State fourth;
                if(fourthOrder()) fourth=fourthOrderPrimitives(grid,index);
                const auto& w=fourthOrder() ? fourth : grid.primitive[index];
                // The containing cell chooses the level. Its prepared halo
                // supplies neighbors across patch, panel and refinement edges.
                for(auto offset : Box(Point::Zeros(),Point::Ones()))
                {
                    double weight=1;
                    for(int d=0;d<DIM;++d) weight*=offset[d] ? fraction[d] : 1-fraction[d];
                    for(int c=0;c<NUMCOMPS;++c) values[c]+=weight*w(lower+offset,c);
                }
                values[iBX]+=testProblemBackgroundBx(ParseInputs::get_init_condition_type());
            }
#ifdef PR_MPI
            MPI_Bcast(values.data(),NUMCOMPS,MPI_DOUBLE,owner,MPI_COMM_WORLD);
#endif
            for(double value : values) if(!std::isfinite(value))
                throw std::runtime_error("Nonfinite AMR probe sample");
            return {l,values};
        }
    }
    return {-1,values};
}

inline void CubedSphereAMR::writeProbes(int iteration,double time,bool giveSpace)
{
    if(!ParseInputs::get_num_probes()) return;
    PR_TIMERS("HelioCubed::AMR::probes");
    prepareStates(time);
    for(int p=0;p<ParseInputs::get_num_probes();++p)
    {
        const auto& trajectory=probeTrajectory(ParseInputs::get_probe_coords()[p]);
        const double epoch=getPhysTime(time);
        const double radius=interpolate(trajectory.time,trajectory.radius,epoch,false)*c_AU;
        const double latitude=interpolate(trajectory.time,trajectory.latitude,epoch,false);
        const double longitude=std::fmod(interpolate(trajectory.time,trajectory.longitude,epoch,false)+181.02+360,360);
        const double theta=(90-latitude)*M_PI/180,phi=longitude*M_PI/180;
        const Array<double,DIM> normal={sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)};
        const Array<double,DIM> position={radius*normal[0],radius*normal[1],radius*normal[2]};
        const auto sample=samplePreparedProbe(position);
        if(sample.first<0)
        {
            pout(0)<<"AMR probe "<<p<<" is outside the shell at iteration "<<iteration<<std::endl;
            continue;
        }
        const auto& w=sample.second;
        // Retain the established probe columns, RTN projection and units.
        auto project=[&](int start,double scale)
        {
            const double x=w[start]*scale,y=w[start+1]*scale,z=w[start+2]*scale;
            return std::array<double,3>{normal[0]*x+normal[1]*y+normal[2]*z,
                -normal[1]*x+normal[0]*y,
                -normal[0]*normal[2]*x-normal[1]*normal[2]*y+(normal[0]*normal[0]+normal[1]*normal[1])*z};
        };
        const auto velocity=project(iVX,1e-5),magnetic=project(iBX,1e6);
        int error=0;
        if(procID()==0)
        {
            std::ofstream file(ParseInputs::get_probe_filenames()[p]+m_outputSuffix,std::ios::app);
            if(giveSpace) file<<'\n';
            file<<std::setw(16)<<std::setprecision(12)<<epoch<<std::setprecision(4);
            for(double value : {radius/c_AU,latitude,longitude,w[iRHO]/c_MP,
                velocity[0],velocity[1],velocity[2],w[iP]*1e12,magnetic[0],magnetic[1],magnetic[2]})
                file<<std::setw(11)<<value;
            file<<'\n';file.close();error=!file;
        }
#ifdef PR_MPI
        MPI_Bcast(&error,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
        if(error) throw std::runtime_error("Could not write AMR probe output");
        pout(0)<<"AMR probe "<<p<<" level = "<<sample.first<<" at iteration "<<iteration<<std::endl;
    }
}
}
