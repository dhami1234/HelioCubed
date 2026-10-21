#pragma once

namespace HelioCubed
{
inline void writeAMRSliceFile(const std::string& filename,const std::string& zones,
    int iteration,double time)
{
    // Stream rank buffers to one writer; root need not hold the entire slice.
    int error=0,ranks=1;
    std::ofstream file;
    if(procID()==0) {file.open(filename);error=!file;}
#ifdef PR_MPI
    MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    MPI_Bcast(&error,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
    if(error) throw std::runtime_error("Could not open AMR slice "+filename);
    const unsigned long long count=zones.size();
    std::vector<unsigned long long> counts(procID()==0 ? ranks : 0);
#ifdef PR_MPI
    MPI_Gather(&count,1,MPI_UNSIGNED_LONG_LONG,counts.data(),1,MPI_UNSIGNED_LONG_LONG,0,MPI_COMM_WORLD);
#endif
    if(procID()==0)
    {
        file<<std::setprecision(17)<<"TITLE = \"t="<<time<<" step="<<iteration<<"\"\n"
            <<"VARIABLES=\"X\" \"Y\" \"Z\"";
        for(const auto& name : outputPrimitiveNames()) file<<" \""<<name<<"\"";
        file<<"\nDATASETAUXDATA PhysTime = \""<<getPhysTime(time)<<"\"\n";
        file.write(zones.data(),zones.size());
    }
#ifdef PR_MPI
    constexpr size_t chunk=65536;
    std::vector<char> receive(procID()==0 ? chunk : 0);
    for(int rank=1;rank<ranks;++rank)
    {
        const auto size=procID()==0 ? counts[rank] : count;
        if(procID()!=0 && int(procID())!=rank) continue;
        for(unsigned long long offset=0;offset<size;)
        {
            const int n=int(std::min<unsigned long long>(chunk,size-offset));
            if(procID()==0)
            {
                MPI_Recv(receive.data(),n,MPI_CHAR,rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                file.write(receive.data(),n);
            }
            else MPI_Send(zones.data()+offset,n,MPI_CHAR,0,0,MPI_COMM_WORLD);
            offset+=n;
        }
    }
#endif
    if(procID()==0) {file.close();error=!file;}
#ifdef PR_MPI
    MPI_Bcast(&error,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
    if(error) throw std::runtime_error("Could not write AMR slice "+filename);
}

inline void CubedSphereAMR::writeSlices(int iteration,double time)
{
    const std::array<bool,4> requested={ParseInputs::get_slice_X(),ParseInputs::get_slice_Y(),
        ParseInputs::get_slice_Z(),ParseInputs::get_slice_R()};
    if(std::none_of(requested.begin(),requested.end(),[](bool value){return value;})) return;
    PR_TIMERS("HelioCubed::AMR::slices");
    prepareStates(time);
    constexpr double radius=1.5e13; // Retain the established spherical slice radius (cm).
    const auto geometry=cubedSphereShellGeometry();
    // Logical direction normal to x/y/z=0 on each cubed-sphere panel.
    // A negative entry means the panel does not intersect that plane.
    constexpr int directions[6][3]={{1,2,-1},{1,2,-1},{1,-1,2},
        {-1,1,2},{1,-1,2},{-1,1,2}};
    const char* labels[] = {"x","y","z","r"};
    const bool scaled=ParseInputs::get_scale_data_with_r()==1;
    const bool cartesian=ParseInputs::get_output_coordinates()==1;
    for(int axis=0;axis<4;++axis)
    {
        if(!requested[axis]) continue;
        double cut=0.5;
        const bool inShell=axis<3 || (radius>=geometry.r0 && radius<=geometry.r1);
        if(axis==3 && inShell)
        {
            Array<double,DIM> xi,position={radius,0,0};
            f_cubedSphereInverseMap(xi,position,3,geometry);
            cut=std::clamp(xi[0],0.0,1.0);
        }
        auto cutCells=[&](int size)
        {
            double coordinate=cut*size;
            const double face=round(coordinate);
            if(abs(coordinate-face)<=32*std::numeric_limits<double>::epsilon()*std::max(1,size))
                coordinate=face;
            const int upper=std::clamp(int(floor(coordinate)),0,size-1);
            const int lower=coordinate==face && face>0 && face<size ? upper-1 : -1;
            return std::array<int,2>{upper,lower};
        };
        std::ostringstream zones;zones<<std::scientific<<std::setprecision(17);
        int invalid=0;
        if(inShell) for(auto* grid : levels())
        {
            for(auto index : grid->layout)
            {
                const unsigned block=grid->layout.block(index);
                const int normal=axis==3 ? 0 : directions[block][axis];
                if(normal<0) continue;
                const auto sizes=grid->layout.domain().getBlock(block).sizes();
                const auto cuts=cutCells(sizes[normal]);
                const Box valid=grid->layout[index];
                if(cuts[0]<valid.low()[normal] || (cuts[1]<0 ? cuts[0] : cuts[1])>valid.high()[normal]) continue;
                std::array<int,2> tangent;int next=0;
                for(int d=0;d<DIM;++d) if(d!=normal) tangent[next++]=d;
                auto exists=[&](const GridLevel& level,Point cell)
                {return level.layout.find(cell/m_boxSizes[level.number],block)!=*level.layout.end();};
                State fourth;
                if(fourthOrder()) fourth=fourthOrderPrimitives(*grid,index);
                const auto& w=fourthOrder() ? fourth : grid->primitive[index];
                std::vector<Array<double,DIM>> nodes;
                std::map<std::array<int,2>,size_t> nodeIndices;
                std::vector<std::array<size_t,4>> elements;
                std::vector<std::array<double,NUMCOMPS>> values;
                for(int cell : cuts)
                {
                    if(cell<valid.low()[normal] || cell>valid.high()[normal]) continue;
                    Point low=valid.low(),high=valid.high();low[normal]=high[normal]=cell;
                    for(auto p : Box(low,high))
                    {
                        // On a cell face choose one side, preferring a finer
                        // level on either side before the upper same-level cell.
                        Point upper=p;upper[normal]=cuts[0];
                        if(cell!=cuts[0] && exists(*grid,upper)) continue;
                        if(grid->number+1<numLevels())
                        {
                            const auto& child=level(grid->number+1);
                            const auto childCuts=cutCells(child.layout.domain().getBlock(block).sizes()[normal]);
                            Point q=p*child.ratio;q[normal]=childCuts[0];
                            bool refined=exists(child,q);
                            if(childCuts[1]>=0) {q[normal]=childCuts[1];refined=refined || exists(child,q);}
                            if(refined) continue;
                        }
                        std::array<size_t,4> element;
                        constexpr int corners[4][2]={{0,0},{1,0},{1,1},{0,1}};
                        for(int corner=0;corner<4;++corner)
                        {
                            const std::array<int,2> key={p[tangent[0]]+corners[corner][0],p[tangent[1]]+corners[corner][1]};
                            auto entry=nodeIndices.emplace(key,nodes.size());
                            if(entry.second)
                            {
                                Array<double,DIM> xi,position;double jacobian;
                                xi[normal]=cut;
                                for(int d=0;d<2;++d) xi[tangent[d]]=double(key[d])/sizes[tangent[d]];
                                f_cubedSphereMap(position,jacobian,xi,block,geometry);
                                if(axis<3) position[axis]=0;
                                nodes.push_back(position);
                            }
                            element[corner]=entry.first->second+1; // Tecplot is one-based.
                        }
                        elements.push_back(element);
                        // Prepared Cartesian halos supply samples across patch,
                        // panel and refinement boundaries. Interpolate only in
                        // the direction normal to this slice, then project.
                        const double coordinate=cut*sizes[normal]-0.5;
                        const int lower=int(floor(coordinate));
                        const double fraction=coordinate-lower;
                        Point a=p;a[normal]=lower;const Point b=a+Point::Basis(normal);
                        std::array<double,NUMCOMPS> value;
                        for(int c=0;c<NUMCOMPS;++c) value[c]=(1-fraction)*w(a,c)+fraction*w(b,c);
                        value[iBX]+=testProblemBackgroundBx(ParseInputs::get_init_condition_type());
                        Array<double,DIM> xi,position;double jacobian;
                        for(int d=0;d<DIM;++d) xi[d]=d==normal ? cut : (p[d]+0.5)/sizes[d];
                        f_cubedSphereMap(position,jacobian,xi,block,geometry);
                        auto presentation=[&](int c)->double& {return value[c];};
                        formatCartesianOutput(presentation,position[0],position[1],position[2],
                            ParseInputs::get_gamma(),cartesian,scaled);
                        for(double v : value) if(!std::isfinite(v)) invalid=1;
                        values.push_back(value);
                    }
                }
                if(elements.empty()) continue;
                zones<<"ZONE T=\"L"<<grid->number<<"_B"<<block<<"_P"<<index.global()
                    <<"\" ZONETYPE=FEQUADRILATERAL, DATAPACKING=BLOCK, N="<<nodes.size()
                    <<", E="<<elements.size()<<", VARLOCATION=([4-11]=CELLCENTERED)\n";
                for(int d=0;d<DIM;++d) {for(const auto& node : nodes) zones<<node[d]<<' ';zones<<'\n';}
                for(int c=0;c<NUMCOMPS;++c) {for(const auto& value : values) zones<<value[c]<<' ';zones<<'\n';}
                for(const auto& element : elements)
                {for(auto node : element) zones<<node<<' ';zones<<'\n';}
            }
        }
#ifdef PR_MPI
        MPI_Allreduce(MPI_IN_PLACE,&invalid,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif
        if(invalid) throw std::runtime_error("Nonfinite AMR slice sample");
        writeAMRSliceFile(ParseInputs::get_data_file_prefix()+m_outputSuffix+"."+labels[axis]+"."
            +std::to_string(iteration)+".dat",zones.str(),iteration,time);
        if(!inShell) pout(0)<<"AMR spherical slice radius "<<radius<<" cm is outside the shell"<<std::endl;
    }
    pout(0)<<"AMR Tecplot slice files written at iteration "<<iteration<<std::endl;
}
}
