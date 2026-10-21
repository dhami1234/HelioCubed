#pragma once

namespace HelioCubed
{
inline void CubedSphereAMR::advanceRungeKuttaSynchronous(double dt,double time,int order)
{
    PR_TIMERS("HelioCubed::AMR::RungeKuttaSynchronous");
    // Match the existing single-level Kutta RK3 and classical RK4 tableaus.
    static constexpr double a3[4][4]={{0,0,0,0},{0.5,0,0,0},{-1,2,0,0},{0,0,0,0}};
    static constexpr double a4[4][4]={{0,0,0,0},{0.5,0,0,0},{0,0.5,0,0},{0,0,1,0}};
    static constexpr double c3[4]={0,0.5,1,0},b3[4]={1.0/6,2.0/3,1.0/6,0};
    static constexpr double c4[4]={0,0.5,0.5,1},b4[4]={1.0/6,1.0/3,1.0/3,1.0/6};
    const auto& a=order==3 ? a3 : a4;
    const auto& c=order==3 ? c3 : c4;
    const auto& b=order==3 ? b3 : b4;
    for(unsigned l=0;l<numLevels();++l) recordLevelStep(l,dt,time);
    for(auto* grid : levels()) for(int stage=0;stage<order;++stage)
        if(!grid->rk[stage]) grid->rk[stage]=std::make_unique<Level>(grid->layout,Point::Zeros());
    for(int stage=0;stage<order;++stage)
    {
        for(auto* grid : levels())
        {
            amrCopyValid(grid->state,grid->stage);
            for(int previous=0;previous<stage;++previous)
                if(a[stage][previous]) grid->stage.increment(*grid->rk[previous],dt*a[stage][previous]);
        }
        for(unsigned l=numLevels()-1;l>0;--l) restrictTo(l,level(l).stage,level(l-1).stage);
        for(unsigned l=0;l<numLevels();++l) prepare(level(l),level(l).stage,time+c[stage]*dt,l ? &level(l-1) : nullptr);
        for(auto* grid : levels()) evaluate(*grid);
        for(auto& reg : m_registers) {reg->initCoarseRegisters();reg->initFineRegisters();reg->initIncrement();}
        for(auto* grid : levels()) accumulateFlux(*grid,1.0);
        for(unsigned l=numLevels()-1;l>0;--l)
        {
            applyReflux(l-1,level(l-1).rhs,1.0,time+c[stage]*dt,"stage",stage+1);
            restrictTo(l,level(l).rhs,level(l-1).rhs);
        }
        for(auto* grid : levels()) amrCopyValid(grid->rhs,*grid->rk[stage]);
    }
    for(auto* grid : levels()) for(int stage=0;stage<order;++stage) grid->state.increment(*grid->rk[stage],dt*b[stage]);
    averageDown();
}

inline void CubedSphereAMR::prepareExtrapolatedRK2Stage(unsigned l,double start,double dt,int stage)
{
    auto& grid=level(l);
    if(grid.timeInterpolated)
    {
        if(stage==0 && grid.preparedChildStage==-1 && grid.preparedTime==start) return;
        if(grid.preparedChildStart==start && grid.preparedChildDt==dt && grid.preparedChildStage==stage) return;
    }
    const double end=grid.startTime+grid.stepDt;
    const double tolerance=32*std::numeric_limits<double>::epsilon()*std::max({1.0,abs(start),abs(end)});
    if(!grid.old || !(grid.stepDt>0) || start<grid.startTime-tolerance || start+dt>end+tolerance)
        throw std::runtime_error("Extrapolated AMR stage lies outside its parent interpolation interval");
    const double chi=std::clamp((start-grid.startTime)/grid.stepDt,0.0,1.0);
    for(auto index : grid.layout) for(auto p : grid.layout[index]) for(int c=0;c<NUMCOMPS;++c)
    {
        double value=(*grid.old)[index](p,c);
        if(chi>0 || stage)
        {
            const double first=(*grid.predictor)[index](p,c);
            const double difference=grid.rhs[index](p,c)-grid.first[index](p,c);
            value+=chi*grid.stepDt*(first+0.5*chi*difference);
            // Match the fine Euler predictor at stage one, rather than the
            // exact solution at its endpoint. This keeps limiter decisions
            // consistent across a coarse-fine boundary in smooth flow.
            if(stage) value+=dt*(first+chi*difference);
        }
        grid.stage[index](p,c)=value;
    }
    if(l) prepareExtrapolatedRK2Stage(l-1,start,dt,stage);
    if(fourthOrder()) prepareFourthOrder(grid,grid.stage,start+stage*dt,l ? &level(l-1) : nullptr,true);
    else prepare(grid,grid.stage,start+stage*dt,l ? &level(l-1) : nullptr);
    grid.timeInterpolated=true;grid.preparedTime=start+stage*dt;
    grid.preparedChildStart=start;grid.preparedChildDt=dt;grid.preparedChildStage=stage;
}

inline void CubedSphereAMR::advanceExtrapolated(double dt,double time,int order)
{
    // Extrapolate the conservative composite RK2 map, whose first derivative
    // includes instantaneous fine-face fluxes and restriction. A coarse-only
    // RK4 dense predictor instead has an O(dt*h^p) interface coupling error.
    // For RK3, one/two RK2 steps use weights (-1/3,4/3) to cancel dt^3.
    // For RK4, one/two/three steps use (1/12,-4/3,9/4) to cancel dt^3
    // and dt^4. Every branch uses the same multirate pattern
    // so those errors have a common expansion as its step size decreases.
    m_extrapolationSubsteps.clear();
    double parentDt=dt;
    for(unsigned l=0;l+1<numLevels();++l)
    {
        const double fraction=parentDt/m_effectiveDt[l+1];
        const double roundoff=32*std::numeric_limits<double>::epsilon()*std::max(1.0,fraction);
        if(!std::isfinite(fraction) || fraction>maxRatio(m_ratios[l])+roundoff)
            throw std::runtime_error("Extrapolated subcycling exceeds its child CFL bound");
        const int count=std::max(1,int(ceil(fraction-roundoff)));
        if(count>maxRatio(m_ratios[l])) throw std::runtime_error("Extrapolated subcycling exceeds its child CFL bound");
        m_extrapolationSubsteps.push_back(count);parentDt/=count;
    }
    struct Reset
    {
        std::vector<int>& steps;bool& capture;bool& reuse;
        ~Reset() {steps.clear();capture=false;reuse=false;}
    } reset{m_extrapolationSubsteps,m_captureInitial,m_reuseInitial};
    std::vector<Level*> initial,increment;
    for(auto* grid : levels())
    {
        for(int k=0;k<2;++k) if(!grid->rk[k])
            grid->rk[k]=std::make_unique<Level>(grid->layout,Point::Zeros());
        initial.push_back(grid->rk[0].get());increment.push_back(grid->rk[1].get());
        amrCopyValid(grid->state,*initial.back());amrSetValLocal(*increment.back(),0);
    }
    static constexpr double weights3[3]={-1.0/3,4.0/3,0},weights4[3]={1.0/12,-4.0/3,9.0/4};
    const auto& weights=order==3 ? weights3 : weights4;
    for(int divisions=1;divisions<order;++divisions)
    {
        m_captureInitial=divisions==1;m_reuseInitial=divisions>1;
        for(auto* grid : levels()) grid->firstStageUsed=false;
        if(divisions>1) for(unsigned l=0;l<numLevels();++l) amrCopyValid(*initial[l],level(l).state);
        const double step=dt/divisions;
        if(ParseInputs::get_amr_verbosity()>=2 && procID()==0)
        {
            std::ostringstream line;
            line<<std::setprecision(17)<<"AMR extrapolation pass = "<<divisions
                <<" coarse_steps = "<<divisions<<" dt = "<<step
                <<" time = "<<time<<" target_time = "<<time+dt<<'\n';
            pout(0)<<line.str();
        }
        for(int part=0;part<divisions;++part)
        {beginStep(0,step,time+part*step);completeStep(0);}
        const double weight=weights[divisions-1];
        for(unsigned l=0;l<numLevels();++l)
            for(auto index : level(l).layout) for(auto p : level(l).layout[index]) for(int c=0;c<NUMCOMPS;++c)
                (*increment[l])[index](p,c)+=weight*(level(l).state[index](p,c)-(*initial[l])[index](p,c));
    }
    for(unsigned l=0;l<numLevels();++l)
    {
        amrCopyValid(*initial[l],level(l).state);level(l).state.increment(*increment[l],1.0);
        level(l).timeInterpolated=false;
    }
    averageDown();
}

}
