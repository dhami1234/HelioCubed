#pragma once

namespace HelioCubed
{
inline bool CubedSphereAMR::insertCME(double time)
{
    if(!m_cmeSchedule.due(time)) return false;
    PR_TIMERS("HelioCubed::AMR::insertCME");
    prepareStates(time);
    long long affected=0;
    for(auto* grid : levels()) for(auto index : grid->layout)
    {
        const Box valid=grid->layout[index],box=fourthOrder() ? valid.grow(1) : valid;
        Scalar radius(box),volume(box);Vector dr(box),area(box);
        auto& op=grid->op[index];const auto block=grid->layout.block(index);
        op.radialMetrics(radius,dr,area,volume,box);
        auto position=forall_p<double,DIM,HOST>(f_cubedSphereMap3,box,radius,op.dx(),0.5,0.5,block);
        State cme(box,0.0),delta(box,0.0),fourth;
        forallInPlace_p(define_CME,cme,position);
        if(fourthOrder()) fourth=fourthOrderPrimitives(*grid,index);
        const auto& w=fourthOrder() ? fourth : grid->primitive[index];
        const double b0=testProblemBackgroundBx(ParseInputs::get_init_condition_type());
        for(auto p : box)
        {
            if(cme(p,iRHO)==0) continue;
            if(valid.containsPoint(p) && !covered(grid->number,p,block)) ++affected;
            // Match the existing CME model: preserve density/thermal pressure,
            // replace velocity and total B inside its mask. Evolve only the
            // conserved increment so cells outside the event remain untouched.
            double energy=0;
            for(int d=0;d<DIM;++d)
            {
                const double oldV=w(p,iVX+d),newV=cme(p,iVX+d);
                const double oldB=w(p,iBX+d)+(d==0 ? b0 : 0),newB=cme(p,iBX+d);
                delta(p,iMOMX+d)=w(p,iRHO)*(newV-oldV);
                delta(p,iBX+d)=newB-oldB;
                energy+=0.5*w(p,iRHO)*(newV-oldV)*(newV+oldV)
                    +(newB-oldB)*(newB+oldB)/(8*M_PI);
            }
            delta(p,iE)=energy;
        }
        State average;
        if(fourthOrder()) average=(1.0*Shift::Zeros()+(1.0/24)*Stencil<double>::Laplacian())(delta,valid);
        const auto& j=grid->jacobian[index];
        for(auto p : valid) for(int c=0;c<NUMCOMPS;++c)
        {
            if(c==iRHO) continue;
            double increment=j(p)*(fourthOrder() ? average(p,c) : delta(p,c));
            if(fourthOrder()) for(int d=0;d<DIM;++d)
            {
                const auto e=Point::Basis(d);
                increment+=(j(p+e)-j(p-e))*(delta(p+e,c)-delta(p-e,c))/48;
            }
            grid->state[index](p,c)+=increment;
        }
    }
    averageDown();
    for(auto* grid : levels()) if(!OP::validConservedState(grid->state))
        throw std::runtime_error("Invalid conserved AMR state after CME insertion");
    m_cmeSchedule.setInserted(true);m_effectiveDt.clear();
#ifdef PR_MPI
    MPI_Allreduce(MPI_IN_PLACE,&affected,1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
#endif
    // Tag the inserted feature now, rather than waiting for the next cadence.
    regrid(time);
    pout(0)<<"CME inserted at simulation time "<<std::setprecision(17)<<time<<" s"<<std::endl;
    pout(0)<<"AMR CME affected composite cells = "<<affected<<std::endl;
    return true;
}
}
