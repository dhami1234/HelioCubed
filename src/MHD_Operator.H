#pragma once

#include "Proto.H"
namespace MHD_Operator {
    using namespace Proto;    
    template<class T, unsigned int C, unsigned char D, unsigned char E, MemType MEM>
    BoxData<T,C,MEM,D,E>
    _cellTensorProduct(
            const BoxData<T,C,MEM,D,E>& a_ql4,
            const BoxData<T,1,MEM,1,1>& a_qr4,
            const BoxData<T,C,MEM,D,E>& a_ql2,
            const BoxData<T,1,MEM,1,1>& a_qr2)
    {
        Box bx = a_ql4.box()&a_qr4.box();
        for (int dir = 0; dir < DIM; dir++)
            {
            bx = bx&a_ql2.box().grow(dir,-1)&a_qr2.box().grow(dir,-1);
            }
        BoxData<T,C,MEM,D,E> quot(bx);
        for (int ee = 0; ee < E;ee++)
            for (int dd = 0; dd < D; dd++)
            for (int cc = 0; cc < C; cc++)
                {
                auto quotcomp = slice(quot,cc,dd,ee);
                auto ql4comp = slice(a_ql4,cc,dd,ee);
                auto ql2comp = slice(a_ql2,cc,dd,ee);
                Operator::_cellProduct(quotcomp,ql4comp,a_qr4,ql2comp,a_qr2);
                }
        return quot;
    }
    template<class T, unsigned int C, unsigned char D, unsigned char E, MemType MEM>
    BoxData<T,C,MEM,D,E>
    _cellTensorProduct2(
                const BoxData<T,C,MEM,D,E>& a_ql2,
                const BoxData<T,1,MEM,1,1>& a_qr2)
    {
        Box bx = a_ql2.box()&a_qr2.box();
        BoxData<T,C,MEM,D,E> quot(bx);
        for (int ee = 0; ee < E;ee++)
            for (int dd = 0; dd < D; dd++)
            for (int cc = 0; cc < C; cc++)
                {
                auto quotcomp = slice(quot,cc,dd,ee);
                auto ql2comp = slice(a_ql2,cc,dd,ee);
                quotcomp=ql2comp*a_qr2;
                }
        return quot;
    }
}