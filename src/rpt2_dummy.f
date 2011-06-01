c
c
c     =====================================================
      subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,
     &                  imp,asdq,bmasdq,bpasdq)
c     =====================================================
      implicit double precision (a-h,o-z)
c
c     # Dummy transverse Riemann solver, for use in dimensionally-split algorithm.
c
       write(*,*) 'Error: Dummy transverse Riemann solver called!'
      return
      end
