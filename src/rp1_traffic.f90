! =========================================================
subroutine rp1(maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =========================================================
!
! solve Riemann problems for the traffic equation.

! waves: 1
! equations: 1

! Conserved quantities:
!       1 q

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell
! On output, wave contains the waves, 
!            s the speeds, 
!            amdq the  left-going flux difference  A^- \Delta q
!            apdq the right-going flux difference  A^+ \Delta q
!
! Note that the i'th Riemann problem has left state qr(i-1,:)
!                                    and right state ql(i,:)
! From the basic clawpack routine step1, rp is called with ql = qr = q.


      implicit double precision (a-h,o-z)
      dimension   ql(meqn, 1-mbc:maxm+mbc)
      dimension   qr(meqn, 1-mbc:maxm+mbc)
      dimension    s(mwaves, 1-mbc:maxm+mbc)
      dimension wave(meqn,  mwaves,1-mbc:maxm+mbc)
      dimension amdq(meqn, 1-mbc:maxm+mbc) 
      dimension apdq(meqn, 1-mbc:maxm+mbc) 
      logical efix
      common /cparam/ umax
!
      efix = .true.   !# Compute correct flux for transonic rarefactions
!
      do 30 i=2-mbc,mx+mbc
!
!        # Compute the wave and speed
!
         wave(1,1,i) = ql(1,i) - qr(1,i-1)
         s(1,i) = umax * (1.d0 - (qr(1,i-1) + ql(1,i)))
!
!
!        # compute left-going and right-going flux differences:
!        ------------------------------------------------------
!
         amdq(1,i) = dmin1(s(1,i), 0.d0) * wave(1,1,i)
         apdq(1,i) = dmax1(s(1,i), 0.d0) * wave(1,1,i)
!
         if (efix) then
!           # entropy fix for transonic rarefactions:
            sim1 = umax*(1.d0 - 2.d0*ql(1,i-1))
            si = umax*(1.d0 - 2.d0*ql(1,i))
            if (sim1.lt.0.d0 .and. si.gt.0.d0) then
               flux0 = 0.25d0*umax
               fluxim1 = qr(1,i-1)*umax*(1.d0 - qr(1,i-1))
               fluxi   = ql(1,i)*umax*(1.d0 - qr(1,i))
               amdq(1,i) = flux0 - fluxim1
               apdq(1,i) = fluxi - flux0
               endif
            endif
   30   continue
!
      return
      end
