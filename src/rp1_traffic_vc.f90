! =========================================================
subroutine rp1(maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! =========================================================
!
! solve Riemann problems for the traffic equation.
! with variable speed limit umax stored in aux(i,1)

! waves: 1
! equations: 1
! aux fields: 1

! Conserved quantities:
!       1 q

! Auxiliary variables:
!       1 u_max

! This is an fwave solver.
!
! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell
! On output, wave contains the waves, 
!            s the speeds, 
!            amdq the  left-going flux difference  A^- \Delta q
!            apdq the right-going flux difference  A^+ \Delta q
!
! Note that the i'th Riemann problem has left state qr(:,i-1)
!                                    and right state ql(:,i)
! From the basic clawpack routine step1, rp is called with ql = qr = q.


      implicit double precision (a-h,o-z)
      dimension   ql(meqn,1-mbc:mx+mbc)
      dimension   qr(meqn,1-mbc:mx+mbc)
      dimension    s(mwaves,1-mbc:mx+mbc)
      dimension fwave(meqn,mwaves,1-mbc:mx+mbc)
      dimension amdq(meqn,1-mbc:mx+mbc)
      dimension apdq(meqn,1-mbc:mx+mbc)
      dimension  auxl(maux,1-mbc:mx+mbc)
      dimension  auxr(maux,1-mbc:mx+mbc)
      common /comlxf/ alxf
!
!
!
      do 30 i=2-mbc,mx+mbc
!
!        # Compute the fwave and speed, and fluctuations
!
!        # compute flux in each cell and flux difference:
         fim1 = auxl(1,i-1)*qr(1,i-1)*(1.d0 - qr(1,i-1))
         fi = auxl(1,i)*ql(1,i)*(1.d0 - ql(1,i))
         fwave(1,1,i) = fi - fim1

!        # compute characteristic speed in each cell:
         sim1 = auxl(1,i-1)*(1.d0 - 2.d0*ql(1,i-1))
         si = auxl(1,i)*(1.d0 - 2.d0*ql(1,i))

         if (sim1 .lt. 0.d0 .and. si .le. 0.d0) then
!             # left-going
              s(1,i) = sim1
              amdq(1,i) = fwave(1,1,i)
              apdq(1,i) = 0.d0
            else if (sim1 .ge. 0.d0 .and. si .gt. 0.d0) then
!             # right-going
              s(1,i) = si
              amdq(1,i) = 0.d0
              apdq(1,i) = fwave(1,1,i)
            else if (sim1 .lt. 0.d0 .and. si .gt. 0.d0) then
!             # transonic rarefaction
!             # split fwave between amdq and apdq:
              s(1,i) = 0.5d0*(sim1 + si)
              dq = ql(1,i) - qr(1,i-1)

!             # entropy fix:  (perhaps doesn't work for all cases!!!)
!             # This assumes the flux in the transonic case should
!             # correspond to q=0.5 on the side with the smaller umax value.
              f0 = dmin1(auxl(1,i-1),auxl(1,i))*0.25d0
              amdq(1,i) = f0 - fim1
              apdq(1,i) = fi - f0

            else
!             # transonic shock
              s(1,i) = 0.5d0*(sim1 + si)
              if (s(1,i) .lt. 0.d0) then 
                   amdq(1,i) = fwave(1,1,i)
                   apdq(1,i) = 0.d0
                else if (s(1,i) .gt. 0.d0) then 
                   amdq(1,i) = 0.d0
                   apdq(1,i) = fwave(1,1,i)
                else
                   amdq(1,i) = 0.5d0 * fwave(1,1,i) 
                   apdq(1,i) = 0.5d0 * fwave(1,1,i)
                endif
            endif
!
   30    continue
!
      return
      end

