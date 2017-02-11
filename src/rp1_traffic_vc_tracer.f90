! =========================================================
subroutine rp1(maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =========================================================
!
! Solve Riemann problems for the variable speed limit LWR traffic model
! with a tracer:

! q_t + (v q(1-q))_x = 0
! p_t + v(1-q) p_x = 0

! Here q is the traffic density, p is the tracer, and 
! v is the speed limit (which may vary with x and t).
! The tracer (p) can be useful for showing vehicle trajectories.

! waves: 2
! equations: 2

! Conserved quantities:
!       1 q
!       2 p

! Auxiliary variables:
!       1 v

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
      dimension auxl(maux,1-mbc:maxm+mbc)
      dimension auxr(maux,1-mbc:maxm+mbc)


      do 30 i=2-mbc,mx+mbc
         q_l = qr(1,i-1)
         q_r = ql(1,i)
         p_l = qr(2,i-1)
         p_r = ql(2,i) 
         v_l = auxr(1,i-1)
         v_r = auxl(1,i)
         wave(1,1,i) = q_r - q_l
         wave(2,1,i) = 0
         wave(1,2,i) = 0
         wave(2,2,i) = p_r - p_l
         
         if ((1.d0 - (q_l + q_r)) .gt. 0.d0) then
            s(1,i) = v_r * (1.d0 - (q_l + q_r))
         else
            s(1,i) = v_l * (1.d0 - (q_l + q_r))
         endif
      
         s(2,i) = v_r*(1.d0 - q_r)

         amdq(1,i) = dmin1(s(1,i), 0.d0) * wave(1,1,i)
         apdq(1,i) = dmax1(s(1,i), 0.d0) * wave(1,1,i)
         amdq(2,i) = 0
         apdq(2,i) = s(2,i)*wave(2,2,i)

         sim1 = v_l*(1.d0 - 2.d0*q_l)
         si   = v_r*(1.d0 - 2.d0*q_r)

         ! entropy fix for transonic rarefactions:
         if (sim1.lt.0.d0 .and. si.gt.0.d0) then
             flux0im1 = 0.25d0 * v_l
             flux0i   = 0.25d0 * v_r
             fluxim1 = v_l*q_l*(1.d0 - q_l)
             fluxi   = v_r*q_r*(1.d0 - q_r)
             amdq(1,i) = flux0im1 - fluxim1
             apdq(1,i) = fluxi - flux0i
         endif
   30   continue

      return
      end
