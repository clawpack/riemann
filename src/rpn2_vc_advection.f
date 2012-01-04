
c
c
c     =====================================================
      subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,
     &			auxl,auxr,wave,s,amdq,apdq)
c     =====================================================
c
c     # Riemann-solver for the advection equation
c     #    q_t  +  u*q_x + v*q_y = 0
c     # where u and v are a given velocity field.
c
c       -----------------------------------------------------------
c     # In advective form, with interface velocities specified in
c     # the auxiliary variable   
c     # aux(i,j,1)  =  u-velocity at left edge of cell (i,j)
c     # aux(i,j,2)  =  v-velocity at bottom edge of cell (i,j)
c       -----------------------------------------------------------
c
c     # solve Riemann problems along one slice of data.
c     # This data is along a slice in the x-direction if ixy=1
c     #                            or the y-direction if ixy=2.
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # On output, wave contains the waves, s the speeds, 
c     # and amdq, apdq the left-going and right-going flux differences,
c     # respectively.  Note that in this advective form, the sum of
c     # amdq and apdq is not equal to a difference of fluxes except in the
c     # case of constant velocities.
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routines, this routine is called with ql = qr
c
c
      implicit real*8(a-h,o-z)
c
      dimension wave(meqn, mwaves, 1-mbc:maxm+mbc)
      dimension    s(mwaves, 1-mbc:maxm+mbc)
      dimension   ql(meqn, 1-mbc:maxm+mbc)
      dimension   qr(meqn, 1-mbc:maxm+mbc)
      dimension apdq(meqn,1-mbc:maxm+mbc)
      dimension amdq(meqn,1-mbc:maxm+mbc)
      dimension auxl(3,1-mbc:maxm+mbc)
      dimension auxr(3,1-mbc:maxm+mbc)
c
c
c     # Set wave, speed, and flux differences:
c     ------------------------------------------
c
      do 30 i = 2-mbc, mx+mbc
         wave(1,1,i) = ql(1,i) - qr(1,i-1)
	 s(1,i) = auxl(ixy,i)
c        # The flux difference df = s*wave  all goes in the downwind direction:
	 amdq(1,i) = dmin1(auxl(ixy,i), 0.d0) * wave(1,1,i)
	 apdq(1,i) = dmax1(auxl(ixy,i), 0.d0) * wave(1,1,i)
   30    continue
c
      return
      end
