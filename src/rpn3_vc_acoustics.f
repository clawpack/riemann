c
c
c
c     ==================================================================
      subroutine rpn3(ixyz,maxm,meqn,mwaves,mbc,mx,ql,qr,
     &			auxl,auxr,maux,wave,s,amdq,apdq)
c     ==================================================================
c
c     # Riemann solver for the acoustics equations in 3d, with varying
c     # material properties.
c
c     # auxl(i,1) holds impedance Z,
c     # auxl(i,2) holds sound speed c, 
c
c     # Note that although there are 4 eigenvectors, two eigenvalues are
c     # always zero and so we only need to compute 2 waves.
c      
c     # Solve Riemann problems along one slice of data.
c     # This data is along a slice in the x-direction if ixyz=1
c     #                               the y-direction if ixyz=2.
c     #                               the z-direction if ixyz=3.
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routines, this routine is called with ql = qr
c
c
      implicit real*8(a-h,o-z)
      dimension wave(meqn,mwaves,1-mbc:maxm+mbc)
      dimension    s(mwaves,1-mbc:maxm+mbc)
      dimension   ql(meqn,1-mbc:maxm+mbc)
      dimension   qr(meqn,1-mbc:maxm+mbc)
      dimension amdq(meqn,1-mbc:maxm+mbc)
      dimension apdq(meqn,1-mbc:maxm+mbc)
      dimension auxl(maux,1-mbc:maxm+mbc)
      dimension auxr(maux,1-mbc:maxm+mbc)
c
c     local arrays
c     ------------
      dimension delta(3)

c
c
c
c     # set mu to point to  the component of the system that corresponds
c     # to velocity in the direction of this slice, mv to the orthogonal
c     # velocity.
c
c
      if (ixyz.eq.1) then
          mu = 2
          mv = 3
          mw = 4
        else if (ixyz.eq.2) then
          mu = 3
          mv = 4
          mw = 2
        else if (ixyz.eq.3) then
          mu = 4
          mv = 2
          mw = 3
        endif

c
c     # split the jump in q at each interface into waves
c     # The jump is split into a leftgoing wave traveling at speed -c
c     # relative to the material properties to the left of the interface,
c     # and a rightgoing wave traveling at speed +c
c     # relative to the material properties to the right of the interface,
c
c     # find a1 and a2, the coefficients of the 2 eigenvectors:
      do 20 i = 2-mbc, mx+mbc
         delta(1) = ql(1,i) - qr(1,i-1)
         delta(2) = ql(mu,i) - qr(mu,i-1)
c        # impedances:
         zi = auxl(1,i)
         zim = auxl(1,i-1)

         a1 = (-delta(1) + zi*delta(2)) / (zim + zi)
         a2 =  (delta(1) + zim*delta(2)) / (zim + zi)


c
c        # Compute the waves.
c
         wave(1,1,i) = -a1*zim
         wave(mu,1,i) = a1
         wave(mv,1,i) = 0.d0
         wave(mw,1,i) = 0.d0
         s(1,i) = -auxl(2,i-1)
c
         wave(1,2,i) = a2*zi
         wave(mu,2,i) = a2
         wave(mv,2,i) = 0.d0
         wave(mw,2,i) = 0.d0
         s(2,i) = auxl(2,i)
c
   20    continue
c
c
c     # compute the leftgoing and rightgoing flux differences:
c     # Note s(i,1) < 0   and   s(i,2) > 0.
c
      do 220 m=1,meqn
         do 220 i = 2-mbc, mx+mbc
            amdq(m,i) = s(1,i)*wave(m,1,i)
            apdq(m,i) = s(2,i)*wave(m,2,i)
  220       continue
c
c
      return
      end

