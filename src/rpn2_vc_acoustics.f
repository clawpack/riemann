c
c
c     =====================================================
      subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,
     &                  auxl,auxr,wave,s,amdq,apdq)
c     =====================================================
c
c     # Riemann solver for the acoustics equations in 2d, with varying
c     # material properties rho and kappa
c
c     # Note that although there are 3 eigenvectors, the second eigenvalue
c     # is always zero and so we only need to compute 2 waves.  
c     # 
c     # solve Riemann problems along one slice of data.
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # auxl(1,i) holds rho, 
c     # auxl(2,i) holds sound speed c, 
c     #   Here it is assumed that auxl=auxr gives the cell values.
c
c
c     # On output, wave contains the waves,
c     #            s the speeds,
c     #            amdq the  left-going flux difference  A^- \Delta q
c     #            apdq the right-going flux difference  A^+ \Delta q
c
c
c     # This data is along a slice in the x-direction if ixy=1 
c     #                            or the y-direction if ixy=2.
c
c     # Note that the i'th Riemann problem has left state qr(:,i-1)
c     #                                    and right state ql(:,i)
c     # From the basic clawpack routines, this routine is called with ql = qr
c
c
      implicit double precision (a-h,o-z)
c
      dimension wave(meqn, mwaves, 1-mbc:maxm+mbc)
      dimension    s(mwaves, 1-mbc:maxm+mbc)
      dimension   ql(meqn, 1-mbc:maxm+mbc)
      dimension   qr(meqn, 1-mbc:maxm+mbc)
      dimension apdq(meqn, 1-mbc:maxm+mbc)
      dimension amdq(meqn, 1-mbc:maxm+mbc)
      dimension auxl(2, 1-mbc:maxm+mbc)  
      dimension auxr(2, 1-mbc:maxm+mbc)  
c
c     local arrays
c     ------------
      dimension delta(3)
c
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
c
c     # set mu to point to  the component of the system that corresponds
c     # to velocity in the direction of this slice, mv to the orthogonal
c     # velocity.
c
c
      if (ixy.eq.1) then
          mu = 2
          mv = 3
        else
          mu = 3
          mv = 2
        endif
c
c     # note that notation for u and v reflects assumption that the 
c     # Riemann problems are in the x-direction with u in the normal
c     # direciton and v in the orthogonal direcion, but with the above
c     # definitions of mu and mv the routine also works with ixy=2
c
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
         zi = auxl(1,i)*auxl(2,i)
         zim = auxl(1,i-1)*auxl(2,i-1)

         a1 = (-delta(1) + zi*delta(2)) / (zim + zi)
         a2 =  (delta(1) + zim*delta(2)) / (zim + zi)

c
c        # Compute the waves.
c
         wave(1,1,i) = -a1*zim
         wave(mu,1,i) = a1
         wave(mv,1,i) = 0.d0
         s(1,i) = -auxl(2,i-1)
c
         wave(1,2,i) = a2*zi
         wave(mu,2,i) = a2
         wave(mv,2,i) = 0.d0
         s(2,i) = auxl(2,i)
c
   20    continue
c
c
c
c     # compute the leftgoing and rightgoing flux differences:
c     # Note s(1,i) < 0   and   s(2,i) > 0.
c
      do 220 m=1,meqn
         do 220 i = 2-mbc, mx+mbc
            amdq(m,i) = s(1,i)*wave(m,1,i)
            apdq(m,i) = s(2,i)*wave(m,2,i)
  220       continue
c
      return
      end
