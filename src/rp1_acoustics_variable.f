c     =====================================================
      subroutine rp(maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &                  wave,s,amdq,apdq)
c     =====================================================
c
c     # Riemann solver for the acoustics equations in 1d, with
c     #  variable coefficients (heterogeneous media)
c
c     # auxl(1,i) should contain the impedance Z in cell i
c     # auxl(2,i) should contain the sound speed c in cell i
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # On output, wave contains the waves,
c     #            s the speeds,
c     #
c     #            amdq = A^- Delta q, 
c     #            apdq = A^+ Delta q,
c     #                   the decomposition of the flux difference
c     #                       f(qr(i-1)) - f(ql(i))
c     #                   into leftgoing and rightgoing parts respectively.
c     #
c
c     # Note that the i'th Riemann problem has left state qr(:,i-1)
c     #                                    and right state ql(:,i)
c     # From the basic clawpack routines, this routine is called with ql = qr
c
c
      implicit double precision (a-h,o-z)
c
      dimension auxl(2, 1-mbc:maxm+mbc)
      dimension auxr(2, 1-mbc:maxm+mbc)
      dimension wave(meqn, mwaves, 1-mbc:maxm+mbc)
      dimension    s(mwaves, 1-mbc:maxm+mbc)
      dimension   ql(meqn, 1-mbc:maxm+mbc)
      dimension   qr(meqn, 1-mbc:maxm+mbc)
      dimension apdq(meqn, 1-mbc:maxm+mbc)
      dimension amdq(meqn, 1-mbc:maxm+mbc)
c
c     local arrays
c     ------------
      dimension delta(2)
c
c
c     # split the jump in q at each interface into waves
c
c     # find a1 and a2, the coefficients of the 2 eigenvectors:
      do 20 i = 2-mbc, mx+mbc
         delta(1) = ql(1,i) - qr(1,i-1)
         delta(2) = ql(2,i) - qr(2,i-1)

c        # impedances:
         zi = auxl(1,i)
         zim = auxr(1,i-1)

         a1 = (-delta(1) + zi*delta(2)) / (zim + zi)
         a2 =  (delta(1) + zim*delta(2)) / (zim + zi)
c
c        # Compute the waves.
c
         wave(1,1,i) = -a1*zim
         wave(2,1,i) = a1
         s(1,i) = -auxr(2,i-1)
c
         wave(1,2,i) = a2*zi
         wave(2,2,i) = a2
         s(2,i) = auxl(2,i)
c
   20    continue
c
c
c     # compute the leftgoing and rightgoing fluctuations:
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
