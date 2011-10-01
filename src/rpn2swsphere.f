
c
c
c     =====================================================
      subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &                wave,s,amdq,apdq)
c     =====================================================
c
c     # Roe-solver for the 2D shallow water equations
c     #  on the sphere, using 3d Cartesian representation of velocities
c
c     # solve Riemann problems along one slice of data.
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # This data is along a slice in the x-direction if ixy=1
c     #                            or the y-direction if ixy=2.
c     # On output, wave contains the waves, s the speeds,
c     # and amdq, apdq the decomposition of the flux difference
c     #   f(qr(i-1)) - f(ql(i))
c     # into leftgoing and rightgoing parts respectively.
c     # With the Roe solver we have
c     #    amdq  =  A^- \Delta q    and    apdq  =  A^+ \Delta q
c     # where A is the Roe matrix.  An entropy fix can also be incorporated
c     # into the flux differences.
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routines, this routine is called with ql = qr
c
c
      implicit double precision (a-h,o-z)
c
      dimension wave(meqn, mwaves,1-mbc:maxm+mbc)
      dimension    s(mwaves, 1-mbc:maxm+mbc)
      dimension   ql(meqn, 1-mbc:maxm+mbc)
      dimension   qr(meqn, 1-mbc:maxm+mbc)
      dimension  apdq(meqn, 1-mbc:maxm+mbc)
      dimension  amdq(meqn, 1-mbc:maxm+mbc)
      dimension auxl(16, 1-mbc:maxm+mbc)
      dimension auxr(16, 1-mbc:maxm+mbc)
c
c     local arrays -- common block comroe is passed to rpt2
c     ------------
c      parameter (1002 = 1002)  !# assumes at most 1000x1000 grid with mbc=2
      dimension delta(3)
      logical efix

      common /sw/  g
      common /comroe/ u(-1:1002),v(-1:1002),a(-1:1002),h(-1:1002)
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
c
      data efix /.true./    !# use entropy fix for transonic rarefactions
c
c      if (-1.gt.1-mbc .or. maxm2 .lt. maxm+mbc) then
c	 write(6,*) 'need to increase maxm2 in rpA'
c	 stop
c	 endif
      
      if(ixy.eq.1) then
        dy = dycom
      else
        dy = dxcom
      endif 

c     The aux array has the following elements:
c         1  kappa = ratio of cell area to dxc*dyc
c         2  enx = x-component of normal vector to left edge in tangent plane
c         3  eny = y-component of normal vector to left edge in tangent plane
c         4  enz = z-component of normal vector to left edge in tangent plane
c         5  etx = x-component of tangent vector to left edge in tangent plane
c         6  ety = y-component of tangent vector to left edge in tangent plane
c         7  etz = z-component of tangent vector to left edge in tangent plane
c         8  enx = x-component of normal vector to bottom edge in tangent plane
c         9  eny = y-component of normal vector to bottom edge in tangent plane
c        10  enz = z-component of normal vector to bottom edge in tangent plane
c        11  etx = x-component of tangent vector to bottom edge in tangent plane
c        12  ety = y-component of tangent vector to bottom edge in tangent plane
c        13  etz = z-component of tangent vector to bottom edge in tangent plane
c        14  erx = x-component of unit vector in radial direction at cell ctr
c        15  ery = y-component of unit vector in radial direction at cell ctr
c        16  erz = z-component of unit vector in radial direction at cell ctr
c
c     # offset to index into aux array for enx, eny, etx, ety, gamma 
c     #    depends on whether ixy=1 (left edge) or ixy=2 (bottom edge).
      ioff = 6*(ixy-1) + 1
c
c
c     # find a1 thru a3, the coefficients of the 3 eigenvectors:

      do 20 i = 2-mbc, mx+mbc

         enx =   auxl(ioff+1,i)
         eny =   auxl(ioff+2,i)
         enz =   auxl(ioff+3,i)
         etx =   auxl(ioff+4,i)
         ety =   auxl(ioff+5,i)
         etz =   auxl(ioff+6,i)
         gamma = dsqrt(etx**2 + ety**2 + etz**2)
         etx =   etx / gamma
         ety =   ety / gamma
         etz =   etz / gamma


c        # compute normal and tangential momentum at cell edge:
         hunl = enx*ql(2,i) + eny*ql(3,i) + enz*ql(4,i)
         hunr = enx*qr(2,i-1) + eny*qr(3,i-1) + enz*qr(4,i-1)

         hutl = etx*ql(2,i) + ety*ql(3,i) + etz*ql(4,i)
         hutr = etx*qr(2,i-1) + ety*qr(3,i-1) + etz*qr(4,i-1)

c        # compute the Roe-averaged variables needed in the Roe solver.
c        # These are stored in the common block comroe since they are
c        # later used in routine rpt2 to do the transverse wave splitting.
c
         hl = ql(1,i)
         hr = qr(1,i-1)
         h(i) = (hl+hr)*0.50d0
         hsqr = dsqrt(hr)
         hsql = dsqrt(hl)
         hsq = hsqr + hsql
         u(i) = (hunr/hsqr + hunl/hsql) / hsq
         v(i) = (hutr/hsqr + hutl/hsql) / hsq
         a(i) = dsqrt(g*h(i))

c        # Split the jump in q at each interface into waves
         delta(1) = hl - hr
         delta(2) = hunl - hunr
         delta(3) = hutl - hutr

         a1 = ((u(i)+a(i))*delta(1) - delta(2))*(0.50d0/a(i))
         a2 = -v(i)*delta(1) + delta(3)
         a3 = (-(u(i)-a(i))*delta(1) + delta(2))*(0.50d0/a(i))

c
c        # Compute the waves.
c
         wave(1,1,i) = a1
         wave(2,1,i) = a1*(u(i)-a(i))*enx + a1*v(i)*etx
         wave(3,1,i) = a1*(u(i)-a(i))*eny + a1*v(i)*ety
         wave(4,1,i) = a1*(u(i)-a(i))*enz + a1*v(i)*etz
         s(1,i) = (u(i)-a(i)) * gamma/dy 
c
         wave(1,2,i) = 0.0d0
         wave(2,2,i) = a2*etx
         wave(3,2,i) = a2*ety
         wave(4,2,i) = a2*etz
         s(2,i) = u(i) * gamma/dy 
c
         wave(1,3,i) = a3
         wave(2,3,i) = a3*(u(i)+a(i))*enx + a3*v(i)*etx
         wave(3,3,i) = a3*(u(i)+a(i))*eny + a3*v(i)*ety
         wave(4,3,i) = a3*(u(i)+a(i))*enz + a3*v(i)*etz
         s(3,i) = (u(i)+a(i)) * gamma/dy 
  281    format(2i4,5d12.4)
  283    format(8x,5d12.4)
   20    continue
c
c
c    # compute flux differences amdq and apdq.
c    ---------------------------------------
c
      if (efix) go to 110
c
c     # no entropy fix
c     ----------------
c
c     # amdq = SUM s*wave   over left-going waves
c     # apdq = SUM s*wave   over right-going waves
c
      do 100 m=1,meqn
         do 100 i=2-mbc, mx+mbc
        amdq(m,i) = 0.d0
        apdq(m,i) = 0.d0
        do 90 mw=1,mwaves
           if (s(mw,i) .lt. 0.d0) then
           amdq(m,i) = amdq(m,i) + s(mw,i)*wave(m,mw,i)
         else
           apdq(m,i) = apdq(m,i) + s(mw,i)*wave(m,mw,i)
         endif
   90          continue
  100       continue
c
c     # project momentum components of amdq and apdq onto tangent plane:

      do i=2-mbc,mx+mbc  
            erx = auxr(14,i-1)
            ery = auxr(15,i-1)
            erz = auxr(16,i-1)
            amn = erx*amdq(2,i) + ery*amdq(3,i) + erz*amdq(4,i)
            amdq(2,i) = amdq(2,i) - amn*erx
            amdq(3,i) = amdq(3,i) - amn*ery
            amdq(4,i) = amdq(4,i) - amn*erz

            erx = auxl(14,i)
            ery = auxl(15,i)
            erz = auxl(16,i)
            apn = erx*apdq(2,i) + ery*apdq(3,i) + erz*apdq(4,i)
            apdq(2,i) = apdq(2,i) - apn*erx
            apdq(3,i) = apdq(3,i) - apn*ery
            apdq(4,i) = apdq(4,i) - apn*erz
            enddo
      go to 900
c
c-----------------------------------------------------
c
  110 continue
c
c     # With entropy fix
c     ------------------
c
c    # compute flux differences amdq and apdq.
c    # First compute amdq as sum of s*wave for left going waves.
c    # Incorporate entropy fix by adding a modified fraction of wave
c    # if s should change sign.
c
         do 200 i=2-mbc,mx+mbc
            do m=1, meqn
                amdq(m,i) = 0.d0
                apdq(m,i) = 0.d0
            enddo
            enx =   auxl(ioff+1,i)
            eny =   auxl(ioff+2,i)
            enz =   auxl(ioff+3,i)
            etx =   auxl(ioff+4,i)
            ety =   auxl(ioff+5,i)
            etz =   auxl(ioff+6,i)
            gamma = dsqrt(etx**2 + ety**2 + etz**2)
            etx =   etx / gamma
            ety =   ety / gamma
            etz =   etz / gamma
c           # compute normal and tangential momentum at cell edge:
            hunl = enx*ql(2,i) + eny*ql(3,i) + enz*ql(4,i)
            hunr = enx*qr(2,i-1) + eny*qr(3,i-1) + enz*qr(4,i-1)

c           check 1-wave
            him1 = qr(1,i-1)
            s0 =  (hunr/him1 - dsqrt(g*him1)) * gamma / dy
c           check for fully supersonic case :
            if (s0.gt.0.0d0.and.s(1,i).gt.0.0d0) then
               do 60 m=1,4
                  amdq(m,i)=0.0d0
   60          continue
               goto 200
            endif
c
            h1 = qr(1,i-1)+wave(1,1,i)
            hu1= hunr + enx*wave(2,1,i) + eny*wave(3,1,i) 
     &                                        + enz*wave(4,1,i)
            s1 = (hu1/h1 - dsqrt(g*h1))*gamma/dy  !speed just to right of 1-wave
            if (s0.lt.0.0d0.and.s1.gt.0.0d0) then
c              transonic rarefaction in 1-wave
               sfract = s0*((s1-s(1,i))/(s1-s0))
            else if (s(1,i).lt.0.0d0) then
c              1-wave is leftgoing
               sfract = s(1,i)
            else
c              1-wave is rightgoing
               sfract = 0.0d0
            endif
            do 120 m=1,4
               amdq(m,i) = sfract*wave(m,1,i)
  120       continue
c           check 2-wave
            if (s(2,i).gt.0.0d0) then
c	       #2 and 3 waves are right-going
           go to 200 
           endif

            do 140 m=1,4
               amdq(m,i) = amdq(m,i) + s(2,i)*wave(m,2,i)
  140       continue
c
c           check 3-wave
c
            hi = ql(1,i)
            s03 = (hunl/hi + dsqrt(g*hi)) * gamma/dy
            h3=ql(1,i)-wave(1,3,i)
            hu3=hunl-(enx*wave(2,3,i)+eny*wave(3,3,i)+enz*wave(4,3,i))
            s3=(hu3/h3 + dsqrt(g*h3)) * gamma/dy
            if (s3.lt.0.0d0.and.s03.gt.0.0d0) then
c              transonic rarefaction in 3-wave
               sfract = s3*((s03-s(3,i))/(s03-s3))
            else if (s(3,i).lt.0.0d0) then
c              3-wave is leftgoing
               sfract = s(3,i)
            else
c              3-wave is rightgoing
               goto 200
            endif
            do 160 m=1,4
               amdq(m,i) = amdq(m,i) + sfract*wave(m,3,i)
  160       continue
  200       continue
c
c           compute rightgoing flux differences :
c
            do 220 i = 2-mbc,mx+mbc  
               do 222 m=1,4
                  df = 0.0d0
                  do 210 mw=1,mwaves
                     df = df + s(mw,i)*wave(m,mw,i)
  210             continue
                  apdq(m,i)=df - amdq(m,i)
 222              continue
c
c                 project momentum components onto tangent plane
c
                  erx = auxr(14,i-1)
                  ery = auxr(15,i-1)
                  erz = auxr(16,i-1)
                  amn = erx*amdq(2,i)+ery*amdq(3,i)+erz*amdq(4,i)
                  amdq(2,i) = amdq(2,i) - amn*erx
                  amdq(3,i) = amdq(3,i) - amn*ery
                  amdq(4,i) = amdq(4,i) - amn*erz

                  erx = auxl(14,i)
                  ery = auxl(15,i)
                  erz = auxl(16,i)
                  apn = erx*apdq(2,i)+ery*apdq(3,i)+erz*apdq(4,i)
                  apdq(2,i) = apdq(2,i) - apn*erx
                  apdq(3,i) = apdq(3,i) - apn*ery
                  apdq(4,i) = apdq(4,i) - apn*erz

  220          continue


  900 continue
      return
      end


