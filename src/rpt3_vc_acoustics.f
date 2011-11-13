c
c
c     ==================================================================
      subroutine rpt3(ixyz,icoor,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,maux,imp,asdq,
     &                  bmasdq,bpasdq)
c     ==================================================================
c
c     # Riemann solver in the transverse direction for the acoustics equations
c     # with varying material properties
c
c     # auxN(i,1,:) holds impedance Z
c     # auxN(i,2,:) holds sound speed c
c     #
c     # On input,
c
c     #    ql,qr is the data along some one-dimensional slice, as in rpn3
c     #         This slice is
c     #             in the x-direction if ixyz=1,
c     #             in the y-direction if ixyz=2, or
c     #             in the z-direction if ixyz=3.
c     #    asdq is an array of flux differences (A^*\Dq).
c     #         asdq(i,:) is the flux difference propagating away from
c     #         the interface between cells i-1 and i.
c     #    Note that asdq represents B^*\Dq if ixyz=2 or C^*\Dq if ixyz=3.
c
c     #    ixyz indicates the direction of the original Riemann solve,
c     #         called the x-like direction in the table below:
c
c     #               x-like direction   y-like direction   z-like direction
c     #      ixyz=1:        x                  y                  z
c     #      ixyz=2:        y                  z                  x
c     #      ixyz=3:        z                  x                  y
c
c     #    icoor indicates direction in which the transverse solve should
c     #         be performed.
c     #      icoor=2: split in the y-like direction.
c     #      icoor=3: split in the z-like direction.
c
c     #    For example,
c     #      ixyz=1, icoor=2 means asdq=A^*\Dq, and should be split in y into
c     #                        bmasdq = B^-A^*\Dq,
c     #                        bpasdq = B^+A^*\Dq.
c     #
c     #      ixyz=2, icoor=2 means asdq=B^*\Dq, and should be split in z into
c     #                        bmasdq = C^-B^*\Dq,
c     #                        bpasdq = C^+B^*\Dq.
c
c     #    The parameter imp is generally needed only if aux
c     #    arrays are being used, in order to access the appropriate
c     #    variable coefficients:

c     #    imp = 1 if asdq = A^- \Dq,  the left-going flux difference
c     #          2 if asdq = A^+ \Dq, the right-going flux difference
c
c     #    aux2(:,:,2) is a 1d slice of the aux array along the row
c     #                 where the data ql, qr lie.
c     #    aux1(:,:,2) and aux3(:,:,2) are neighboring rows in the
c     #                 y-like direction
c     #    aux2(:,:,1) and aux2(:,:,3) are neighboring rows in the
c     #                z-like direction
c
c
      implicit real*8(a-h,o-z)
      dimension     ql(meqn,1-mbc:maxm+mbc)
      dimension     qr(meqn,1-mbc:maxm+mbc)
      dimension   asdq(meqn,1-mbc:maxm+mbc)
      dimension bmasdq(meqn,1-mbc:maxm+mbc)
      dimension bpasdq(meqn,1-mbc:maxm+mbc)
      dimension   aux1(maux,1-mbc:maxm+mbc,3)
      dimension   aux2(maux,1-mbc:maxm+mbc,3)
      dimension   aux3(maux,1-mbc:maxm+mbc,3)
c
c
c     # set iuvw = 2,3,4, depending on which component of q represents
c     # velocity in the transverse direction in which splitting is to 
c     # be performed:
      iuvw = ixyz + icoor 
      if (iuvw.gt.4) iuvw = iuvw-3

c
      do 10 i=2-mbc,mx+mbc

c        # The flux difference asdq is split into downward moving part
c        # traveling at speed -c relative to the medium below and
c        # an upward moving part traveling
c        # at speed +c relative to the medium above.
c
c        # Note that the sum of these parts does not give all of asdq
c        # since there is also reflection at the interfaces which decreases
c        # the flux.
c
c        # set impendance and sound speed in each row of cells,
c        # by selecting appropriate values from aux arrays depending on
c        # values of icoor, imp:

c        # imp is used to flag whether wave is going to left or right.
         i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq
c
         if (icoor .eq. 2) then  
c           # transverse direction is y-like direction so
c           # auxN(:,:,2) holds data in appropriate plane and N=(1,2,3)
c           # for row (below,at,above) the slice of q data
            zm = aux1(1,i1,2)
            zz = aux2(1,i1,2)
            zp = aux3(1,i1,2)
            cm = aux1(2,i1,2)
            c  = aux2(2,i1,2)
            cp = aux3(2,i1,2)
         else !! (icoor .eq. 3) 
c           # transverse direction is z-like direction so
c           # aux2(:,:,N) holds data in appropriate plane and N=(1,2,3)
c           # for row (below,at,above) the slice of q data
            zm = aux2(1,i1,1)
            zz = aux2(1,i1,2)
            zp = aux2(1,i1,3)
            cm = aux2(2,i1,1)
            c  = aux2(2,i1,2)
            cp = aux2(2,i1,3)
         endif
c

c        # transmitted part of down-going wave:
         a1 = (-asdq(1,i) + asdq(iuvw,i)*zz) /
     &         (zm + zz)

c        # transmitted part of up-going wave:
         a2 = (asdq(1,i) + asdq(iuvw,i)*zz) /
     &         (zz + zp)
c
c        # The down-going flux difference bmasdq is the product  -c * wave
c
         bmasdq(1,i) = cm * a1*zm
         bmasdq(2,i) = 0.d0
         bmasdq(3,i) = 0.d0
         bmasdq(4,i) = 0.d0
         bmasdq(iuvw,i) = -cm * a1
c
c        # The up-going flux difference bpasdq is the product  c * wave
c
         bpasdq(1,i) = cp * a2*zp
         bpasdq(2,i) = 0.d0
         bpasdq(3,i) = 0.d0
         bpasdq(4,i) = 0.d0
         bpasdq(iuvw,i) = cp * a2
c
   10    continue

c
      return
      end


