! =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! =====================================================
!
!     # Riemann solver in the transverse direction for the elastic equations
!     # with varying material properties 
!
!
!     # Contents of ql and qr:
!     # 
!     # q(1,:) = sigma^{11} if ixy=1   or   sigma^{22} if ixy=2
!     # q(2,:) = sigma^{22} if ixy=1   or   sigma^{11} if ixy=2
!     # q(3,:) = sigma^{12} = sigma^{21}
!     # q(4,:) = u          if ixy=1   or   v          if ixy=2
!     # q(5,:) = v          if ixy=1   or   u          if ixy=2
!     # 
!     # auxN holds corresponding slices of the aux array:
!     #  N = 1 for row below
!     #      2 for this row
!     #      3 for row above
!     # 
!     #  auxN(1,i) = rho 
!     #  auxN(2,i) = lambda 
!     #  auxN(3,i) = mu
!     #  auxN(4,i) = cp 
!     #  auxN(5,i) = cs
!
!
!
!     # Split asdq into down-going flux bmasdq and up-going flux bpasdq.
!
!     # imp=1  means  asdq=amdq,    imp=2 means asdq=apdq
!
      implicit none
      !Input
      integer, intent(in)  :: ixy,imp,maxm,meqn,mwaves,maux,mbc,mx
      double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: ql,qr
      double precision, dimension(maux,1-mbc:maxm+mbc), intent(in) :: aux1,aux2,aux3
      double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: asdq

      !Output
      double precision, dimension(meqn,1-mbc:maxm+mbc), intent(out) :: bmasdq,bpasdq

      !Local
      integer :: i,i1,ksig11,ksig22,ku,kv
      double precision :: a1,a2,a3,a4
      double precision :: alamm,alam,alamp,amum,amu,amup,bulkm,bulk,bulkp
      double precision :: cpm,cp,cpp,csm,cs,csp
      double precision :: dsig11,dsig22,dsig12,du,dv
      double precision :: det

!
!
!
!     # set ku to point to  the component of the system that corresponds
!     # to velocity in the direction of this slice, kv to the orthogonal
!     # velocity.  Similarly ksig11 and ksig22 point to normal stresses.
!     # 3rd component is always shear stress sig12.
!
!
      if (ixy.eq.1) then
         ksig11 = 1
         ksig22 = 2
         ku = 4
         kv = 5
      else
         ksig11 = 2
         ksig22 = 1
         ku = 5
         kv = 4
      endif
!
!
      do i = 2-mbc, mx+mbc
!
!        # imp is used to flag whether wave is going to left or right,
!        # since material properties are different on the two sides
!
         if (imp.eq.1) then 
!            # asdq = amdq, moving to left
             i1 = i-1
         else
!            # asdq = apdq, moving to right
             i1 = i
         endif
!
!        # The flux difference asdq is split into downward moving parts
!        # traveling at speeds -cp and -cs relative to the medium below and
!        # upward moving parts traveling
!        # at speeds +cp and +cs relative to the medium above.
!
!        # Note that the sum of these parts does not give all of asdq
!        # since there is also reflection at the interfaces which decreases
!        # the flux.
!
!        # jumps in asdq:
         dsig11 = asdq(ksig11,i)
         dsig22 = asdq(ksig22,i)
         dsig12 = asdq(3,i)
         du     = asdq(ku,i)
         dv     = asdq(kv,i)
!
!
!        # Material parameters in each row of cells:
         alamm = aux1(2,i1)
         alam  = aux2(2,i1)
         alamp = aux3(2,i1)
         amum  = aux1(3,i1)
         amu   = aux2(3,i1)
         amup  = aux3(3,i1)
         bulkm = alamm + 2.d0*amum
         bulk  = alam  + 2.d0*amu 
         bulkp = alamp + 2.d0*amup
     
!        # P-wave and S-wave speeds in each row of cells:
         cpm = aux1(4,i1)
         cp  = aux2(4,i1)
         cpp = aux3(4,i1)
         csm = aux1(5,i1)
         cs  = aux2(5,i1)
         csp = aux3(5,i1)
!

!        # transmitted part of down-going P-wave:
         det = bulkm*cp + bulk*cpm
         if (det .eq. 0.d0) then
            write(6,*) 'det1 = 0 in rpt2'
            stop
            endif
         a1 = (cp*dsig22 + bulk*dv) / det

!        # transmitted part of up-going P-wave:
         det = bulk*cpp + bulkp*cp
         if (det .eq. 0.d0) then
            write(6,*) 'det2 = 0 in rpt2'
            stop
            endif
         a2 = (cp*dsig22 - bulk*dv) / det
!
!        # transmitted part of down-going S-wave:
         det = -(amum*cs + amu*csm)
         if (det .eq. 0.d0) then
             a3 = 0.d0
         else
             a3 = (cs*dsig12 + amu*du) / det
         endif

!        # transmitted part of up-going S-wave:
         det = -(amu*csp + amup*cs)
         if (det .eq. 0.d0) then
             a4 = 0.d0
         else
             a4 = (cs*dsig12 - amu*du) / det
         endif
!
!        # The down-going flux difference bmasdq is the product  -c * wave
!        # summed over down-going P-wave and S-wave:
!
         bmasdq(ksig11,i) = -cpm*a1*alamm
         bmasdq(ksig22,i) = -cpm*a1*bulkm
         bmasdq(3,i) =      -csm*a3*(-amum)
         bmasdq(ku,i) =     -csm*a3*(-csm)
         bmasdq(kv,i) =     -cpm*a1*(cpm)
!
!        # The up-going flux difference bpasdq is the product  c * wave
!        # summed over up-going P-wave and S-wave:
!
         bpasdq(ksig11,i) =  cpp*a2*alamp
         bpasdq(ksig22,i) =  cpp*a2*bulkp
         bpasdq(3,i) =       csp*a4*(-amup)
         bpasdq(ku,i) =      csp*a4*(csp)
         bpasdq(kv,i) =      cpp*a2*(-cpp)
!
      enddo
!
      return
      end
