! =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! =====================================================
!
!     # Riemann solver in the transverse direction for the elastic equations
!     # with varying material properties 
!
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


    integer, intent(in) :: ixy, imp, maxm, meqn, mwaves, maux, mbc, mx
    double precision, intent(in) :: ql, qr, aux1, aux2, aux3, asdq
    double precision, intent(out) :: bmasdq, bpasdq

    double precision :: wave, s

    dimension     ql(meqn, 1-mbc:maxm+mbc)
    dimension     qr(meqn, 1-mbc:maxm+mbc)
    dimension   asdq(meqn, 1-mbc:maxm+mbc)
    dimension bmasdq(meqn, 1-mbc:maxm+mbc)
    dimension bpasdq(meqn, 1-mbc:maxm+mbc)
    dimension   aux1(maux, 1-mbc:maxm+mbc)
    dimension   aux2(maux, 1-mbc:maxm+mbc)
    dimension   aux3(maux, 1-mbc:maxm+mbc)
    dimension wave( meqn, mwaves)
    dimension    s(mwaves)
    integer :: sig_ii, sig_kk, sig_ij, u_i, u_j, m, i, j, i1
    double precision :: dsig_ii, dsig_ij, du_i, du_j
    double precision :: lamp, mup, bulkp, cpp, csp
    double precision :: lam, mu, bulk, cp, cs
    double precision :: lamm, mum, bulkm, cpm, csm
    double precision :: det, a1, a2, a3, a4


    ! set sig_ii and u_i to correspond to the split direction
    ! set sig_kk to the diagonal element of stress in the orthogonal direction
    ! set sig_ij to the off-diagonal element of stress
    ! set u_j to the orthoginal direction

!
    if (ixy.eq.1) then
        sig_ii = 2
        sig_kk = 1
        sig_ij = 3
        u_i = 6
        u_j = 5
    else
        sig_ii = 1
        sig_kk = 2
        sig_ij = 3
        u_i = 5
        u_j = 6
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
         dsig_ii = asdq(sig_ii,i)
         dsig_ij = asdq(sig_ij,i)
         du_i     = asdq(u_i,i)
         du_j     = asdq(u_j,i)

!        # Material parameters in each row of cells:
         lamm = aux1(2,i1)
         lam  = aux2(2,i1)
         lamp = aux3(2,i1)
         mum  = aux1(3,i1)
         mu   = aux2(3,i1)
         mup  = aux3(3,i1)
         bulkm = lamm + 2.d0*mum
         bulk  = lam  + 2.d0*mu
         bulkp = lamp + 2.d0*mup
     
!        # P-wave and S-wave speeds in each row of cells:
         cpm = aux1(4,i1)
         cp  = aux2(4,i1)
         cpp = aux3(4,i1)
         csm = aux1(5,i1)
         cs  = aux2(5,i1)
         csp = aux3(5,i1)
!

!        # transmitted part of up-going P-wave:
         det = bulk*cpp + bulkp*cp
         if (det .eq. 0.d0) then
            write(6,*) 'det1 = 0 in rpt2'
            stop
            endif
         a1 = (cp*dsig_ii - bulk*du_i) / det

!        # transmitted part of down-going P-wave:
         det = bulkm*cp + bulk*cpm
         if (det .eq. 0.d0) then
            write(6,*) 'det2 = 0 in rpt2'
            stop
            endif
         a2 = (cp*dsig_ii + bulk*du_i) / det
!
!        # transmitted part of up-going S-wave:
         det = mu*csp + mup*cs
         if (det .eq. 0.d0) then
             a3 = 0.d0
         else
             a3 = (cs*dsig_ij - mu*du_j) / det
         endif

!        # transmitted part of down-going S-wave:
         det = mum*cs + mu*csm
         if (det .eq. 0.d0) then
             a4 = 0.d0
         else
             a4 = (cs*dsig_ij + mu*du_j) / det
         endif

        ! Compute the P-waves parts (a1 upward, a2 downward)
        wave(sig_ii,1) = a1 * bulkp
        wave(sig_kk,1) = a1 * lamp
        wave(sig_ij,1) = 0.d0
        wave(4,1) = a1 * lamp
        wave(u_i,1) = -a1 * cpp
        wave(u_j,1) = 0.d0
        s(1) = cpp

        wave(sig_ii,2) = a2 * bulkm
        wave(sig_kk,2) = a2 * lamm
        wave(sig_ij,2) = 0.d0
        wave(4,2) = a2 * lamm
        wave(u_i,2) = a2 * cpm
        wave(u_j,2) = 0.d0
        s(2) = -cpm

        ! Compute the S-waves parts (a3 upward, a4 downward
        wave(sig_ii,3) = 0.d0
        wave(sig_kk,3) = 0.d0
        wave(sig_ij,3) = a3 * mup
        wave(4,3) = 0.d0
        wave(u_i,3) = 0.d0
        wave(u_j,3) = -a3 * csp
        s(3) = csp

        wave(sig_ii,4) = 0.d0
        wave(sig_kk,4) = 0.d0
        wave(sig_ij,4) = a4 * mum
        wave(4,4) = 0.d0
        wave(u_i,4) = 0.d0
        wave(u_j,4) = a4 * csm
        s(4) = -csm

!        # Compute downward and upward flux difference:
!       Remember that s1,s3 > 0 and s2,s4 < 0

        do j=1,meqn
            bmasdq(j,i) = s(2)*wave(j,2) + s(4)*wave(j,4)
            bpasdq(j,i) = s(1)*wave(j,1) + s(3)*wave(j,3)
        end do

!
      enddo
!
      return
      end
