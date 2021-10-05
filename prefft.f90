module mod_prefft
use mod_var
use mod_cons
use decomp_2d
implicit none
private
public init_tri,initsolver
contains
!
subroutine init_tri
integer :: k
!
! Generate tridiagonal matrix
! l is k-1, m is k, r is k+1
do k=1,kmax
  tri1_l(k) = dzi*dzi
  tri2_l(k) = dzi*dzi
  tri1_r(k) = dzi*dzi
  tri2_r(k) = dzi*dzi
  tri1_m(k) = -(tri1_l(k) + tri1_r(k))
  tri2_m(k) = -(tri2_l(k) + tri2_r(k))
enddo

!
! Neumann boundary condition, wall at both top and bottom
! a is k-1, b is k, c is k+1
tri1_m(1)    = tri1_m(1)    + tri1_l(1)
tri1_m(kmax) = tri1_m(kmax) + tri1_r(kmax)
tri1_l(1) = 0.
tri1_r(kmax) = 0.

! outlet at top and wall at bottom

tri2_m(1)    = tri2_m(1) + tri2_l(1)
tri2_m(kmax) = tri2_m(kmax)
tri2_l(1)    = 0.
tri2_r(kmax) = 0.
!end select


end subroutine init_tri

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

subroutine initsolver
integer :: i,j,iv,jv
real :: xrt(itot), yrt(jtot)

! set lookup tables.
!
call vrffti(itot,wi)
call vrffti(jtot,wj)
!
! generate eigenvalues ( xrt and yrt ).
!
!
! x direction
!
do i=3,itot,2
  xrt(i-1) = -4.*dxi*dxi*(sin(float((i-1))*pi/(2.*itot)))**2
  xrt(i) = xrt(i-1)
enddo
xrt(1   ) = 0.
xrt(itot) = -4.*dxi*dxi
!
! y direction
!
do j=3,jtot,2
  yrt(j-1) = -4.*dyi*dyi*(sin(float((j-1))*pi/(2.*jtot)))**2
  yrt(j  ) = yrt(j-1)
enddo
yrt(1   ) = 0.
yrt(jtot) = -4.*dyi*dyi
!
do j=1,jmax
  jv = j + zstart(2) - 1
  do i=1,imax
    iv = i + zstart(1) - 1
    xyrt(i,j) = xrt(iv)+yrt(jv)
  enddo
enddo
!
return
end subroutine initsolver
!
end module mod_prefft
