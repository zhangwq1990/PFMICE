module mod_solver
use mod_cons
use mod_var
implicit none
private
public solver2d
contains
subroutine solver2d(pz,tri0_l,tri0_m,tri0_r)
use decomp_2d
implicit none
real, intent(inout), dimension(1:,1:,1:) :: pz
real, intent(in), dimension(1:kmax) :: tri0_l,tri0_m,tri0_r
real, dimension(itot/dims(1),jtot,ktot/dims(2)) :: py
real, dimension(itot,jtot/dims(1),ktot/dims(2)) :: px
real :: bb
real :: z,d(imax,jmax,kmax)
real :: di(itot),dj(jtot)
integer :: i,j,k
!
call transpose_z_to_y(pz,py)
call transpose_y_to_x(py,px)
!
do k=1,xsize(3)
  do j=1,xsize(2)
    call vrfftf(1,itot,px(1:itot,j,k),di,1,wi)
  enddo
enddo
!
call transpose_x_to_y(px,py)
!
do k=1,ysize(3)
  do i=1,ysize(1)
    call vrfftf(1,jtot,py(i,1:jtot,k),dj,1,wj)
  enddo
enddo
!
call transpose_y_to_z(py,pz)
!
do j=1,jmax
  do i=1,imax
    z        = 1./(tri0_m(1)+xyrt(i,j))
    d(i,j,1) = tri0_r(1)*z
    pz(i,j,1) = pz(i,j,1)*z
  enddo
enddo
do k=2,kmax-1
   do j=1,jmax
     do i=1,imax
       bb       = tri0_m(k)+xyrt(i,j)
       z        = 1./(bb-tri0_l(k)*d(i,j,k-1))
       d(i,j,k) = tri0_r(k)*z
       pz(i,j,k) = (pz(i,j,k)-tri0_l(k)*pz(i,j,k-1))*z
     enddo
  enddo
enddo
do j=1,jmax
  do i=1,imax
    bb       = tri0_m(kmax)+xyrt(i,j)
    z        = bb-tri0_l(kmax)*d(i,j,kmax-1)
    if(z.ne.0.) then
      pz(i,j,kmax) = (pz(i,j,kmax)-tri0_l(kmax)*pz(i,j,kmax-1))/z
    else
      pz(i,j,kmax) =0.
    endif
  enddo
enddo

do k=kmax-1,1,-1
  do j=1,jmax
    do i=1,imax
      pz(i,j,k) = pz(i,j,k)-d(i,j,k)*pz(i,j,k+1)
    enddo
  enddo
enddo

call transpose_z_to_y(pz,py)
do k=1,ysize(3)
  do i=1,ysize(1)
    call vrfftb(1,jtot,py(i,1:jtot,k),dj,1,wj)
  enddo
enddo
!
call transpose_y_to_x(py,px)
!
do k=1,xsize(3)
  do j=1,xsize(2)
    call vrfftb(1,itot,px(1:itot,j,k),di,1,wi)
  enddo
enddo
!
call transpose_x_to_y(px,py)
call transpose_y_to_z(py,pz)
!
return
end subroutine solver2d
!
end module mod_solver
