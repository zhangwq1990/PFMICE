module mod_correc
use mod_cons
use mod_var
implicit none
private
public correc
contains
!
! Corrects the velocity so that it is divergence free
!
subroutine correc(p)
real, intent(in), dimension(0:,0:,0:) :: p
real :: factori,factorj,factork,rho0,corrx,corry,corrz
integer :: i,j,k
!
rho0=MIN(rho_1,rho_2)
factori = dt*dxi
factorj = dt*dyi
factork = dt*dzi
!

do k=1,kmax
  do j=1,jmax
    do i=1,imax
      corrx = - factori * ( p(i+1,j,k)-p(i,j,k) )
      corry = - factorj * ( p(i,j+1,k)-p(i,j,k) )
      corrz = - factork * ( p(i,j,k+1)-p(i,j,k) )

        unew(i,j,k) = ustar(i,j,k) + corrx/rho0 - &
  factori * (1.0/(rhol(i,j,k)+rhol(i+1,j,k))*2.0-1.0/rho0) * ( phat(i+1,j,k)-phat(i,j,k) )

        vnew(i,j,k) = vstar(i,j,k) + corry/rho0 - &
  factorj * (1.0/(rhol(i,j,k)+rhol(i,j+1,k))*2.0-1.0/rho0) * ( phat(i,j+1,k)-phat(i,j,k) )


        wnew(i,j,k) = wstar(i,j,k) + corrz/rho0 - &
  factork * (1.0/(rhol(i,j,k)+rhol(i,j,k+1))*2.0-1.0/rho0) * ( phat(i,j,k+1)-phat(i,j,k) )

    enddo
  enddo
enddo
!
return
end subroutine correc
!
end module mod_correc
