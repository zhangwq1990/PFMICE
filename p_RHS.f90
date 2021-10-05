module mod_p_RHS
use mod_cons
use mod_var
implicit none
private
public p_rhs
contains
subroutine p_rhs(p)
real, intent(out), dimension(0:,0:,0:) :: p
real :: dti,dtidxi,dtidyi,dtidzi,rho0
real :: as,cd
integer :: i,j,k,im,jm,km
!
rho0=MIN(rho_1,rho_2)
dti = 1./dt
dtidxi = dti*dxi
dtidyi = dti*dyi
dtidzi = dti*dzi





do k=1,kmax
  do j=1,jmax
    do i=1,imax
    
      p(i,j,k) = ( &
                  ( wstar(i,j,k)-wstar(i,j,k-1))*dtidzi+ &
                  ( vstar(i,j,k)-vstar(i,j-1,k))*dtidyi+ &
                  ( ustar(i,j,k)-ustar(i-1,j,k))*dtidxi &
                  -delta_u(i,j,k)*(rho_2-rho_3)/rhol(i,j,k)/dt ) !delta_u is not right, not this time step!



      p(i,j,k) = p(i,j,k)*rho0 + (&
               ((1.0-rho0/(rhol(i,j,k)+rhol(i+1,j,k))*2.0)*(phat(i+1,j,k)-phat(i,j,k))*dxi- &
                (1.0-rho0/(rhol(i,j,k)+rhol(i-1,j,k))*2.0)*(phat(i,j,k)-phat(i-1,j,k))*dxi)*dxi+ &
               ((1.0-rho0/(rhol(i,j,k)+rhol(i,j+1,k))*2.0)*(phat(i,j+1,k)-phat(i,j,k))*dyi- &
                (1.0-rho0/(rhol(i,j,k)+rhol(i,j-1,k))*2.0)*(phat(i,j,k)-phat(i,j-1,k))*dyi)*dyi+ &
               ((1.0-rho0/(rhol(i,j,k)+rhol(i,j,k+1))*2.0)*(phat(i,j,k+1)-phat(i,j,k))*dzi- &
                (1.0-rho0/(rhol(i,j,k)+rhol(i,j,k-1))*2.0)*(phat(i,j,k)-phat(i,j,k-1))*dzi)*dzi) 

    enddo
  enddo
enddo






!
return
end subroutine p_rhs
!
end module mod_p_rhs
