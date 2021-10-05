module mod_phi
use mod_var
use mod_global
use mod_bou
use mod_conserve
implicit none
private
public updatePfm
contains
!
subroutine updatePfm(dPFM)
implicit none
integer i,j,k
real :: factor1,factor2
real ::dPFMnew(0:i1,0:j1,0:k1)
real, intent(inout) :: dPFM(0:,0:,0:)
real:: Phi_wall_n(0:i1,0:j1),phi_wall_n_top(0:i1,0:j1)
real:: dPFM_bound_new,Phi_wall_n1,Phi_wall_n1_top
real, dimension(0:i1,0:j1) ::dPFM_bound_n,dPFM_bound_n_top
real:: RHS_all
real::advective
real::advectivex,advectivey, advectivez
real:: Cahn_Hilliard_RHS,d2potdx2,d2potdy2,d2potdz2,source
!
Phi_wall_n=0.
phi_wall_n_top=0.
!
factor1 = dt*( 1.5)
factor2 = dt*(-0.5)
!factor1 = dt*( 2.0)
!factor2 = dt*(-1.0)  !according to what Ziyang said
!
!!at time n:
do j=1,jmax
  do i=1,imax
    phi_wall_n(i,j)=0.5*(PFM_phi(i,j,1)+PFM_phi(i,j,0))
    phi_wall_n_top(i,j)=0.5*(PFM_phi(i,j,kmax)+PFM_phi(i,j,k1))

  enddo
enddo
call boundPFM(PFM_phi,dPFM_bound_n,dPFM_bound_n_top,unew,vnew,wnew)


!call m_phi

!update phi to time n+1
do k=1,kmax
  do j=1,jmax
    do i=1,imax
       source = PFM_phi(i,j,k)*( (unew(i,j,k)-unew(i-1,j,k))/dx + (vnew(i,j,k)-vnew(i,j-1,k))/dy + (wnew(i,j,k)-wnew(i,j,k-1))/dz )
       !dPFMnew(i,j,k) = advective + Cahn_Hilliard_RHS +source
       dPFMnew(i,j,k) = -(m_x_phi(i,j,k)-m_x_phi(i-1,j,k))/dx - &
                         (m_y_phi(i,j,k)-m_y_phi(i,j-1,k))/dy - &
                         (m_z_phi(i,j,k)-m_z_phi(i,j,k-1))/dz + source

    enddo
  enddo
enddo


do k=1,kmax
  do j=1,jmax
    do i=1,imax
      PFM_phi(i,j,k) =PFM_phi_old(i,j,k) + dt*dPFMnew(i,j,k)
      !PFM_phi(i,j,k) =PFM_phi(i,j,k) + dt*dPFMnew(i,j,k)
      !PFM_phi(i,j,k) =( 2*PFM_phi_old(i,j,k)-0.5*PFM_phi_old2(i,j,k)+factor1*dPFMnew(i,j,k) + factor2*dPFM(i,j,k) )/1.5
      !dPFM(i,j,k)=dPFMnew(i,j,k)
      !PFM_phi(i,j,k) =( PFM_phi_old(i,j,k)*rholold(i,j,k) + dt*dPFMnew(i,j,k) )/rhol(i,j,k)
    enddo
  enddo
enddo




!update the boundary to n+1
do j=1,jmax
  do i=1,imax
    Phi_wall_n1=phi_wall_n(i,j)+factor1*dPFM_bound_n(i,j)+factor2*dPFM_bound_old(i,j)
    Phi_wall_n1_top=phi_wall_n_top(i,j)+factor1*dPFM_bound_n_top(i,j)+factor2*dPFM_bound_old_top(i,j)
    if (abs(PFM_mu_f).gt.1e-8) then
       PFM_phi(i,j,0) = 2*Phi_wall_n1-PFM_phi(i,j,1)
       dPFM_bound_old(i,j)=dPFM_bound_n(i,j)
    endif

  enddo
enddo
!



end subroutine updatePfm

end module mod_phi
