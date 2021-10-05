module mod_c2
use mod_cons
use mod_var
use mod_global
use mod_bou
use mod_conserve
implicit none
private
public updatec_t, update_delta
contains
!

subroutine updatec_t
implicit none
integer i,j,k
do k=1,kmax
   do j=1,jmax
      do i=1,imax  !here use new phi value
         ct_change(i,j,k)=(Phi_c(i,j,k)-Phi_c_old(i,j,k))/dt+&
                          (m_x_phi_c(i,j,k)-m_x_phi_c(i-1,j,k))/dx+&
                          (m_y_phi_c(i,j,k)-m_y_phi_c(i,j-1,k))/dy+&
                          (m_z_phi_c(i,j,k)-m_z_phi_c(i,j,k-1))/dz
      enddo
   enddo
enddo
end subroutine updatec_t





subroutine update_delta
implicit none
integer i,j,k
real:: lap_c,c_star,temp,source,pterm


!this is used for momentum equation, now c has been updated
! and here use new phi value
do k=1,kmax
   do j=1,jmax
      do i=1,imax
         lap_c= -lamda_c  * ( ( (PFM_c(i+1,j,k)-PFM_c(i,j,k))/dx*(PFM_phi(i+1,j,k)+PFM_phi(i,j,k))/2.0 - &
                               (PFM_c(i,j,k)-PFM_c(i-1,j,k))/dx*(PFM_phi(i,j,k)+PFM_phi(i-1,j,k))/2.0  )/dx +&
                             ( (PFM_c(i,j+1,k)-PFM_c(i,j,k))/dy*(PFM_phi(i,j+1,k)+PFM_phi(i,j,k))/2.0 - &
                               (PFM_c(i,j,k)-PFM_c(i,j-1,k))/dy*(PFM_phi(i,j,k)+PFM_phi(i,j-1,k))/2.0  )/dy +& 
                             ( (PFM_c(i,j,k+1)-PFM_c(i,j,k))/dz*(PFM_phi(i,j,k+1)+PFM_phi(i,j,k))/2.0 - &
                               (PFM_c(i,j,k)-PFM_c(i,j,k-1))/dz*(PFM_phi(i,j,k)+PFM_phi(i,j,k-1))/2.0  )/dz & 
                           )

         c_star=2.0*PFM_c(i,j,k)*(PFM_c(i,j,k)-1)*(2*PFM_c(i,j,k)-1)*lamda_c/PFM_l_c**2*PFM_phi(i,j,k)

         if (Tnew(i,j,k)>t_melt .and. PFM_c(i,j,k)==0.) then !here the temperature is in degree, not Kelvin
            pterm=1
         else if (Tnew(i,j,k)<t_melt .and. PFM_c(i,j,k)==1 ) then
            pterm=1
         else
           pterm=30*(PFM_c(i,j,k)**2*(PFM_c(i,j,k)-1)**2)
         endif
         temp=rho_2*latent/t_melt*PFM_phi(i,j,k)*(t_melt-Tnew(i,j,k)) * pterm


         delta_u(i,j,k)=mobility_c*(lap_c+c_star+temp) !this is used in momentum equation, not here
           
      enddo
   enddo
enddo

end subroutine update_delta



end module mod_c2
