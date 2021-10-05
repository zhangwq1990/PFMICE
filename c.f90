module mod_c
use mod_cons
use mod_var
use mod_global
use mod_bou
use mod_conserve
implicit none
private
public updatephi_c, updatec
contains
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! c range now is 0 to 1, 1 means full water and 0 means full ice
!subroutine updatec(dPFM_c_old)
subroutine updatephi_c
implicit none
integer i,j,k
real ::dPFM_c(0:i1,0:j1,0:k1)
!real, intent(inout) :: dPFM_c_old(0:,0:,0:)
real:: advectivex,advectivey,advectivez,advective,lap_c,c_star,temp,source,pterm



do k=1,kmax
   do j=1,jmax
      do i=1,imax

         advective = (m_x_phi_c(i,j,k)-m_x_phi_c(i-1,j,k))/dx + &
                     (m_y_phi_c(i,j,k)-m_y_phi_c(i,j-1,k))/dy + &
                     (m_z_phi_c(i,j,k)-m_z_phi_c(i,j,k-1))/dz 



         !use phi_old, because phi has already been updated
         lap_c= -lamda_c  *( ( (PFM_c(i+1,j,k)-PFM_c(i,j,k))/dx*(PFM_phi_old(i+1,j,k)+PFM_phi_old(i,j,k))/2.0 - &
                               (PFM_c(i,j,k)-PFM_c(i-1,j,k))/dx*(PFM_phi_old(i,j,k)+PFM_phi_old(i-1,j,k))/2.0  )/dx +&
                             ( (PFM_c(i,j+1,k)-PFM_c(i,j,k))/dy*(PFM_phi_old(i,j+1,k)+PFM_phi_old(i,j,k))/2.0 - &
                               (PFM_c(i,j,k)-PFM_c(i,j-1,k))/dy*(PFM_phi_old(i,j,k)+PFM_phi_old(i,j-1,k))/2.0  )/dy +& 
                             ( (PFM_c(i,j,k+1)-PFM_c(i,j,k))/dz*(PFM_phi_old(i,j,k+1)+PFM_phi_old(i,j,k))/2.0 - &
                               (PFM_c(i,j,k)-PFM_c(i,j,k-1))/dz*(PFM_phi_old(i,j,k)+PFM_phi_old(i,j,k-1))/2.0  )/dz & 
                           )

         c_star=2.0*PFM_c(i,j,k)*(PFM_c(i,j,k)-1)*(2*PFM_c(i,j,k)-1)*lamda_c/PFM_l_c**2*PFM_phi_old(i,j,k)

         if (Tnew(i,j,k)>t_melt .and. PFM_c(i,j,k)==0.) then !here the temperature is in degree, not Kelvin
            pterm=1.0
         else if (Tnew(i,j,k)<t_melt .and. PFM_c(i,j,k)==1 ) then
            pterm=1.0
         else
           pterm=30*(PFM_c(i,j,k)**2*(PFM_c(i,j,k)-1)**2)
         endif
         temp=rho_2*latent/t_melt*PFM_phi_old(i,j,k)*(t_melt-Tnew(i,j,k)) * pterm


         source = Phi_c_old(i,j,k)*( (unew(i,j,k)-unew(i-1,j,k))/dx +&
                                     (vnew(i,j,k)-vnew(i,j-1,k))/dy +&
                                     (wnew(i,j,k)-wnew(i,j,k-1))/dz )

         dPFM_c(i,j,k)=-advective-mobility_c*(lap_c+c_star+temp)+source
      enddo
   enddo
enddo




do k=1,kmax
   do j=1,jmax
      do i=1,imax
         Phi_c(i,j,k)=Phi_c_old(i,j,k)+dt*dPFM_c(i,j,k)
      enddo
   enddo
enddo


end subroutine updatephi_c



subroutine updatec
implicit none
integer i,j,k


do k=1,kmax
   do j=1,jmax
      do i=1,imax  !here use new phi value
         if ( (abs(PFM_phi(i,j,k)-Phi_c(i,j,k))<1e-15) .or. (PFM_phi(i,j,k)<1e-15) ) then
            PFM_c(i,j,k)=1.0
         else
         PFM_c(i,j,k)=Phi_c(i,j,k)/(PFM_phi(i,j,k)+1e-30)   !this very small number can be selected
         PFM_c(i,j,k)=MIN(MAX(PFM_c(i,j,k), 0.0),1.0)
         endif
      enddo
   enddo
enddo

end subroutine updatec


end module mod_c
