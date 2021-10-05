module mod_temp
use mod_var
use mod_global
use mod_bou
use mod_conserve
implicit none
private
public temp
contains
!
subroutine temp
implicit none
integer i,j,k,info
!real, intent(inout),dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: Tnew
real ,dimension(1:imax,1:jmax,1:kmax) :: RHS
real:: advectivex,advectivey, advectivez,advective
!real:: pressure
real:: dissipate1, dissipate2, dissipate3,dissipate
real:: conductionx, conductiony, conductionz, conduction,gravity
real:: melt




!call m_rho_cp


do k=1,kmax
   do j=1,jmax
      do i=1,imax

         !&&&&&&&&&&&&&&&&&&&&&&&&&&&
         ! dissipation term
         !&&&&&&&&&&&&&&&&&&&&&&&&&&&
         dissipate1= -2.0/3.0*vislold(i,j,k)*( (unew(i,j,k)-unew(i-1,j,k))/dx + &
                                             (vnew(i,j,k)-vnew(i,j-1,k))/dy + &
                                             (wnew(i,j,k)-wnew(i,j,k-1))/dz )**2
         dissipate2= vislold(i,j,k)*2.0*&
                     (((unew(i,j,k)-unew(i-1,j,k))/dx)**2+&
                      ((vnew(i,j,k)-vnew(i,j-1,k))/dy)**2+&
                      ((wnew(i,j,k)-wnew(i,j,k-1))/dz)**2)
         dissipate3= vislold(i,j,k)*&
                     ( ( (unew(i,j+1,k)-unew(i,j-1,k))/dy/2.0 + (unew(i-1,j+1,k)-unew(i-1,j-1,k))/dy/2.0 )/2.0 + &
                       ( (vnew(i+1,j,k)-vnew(i-1,j,k))/dx/2.0 + (vnew(i+1,j-1,k)-vnew(i-1,j-1,k))/dx/2.0 )/2.0 )**2
         dissipate3= dissipate3 + vislold(i,j,k)*&
                     ( ( (wnew(i+1,j,k)-wnew(i-1,j,k))/dx/2.0 + (wnew(i+1,j,k-1)-wnew(i-1,j,k-1))/dx/2.0 )/2.0 + &
                       ( (unew(i,j,k+1)-unew(i,j,k-1))/dz/2.0 + (unew(i-1,j,k+1)-unew(i-1,j,k-1))/dz/2.0 )/2.0 )**2
         dissipate3= dissipate3 + vislold(i,j,k)*&
                     ( ( (vnew(i,j,k+1)-vnew(i,j,k-1))/dz/2.0 + (vnew(i,j-1,k+1)-vnew(i,j-1,k-1))/dz/2.0 )/2.0 + &
                       ( (wnew(i,j+1,k)-wnew(i,j-1,k))/dy/2.0 + (wnew(i,j+1,k-1)-wnew(i,j-1,k-1))/dy/2.0 )/2.0 )**2
         dissipate = dissipate1 + dissipate2 + dissipate3
         !&&&&&&&&&&&&&&&&&&&&&&&&&&&
         ! heat conduction
         !&&&&&&&&&&&&&&&&&&&&&&&&&&&
         conductionx= ( (kkold(i+1,j,k)+kkold(i,j,k))/2.0*(Tnew(i+1,j,k)-Tnew(i,j,k))/dx - &
                        (kkold(i,j,k)+kkold(i-1,j,k))/2.0*(Tnew(i,j,k)-Tnew(i-1,j,k))/dx )/dx
         conductiony= ( (kkold(i,j+1,k)+kkold(i,j,k))/2.0*(Tnew(i,j+1,k)-Tnew(i,j,k))/dy - &
                        (kkold(i,j,k)+kkold(i,j-1,k))/2.0*(Tnew(i,j,k)-Tnew(i,j-1,k))/dy )/dy
         conductionz= ( (kkold(i,j,k+1)+kkold(i,j,k))/2.0*(Tnew(i,j,k+1)-Tnew(i,j,k))/dz - &
                        (kkold(i,j,k)+kkold(i,j,k-1))/2.0*(Tnew(i,j,k)-Tnew(i,j,k-1))/dz )/dz
         conduction=conductionx+conductiony+conductionz
         !&&&&&&&&&&&&&&&&&&&&&&&&&&&
         ! heat convection
         !&&&&&&&&&&&&&&&&&&&&&&&&&&&

         advectivex= (m_x_rho_cp(i,j,k)*(Tnew(i+1,j,k)+Tnew(i,j,k))/2.0-m_x_rho_cp(i-1,j,k)*(Tnew(i,j,k)+Tnew(i-1,j,k))/2.0 )/dx
         advectivey= (m_y_rho_cp(i,j,k)*(Tnew(i,j+1,k)+Tnew(i,j,k))/2.0-m_y_rho_cp(i,j-1,k)*(Tnew(i,j,k)+Tnew(i,j-1,k))/2.0 )/dy
         advectivez= (m_z_rho_cp(i,j,k)*(Tnew(i,j,k+1)+Tnew(i,j,k))/2.0-m_z_rho_cp(i,j,k-1)*(Tnew(i,j,k)+Tnew(i,j,k-1))/2.0 )/dz
         advective = advectivex+advectivey+advectivez

         !melt= ( (PFM_phi(i,j,k)*PFM_c(i,j,k) - PFM_phi_old(i,j,k)*PFM_c_old(i,j,k))/dt+ &
         !        (m_x_phi_c(i,j,k)-m_x_phi_c(i-1,j,k))/dx+(m_y_phi_c(i,j,k)-m_y_phi_c(i,j-1,k))/dy+&
         !        (m_z_phi_c(i,j,k)-m_z_phi_c(i,j,k-1))/dz )*&
         !      rho_3*latent
         !melt= ( (Phi_c(i,j,k) - Phi_c_old(i,j,k))/dt+ &
         !        (m_x_phi_c(i,j,k)-m_x_phi_c(i-1,j,k))/dx+(m_y_phi_c(i,j,k)-m_y_phi_c(i,j-1,k))/dy+&
         !        (m_z_phi_c(i,j,k)-m_z_phi_c(i,j,k-1))/dz )*&
         !      rho_3*latent
         melt= ct_change(i,j,k)*rho_2*latent
         !gravity = rhol(i,j,k)*(g_x+g_y+g_z) 
         gravity = 0

         RHS(i,j,k) = gravity + conduction + dissipate - advective - melt


      enddo
   enddo
enddo




do k=1,kmax
   do j=1,jmax
      do i=1,imax
         Tnew(i,j,k)=( Told(i,j,k)*rholcp_old(i,j,k) + dt*RHS(i,j,k) )/(rholcp(i,j,k))
      enddo
   enddo
enddo



return
end subroutine temp
!
end module mod_temp
