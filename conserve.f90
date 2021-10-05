module mod_conserve
use mod_var
use mod_cons
!use mod_bou
implicit none
private
public m_phi,m_imphi,m_phi_c, m_rho_cp, m_rho
contains
!
subroutine m_phi
implicit none
integer i,j,k

do k=0,k1   !I am not sure if I need +1
   do j=0,j1
      do i=0,i1
         !x direction,defined at surface
         m_x_phi(i,j,k)=unew(i,j,k)*(PFM_phi(i+1,j,k)+PFM_phi(i,j,k))/2.0-mobility*&
                      (chem_pot(i+1,j,k)-chem_pot(i,j,k))/dx

         !y direction, defined at surface
         m_y_phi(i,j,k)=vnew(i,j,k)*(PFM_phi(i,j+1,k)+PFM_phi(i,j,k))/2.0-mobility*&
                      (chem_pot(i,j+1,k)-chem_pot(i,j,k))/dy


         !z direction, defined at surface
         m_z_phi(i,j,k)=wnew(i,j,k)*(PFM_phi(i,j,k+1)+PFM_phi(i,j,k))/2.0-mobility*&
                      (chem_pot(i,j,k+1)-chem_pot(i,j,k))/dz
      enddo
   enddo
enddo

return
end subroutine m_phi
!



subroutine m_imphi
implicit none
integer i,j,k

do k=0,k1   !I am not sure if I need +1
   do j=0,j1
      do i=0,i1
         !x direction,defined at surface
         m_x_phi(i,j,k)=unew(i,j,k)*(PFM_phi_old(i+1,j,k)+PFM_phi_old(i,j,k))/2.0-mobility*&
                      (chem_pot(i+1,j,k)-chem_pot(i,j,k))/dx

         !y direction, defined at surface
         m_y_phi(i,j,k)=vnew(i,j,k)*(PFM_phi_old(i,j+1,k)+PFM_phi_old(i,j,k))/2.0-mobility*&
                      (chem_pot(i,j+1,k)-chem_pot(i,j,k))/dy


         !z direction, defined at surface
         m_z_phi(i,j,k)=wnew(i,j,k)*(PFM_phi_old(i,j,k+1)+PFM_phi_old(i,j,k))/2.0-mobility*&
                      (chem_pot(i,j,k+1)-chem_pot(i,j,k))/dz 
      enddo
   enddo
enddo

return
end subroutine m_imphi




subroutine m_phi_c
implicit none
integer i,j,k

do k=0,k1   !I am not sure if I need +1
   do j=0,j1
      do i=0,i1
         !x direction,defined at surface
         m_x_phi_c(i,j,k)=m_x_phi(i,j,k)*(PFM_c(i+1,j,k)+PFM_c(i,j,k))/2.0

         !y direction, defined at surface
         m_y_phi_c(i,j,k)=m_y_phi(i,j,k)*(PFM_c(i,j+1,k)+PFM_c(i,j,k))/2.0

         !z direction, defined at surface
         m_z_phi_c(i,j,k)=m_z_phi(i,j,k)*(PFM_c(i,j,k+1)+PFM_c(i,j,k))/2.0
      enddo
   enddo
enddo

return
end subroutine m_phi_c




subroutine m_rho_cp
implicit none
integer i,j,k

do k=0,k1   !I am not sure if I need +1
   do j=0,j1
      do i=0,i1
         !x direction,defined at surface
         m_x_rho_cp(i,j,k)=rho_1*cp_1*unew(i,j,k) + &
                      (rho_3*cp_3-rho_1*cp_1)*m_x_phi(i,j,k) + &
                      (rho_2*cp_2-rho_3*cp_3)*m_x_phi_c(i,j,k)

         !y direction, defined at surface
         m_y_rho_cp(i,j,k)=rho_1*cp_1*vnew(i,j,k) + &
                      (rho_3*cp_3-rho_1*cp_1)*m_y_phi(i,j,k) + &
                      (rho_2*cp_2-rho_3*cp_3)*m_y_phi_c(i,j,k)

         !z direction, defined at surface
         m_z_rho_cp(i,j,k)=rho_1*cp_1*wnew(i,j,k) + &
                      (rho_3*cp_3-rho_1*cp_1)*m_z_phi(i,j,k) + &
                      (rho_2*cp_2-rho_3*cp_3)*m_z_phi_c(i,j,k)
      enddo
   enddo
enddo

return
end subroutine m_rho_cp



subroutine m_rho
implicit none
integer i,j,k
!real mxphi,myphi,mzphi,mxphic,myphic,mzphic

do k=0,k1   !I am not sure if I need +1
   do j=0,j1
      do i=0,i1
         !x direction,defined at surface
         m_x_rho(i,j,k)=unew(i,j,k)*rho_1 + (rho_3-rho_1)*m_x_phi(i,j,k) + &
                        (rho_2-rho_3)*m_x_phi_c(i,j,k)

         !y direction, defined at surface
         m_y_rho(i,j,k)=vnew(i,j,k)*rho_1 + (rho_3-rho_1)*m_y_phi(i,j,k) + &
                        (rho_2-rho_3)*m_y_phi_c(i,j,k)


         !z direction, defined at surface
         m_z_rho(i,j,k)=wnew(i,j,k)*rho_1 + (rho_3-rho_1)*m_z_phi(i,j,k) + &
                        (rho_2-rho_3)*m_z_phi_c(i,j,k)

      enddo
   enddo
enddo

return
end subroutine m_rho




end module mod_conserve
