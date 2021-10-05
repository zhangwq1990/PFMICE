module mod_field
use mod_var
implicit none
private
public updatefield
contains
!
subroutine updatefield(choose)
implicit none
integer i,j,k
integer, intent(in) ::  choose

select case(choose)

case(0)
  do k=-2,k1+2
   do j=-2,j1+2
    do i=-2,i1+2
     rhol(i,j,k)= rho_1+(rho_3-rho_1)*PFM_phi(i,j,k)+(rho_2-rho_3)*Phi_c(i,j,k)
     visl(i,j,k)= vis_1+(vis_3-vis_1)*PFM_phi(i,j,k)+(vis_2-vis_3)*Phi_c(i,j,k)
     kk(i,j,k)= kk_1+(kk_3-kk_1)*PFM_phi(i,j,k)+(kk_2-kk_3)*Phi_c(i,j,k)
     rholcp(i,j,k)=rho_1*cp_1+(rho_3*cp_3-rho_1*cp_1)*PFM_phi(i,j,k)+&
                   (rho_2*cp_2-rho_3*cp_3)*Phi_c(i,j,k)
    enddo
   enddo
  enddo
case(1)
  do k=-2,k1+2
   do j=-2,j1+2
    do i=-2,i1+2
     rhol(i,j,k)= rho_1+(rho_3-rho_1)*PFM_phi(i,j,k)+(rho_2-rho_3)*Phi_c(i,j,k)
    enddo
   enddo
  enddo
case(2)
  do k=-2,k1+2
   do j=-2,j1+2
    do i=-2,i1+2
     visl(i,j,k)= vis_1+(vis_3-vis_1)*PFM_phi(i,j,k)+(vis_2-vis_3)*Phi_c(i,j,k)
    enddo
   enddo
  enddo
case(3)
  do k=-2,k1+2
   do j=-2,j1+2
    do i=-2,i1+2
     kk(i,j,k)= kk_1+(kk_3-kk_1)*PFM_phi(i,j,k)+(kk_2-kk_3)*Phi_c(i,j,k)
    enddo
   enddo
  enddo
case(4)
  do k=-2,k1+2
   do j=-2,j1+2
    do i=-2,i1+2
     rholcp(i,j,k)=rho_1*cp_1+(rho_3*cp_3-rho_1*cp_1)*PFM_phi(i,j,k)+&
                   (rho_2*cp_2-rho_3*cp_3)*Phi_c(i,j,k)
    enddo
   enddo
  enddo
end select
end subroutine updatefield

end module mod_field
