module mod_cim
use mod_cons
use mod_var
use mod_global
use mod_bou
use mod_conserve
use mod_thomas
use mod_triperiodic
implicit none
private
public updatecim
contains
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!I have a question: how to consider the boundary condition regariding the icing near wall
! c range now is 0 to 1, 1 means full water and 0 means full ice
subroutine updatecim
!implicit none
integer i,j,k
real:: advective,lap_c,c_star,temp,source,pterm
real:: t_step,num_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!for 2D pencil transfer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real, dimension(itot/dims(1),jtot/dims(2),ktot) :: AA_z,BB_z,CC_z,RHS_z  !I think maybe it is not essential to be pointers
real, dimension(itot/dims(1),jtot,ktot/dims(2)) :: AA_y,BB_y,CC_y,RHS_y
!real, dimension(itot,jtot/dims(1),ktot/dims(2)) :: cox1,cox2,cox3,RHS_x 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


t_step=1.0   !we use two sub steps
num_s=1e-30

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!first, let us calculate the z direction
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!calculate RHS 
do k=1,kmax
   do j=1,jmax
      do i=1,imax

         advective = (m_x_phi_c(i,j,k)-m_x_phi_c(i-1,j,k))/dx + &
                     (m_y_phi_c(i,j,k)-m_y_phi_c(i,j-1,k))/dy + &
                     (m_z_phi_c(i,j,k)-m_z_phi_c(i,j,k-1))/dz 
         !if (k==kmax) then
         !   advective = (m_x_phi_c(i,j,k)-m_x_phi_c(i-1,j,k))/dx + &
         !               (m_y_phi_c(i,j,k)-m_y_phi_c(i,j-1,k))/dy + &
         !               (m_z_phi_c(i,j,k-1)-m_z_phi_c(i,j,k-2))/dz 
         !endif

         !use phi_old, because phi has already been updated
         lap_c= -lamda_c  *( ( (PFM_c(i+1,j,k)-PFM_c(i,j,k))/dx*(PFM_phi(i+1,j,k)+PFM_phi(i,j,k))/2.0 - &
                               (PFM_c(i,j,k)-PFM_c(i-1,j,k))/dx*(PFM_phi(i,j,k)+PFM_phi(i-1,j,k))/2.0  )/dx +&
                             ( (PFM_c(i,j+1,k)-PFM_c(i,j,k))/dy*(PFM_phi(i,j+1,k)+PFM_phi(i,j,k))/2.0 - &
                               (PFM_c(i,j,k)-PFM_c(i,j-1,k))/dy*(PFM_phi(i,j,k)+PFM_phi(i,j-1,k))/2.0  )/dy &
                           )

         c_star=2.0*PFM_c(i,j,k)**2*(3-4.0*PFM_c(i,j,k))*lamda_c/PFM_l_c**2*PFM_phi(i,j,k)

         if (Tnew(i,j,k)>t_melt .and. PFM_c(i,j,k)==0.) then !here the temperature is in degree, not Kelvin
            pterm=1.0
         else if (Tnew(i,j,k)<t_melt .and. PFM_c(i,j,k)==1 ) then
            pterm=1.0
         else
           !pterm=-30.0*PFM_c(i,j,k)**2*(3.0*PFM_c(i,j,k)**2-4.0*PFM_c(i,j,k)+1.0)
           pterm=30.0*PFM_c(i,j,k)**2*(3.0*PFM_c(i,j,k)**2-4.0*PFM_c(i,j,k)+1.0)
         endif
         temp=rho_2*latent/t_melt*PFM_phi(i,j,k)*(t_melt-Tnew(i,j,k)) * pterm


         source = Phi_c_old(i,j,k)*( (unew(i,j,k)-unew(i-1,j,k))/dx +&
                                     (vnew(i,j,k)-vnew(i,j-1,k))/dy +&
                                     (wnew(i,j,k)-wnew(i,j,k-1))/dz )

         RHS_z(i,j,k)=-advective-mobility_c*(lap_c+c_star+temp)+source+Phi_c_old(i,j,k)/dt*t_step
         debug1(i,j,k)=RHS_z(i,j,k)
      enddo
   enddo
enddo


!calculate the coefficient of the matrix
!aa is coefficient for k-1, bb for k, cc for k+1
do k=1,kmax
   do j=1,jmax
      do i=1,imax          !aa is coefficient for k-1, bb for k, cc for k+1
         AA_z(i,j,k)=-mobility_c*lamda_c*(PFM_phi(i,j,k)+PFM_phi(i,j,k-1)+2.0*num_s)/2.0/dz2

         BB_z(i,j,k)=(PFM_phi(i,j,k)+num_s)/dt*t_step  &
                     +mobility_c*lamda_c*(PFM_phi(i,j,k+1)+2.0*PFM_phi(i,j,k)+PFM_phi(i,j,k-1)+4.0*num_s)/2.0/dz2 &
                     +mobility_c*lamda_c/PFM_l_c**2*(PFM_phi(i,j,k)+num_s)*(12.0*PFM_c(i,j,k)**2-12.0*PFM_c(i,j,k)+2.0) &
                     +mobility_c*rho_2*latent/t_melt*(PFM_phi(i,j,k)+num_s)*(t_melt-Tnew(i,j,k)) &
                     *(-60.0)*PFM_c(i,j,k)*(2.0*PFM_c(i,j,k)-1.0)*(PFM_c(i,j,k)-1.0)

                     
         CC_z(i,j,k)=-mobility_c*lamda_c*(PFM_phi(i,j,k+1)+PFM_phi(i,j,k)+2.0*num_s)/2.0/dz2
      enddo
   enddo
enddo


!revise the coefficient of the matrix according to the boundary condition
!aa is coefficient for k-1, bb for k, cc for k+1

do j=1,jmax
   do i=1,imax 

      !for Neumann at top, just correct bb, as aa and cc are not used in thomas
      BB_z(i,j,kmax)=BB_z(i,j,kmax)+CC_z(i,j,kmax)

      !for Neumann at bottom, the RHS need to add the bottom value
      BB_z(i,j,1)   =BB_z(i,j,1)   +AA_z(i,j,1)
   enddo
enddo


!then calculate the C value
do j=1,jmax
   do i=1,imax
      call thomas(RHS_z(i,j,1:ktot),AA_z(i,j,1:ktot),BB_z(i,j,1:ktot),CC_z(i,j,1:ktot),ktot)
   enddo
enddo


! then transfer the value from RHS to C
do k=1,kmax
   do j=1,jmax
      do i=1,imax
         PFM_c(i,j,k)=RHS_z(i,j,k)
         debug2(i,j,k)=RHS_z(i,j,k)
      enddo
   enddo
enddo

call boundc(PFM_c)


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!then, let us calculate the y direction
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!calculate RHS 
do k=1,kmax
   do j=1,jmax
      do i=1,imax
         ! in fact, this advective does not need to computed again, should save it
         advective = (m_x_phi_c(i,j,k)-m_x_phi_c(i-1,j,k))/dx + &
                     (m_y_phi_c(i,j,k)-m_y_phi_c(i,j-1,k))/dy + &
                     (m_z_phi_c(i,j,k)-m_z_phi_c(i,j,k-1))/dz 
         !use phi_old, because phi has already been updated
         lap_c= -lamda_c  *( ( (PFM_c(i+1,j,k)-PFM_c(i,j,k))/dx*(PFM_phi(i+1,j,k)+PFM_phi(i,j,k))/2.0 - &
                               (PFM_c(i,j,k)-PFM_c(i-1,j,k))/dx*(PFM_phi(i,j,k)+PFM_phi(i-1,j,k))/2.0  )/dx +&
                             ( (PFM_c(i,j,k+1)-PFM_c(i,j,k))/dz*(PFM_phi(i,j,k+1)+PFM_phi(i,j,k))/2.0 - &
                               (PFM_c(i,j,k)-PFM_c(i,j,k-1))/dz*(PFM_phi(i,j,k)+PFM_phi(i,j,k-1))/2.0  )/dz &
                           )

         c_star=2.0*PFM_c(i,j,k)**2*(3-4.0*PFM_c(i,j,k))*lamda_c/PFM_l_c**2*PFM_phi(i,j,k)

         if (Tnew(i,j,k)>t_melt .and. PFM_c(i,j,k)==0.) then !here the temperature is in degree, not Kelvin
            pterm=1.0
         else if (Tnew(i,j,k)<t_melt .and. PFM_c(i,j,k)==1 ) then
            pterm=1.0
         else
           pterm=30.0*PFM_c(i,j,k)**2*(3.0*PFM_c(i,j,k)**2-4.0*PFM_c(i,j,k)+1.0)
         endif
         temp=rho_2*latent/t_melt*PFM_phi(i,j,k)*(t_melt-Tnew(i,j,k)) * pterm


         source = Phi_c_old(i,j,k)*( (unew(i,j,k)-unew(i-1,j,k))/dx +&
                                     (vnew(i,j,k)-vnew(i,j-1,k))/dy +&
                                     (wnew(i,j,k)-wnew(i,j,k-1))/dz )

         RHS_z(i,j,k)=-advective-mobility_c*(lap_c+c_star+temp)+source+ Phi_c_old(i,j,k)/dt*t_step
      enddo
   enddo
enddo


!calculate the coefficient of the matrix
!aa is coefficient for j-1, bb for j, cc for j+1
do k=1,kmax
   do j=1,jtot/dims(2)
      do i=1,imax
         AA_z(i,j,k)=-mobility_c*lamda_c*(PFM_phi(i,j,k)+PFM_phi(i,j-1,k)+2.0*num_s)/2.0/dy2

         BB_z(i,j,k)=(PFM_phi(i,j,k)+num_s)/dt*t_step &
                     +mobility_c*lamda_c*(PFM_phi(i,j+1,k)+2.0*PFM_phi(i,j,k)+PFM_phi(i,j-1,k)+4.0*num_s)/2.0/dy2 &
                     +mobility_c*lamda_c/PFM_l_c**2*(PFM_phi(i,j,k)+num_s)*(12.0*PFM_c(i,j,k)**2-12.0*PFM_c(i,j,k)+2.0) &
                     +mobility_c*rho_2*latent/t_melt*(PFM_phi(i,j,k)+num_s)*(t_melt-Tnew(i,j,k)) &
                     *(-60.0)*PFM_c(i,j,k)*(2.0*PFM_c(i,j,k)-1.0)*(PFM_c(i,j,k)-1.0)
                     
         CC_z(i,j,k)=-mobility_c*lamda_c*(PFM_phi(i,j+1,k)+PFM_phi(i,j,k)+2.0*num_s)/2.0/dy2
      enddo
   enddo
enddo



!for periodic boundary condition, no additional revision is needed for the coefficient, I will do it with triperiodic


! then I need to do 2D pencil
! then transfer the y and z direction
AA_y=0.
BB_y=0.
CC_y=0.
RHS_y=0.
call transpose_z_to_y(AA_z, AA_y)
call transpose_z_to_y(BB_z, BB_y)
call transpose_z_to_y(CC_z, CC_y)
call transpose_z_to_y(RHS_z, RHS_y)

! then calculate along y direction as implicit, Z as explicit


do i=1,imax   !here the sequnce can be changed
   do k=1,ktot/dims(2)
      call triperiodic(RHS_y(i,1:jtot,k),AA_y(i,1:jtot,k),BB_y(i,1:jtot,k),CC_y(i,1:jtot,k),jtot)
   enddo
enddo

!then transpose back
call transpose_y_to_z(RHS_y, RHS_z)


! then transfer the value from RHS to C
do k=1,kmax
   do j=1,jmax
      do i=1,imax
         PFM_c(i,j,k)=RHS_z(i,j,k)
         !debug2(i,j,k)=RHS_z(i,j,k)
         Phi_c(i,j,k)=PFM_c(i,j,k)*(PFM_phi(i,j,k)+1e-30)   !this is right, or temperature will NAN
      enddo
   enddo
enddo



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! at last, let us do the mapping
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

do k=1,kmax
   do j=1,jmax
      do i=1,imax  !here use new phi value
         if ( (abs(PFM_phi(i,j,k)-Phi_c(i,j,k))<1e-12) .or. (PFM_phi(i,j,k)<1e-15) ) then
            PFM_c(i,j,k)=1.0
         else
            PFM_c(i,j,k)=MIN(MAX(PFM_c(i,j,k), 0.0),1.0)
         endif
      enddo
   enddo
enddo



return
end subroutine updatecim

end module mod_cim
