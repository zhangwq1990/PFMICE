module mod_energy
use mod_var
use mod_global
use mod_bou
use mod_thomas
use mod_triperiodic
use mod_field
implicit none
private
!public energy,oldenergy
public energy
contains
!
subroutine energy
implicit none
integer i,j,k
real:: advectivex,advectivey, advectivez,advective
real:: dissipate1, dissipate2, dissipate3,dissipate
real:: conductionx, conductiony, conductionz, conduction,gravity
real:: melt
real:: t_step

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!for 2D pencil transfer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


real, dimension(itot/dims(1),jtot/dims(2),ktot) :: AA_z,BB_z,CC_z,RHS_z  !I think maybe it is not essential to be pointers
real, dimension(itot/dims(1),jtot,ktot/dims(2)) :: AA_y,BB_y,CC_y,RHS_y
!real, dimension(itot,jtot/dims(1),ktot/dims(2)) :: cox1,cox2,cox3,RHS_x 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

t_step=1.0   !we use two sub steps

!first, kk need to be updated
call updatefield(3)


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! first, z direction                               ^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
         conduction = conductionx+conductiony

         !&&&&&&&&&&&&&&&&&&&&&&&&&&&
         ! heat convection
         !&&&&&&&&&&&&&&&&&&&&&&&&&&&

         advectivex= (m_x_rho_cp(i,j,k)*(Tnew(i+1,j,k)+Tnew(i,j,k))/2.0-m_x_rho_cp(i-1,j,k)*(Tnew(i,j,k)+Tnew(i-1,j,k))/2.0 )/dx
         advectivey= (m_y_rho_cp(i,j,k)*(Tnew(i,j+1,k)+Tnew(i,j,k))/2.0-m_y_rho_cp(i,j-1,k)*(Tnew(i,j,k)+Tnew(i,j-1,k))/2.0 )/dy
         advective = advectivex+advectivey

         melt =ct_change(i,j,k)*rho_2*latent

         !gravity = rhol(i,j,k)*(g_x+g_y+g_z) 
         gravity = 0.0

         RHS_z(i,j,k) = gravity + conduction + dissipate - advective - melt + Told(i,j,k)*rholcp_old(i,j,k)/dt*t_step !this 2.000 can be changed

      enddo
   enddo
enddo

do k=1,kmax
   do j=1,jmax
      do i=1,imax          !aa is coefficient for k-1, bb for k, cc for k+1
         AA_z(i,j,k)=-(kk(i,j,k-1)+kk(i,j,k))/2.0/dz2 - m_z_rho_cp(i,j,k-1)/2.0/dz
         BB_z(i,j,k)=rholcp(i,j,k)/dt*t_step+(kk(i,j,k+1)+2.0*kk(i,j,k)+kk(i,j,k-1))/2.0/dz2 + &
                     (m_z_rho_cp(i,j,k)-m_z_rho_cp(i,j,k-1))/2.0/dz !this 2.000 can be changed
         CC_z(i,j,k)=-(kk(i,j,k+1)+kk(i,j,k))/2.0/dz2 + m_z_rho_cp(i,j,k)/2.0/dz
      enddo
   enddo
enddo


   do j=1,jmax
      do i=1,imax          !aa is coefficient for k-1, bb for k, cc for k+1

         !for Neumann at top, just correct bb, as aa and cc are not used in thomas
         BB_z(i,j,kmax)=BB_Z(i,j,kmax)+CC_Z(i,j,kmax)

         !for Dirichlet at bottom, the RHS need to add the bottom value
         RHS_z(i,j,1)=RHS_z(i,j,1)+(kk(i,j,0)+kk(i,j,1))/2.0/dz2*Tnew(i,j,0) + m_z_rho_cp(i,j,0)/2.0/dz*Tnew(i,j,0)
      enddo
   enddo



do j=1,jmax
   do i=1,imax
      call thomas(RHS_z(i,j,1:ktot),AA_z(i,j,1:ktot),BB_z(i,j,1:ktot),CC_z(i,j,1:ktot),ktot)
   enddo
enddo

do k=1,kmax                                   !all these kmax, jmax and imax need to be changed!
   do j=1,jmax
      do i=1,imax
         Tnew(i,j,k)=RHS_z(i,j,k)
      enddo
   enddo
enddo

!then give the boundary condition
call boundT(Tnew)


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! then, y direction                               ^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



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
         conductionz= ( (kkold(i,j,k+1)+kkold(i,j,k))/2.0*(Tnew(i,j,k+1)-Tnew(i,j,k))/dz - &
                        (kkold(i,j,k)+kkold(i,j,k-1))/2.0*(Tnew(i,j,k)-Tnew(i,j,k-1))/dz )/dz
         conduction = conductionx+conductionz

         !&&&&&&&&&&&&&&&&&&&&&&&&&&&
         ! heat convection
         !&&&&&&&&&&&&&&&&&&&&&&&&&&&

         advectivex= (m_x_rho_cp(i,j,k)*(Tnew(i+1,j,k)+Tnew(i,j,k))/2.0-m_x_rho_cp(i-1,j,k)*(Tnew(i,j,k)+Tnew(i-1,j,k))/2.0 )/dx
         advectivez= (m_z_rho_cp(i,j,k)*(Tnew(i,j,k+1)+Tnew(i,j,k))/2.0-m_z_rho_cp(i,j,k-1)*(Tnew(i,j,k)+Tnew(i,j,k-1))/2.0 )/dz
         advective = advectivex+advectivez


         !this and dissipation can be saved but not computed again
         melt= ct_change(i,j,k)*rho_2*latent

         !gravity = rhol(i,j,k)*(g_x+g_y+g_z) 
         gravity = 0.0

         RHS_z(i,j,k) = gravity + conduction + dissipate - advective - melt + Told(i,j,k)*rholcp_old(i,j,k)/dt*t_step !this 2.0 can be changed, and here is Tnew
      enddo
   enddo
enddo



do k=1,ktot
   do j=1,jtot/dims(2)
      do i=1,imax          !aa is coefficient for k-1, bb for k, cc for k+1
         AA_z(i,j,k)=-(kk(i,j-1,k)+kk(i,j,k))/2.0/dy2 - m_y_rho_cp(i,j-1,k)/2.0/dy
         BB_z(i,j,k)=rholcp(i,j,k)/dt*t_step+(kk(i,j+1,k)+2.0*kk(i,j,k)+kk(i,j-1,k))/2.0/dy2 +&
                     (m_y_rho_cp(i,j,k)-m_y_rho_cp(i,j-1,k))/2.0/dy  !this 2.0 can be changed
         CC_z(i,j,k)=-(kk(i,j+1,k)+kk(i,j,k))/2.0/dy2 + m_y_rho_cp(i,j,k)/2.0/dy
      enddo
   enddo
enddo

!for periodic at right, no additional correct is needed
! because I use triperiodic

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



do k=1,kmax
   do j=1,jmax
      do i=1,imax
         Tnew(i,j,k)=RHS_z(i,j,k)
      enddo
   enddo
enddo


return
end subroutine energy
!

end module mod_energy
