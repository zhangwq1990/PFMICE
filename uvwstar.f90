module mod_uvwstar
use mod_cons
use mod_var
use mod_global
use mod_conserve
!use mod_bou
implicit none
private
public uvwstar
contains
!
!
subroutine uvwstar
!
!implicit none
integer i,j,k
real :: convection_x,convection_y,convection_z
real :: diffusion_x,diffusion_y,diffusion_z
real :: surface_phi_x,surface_phi_y,surface_phi_z
real :: as,cd,sm
real :: body_x,body_y,body_z



do k=0,k1
  do j=0,j1
    do i=0,i1
       !as= abs(PFM_phi(i,j,k)- Phi_c(i,j,k))
       as= abs(PFM_phi(i,j,k)*(1.0- PFM_c(i,j,k)))
       !debug5(i,j,k)=as
       cd=rhol(i,j,k)/dt+visl(i,j,k)/dz2  !here use ziyang's creterion
       au(i,j,k)=100000.0*cd*as**0.1/( (1-as)**8+1e-3 )
       !debug6(i,j,k)=au(i,j,k)
    enddo
  enddo
enddo


!call m_rho


do k=1,kmax
   do j=1,jmax
      do i=1,imax                       
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! X direction
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Y direction, convection
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      ! x component
      convection_y= ( (m_x_rho(i,j,k)+m_x_rho(i,j+1,k))/2.0*(vnew(i,j,k)+vnew(i+1,j,k))/2.0 - &
                      (m_x_rho(i-1,j,k)+m_x_rho(i-1,j+1,k))/2.0*(vnew(i-1,j,k)+vnew(i,j,k))/2.0 &
                    )/dx
      ! y component
      convection_y= convection_y + &
                    (m_y_rho(i,j+1,k)*vnew(i,j+1,k)-m_y_rho(i,j-1,k)*vnew(i,j-1,k))/2.0/dy
      ! z component
      convection_y= convection_y + &
                    ( (m_z_rho(i,j,k)+m_z_rho(i,j+1,k))/2.0*(vnew(i,j,k+1)+vnew(i,j,k))/2.0 - &
                      (m_z_rho(i,j,k-1)+m_z_rho(i,j+1,k-1))/2.0*(vnew(i,j,k)+vnew(i,j,k-1))/2.0 &
                    )/dz
      !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&! Y dirction, diffusion
      ! x component
      diffusion_y= ( (visl(i,j,k)+visl(i+1,j,k)+visl(i,j+1,k)+visl(i+1,j+1,k))/4.0*&
                     ( (vnew(i+1,j,k)-vnew(i,j,k))/dx+(unew(i,j+1,k)-unew(i,j,k))/dy ) -&
                     (visl(i,j,k)+visl(i-1,j,k)+visl(i,j+1,k)+visl(i-1,j+1,k))/4.0*&
                     ( (vnew(i,j,k)-vnew(i-1,j,k))/dx+(unew(i-1,j+1,k)-unew(i-1,j,k))/dy ) )/dx
      ! y component
      diffusion_y= diffusion_y+&
                   ( -2.0/3.0*( visl(i,j+1,k)*((unew(i,j+1,k)-unew(i-1,j+1,k))/dx + &
                              (vnew(i,j+1,k)-vnew(i,j,k))/dy+(wnew(i,j+1,k)-wnew(i,j+1,k-1))/dz) -& 
                                visl(i,j,k)*((unew(i,j,k)-unew(i-1,j,k))/dx+(vnew(i,j,k)-vnew(i,j-1,k))/dy+&
                                (wnew(i,j,k)-wnew(i,j,k-1))/dz) )+&
                     2.0*visl(i,j+1,k)*(vnew(i,j+1,k)-vnew(i,j,k))/dy - 2.0*visl(i,j,k)*(vnew(i,j,k)-vnew(i,j-1,k))/dy  )/dy
      ! z component
      diffusion_y=diffusion_y+&
                   ( (visl(i,j,k)+visl(i,j+1,k)+visl(i,j,k+1)+visl(i,j+1,k+1))/4.0*&
                     ( (wnew(i,j+1,k)-wnew(i,j,k))/dy+(vnew(i,j,k+1)-vnew(i,j,k))/dz )-&  
                     (visl(i,j,k)+visl(i,j+1,k)+visl(i,j,k-1)+visl(i,j+1,k-1))/4.0*&
                     ( (wnew(i,j+1,k-1)-wnew(i,j,k-1))/dy+(vnew(i,j,k)-vnew(i,j,k-1))/dz )  )/dz

      surface_phi_y=(PFM_c(i,j,k)+PFM_c(i,j+1,k))/2.0*&
                    !(chem_pot_old(i,j+1,k)+chem_pot_old(i,j,k))/2.0*(PFM_phi_old(i,j+1,k)-PFM_phi_old(i,j,k))/dy
                    (chem_pot(i,j+1,k)+chem_pot(i,j,k))/2.0*(PFM_phi(i,j+1,k)-PFM_phi(i,j,k))/dy


      body_y=0.
      
      RHS_y(i,j,k)=-convection_y+diffusion_y+body_y+surface_phi_y


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Z direction
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo convection term
      ! x component
      convection_z= ( (m_x_rho(i,j,k)+m_x_rho(i,j,k))/2.0*(wnew(i,j,k)+wnew(i+1,j,k))/2.0 - &
                      (m_x_rho(i-1,j,k)+m_x_rho(i-1,j,k))/2.0*(wnew(i-1,j,k)+wnew(i,j,k))/2.0 &
                    )/dx
      ! y component
      convection_z= convection_z + &
                    ( (m_y_rho(i,j,k)+m_y_rho(i,j,k+1))/2.0*(wnew(i,j,k)+wnew(i,j+1,k))/2.0 - &
                      (m_y_rho(i,j-1,k)+m_y_rho(i,j-1,k+1))/2.0*(wnew(i,j-1,k)+wnew(i,j,k))/2.0 &
                    )/dy
      ! z component
      convection_z= convection_z + &
                    ( m_z_rho(i,j,k+1)*wnew(i,j,k+1)-m_z_rho(i,j,k-1)*wnew(i,j,k-1) )/2.0/dz
      !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&! diffusion term
      ! x component
      diffusion_z=( (visl(i,j,k)+visl(i+1,j,k)+visl(i,j,k+1)+visl(i+1,j,k+1))/4.0*&
                    ( (unew(i,j,k+1)-unew(i,j,k))/dz + (wnew(i+1,j,k)-wnew(i,j,k))/dx )-&
                    (visl(i,j,k)+visl(i-1,j,k)+visl(i,j,k+1)+visl(i-1,j,k+1))/4.0*&
                    ( (unew(i-1,j,k+1)-unew(i-1,j,k))/dz + (wnew(i,j,k)-wnew(i-1,j,k))/dx )   )/dx
      ! y component
      diffusion_z=diffusion_z+&
                  ( (visl(i,j,k)+visl(i,j+1,k)+visl(i,j,k+1)+visl(i,j+1,k+1))/4.0*&
                    ( (wnew(i,j+1,k)-wnew(i,j,k))/dy + (vnew(i,j,k+1)-vnew(i,j,k))/dz )-&
                    (visl(i,j,k)+visl(i,j-1,k)+visl(i,j,k+1)+visl(i,j-1,k+1))/4.0*&
                    ( (wnew(i,j,k)-wnew(i,j-1,k))/dy + (vnew(i,j-1,k+1)-vnew(i,j-1,k))/dz )   )/dy
      ! z component
      diffusion_z=diffusion_z&
                   -2.0/3.0*( visl(i,j,k+1)*((unew(i,j,k+1)-unew(i-1,j,k+1))/dx+ &
                             (vnew(i,j,k+1)-vnew(i,j-1,k+1))/dy+(wnew(i,j,k+1)-wnew(i,j,k))/dz) -&
                               visl(i,j,k)*((unew(i,j,k)-unew(i-1,j,k))/dx+&
                             (vnew(i,j,k)-vnew(i,j-1,k))/dy+(wnew(i,j,k)-wnew(i,j,k-1))/dz) )/dz +&
                   ( 2.0*visl(i,j,k+1)*(wnew(i,j,k+1)-wnew(i,j,k))/dz - 2.0*visl(i,j,k)*(wnew(i,j,k)-wnew(i,j,k-1))/dz  )/dz


      surface_phi_z=(PFM_c(i,j,k)+PFM_c(i,j,k+1))/2.0*&
                    !(chem_pot_old(i,j,k+1)+chem_pot_old(i,j,k))/2.0*(PFM_phi_old(i,j,k+1)-PFM_phi_old(i,j,k))/dz
                    (chem_pot(i,j,k+1)+chem_pot(i,j,k))/2.0*(PFM_phi(i,j,k+1)-PFM_phi(i,j,k))/dz

      !body_z=g_z*(rhol(i,j,k+1)+rhol(i,j,k))/2.0
      body_z=0.
      
      RHS_z(i,j,k)=-convection_z+diffusion_z+body_z+surface_phi_z



      

      enddo
   enddo
enddo





! this can be incorporate into last loop

do k=1,kmax
   do j=1,jmax
      do i=1,imax
         ustar(i,j,k)= 0.
         vstar(i,j,k)= (RHS_y(i,j,k)*dt+(rholold(i,j,k)+rholold(i,j+1,k))/2.0*vold(i,j,k))/( (rhol(i,j,k)+rhol(i,j+1,k))/2.0 -&
                                              (residual1(i,j,k)+residual1(i,j+1,k))/2.0*dt +&
                                              (au(i,j,k)+au(i,j+1,k))/2.0*dt)
         wstar(i,j,k)= (RHS_z(i,j,k)*dt+(rholold(i,j,k)+rholold(i,j,k+1))/2.0*wold(i,j,k))/( (rhol(i,j,k)+rhol(i,j,k+1))/2.0 -&
                                              (residual1(i,j,k)+residual1(i,j,k+1))/2.0*dt +&
                                              (au(i,j,k)+au(i,j,k+1))/2.0*dt)
      enddo
   enddo
enddo







return
end subroutine uvwstar

!
end module mod_uvwstar
