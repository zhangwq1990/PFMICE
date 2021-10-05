module mod_init
use mod_var
use decomp_2d
use mod_cons
use mod_global
implicit none
private
public init
contains
subroutine init
implicit none
real :: distance,xxx,yyy,zzz,DropletsCenty,DropletsCentz
integer::i,j,k
unew(:,:,:)    = 0.
vnew(:,:,:)    = 0.
wnew(:,:,:)    = 0.
pnew(:,:,:)    = 0.
PFM_phi(:,:,:) = 0.
PFM_c(:,:,:) = 1.0 !all liquid
Tnew(:,:,:)= 274  !higher than freezing temperature

! set variable names for output
variables(1)='PFM_phi'
variables(2)='PFM_c'
variables(3)='rhol'
variables(4)='Tnew'
variables(5)='phi_c'
variables(6)='residual'



DropletsCenty= ly/2. !ly/2 !(droplet_radius+0.5*(Droplets_separation-2*droplet_radius))
DropletsCentz= 0.85*droplet_radius   !this is for quasi-2D
!DropletsCentz= 3.85*droplet_radius   !this is for quasi-2D
do k=0,k1
  do j=0,j1
   do i=0,i1
      xxx=(i+coords(1)*imax-1)*dx+0.5*dx
      yyy=(j+coords(2)*jmax-1)*dy+0.5*dy
      zzz=(k               -1)*dz+0.5*dz
      distance= sqrt(( yyy-DropletsCenty )**2+(zzz-DropletsCentz )**2)
      PFM_phi(i,j,k) = -(tanh((distance-droplet_radius)/(PFM_l*sqrt(2.)))-1)/2.0
      !PFM_phi(i,j,k) = -(tanh((zzz-droplet_radius*4.0)/(PFM_l*sqrt(2.)))-1)/2.0

      Phi_c(i,j,k)=PFM_phi(i,j,k)*PFM_c(i,j,k)
   enddo
  enddo
enddo

return
end subroutine init
!
end module mod_init
