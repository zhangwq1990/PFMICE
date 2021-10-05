module mod_phi_im
use mod_var
use mod_global
use mod_bou
use mod_cons
use mod_chem
use decomp_2d
implicit none
private init_bigphi, solver2ddd
public updatephi_im
contains
!
subroutine updatephi_im
implicit none
integer i,j,k
!real :: factor1,factor2
real::advectivex,advectivey, advectivez,advective,source,phi_lap
!real :: t11,t22,t33,t44,t55
!
real :: middle, SSSS,aaaa,xishu
real, target, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: Q,dPFM
real, pointer, dimension(:,:,:) :: qqq
real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: g_prime


real :: xst(itot), yst(jtot)
real, dimension(imax,jmax) :: xyst
!real :: bb
real(mytype) ::  si(itot+15), sj(jtot+15)
real, dimension(kmax) :: a0,b0,c0



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
middle=3.0*sqrt(2.)*PFM_Sigma_ref*PFM_l  !this is the gamma1 in the paper
SSSS=1.0*PFM_l**2*sqrt(4.0*1.0/mobility/middle/dt)  !this is the S in the paper
aaaa=-SSSS/2.0/PFM_l**2*( 1+sqrt(1-4.0*1.0/middle/mobility/dt*PFM_l**4/SSSS**2) )  !this is the alpha in the paper
xst=0
yst=0
xyst=0
a0=0.
b0=0.
c0=0.
si=0
sj=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! here start my own code
xishu= -(aaaa+SSSS/PFM_l**2)    !it is negative
call init_bigphi(xishu,xst,yst,xyst,si,sj,a0,b0,c0)
qqq => Q(1:imax,1:jmax,1:kmax)



do k=0,k1
  do j=0,j1
    do i=0,i1
!do k=-2,k1+2
!  do j=-2,j1+2
!    do i=-2,i1+2
       g_prime(i,j,k)=(2.0*PFM_phi(i,j,k)*(PFM_phi(i,j,k)-1)*(2.0*PFM_phi(i,j,k)-1)-SSSS*PFM_phi(i,j,k))/PFM_l**2  !this has an additional term
    enddo
  enddo
enddo


do k=1,kmax
  do j=1,jmax
    do i=1,imax
       source    = PFM_phi(i,j,k)*( (unew(i,j,k)-unew(i-1,j,k))/dx + &
                                    (vnew(i,j,k)-vnew(i,j-1,k))/dy + (wnew(i,j,k)-wnew(i,j,k-1))/dz )

       advectivex=( unew(i,j,k)*(PFM_phi(i+1,j,k)+PFM_phi(i,j,k))/2.0-unew(i-1,j,k)*(PFM_phi(i,j,k)+PFM_phi(i-1,j,k))/2.0 )/dx !conservative
       advectivey=( vnew(i,j,k)*(PFM_phi(i,j+1,k)+PFM_phi(i,j,k))/2.0-vnew(i,j-1,k)*(PFM_phi(i,j,k)+PFM_phi(i,j-1,k))/2.0 )/dy
       advectivez=( wnew(i,j,k)*(PFM_phi(i,j,k+1)+PFM_phi(i,j,k))/2.0-wnew(i,j,k-1)*(PFM_phi(i,j,k)+PFM_phi(i,j,k-1))/2.0 )/dz
       advective=advectivex+advectivey+advectivez

       phi_lap= (g_prime(i+1,j,k)+g_prime(i-1,j,k)-2.0*g_prime(i,j,k))/dx2 + &
                (g_prime(i,j+1,k)+g_prime(i,j-1,k)-2.0*g_prime(i,j,k))/dy2 + &
                (g_prime(i,j,k+1)+g_prime(i,j,k-1)-2.0*g_prime(i,j,k))/dz2 
               
       Q(i,j,k)=( PFM_phi(i,j,k)/dt-advective+source )/middle/mobility + phi_lap
       !debug1(i,j,k)=Q(i,j,k)
    enddo
  enddo
enddo


call solver2ddd(qqq,xyst,si,sj,a0,b0,c0)   !then here we get the value of bigphi
!call boundChem(Q)    !here Q has the value of bigphi, this is not essential


!do k=1,kmax
!  do j=1,jmax
!    do i=1,imax
!       debug2(i,j,k)=(Q(i-1,j,k)+Q(i+1,j,k)-2.0*Q(i,j,k))/dx2+&
!                     (Q(i,j-1,k)+Q(i,j+1,k)-2.0*Q(i,j,k))/dy2+&
!                     (Q(i,j,k-1)+Q(i,j,k+1)-2.0*Q(i,j,k))/dz2+&
!                     xishu*Q(i,j,k)-debug1(i,j,k)
!       debug3(i,j,k)=Q(i,j,k)
!    enddo
!  enddo
!enddo




xishu= aaaa    !it is negative
call init_bigphi(xishu,xst,yst,xyst,si,sj,a0,b0,c0)
call solver2ddd(qqq,xyst,si,sj,a0,b0,c0)   !then here we get the value of bigphi
!call boundChem(Q)    !here Q has the value of bigphi, this is not essential
!unless different boundary condition is used, not the Neumann scheme



!do k=1,kmax
!  do j=1,jmax
!    do i=1,imax
!       debug4(i,j,k)=(Q(i-1,j,k)+Q(i+1,j,k)-2.0*Q(i,j,k))/dx2+&
!                     (Q(i,j-1,k)+Q(i,j+1,k)-2.0*Q(i,j,k))/dy2+&
!                     (Q(i,j,k-1)+Q(i,j,k+1)-2.0*Q(i,j,k))/dz2+&
!                     xishu*Q(i,j,k)-debug3(i,j,k)
!       debug5(i,j,k)=Q(i,j,k)
!    enddo
!  enddo
!enddo




do k=1,kmax
  do j=1,jmax
    do i=1,imax
       PFM_phi(i,j,k)=Q(i,j,k)
       !debug1(i,j,k)=PFM_phi(i,j,k)
    enddo
  enddo
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calculate the overall residual
!call boundPFM(PFM_phi,dPFM_bound_old,dPFM_bound_old_top,unew,vnew,wnew)
!call ChemicalPotential_im
!do k=1,kmax
!  do j=1,jmax
!    do i=1,imax
!        dPFM(i,j,k)=mobility*middle*SSSS/PFM_l**2*(PFM_phi(i,j,k)-PFM_phi_old(i,j,k))
!!       t11 = 2.0*PFM_phi_old(i,j,k)*(PFM_phi_old(i,j,k)-1)*(2.0*PFM_phi_old(i,j,k)-1)
!!       t22 = (PFM_phi(i+1,j,k)-2.*PFM_phi(i,j,k)+PFM_phi(i-1,j,k))/dx2
!!       t33 = (PFM_phi(i,j+1,k)-2.*PFM_phi(i,j,k)+PFM_phi(i,j-1,k))/dy2
!!       t44 = (PFM_phi(i,j,k+1)-2.*PFM_phi(i,j,k)+PFM_phi(i,j,k-1))/dz2
!!       t55= t22+t33+t44
!!       Chem_Pot(i,j,k) =3.0*sqrt(2.0)*(PFM_Sigma_ref/PFM_l*t11-PFM_Sigma_ref*PFM_l*t55)
!    enddo
!  enddo
!enddo
!call boundChem(dPFM)
!call boundChem(chem_pot)
!do k=1,kmax
!  do j=1,jmax
!    do i=1,imax
!       source    = PFM_phi_old(i,j,k)*( (unew(i,j,k)-unew(i-1,j,k))/dx + &
!                                    (vnew(i,j,k)-vnew(i,j-1,k))/dy + (wnew(i,j,k)-wnew(i,j,k-1))/dz )
!
!       advectivex=( unew(i,j,k)*(PFM_phi_old(i+1,j,k)+PFM_phi_old(i,j,k))/2.0-&
!                    unew(i-1,j,k)*(PFM_phi_old(i,j,k)+PFM_phi_old(i-1,j,k))/2.0 )/dx !conservative
!       advectivey=( vnew(i,j,k)*(PFM_phi_old(i,j+1,k)+PFM_phi_old(i,j,k))/2.0-&
!                    vnew(i,j-1,k)*(PFM_phi_old(i,j,k)+PFM_phi_old(i,j-1,k))/2.0 )/dy
!       advectivez=( wnew(i,j,k)*(PFM_phi_old(i,j,k+1)+PFM_phi_old(i,j,k))/2.0-&
!                    wnew(i,j,k-1)*(PFM_phi_old(i,j,k)+PFM_phi_old(i,j,k-1))/2.0 )/dz
!
!       advective=advectivex+advectivey+advectivez
!
!       debug6(i,j,k)=(PFM_phi(i,j,k)-PFM_phi_old(i,j,k))/dt+&
!                     advective-source- mobility*&                    
!                     ((Chem_pot(i-1,j,k)+Q(i+1,j,k)-2.0*Chem_pot(i,j,k))/dx2+&
!                     (Chem_pot(i,j-1,k)+Q(i,j+1,k)-2.0*Chem_pot(i,j,k))/dy2+&
!                     (Chem_pot(i,j,k-1)+Q(i,j,k+1)-2.0*Chem_pot(i,j,k))/dz2 )
!                     !((dPFM(i-1,j,k)+dPFM(i+1,j,k)-2.0*dPFM(i,j,k))/dx2+&
!                     !(dPFM(i,j-1,k)+dPFM(i,j+1,k)-2.0*dPFM(i,j,k))/dy2+&
!                     !(dPFM(i,j,k-1)+dPFM(i,j,k+1)-2.0*dPFM(i,j,k))/dz2 )
!    enddo
!  enddo
!enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






!!update the boundary to n+1
!do j=1,jmax
!  do i=1,imax
!    Phi_wall_n1=phi_wall_n(i,j)+factor1*dPFM_bound_n(i,j)+factor2*dPFM_bound_old(i,j)
!    Phi_wall_n1_top=phi_wall_n_top(i,j)+factor1*dPFM_bound_n_top(i,j)+factor2*dPFM_bound_old_top(i,j)
!    if (abs(PFM_mu_f).gt.1e-8) then
!       PFM_phi(i,j,0) = 2*Phi_wall_n1-PFM_phi(i,j,1)
!       dPFM_bound_old(i,j)=dPFM_bound_n(i,j)
!    endif
!
!  enddo
!enddo
!


return
end subroutine updatephi_im


!!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
subroutine init_bigphi(xishu,xst,yst,xyst,si,sj,a0,b0,c0)
implicit none
integer :: i,j,k,iv,jv
real, dimension(itot),intent(inout) :: xst
real, dimension(jtot),intent(inout) :: yst
real, dimension(imax,jmax),intent(inout) :: xyst
real(mytype),dimension(itot+15),intent(inout) ::  si
real(mytype),dimension(jtot+15),intent(inout) ::  sj
real, dimension(kmax),intent(inout) :: a0,b0,c0
real :: xishu
!
! Generate tridiagonal matrix
!
!print*, xishu
do k=1,kmax
  a0(k) = dzi*dzi
  c0(k) = dzi*dzi
  b0(k) = -(a0(k) + c0(k))
enddo

! Neumann boundary condition for correction big_phi in z-direction;
! For big_phi, I only need Neumann, but Q need to be changed as well, because this is no a 
! simple Neumann
!
b0(1) = b0(1) + a0(1)
b0(kmax) = b0(kmax) + c0(kmax)
a0(1) = 0.
c0(kmax) = 0.
!


! set lookup tables.
!
call vrffti(itot,si)
call vrffti(jtot,sj)
!
! generate eigenvalues ( xrt and yrt ).
!
!
! x direction
!
do i=3,itot,2
  xst(i-1) = -4.*dxi*dxi*(sin(float((i-1))*pi/(2.*itot)))**2
  xst(i) = xst(i-1)
enddo
xst(1   ) = 0.
xst(itot) = -4.*dxi*dxi
!
! y direction
!
do j=3,jtot,2
  yst(j-1) = -4.*dyi*dyi*(sin(float((j-1))*pi/(2.*jtot)))**2
  yst(j) = yst(j-1)
enddo
yst(1   ) = 0.
yst(jtot) = -4.*dyi*dyi

!!Neumann-Neumann
!do j=1,jtot
!  yst(j) = -4.*dyi*dyi*(sin(float((j-1))*pi/(2.*jtot)))**2
!enddo


do j=1,jmax
  jv = j + zstart(2) - 1
  do i=1,imax
    iv = i + zstart(1) - 1
    xyst(i,j) = xst(iv)+yst(jv)+xishu   !this is the difference between Poisson and Helmholtz
  enddo
enddo
!
return
end subroutine init_bigphi
!


subroutine solver2ddd(phiz,xyst,si,sj,a0,b0,c0)
use decomp_2d
implicit none
real, intent(inout), dimension(1:,1:,1:) :: phiz
real, dimension(itot/dims(1),jtot,ktot/dims(2)) :: phiy
real, dimension(itot,jtot/dims(1),ktot/dims(2)) :: phix
real :: bb
real :: z,d(imax,jmax,kmax)
real :: di(itot),dj(jtot)
integer :: i,j,k
!real, intent(in) :: xst(itot), yst(jtot)
real, dimension(imax,jmax),intent(in) :: xyst
real(mytype),intent(in) ::  si(itot+15), sj(jtot+15)
real, dimension(kmax),intent(in) :: a0,b0,c0




!
call transpose_z_to_y(phiz,phiy)
call transpose_y_to_x(phiy,phix)
!
do k=1,xsize(3)
  do j=1,xsize(2)
    call vrfftf(1,itot,phix(1:itot,j,k),di,1,si)
  enddo
enddo
!

!!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!do i=1,itot
!   phix(1:itot,j,k)=phix(1:itot,j,k)*0.5*(1+cos(pi*float(i-1)/itot))
!enddo
!!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



call transpose_x_to_y(phix,phiy)
!
do k=1,ysize(3)
  do i=1,ysize(1)
    call vrfftf(1,jtot,phiy(i,1:jtot,k),dj,1,sj)
  enddo
enddo
!

!!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!do i=1,itot
!   phiy(i,1:jtot,k)=phiy(i,1:jtot,k)*0.5*(1+cos(pi*float(i-1)/jtot))
!enddo
!!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



call transpose_y_to_z(phiy,phiz)
!
do j=1,jmax
  do i=1,imax
    z        = 1./(b0(1)+xyst(i,j))
    d(i,j,1) = c0(1)*z
    phiz(i,j,1) = phiz(i,j,1)*z
  enddo
enddo
do k=2,kmax-1
   do j=1,jmax
     do i=1,imax
       bb       = b0(k)+xyst(i,j)
       z        = 1./(bb-a0(k)*d(i,j,k-1))
       d(i,j,k) = c0(k)*z
       phiz(i,j,k) = (phiz(i,j,k)-a0(k)*phiz(i,j,k-1))*z
     enddo
  enddo
enddo
do j=1,jmax
  do i=1,imax
    bb       = b0(kmax)+xyst(i,j)
    z        = bb-a0(kmax)*d(i,j,kmax-1)
    if(z.ne.0.) then
      phiz(i,j,kmax) = (phiz(i,j,kmax)-a0(kmax)*phiz(i,j,kmax-1))/z
    else
      phiz(i,j,kmax) =0.
    endif
  enddo
enddo

do k=kmax-1,1,-1
  do j=1,jmax
    do i=1,imax
      phiz(i,j,k) = phiz(i,j,k)-d(i,j,k)*phiz(i,j,k+1)
    enddo
  enddo
enddo

call transpose_z_to_y(phiz,phiy)
do k=1,ysize(3)
  do i=1,ysize(1)
    call vrfftb(1,jtot,phiy(i,1:jtot,k),dj,1,sj)
  enddo
enddo
!
call transpose_y_to_x(phiy,phix)
!
do k=1,xsize(3)
  do j=1,xsize(2)
    call vrfftb(1,itot,phix(1:itot,j,k),di,1,si)
  enddo
enddo
!
call transpose_x_to_y(phix,phiy)
call transpose_y_to_z(phiy,phiz)
!
return
end subroutine solver2ddd





end module mod_phi_im
