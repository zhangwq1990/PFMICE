module mod_bou
use mod_global
implicit none
private

public bounduvw, boundp, boundc,boundT,updthalos, updthalosBig,boundPFM,boundChem,boundloadd,&
       boundBigQ,boundphib
contains
!
subroutine bounduvw(u,v,w,PFM_phiB,rholB,vislB)
use mod_cons
implicit none
integer :: i,j,k

real, dimension(0:i1,0:j1,0:k1), intent(inout) :: u,v,w
real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in):: PFM_phiB
real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(in):: rholB,vislB
real:: DphiDx,DphiDy,DphiDz,Phi_grad_abs
real ::  first_term ,second_term_a,second_term_b,second_term
real :: mu_bound,dvdz_wall
real::g_prime_phi
real :: phi_wall 
real:: phi_jp, phi_jm,phi_kp,phi_km

select case(top)
case ('wall')

  !this is the wall boundary
  do j=0,j1
    do i=0,i1
      u(i,j,0)    = -u(i,j,1)
      w(i,j,0)    =  0.0
      u(i,j,k1)   = -u(i,j,kmax)
      w(i,j,kmax) =  0.0
      w(i,j,k1)   =  w(i,j,kmax-1)
      select case (PhaseField)
        case ('Droplet')
        v(i,j,k1)   = -v(i,j,kmax)
     end select
    enddo
  enddo


case ('outlet')

  !here I use outlet boundary condition
  do j=0,j1
    do i=0,i1
      u(i,j,0)    = -u(i,j,1)
      w(i,j,0)    =  0.0
      u(i,j,k1)   =  u(i,j,kmax)
      w(i,j,k1)   =  w(i,j,kmax)
      v(i,j,k1)   =  v(i,j,kmax)
    enddo
  enddo

end select


    call updthalos(u,1)
    call updthalos(v,1)
    call updthalos(w,1)

    call updthalos(u,2)
    call updthalos(v,2)
    call updthalos(w,2)

!bottom wall
    do j=0,j1
      do i=0,i1
        v(i,j,0)    = 0
      enddo
    enddo



! communicate data in x direction (periodic b.c.'s incorporated)

call updthalos(u,1)
call updthalos(v,1)
call updthalos(w,1)

! communicate data in y direction (periodic b.c.'s incorporated)

call updthalos(u,2)
call updthalos(v,2)
call updthalos(w,2)
return
end subroutine bounduvw
!
subroutine boundp(p)
use mod_cons
implicit none
integer :: i,j
real, dimension(0:,0:,0:) :: p
!

select case(top)
case ('wall')
     do j=0,j1
       do i=0,i1
         p(i,j,0) = p(i,j,1)     ! Newmann (consistent with no/free-slip)
         p(i,j,k1) = p(i,j,kmax) ! Newmann (consistent with no/free-slip)
       enddo
     enddo
case ('outlet')
     do j=0,j1
       do i=0,i1
         p(i,j,0) = p(i,j,1)     ! Newmann (consistent with no/free-slip)
         p(i,j,k1) = 0.0 ! top wall
       enddo
     enddo
end select


call updthalos(p,1)
!
call updthalos(p,2)
!
return
end subroutine boundp





subroutine boundBigQ(BigQ)
use mod_cons
implicit none
integer :: i,j
real, dimension(0:,0:,0:) :: BigQ
!

do j=0,j1
  do i=0,i1
    BigQ(i,j,0) = BigQ(i,j,1)     ! Newmann (consistent with no/free-slip)
    BigQ(i,j,k1) = BigQ(i,j,kmax) ! Newmann (consistent with no/free-slip)
  enddo
enddo


call updthalos(BigQ,1)
!
call updthalos(BigQ,2)
!
return
end subroutine boundBigQ








subroutine boundc(c)
use mod_cons
implicit none
integer :: i,j
real, dimension(-2:i1+2,-2:j1+2,-2:k1+2),intent(inout)  :: c
!
do j=-2,j1+2
  do i=-2,i1+2
    ! I am not sure about the wall condition
    c(i,j,0) = c(i,j,1)     ! Newmann (consistent with no/free-slip)
    c(i,j,k1) = c(i,j,kmax) ! Newmann (consistent with no/free-slip)
  enddo
enddo

do j=0,j1
  do i=0,i1
    c(i,j,-1)    = c(i,j,0)
    c(i,j,-2)    = c(i,j,0)
    c(i,j,k1+1)  = c(i,j,k1)
    c(i,j,k1+2)  = c(i,j,k1)
   enddo
enddo



!
call updthalosBig(c,1)
!
call updthalosBig(c,2)
!
return
end subroutine boundc



subroutine boundT(T)
use mod_cons
implicit none
integer :: i,j
real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(inout)  :: T
!
do j=-2,j1+2
  do i=-2,i1+2
    T(i,j,0) = 270     ! fix value
    T(i,j,k1) = T(i,j,kmax) ! Newmann, adiabetic
  enddo
enddo

!do j=0,j1
!  do i=0,i1
do j=-2,j1+2
  do i=-2,i1+2
    T(i,j,-1)    = T(i,j,0)
    T(i,j,-2)    = T(i,j,0)
    T(i,j,k1+1)  = T(i,j,k1)
    T(i,j,k1+2)  = T(i,j,k1)
   enddo
enddo

call updthalosBig(T,1)
call updthalosBig(T,2)
!
return
!return
end subroutine boundT



subroutine boundPFM(PFM_phiB,dPFM_boundd,dPFM_boundd_top,u,v,w)
use mod_cons
use mod_chem
implicit none
integer :: i,j,k
real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(inout) :: PFM_phiB
real, dimension(0:i1,0:j1,0:k1), intent(in) :: u,v,w
real, dimension(0:i1,0:j1), intent(out) ::dPFM_boundd,dPFM_boundd_top
real:: first_term,second_term,third_term,DphiDx,DphiDy,DphiDz,Phi_grad_abs
real:: g_prime_phi,v_wall
real::phi_jp,phi_jm,phi_wall
real:: DphiDzii
real::PFM_mu_f_1
select case (PhaseField)
 case ('Droplet')
   do j=-2,j1+2
     do i=-2,i1+2
       PFM_phiB(i,j,k1:k1+2) = PFM_phiB(i,j,kmax) 
     enddo
   enddo
   dPFM_boundd_top(:,:) = 0.
end select
!Bottom wall
   do j=1,jmax
    do i=1,imax
       PFM_phiB(i,j,0)=PFM_phiB(i,j,1)
     enddo
   enddo

   do j=0,j1
     do i=0,i1
       PFM_phiB(i,j,-1)    = PFM_phiB(i,j,0)
       PFM_phiB(i,j,-2)    = PFM_phiB(i,j,0)
       PFM_phiB(i,j,k1+1)  = PFM_phiB(i,j,k1)
       PFM_phiB(i,j,k1+2)  = PFM_phiB(i,j,k1)
      enddo
   enddo
call updthalosBig(PFM_phiB,1)
call updthalosBig(PFM_phiB,2)
!
return
return
end subroutine boundPFM







subroutine boundphib(phi_b)
use mod_cons
implicit none
integer :: i,j,k
real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(inout) :: phi_b
real:: g_prime_phi
real::phi_wall
!top surface
   do j=-2,j1+2
     do i=-2,i1+2
       phi_b(i,j,k1:k1+2) = phi_b(i,j,kmax) 
     enddo
   enddo
!Bottom wall
   do j=1,jmax
    do i=1,imax
      phi_wall = 0.5*(phi_b(i,j,0)+phi_b(i,j,1))
      g_prime_phi = 1.5*phi_wall*(1-phi_wall) 
      phi_b(i,j,0) = phi_b(i,j,1)+ (dz/PFM_l)*g_prime_phi*cos(PFM_thetta)*(2.*sqrt(2.)/3)
     enddo
   enddo

   do j=0,j1
     do i=0,i1
       phi_b(i,j,-1)    = phi_b(i,j,0)
       phi_b(i,j,-2)    = phi_b(i,j,0)
       phi_b(i,j,k1+1)  = phi_b(i,j,k1)
       phi_b(i,j,k1+2)  = phi_b(i,j,k1)
      enddo
   enddo
call updthalosBig(phi_b,1)
call updthalosBig(phi_b,2)
!
return
return
end subroutine boundphib






subroutine boundChem(chem_potB)
use mod_cons
implicit none
real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(inout):: chem_potB
integer:: i,j,k
do j=-2,j1+2
  do i=-2,i1+2
    chem_PotB(i,j,0)    = Chem_potB(i,j,1)
    chem_PotB(i,j,-1)   = Chem_potB(i,j,1)
    chem_PotB(i,j,-2)   = Chem_potB(i,j,1)
    chem_PotB(i,j,k1+2) = Chem_potB(i,j,kmax)
    chem_PotB(i,j,k1+1) = Chem_potB(i,j,kmax)
    chem_PotB(i,j,k1)   = Chem_potB(i,j,kmax)
 enddo
enddo
!
call updthalosBig(chem_PotB,1)
call updthalosBig(chem_PotB,2)
end subroutine boundChem

subroutine boundloadd(u,v,w,p,PFM_phiB)
use mod_cons
implicit none
real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(inout) :: PFM_phiB
real, dimension(0:i1,0:j1,0:k1), intent(inout) :: u,v,w,p
integer::i,j

 do j=0,j1
     do i=0,i1   !I changed here!
       PFM_phiB(i,j,-1)    = PFM_phiB(i,j,0)
       PFM_phiB(i,j,-2)    = PFM_phiB(i,j,0)
       PFM_phiB(i,j,k1+1)  = PFM_phiB(i,j,k1)
       PFM_phiB(i,j,k1+2)  = PFM_phiB(i,j,k1)
     enddo
 enddo
call updthalos(u,1)
call updthalos(v,1)
call updthalos(w,1)
call updthalos(u,2)
call updthalos(v,2)
call updthalos(w,2)
call updthalos(p,1)
call updthalos(p,2)
call updthalosBig(PFM_phiB,1)
call updthalosBig(PFM_phiB,2)



end subroutine boundloadd



subroutine updthalos(var,dir)
use mpi
use mod_cons
use mod_global
implicit none
real, dimension(0:,0:,0:), intent(inout) :: var
integer, intent(in) :: dir
integer :: requests(4), statuses(MPI_STATUS_SIZE,4)
!
!
select case(dir)
case(1) ! x direction
call MPI_SENDRECV(var(1,0,0),1,xhalo,left,0,   &
                    var(i1,0,0),1,xhalo,right,0, &
                    comm_cart,status,error)
call MPI_SENDRECV(var(imax,0,0),1,xhalo,right,0, &
                    var(0,0,0),1,xhalo,left,0,     &
                    comm_cart,status,error)
case(2) ! y direction
  call MPI_SENDRECV(var(0,1,0),1,yhalo,front,0, &       !jmax+1=1
                    var(0,j1,0),1,yhalo,back,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(0,jmax,0),1,yhalo,back,0, &     !jmax=0
                    var(0,0,0),1,yhalo,front,0,   &
                    comm_cart,status,error)
end select
!
return
end subroutine updthalos
!
!
subroutine updthalosBig(var,dir)
use mpi
use mod_cons
use mod_global
implicit none
real, dimension(-2:i1+2,-2:j1+2,-2:k1+2), intent(inout) :: var
integer, intent(in) :: dir
!
!
select case(dir)
case(1) ! x direction

  call MPI_SENDRECV(var(1,-2,-2),1,xhalo3,left,0,   &
                    var(i1,-2,-2),1,xhalo3,right,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(imax,-2,-2),1,xhalo3,right,0, &
                    var(0,-2,-2),1,xhalo3,left,0,     &
                    comm_cart,status,error)

  call MPI_SENDRECV(var(1+1,-2,-2),1,xhalo3,left,0,   &
                    var(i1+1,-2,-2),1,xhalo3,right,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(imax-1,-2,-2),1,xhalo3,right,0, &
                    var(0-1,-2,-2),1,xhalo3,left,0,     &
                    comm_cart,status,error)

  call MPI_SENDRECV(var(1+2,-2,-2),1,xhalo3,left,0,   &
                    var(i1+2,-2,-2),1,xhalo3,right,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(imax-2,-2,-2),1,xhalo3,right,0, &
                    var(0-2,-2,-2),1,xhalo3,left,0,     &
                    comm_cart,status,error)

case(2) ! y direction

  call MPI_SENDRECV(var(-2,1,-2),1,yhalo3,front,0, &
                    var(-2,j1,-2),1,yhalo3,back,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(-2,jmax,-2),1,yhalo3,back,0, &
                    var(-2,0,-2),1,yhalo3,front,0,   &
                    comm_cart,status,error)

  call MPI_SENDRECV(var(-2,1+1,-2),1,yhalo3,front,0, &
                    var(-2,j1+1,-2),1,yhalo3,back,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(-2,jmax-1,-2),1,yhalo3,back,0, &
                    var(-2,0-1,-2),1,yhalo3,front,0,   &
                    comm_cart,status,error)

  call MPI_SENDRECV(var(-2,1+2,-2),1,yhalo3,front,0, &
                    var(-2,j1+2,-2),1,yhalo3,back,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(-2,jmax-2,-2),1,yhalo3,back,0, &
                    var(-2,0-2,-2),1,yhalo3,front,0,   &
                    comm_cart,status,error)
end select
!
return
end subroutine updthalosBig
!
end module mod_bou
