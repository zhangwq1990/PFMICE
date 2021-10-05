module mod_phi_RHS
use mod_cons
use mod_var
implicit none
private
public phi_rhs,c_rhs
contains
subroutine phi_rhs(BigQ)
real, intent(out), dimension(0:,0:,0:) :: BigQ
integer :: i,j,k
!


do k=1,kmax
  do j=1,jmax
    do i=1,imax
      BigQ(i,j,k) = (phi_b(i,j,k)-PFM_phi(i,j,k))/dt 
    enddo
  enddo
enddo



!
return
end subroutine phi_rhs


subroutine c_rhs(BigQ)
real, intent(out), dimension(0:,0:,0:) :: BigQ
integer :: i,j,k
!

do k=1,kmax
  do j=1,jmax
    do i=1,imax
      BigQ(i,j,k) = (phi_c_b(i,j,k)-Phi_c(i,j,k))/dt 
    enddo
  enddo
enddo



!
return
end subroutine c_rhs


!
end module mod_phi_RHS
