module mod_triperiodic
use mod_thomas
implicit none
private
public triperiodic
contains
subroutine triperiodic(q,a,b,c,k)
!use decomp_2d

!here default is z direction, but it can be other direction as well

implicit none
real, intent(inout), dimension(1:) :: q
real, intent(in), dimension(1:) :: a
real, intent(in), dimension(1:) :: b
real, intent(in), dimension(1:) :: c
integer,intent(in) :: k
integer :: i
real, allocatable :: q0(:),xx(:),xx1(:),xx2(:)
real :: xk

!these allocate and deallocate part can be modified to make them faster
allocate(q0(k-1))  !this is to solve the second part solution
allocate(xx(k))  !this is the final solution
allocate(xx1(k-1))  !first part solution
allocate(xx2(k-1))  !second part solution

q0=0.
q0(1)=-a(1)
q0(k-1)=-c(k-1)

call thomas(q(1:k-1), a(1:k-1),b(1:k-1),c(1:k-1),k-1)
call thomas(q0(1:k-1),a(1:k-1),b(1:k-1),c(1:k-1),k-1)


xk= ( q(k)-c(k)*q(1)-a(k)*q(k-1) )/( c(k)*q0(1)+a(k)*q0(k-1)+b(k) )

do i=1,k-1
   xx(i)=q(i)+q0(i)*xk
enddo
xx(k)=xk

!give the true value back to q
q=xx

deallocate(xx)
deallocate(q0)
deallocate(xx1)
deallocate(xx2)



!!
return
end subroutine triperiodic
!
end module mod_triperiodic
