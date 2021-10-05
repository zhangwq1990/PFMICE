module mod_revise
use mod_cons
use mod_var
implicit none
private
public revise,revisec
contains
!
! Corrects the m_phi to make it consistent
!
subroutine revise(BigQ)
real, intent(in), dimension(0:i1,0:j1,0:k1) :: BigQ
integer :: i,j,k
!


do k=1,kmax
  do j=1,jmax
    do i=1,imax
       !m_x_phi(i,j,k)=m_x_phi(i,j,k)-(1.0-(phi_b(i+1,j,k)+phi_b(i,j,k)-1.0)**2)*&
       !                              (BigQ(i+1,j,k)-BigQ(i,j,k))/dx
       !m_y_phi(i,j,k)=m_y_phi(i,j,k)-(1.0-(phi_b(i,j+1,k)+phi_b(i,j,k)-1.0)**2)*&
       !                              (BigQ(i,j+1,k)-BigQ(i,j,k))/dy
       !m_z_phi(i,j,k)=m_z_phi(i,j,k)-(1.0-(phi_b(i,j,k+1)+phi_b(i,j,k)-1.0)**2)*&
       !                              (BigQ(i,j,k+1)-BigQ(i,j,k))/dz

       m_x_phi(i,j,k)=m_x_phi(i,j,k)-(BigQ(i+1,j,k)-BigQ(i,j,k))/dx
       m_y_phi(i,j,k)=m_y_phi(i,j,k)-(BigQ(i,j+1,k)-BigQ(i,j,k))/dy
       m_z_phi(i,j,k)=m_z_phi(i,j,k)-(BigQ(i,j,k+1)-BigQ(i,j,k))/dz
    enddo
  enddo
enddo
!
return
end subroutine revise
!


subroutine revisec(BigQ)
real, intent(in), dimension(0:i1,0:j1,0:k1) :: BigQ
integer :: i,j,k
!


do k=1,kmax
  do j=1,jmax
    do i=1,imax
       m_x_phi_c(i,j,k)=m_x_phi_c(i,j,k)-(BigQ(i+1,j,k)-BigQ(i,j,k))/dx
       m_y_phi_c(i,j,k)=m_y_phi_c(i,j,k)-(BigQ(i,j+1,k)-BigQ(i,j,k))/dy
       m_z_phi_c(i,j,k)=m_z_phi_c(i,j,k)-(BigQ(i,j,k+1)-BigQ(i,j,k))/dz
    enddo
  enddo
enddo
!
return
end subroutine revisec


end module mod_revise
