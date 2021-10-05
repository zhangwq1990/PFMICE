module mod_map
use mod_cons
use mod_var
use mod_global
use mod_bou
implicit none
private
!public check_map_c,mapc,check_map_phi, mapphi
public check_map_phi, mapphi,check_map_c,mapc
contains
!




subroutine check_map_c
implicit none
integer i,j,k

! check the range of c to see if we need to use mapping
use_c_map=0.
phi_c_wb=0.  !give them zero value
phi_c_diff=0.

outer2: do k=1,kmax
   do j=1,jmax
      do i=1,imax
          if (Phi_c(i,j,k)<0.0 .or. Phi_c(i,j,k)>1.0) then
             use_c_map=1
             exit outer2
          endif
      enddo
   enddo
enddo outer2

! then begin to use mapping
   do k=1,kmax
      do j=1,jmax
         do i=1,imax
            Phi_c_bstar(i,j,k)=PFM_phi(i,j,k)*PFM_c(i,j,k)
            phi_c_wb = phi_c_wb + 1.0-(2.0*Phi_c_bstar(i,j,k)-1.0)**2
            phi_c_diff = phi_c_diff + (Phi_c(i,j,k)-Phi_c_bstar(i,j,k))
          enddo
      enddo
   enddo
!endif

return
end subroutine check_map_c
!
!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! to force phi_c back to the range of [-1,1]
subroutine mapc
implicit none
integer i,j,k
real coefficient

coefficient=Phi_c_diff/Phi_c_wb

do k=1,kmax
   do j=1,jmax
      do i=1,imax
         Phi_c_b(i,j,k)=Phi_c_bstar(i,j,k) + (1.0-(2.0*Phi_c_bstar(i,j,k)-1.0)**2)*coefficient
      enddo
   enddo
enddo

return

end subroutine mapc





subroutine check_map_phi
implicit none
integer i,j,k

! check the range of c to see if we need to use mapping
use_phi_map=0.
phi_wb=0.  !give them zero value
phi_diff=0.

outer1: do k=1,kmax
   do j=1,jmax
      do i=1,imax
          if (PFM_phi(i,j,k)<0.0 .or. PFM_phi(i,j,k)>1.0) then
             use_phi_map=1
             exit outer1
          endif
      enddo
   enddo
enddo outer1

! then begin to use mapping
   do k=1,kmax
      do j=1,jmax
         do i=1,imax
            ! to get the phi_bstar value
            phi_bstar(i,j,k)=MIN(MAX(PFM_phi(i,j,k), 0.0),1.0)
            phi_wb = phi_wb + 1.0-(2.0*phi_bstar(i,j,k)-1.0)**2
            phi_diff = phi_diff + (PFM_phi(i,j,k)-phi_bstar(i,j,k))
          enddo
      enddo
   enddo
!endif

return
end subroutine check_map_phi





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!to force c back to the range of [-1,1]
subroutine mapphi
implicit none
integer i,j,k
real coefficient

coefficient=phi_diff/phi_wb

do k=1,kmax
   do j=1,jmax
      do i=1,imax
         phi_b(i,j,k)=phi_bstar(i,j,k) + (1.0-(2.0*phi_bstar(i,j,k)-1.0)**2)*coefficient
      enddo
   enddo
enddo


return
end subroutine mapphi



end module mod_map
