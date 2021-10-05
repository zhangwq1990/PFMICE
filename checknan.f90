module mod_chknan
implicit none
private
public chk
contains
subroutine chk(var)
use mpi
use mod_cons
use mod_var
use mod_global
implicit none
!
real, intent(in) :: var(0:,0:,0:)
real  :: infinity,dbl_prec_var
integer :: i,j,k
!

! check if infinity or nan
infinity = HUGE(dbl_prec_var)
outer: do k=0,k1
          do j=0,j1
             do i=0,i1
                if((var(i,j,k) > infinity) .or.(var(i,j,k) /= var(i,j,k)) ) then
                  kill=0
                  exit outer
                endif
             enddo
          enddo
       enddo outer


call mpi_allreduce(MPI_IN_PLACE,kill,1,mpi_real8,mpi_min,comm_cart,error)

return
end subroutine chk
!
end module mod_chknan
