module mod_restart
use mod_global
implicit none
private
public ck_restart
contains
subroutine ck_restart(filename)
implicit none
character(len=*), intent(in) :: filename
logical :: exist

!print*, "here is the name!", filename
inquire(file='data/fld'//filename, exist=exist)
if (exist) then
   if (myid .eq. 0) then
          write(6,*) "restart file is ready...."
      endif
   else
      if (myid .eq. 0) then
         write(6,*) "restart file is missing, stop programme..."
      endif
   call decomp_2d_finalize
   call MPI_FINALIZE(error)
endif
                                   

return
end subroutine ck_restart
!
end module mod_restart
