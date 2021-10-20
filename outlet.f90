module mod_outlet
implicit none
private
public outmass
contains
subroutine outmass(totalmass)
use mpi
use mod_cons
use mod_var
use mod_global, only: myid, error,comm_cart
implicit none
!
real  :: mass
real, intent(inout) :: totalmass
integer :: i,j,k
!
mass=0.
! to compute the mass flow across the outlet, here the top surface is outlet
do j=1,jmax
   do i=1,imax
      mass=mass+wnew(i,j,kmax)*(rhol(i,j,kmax)+rhol(i,j,kmax-1))/2.0*dt/dz
   enddo
enddo

call mpi_allreduce(MPI_IN_PLACE,mass,1,mpi_real8,mpi_sum,comm_cart,error)

totalmass=totalmass+mass

if (myid.eq.0) write(6,*) "out mass is", mass, "total mass is", totalmass
return
end subroutine outmass
!
end module mod_outlet
