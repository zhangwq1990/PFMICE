module mod_parallel
use decomp_2d
use mod_cons
use mod_global
implicit none
private
public initmpi
contains
!
subroutine initmpi
implicit none
integer :: coordsleft(1:ndims) ,coordsright(1:ndims)
integer :: coordsfront(1:ndims),coordsback(1:ndims)
integer :: coordsneighbor(1:ndims)
!
call MPI_INIT(error)
!
periods(1) = .true.
periods(2) = .true.
periods(3) = .false.
!
call decomp_2d_init(itot, jtot, ktot, dims(1), dims(2), periods)
myid = nrank
!
coords(1) = (zstart(1)-1)*dims(1)/itot
coords(2) = (zstart(2)-1)*dims(2)/jtot
boundleftmyid  = coords(1)*lx/(1.*dims(1)) ! left  boundary
boundfrontmyid = coords(2)*ly/(1.*dims(2)) ! front boundary
!
if ( jtot .lt. itot) then
  if (myid.eq.0) then
    write(6,*) 'jtot smaller than itot: this is computationally less efficient!'
    write(6,*) 'Interchange the coordinate directions.'
    write(6,*) 'Program aborted...'
  endif
  call mpi_finalize(error)
  stop
endif
if ( mod(itot,dims(1)).ne.0 ) then
  if (myid.eq.0) then
    write(6,*) 'itot not divisable by the number of threads in its direction'
    write(6,*) 'Change the value of itot or dims(1) in "param.f90".'
    write(6,*) 'Program aborted...'
  endif
  call mpi_finalize(error)
  stop
endif
if ( mod(jtot,dims(2)).ne.0 ) then
  if (myid.eq.0) then
    write(6,*) 'jtot not divisable by the number of threads in its direction'
    write(6,*) 'Change the value of jtot or dims(2) in "param.f90".'
    write(6,*) 'Program aborted...'
  endif
  call mpi_finalize(error)
  stop
endif
!
comm_cart = DECOMP_2D_COMM_CART_Z
!
call MPI_CART_SHIFT(comm_cart,0,1,left,right,error)
call MPI_CART_SHIFT(comm_cart,1,1,front,back,error)
call MPI_CART_COORDS(comm_cart,right,ndims,coordsright,error)
call MPI_CART_COORDS(comm_cart,front,ndims,coordsfront,error)
coordsneighbor(1) = coordsright(1)
coordsneighbor(2) = coordsfront(2)
call MPI_CART_RANK(comm_cart,coordsneighbor,rightfront,error)
call MPI_CART_COORDS(comm_cart,back ,ndims,coordsback ,error)
coordsneighbor(1) = coordsright(1)
coordsneighbor(2) = coordsback(2)
call MPI_CART_RANK(comm_cart,coordsneighbor,rightback,error)
call MPI_CART_COORDS(comm_cart,left,ndims,coordsleft,error)
call MPI_CART_COORDS(comm_cart,front,ndims,coordsfront,error)
coordsneighbor(1) = coordsleft(1)
coordsneighbor(2) = coordsfront(2)
call MPI_CART_RANK(comm_cart,coordsneighbor,leftfront,error)
call MPI_CART_COORDS(comm_cart,back ,ndims,coordsback ,error)
coordsneighbor(1) = coordsleft(1)
coordsneighbor(2) = coordsback(2)
call MPI_CART_RANK(comm_cart,coordsneighbor,leftback,error)
!
neighbor(0) = myid
neighbor(1) = right
neighbor(2) = rightfront
neighbor(3) = front
neighbor(4) = leftfront
neighbor(5) = left
neighbor(6) = leftback
neighbor(7) = back
neighbor(8) = rightback
!
!
call MPI_TYPE_VECTOR((j1+1)*(k1+1),1,(i1+1),MPI_REAL8,xhalo,error) 
call MPI_TYPE_COMMIT(xhalo,error)
call MPI_TYPE_VECTOR((k1+1),(i1+1),(i1+1)*(j1+1),MPI_REAL8,yhalo,error) 
call MPI_TYPE_COMMIT(yhalo,error)

call MPI_TYPE_VECTOR((j1+1+4)*(k1+1),1,(i1+1+4),MPI_REAL8,xhalo2,error)
call MPI_TYPE_COMMIT(xhalo2,error)
call MPI_TYPE_VECTOR((k1+1),(i1+1+4),(i1+1+4)*(j1+1+4),MPI_REAL8,yhalo2,error)
call MPI_TYPE_COMMIT(yhalo2,error)

call MPI_TYPE_VECTOR((j1+1+4)*(k1+1+4),1,(i1+1+4),MPI_REAL8,xhalo3,error)
call MPI_TYPE_COMMIT(xhalo3,error)
call MPI_TYPE_VECTOR((k1+1+4),(i1+1+4),(i1+1+4)*(j1+1+4),MPI_REAL8,yhalo3,error)
call MPI_TYPE_COMMIT(yhalo3,error)

!
return
end subroutine initmpi
!
end module mod_parallel
