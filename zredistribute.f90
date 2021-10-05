module mod_zredistribute
use mod_cons
use mod_global
use mpi
implicit none
contains
subroutine zredistribute(p,p1,in)
!
!     in=0 : redistribution in the wall-normal direction (bottom corr. to
!            process 0); p not changed
!     in=1 : redistribution from bottom/top like structure to
!            cartesian structure; p1 not changed
!
implicit none
integer nprocs,ksol
parameter(nprocs=dims(1)*dims(2),ksol=kmax/nprocs)
real p (0:,0:,0:)
real p1(:,:,:)
real,dimension(imax,jmax,ksol) :: bufin,bufout
integer in,ib,jb,kb,sub_rank,sub_coords(2),i,j

integer(MPI_ADDRESS_KIND) ereal8, lb
integer mytype1, mytype2, mytype3, montype1, montype2, montype3

! construct a derived type (mytype3) for handling 
!            p(1:imax,1:jmax,kb:kb+ksol-1)

call MPI_Type_get_extent(MPI_REAL8, lb, ereal8, error)

call MPI_Type_contiguous(imax, MPI_REAL8, mytype1, error)
call MPI_Type_commit(mytype1, error)

call MPI_Type_create_hvector(jmax, 1, size(p,1)*ereal8, mytype1, mytype2, error)
call MPI_Type_commit(mytype2, error)

call MPI_Type_create_hvector(ksol, 1, size(p,1)*size(p,2)*ereal8, mytype2, mytype3, error)
call MPI_Type_commit(mytype3, error)

! construct a derived type (montype3) for handling 
!            p1(ib:ib+imax-1,jb:jb+jmax-1,:)

call MPI_Type_contiguous(imax, MPI_REAL8, montype1, error)
call MPI_Type_commit(montype1, error)

call MPI_Type_create_hvector(jmax, 1, size(p1,1)*ereal8, montype1, montype2, error)
call MPI_Type_commit(montype2, error)

call MPI_Type_create_hvector(ksol, 1, size(p1,1)*size(p1,2)*ereal8, montype2, montype3, error)
call MPI_Type_commit(montype3, error)

! make here your choice between the different methods:

call usealltoallw
! call usemytype
! call orig

call MPI_Type_free(mytype1, error)
call MPI_Type_free(mytype2, error)
call MPI_Type_free(mytype3, error)
call MPI_Type_free(montype1, error)
call MPI_Type_free(montype2, error)
call MPI_Type_free(montype3, error)

return

contains
subroutine usealltoallw

! like usemytype but now using alltoallw 

integer(MPI_ADDRESS_KIND) :: base,base1,address
integer, dimension(dims(1)*dims(2)) :: disp,disp1

integer, dimension(dims(1)*dims(2)) :: counts, types ,types1

! compute displacements for the mpi_alltoallw call

call MPI_Get_address( p(1,1,1), base,error)
call MPI_Get_address(p1(1,1,1),base1,error)

do j = 0, dims(2) - 1
  do i = 0, dims(1) - 1
    call MPI_Cart_rank(comm_cart,(/i,j/),sub_rank,error)
    kb     = (i*dims(2)+j)*ksol + 1
    call MPI_Get_address(p(1,1,kb), address, error)
    disp(sub_rank+1) = address - base
    ib = i*imax + 1
    jb = j*jmax + 1
    call MPI_Get_address(p1(ib,jb,1), address, error)
    disp1(sub_rank+1) = address - base1
  enddo
enddo

counts(:) = 1
types(:)  = mytype3
types1(:) = montype3

if(in .eq. 0) then
  call MPI_Alltoallw(p(1,1,1),counts,disp, types, &
                    p1(1,1,1),counts,disp1,types1,comm_cart,error)
endif

if(in .eq. 1) then
  call MPI_Alltoallw(p1(1,1,1),counts,disp1,types1, &
                      p(1,1,1),counts,disp, types,comm_cart,error)
endif

end subroutine usealltoallw

subroutine usemytype

! like original code, but now using derived types

if (in.eq.0) then
  do j = 0,dims(2)-1
    do i = 0,dims(1)-1
      call MPI_Cart_rank(comm_cart,(/i,j/),sub_rank,error)
      kb     = (i*dims(2)+j)*ksol + 1

      call MPI_Issend(p(1,1,kb),1,mytype3,sub_rank, &
                      0,comm_cart,request,error)

      ib = i*imax + 1
      jb = j*jmax + 1

      call MPI_Recv(p1(ib,jb,lbound(p1,3)),1,montype3,sub_rank, &
                    0,comm_cart,status,error)

      call MPI_Wait(request,status,error)
    enddo
  enddo
endif

if (in.eq.1) then
  do j = 0,dims(2)-1
    do i = 0,dims(1)-1
      call MPI_Cart_rank(comm_cart,(/i,j/),sub_rank,error)
      ib = i*imax + 1
      jb = j*jmax + 1
      call MPI_Issend(p1(ib,jb,lbound(p1,3)),1,montype3,sub_rank, &
                      0,comm_cart,request,error)

      kb = (i*dims(2)+j)*ksol + 1

      call MPI_Recv(p(1,1,kb),1,mytype3,sub_rank, &
                    0,comm_cart,status,error)
      call MPI_Wait(request,status,error)
    enddo
  enddo
endif
end subroutine usemytype

subroutine orig
! original code
if (in.eq.0) then
  do j = 0,dims(2)-1
    do i = 0,dims(1)-1
      sub_coords(1) = i
      sub_coords(2) = j
      call MPI_Cart_rank(comm_cart,sub_coords,sub_rank,error)
      kb     = (i*dims(2)+j)*ksol + 1
      bufout = p(1:imax,1:jmax,kb:kb+ksol-1)
      call MPI_Issend(bufout,imax*jmax*ksol,MPI_REAL8,sub_rank, &
                      0,comm_cart,request,error)
      call MPI_Recv(bufin,imax*jmax*ksol,MPI_REAL8,sub_rank, &
                    0,comm_cart,status,error)
      ib = i*imax + 1
      jb = j*jmax + 1
      p1(ib:ib+imax-1,jb:jb+jmax-1,:)=bufin
      call MPI_Wait(request,status,error)
    enddo
  enddo
endif

if (in.eq.1) then
  do j = 0,dims(2)-1
    do i = 0,dims(1)-1
      sub_coords(1) = i
      sub_coords(2) = j
      call MPI_Cart_rank(comm_cart,sub_coords,sub_rank,error)
      ib = i*imax + 1
      jb = j*jmax + 1
      bufout = p1(ib:ib+imax-1,jb:jb+jmax-1,:)
      call MPI_Issend(bufout,imax*jmax*ksol,MPI_REAL8,sub_rank, &
                      0,comm_cart,request,error)
      call MPI_Recv(bufin,imax*jmax*ksol,MPI_REAL8,sub_rank, &
                    0,comm_cart,status,error)
      kb = (i*dims(2)+j)*ksol + 1
      p(1:imax,1:jmax,kb:kb+ksol-1)=bufin
      call MPI_Wait(request,status,error)
    enddo
  enddo
endif

return
end subroutine orig
end subroutine zredistribute
end module mod_zredistribute
