module mod_debug_write
  use mod_var
  use mod_global
  use decomp_2d
  use decomp_2d_io
  use mod_zredistribute
  implicit none
  private
  public ddebug,todebug
contains
  !


subroutine ddebug(istep,type_plane,pl_num,cfld)
implicit none
integer, intent(in) :: istep,type_plane,pl_num
real, dimension(0:i1,0:j1,0:k1,1:6), intent(in) :: cfld
real    :: c1_11(1:jmax,1:kmax)
real    :: c1_12(1:jmax,1:kmax)
real    :: c1_13(1:jmax,1:kmax)
real    :: c1_22(1:jmax,1:kmax)
real    :: c1_23(1:jmax,1:kmax)
real    :: c1_33(1:jmax,1:kmax)
real    :: c1_11r(0:(jtot+1),0:(kmax+1))
real    :: c1_12r(0:(jtot+1),0:(kmax+1))
real    :: c1_13r(0:(jtot+1),0:(kmax+1))
real    :: c1_22r(0:(jtot+1),0:(kmax+1))
real    :: c1_23r(0:(jtot+1),0:(kmax+1))
real    :: c1_33r(0:(jtot+1),0:(kmax+1))
!.............................................................
real    :: c2_11(1:imax,1:kmax)
real    :: c2_12(1:imax,1:kmax)
real    :: c2_13(1:imax,1:kmax)
real    :: c2_22(1:imax,1:kmax)
real    :: c2_23(1:imax,1:kmax)
real    :: c2_33(1:imax,1:kmax)
real    :: c2_11r(0:(itot+1),0:(kmax+1))
real    :: c2_12r(0:(itot+1),0:(kmax+1))
real    :: c2_13r(0:(itot+1),0:(kmax+1))
real    :: c2_22r(0:(itot+1),0:(kmax+1))
real    :: c2_23r(0:(itot+1),0:(kmax+1))
real    :: c2_33r(0:(itot+1),0:(kmax+1))
!..............................................................
real    :: c3_11r(0:(itot+1),0:(jtot+1))
real    :: c3_12r(0:(itot+1),0:(jtot+1))
real    :: c3_13r(0:(itot+1),0:(jtot+1))
real    :: c3_22r(0:(itot+1),0:(jtot+1))
real    :: c3_23r(0:(itot+1),0:(jtot+1))
real    :: c3_33r(0:(itot+1),0:(jtot+1))
integer :: ksol,nprocs
parameter(nprocs=dims(1)*dims(2),ksol=kmax/nprocs)
real ::   w3(0:i1,0:j1,0:k1)
real ::   c11new3(itot,jtot,ksol)
real ::   c12new3(itot,jtot,ksol)
real ::   c13new3(itot,jtot,ksol)
real ::   c22new3(itot,jtot,ksol)
real ::   c23new3(itot,jtot,ksol)
real ::   c33new3(itot,jtot,ksol)
!..............................................................
integer :: l,q,rank,sub_rank,sub_coords(ndims)
integer :: npoints
integer :: i,j,k,im,jm,km,it,jt,iloc,jloc
 character*3 number
 character*7 filenumber

if (type_plane.eq.1) then

 !!!!!!!!!!!!!!!!YZ Plane
      i = pl_num
      q=-1
  114 q=q+1
      if ( (1+q*imax) .lt. i ) go to 114
      q=q-1
      iloc = i - q*imax
      sub_coords(1) = q
      sub_coords(2) = 0
      call MPI_CART_RANK(comm_cart,sub_coords,rank,error) !rank of process to which data has to be sent
      do l=1,dims(2)-1
        sub_coords(1) = q
        sub_coords(2) = l
        call MPI_CART_RANK(comm_cart,sub_coords,sub_rank,error)
        if (myid .eq. sub_rank) then
          do k=1,kmax
            do j=1,jmax
              c1_11(j,k)=cfld(iloc,j,k,1)
              c1_12(j,k)=cfld(iloc,j,k,2)
              c1_13(j,k)=cfld(iloc,j,k,3)
              c1_22(j,k)=cfld(iloc,j,k,4)
              c1_23(j,k)=cfld(iloc,j,k,5)
              c1_33(j,k)=cfld(iloc,j,k,6)
            enddo
          enddo
          call MPI_SSEND(c1_11,jmax*kmax,MPI_REAL8,rank,12,comm_cart,error)
          call MPI_SSEND(c1_12,jmax*kmax,MPI_REAL8,rank,13,comm_cart,error)
          call MPI_SSEND(c1_13,jmax*kmax,MPI_REAL8,rank,14,comm_cart,error)
          call MPI_SSEND(c1_22,jmax*kmax,MPI_REAL8,rank,15,comm_cart,error)
          call MPI_SSEND(c1_23,jmax*kmax,MPI_REAL8,rank,16,comm_cart,error)
          call MPI_SSEND(c1_33,jmax*kmax,MPI_REAL8,rank,17,comm_cart,error)
        endif
        if (myid .eq. rank) then
          call MPI_RECV(c1_11,jmax*kmax,MPI_REAL8,sub_rank,12, &
                        comm_cart,status,error)
          call MPI_RECV(c1_12,jmax*kmax,MPI_REAL8,sub_rank,13, &
                        comm_cart,status,error)
          call MPI_RECV(c1_13,jmax*kmax,MPI_REAL8,sub_rank,14, &
                        comm_cart,status,error)
          call MPI_RECV(c1_22,jmax*kmax,MPI_REAL8,sub_rank,15, &
                        comm_cart,status,error)
          call MPI_RECV(c1_23,jmax*kmax,MPI_REAL8,sub_rank,16, &
                        comm_cart,status,error)
          call MPI_RECV(c1_33,jmax*kmax,MPI_REAL8,sub_rank,17, &
                        comm_cart,status,error)

          do k=1,kmax
            do j=1,jmax
              c1_11r(j+l*jmax,k)=c1_11(j,k)
              c1_12r(j+l*jmax,k)=c1_12(j,k)
              c1_13r(j+l*jmax,k)=c1_13(j,k)
              c1_22r(j+l*jmax,k)=c1_22(j,k)
              c1_23r(j+l*jmax,k)=c1_23(j,k)
              c1_33r(j+l*jmax,k)=c1_33(j,k)
            enddo
          enddo
        endif
      enddo
      if (myid .eq. rank) then
        do k=1,kmax
          do j=1,jmax
            c1_11r(j,k)=cfld(iloc,j,k,1)
            c1_12r(j,k)=cfld(iloc,j,k,2)
            c1_13r(j,k)=cfld(iloc,j,k,3)
            c1_22r(j,k)=cfld(iloc,j,k,4)
            c1_23r(j,k)=cfld(iloc,j,k,5)
            c1_33r(j,k)=cfld(iloc,j,k,6)
          enddo
        enddo


        call todebug(c1_11r,c1_12r,c1_13r,c1_22r,c1_23r,c1_33r,istep,1,pl_num)


      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

elseif (type_plane.eq.2) then

!!!!!XZ Plane
      j = pl_num
      q=-1
  112 q=q+1
      if ( (1+q*jmax) .lt. j ) go to 112
      q=q-1
      jloc = j-q*jmax
      sub_coords(1) = 0
      sub_coords(2) = q
      call MPI_CART_RANK(comm_cart,sub_coords,rank,error) !rank of process to which data has to be sent
      do l=1,dims(1)-1
        sub_coords(1) = l
        sub_coords(2) = q
        call MPI_CART_RANK(comm_cart,sub_coords,sub_rank,error)
        if (myid .eq. sub_rank) then
          do k=1,kmax
            do i=1,imax
              c2_11(i,k)=cfld(i,jloc,k,1)
              c2_12(i,k)=cfld(i,jloc,k,2)
              c2_13(i,k)=cfld(i,jloc,k,3)
              c2_22(i,k)=cfld(i,jloc,k,4)
              c2_23(i,k)=cfld(i,jloc,k,5)
              c2_33(i,k)=cfld(i,jloc,k,6)
            enddo
          enddo
          call MPI_SSEND(c2_11,imax*kmax,MPI_REAL8,rank,18,comm_cart,error)
          call MPI_SSEND(c2_12,imax*kmax,MPI_REAL8,rank,19,comm_cart,error)
          call MPI_SSEND(c2_13,imax*kmax,MPI_REAL8,rank,20,comm_cart,error)
          call MPI_SSEND(c2_22,imax*kmax,MPI_REAL8,rank,21,comm_cart,error)
          call MPI_SSEND(c2_23,imax*kmax,MPI_REAL8,rank,22,comm_cart,error)
          call MPI_SSEND(c2_33,imax*kmax,MPI_REAL8,rank,23,comm_cart,error)
        endif
        if (myid .eq. rank) then
          call MPI_RECV(c2_11,imax*kmax,MPI_REAL8,sub_rank,18, &
                        comm_cart,status,error)
          call MPI_RECV(c2_12,imax*kmax,MPI_REAL8,sub_rank,19, &
                        comm_cart,status,error)
          call MPI_RECV(c2_13,imax*kmax,MPI_REAL8,sub_rank,20, &
                        comm_cart,status,error)
          call MPI_RECV(c2_22,imax*kmax,MPI_REAL8,sub_rank,21, &
                        comm_cart,status,error)
          call MPI_RECV(c2_23,imax*kmax,MPI_REAL8,sub_rank,22, &
                        comm_cart,status,error)
          call MPI_RECV(c2_33,imax*kmax,MPI_REAL8,sub_rank,23, &
                        comm_cart,status,error)
          do k=1,kmax
            do i=1,imax
              c2_11r(i+l*imax,k)=c2_11(i,k)
              c2_12r(i+l*imax,k)=c2_12(i,k)
              c2_13r(i+l*imax,k)=c2_13(i,k)
              c2_22r(i+l*imax,k)=c2_22(i,k)
              c2_23r(i+l*imax,k)=c2_23(i,k)
              c2_33r(i+l*imax,k)=c2_33(i,k)
            enddo
          enddo
        endif
      enddo
      if (myid .eq. rank) then
        do k=1,kmax
          do i=1,imax
            c2_11r(i,k)=cfld(i,jloc,k,1)
            c2_12r(i,k)=cfld(i,jloc,k,2)
            c2_13r(i,k)=cfld(i,jloc,k,3)
            c2_22r(i,k)=cfld(i,jloc,k,4)
            c2_23r(i,k)=cfld(i,jloc,k,5)
            c2_33r(i,k)=cfld(i,jloc,k,6)
          enddo
        enddo
        call todebug(c2_11r,c2_12r,c2_13r,c2_22r,c2_23r,c2_33r,istep,2,pl_num)

      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
elseif (type_plane.eq.3) then
!!!!!!XY
      call zredistribute(cfld(:,:,:,1),c11new3,0)
      call zredistribute(cfld(:,:,:,2),c12new3,0)
      call zredistribute(cfld(:,:,:,3),c13new3,0)
      call zredistribute(cfld(:,:,:,4),c22new3,0)
      call zredistribute(cfld(:,:,:,5),c23new3,0)
      call zredistribute(cfld(:,:,:,6),c33new3,0)
      k = pl_num
      l = -1
  199 l = l+1
      if ( (1+l*ksol) .lt. k ) go to 199
      l = l-1
      rank = l
      if (myid .eq. rank) then
        write(number,'(i3.3)') k
        k    = k - rank*ksol
        do j=0,jtot
          do i=0,itot

             c3_11r(i,j) = c11new3(i,j,k) 
             c3_12r(i,j) = c12new3(i,j,k) 
             c3_13r(i,j) = c13new3(i,j,k) 
             c3_22r(i,j) = c22new3(i,j,k) 
             c3_23r(i,j) = c23new3(i,j,k) 
             c3_33r(i,j) = c33new3(i,j,k) 

          enddo
        enddo

        call todebug(c3_11r,c3_12r,c3_13r,c3_22r,c3_23r,c3_33r,istep,3,pl_num)

      endif


endif
!
return
end subroutine ddebug



  subroutine todebug(c11r,c12r,c13r,c22r,c23r,c33r,istep,type2d,origin)


    implicit none

    !----- Inputs (global variables) -----!

    integer, intent(in) :: istep
    real, dimension(0:, 0: ), intent(in) :: c11r, c12r, c13r, c22r, c23r, c33r
    integer, intent(in) :: type2d,origin

    !----- Miscellaneous -----!

    integer :: i, j, k, ifich,nout,nin
    integer :: iostatus
    character :: q
    character(len=1), parameter :: newline = char(10)
    character(len=100) :: s_buffer
    real :: coorz,valme

    integer :: output_unit, input_unit
    parameter ( input_unit = 8, output_unit = 9 )

    character(len=8) :: istepchar
    character(len=48) :: filename

    !----- Create file -----!

    ifich = 10
    q = char(34)

    write(istepchar,'(i8.8)') istep
    if (type2d.eq.1) then  
      filename = datadir//'debugYZ_'//istepchar//'.vtk'
    elseif (type2d.eq.2) then
      filename = datadir//'debugXZ_'//istepchar//'.vtk'
    elseif (type2d.eq.3) then
      filename = datadir//'debugXY_'//istepchar//'.vtk'
    endif

    open( unit = ifich , file = filename , form = 'unformatted' , access = 'stream' , &
          action = 'write' , convert = 'BIG_ENDIAN' , iostat = iostatus )

    write(unit = ifich, iostat = iostatus) '# vtk DataFile Version 3.0' // newline
    write(unit = ifich, iostat = iostatus) 'test file' // newline
    write(unit = ifich, iostat = iostatus) 'BINARY' // newline
    write(unit = ifich, iostat = iostatus) newline
    write(unit = ifich, iostat = iostatus) 'DATASET ' // 'STRUCTURED_POINTS' // newline


   
    if (type2d.eq.1) then !!YZ plane 

    !---- Define the domain ----!

      write(s_buffer, FMT = '(A,3I12)', iostat = iostatus) 'DIMENSIONS', 1, jtot, kmax
      write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
      write(s_buffer, FMT = '(A,3E14.6E2)', iostat = iostatus) 'ORIGIN', (dx/2.)*(2*origin-1), dy/2., dz/2.
      write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
      write(s_buffer, FMT = '(A,3E14.6E2)', iostat = iostatus) 'SPACING', dx, dy, dz
      write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline

    !---- Write data ----!

      write(unit = ifich, iostat = iostatus) newline
      write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'POINT_DATA', 1*jtot*kmax
      write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline

      nout = kmax
      nin  = jtot


    elseif (type2d.eq.2) then !! XZ plane

    !---- Define the domain ----!

      write(s_buffer, FMT = '(A,3I12)', iostat = iostatus) 'DIMENSIONS', itot, 1, kmax
      write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
      write(s_buffer, FMT = '(A,3E14.6E2)', iostat = iostatus) 'ORIGIN', dx/2., (dy/2.)*(2*origin-1), dz/2.
      write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
      write(s_buffer, FMT = '(A,3E14.6E2)', iostat = iostatus) 'SPACING', dx, dy, dz
      write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline

    !---- Write data ----!

      write(unit = ifich, iostat = iostatus) newline
      write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'POINT_DATA', itot*1*kmax
      write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline

      nout = kmax
      nin  = itot



    elseif (type2d.eq.3) then !! XY plane


    !---- Define the domain ----!

      write(s_buffer, FMT = '(A,3I12)', iostat = iostatus) 'DIMENSIONS', itot, jtot, 1
      write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
      write(s_buffer, FMT = '(A,3E14.6E2)', iostat = iostatus) 'ORIGIN', dx/2., dy/2., (dz/2.)*(2*origin-1)
      write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
      write(s_buffer, FMT = '(A,3E14.6E2)', iostat = iostatus) 'SPACING', dx, dy, dz
      write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline

    !---- Write data ----!

      write(unit = ifich, iostat = iostatus) newline
      write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'POINT_DATA', itot*jtot*1
      write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline

      nout = jtot
      nin  = itot

    endif

    !---- 1. C11 ----!

    write(s_buffer, FMT = '(A,A,A,I12)', iostat = iostatus) 'SCALARS ', 'debug1', ' double', 1
    write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
    write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

!    do k = 1,ktot
       do j = 1,nout
          do i = 1,nin

             write(unit = ifich, iostat = iostatus) c11r(i,j)

          enddo
       enddo
!    enddo

    !---- 2. C12 ----!

    write(s_buffer, FMT = '(A,A,A,I12)', iostat = iostatus) 'SCALARS ', 'debug2', ' double', 1
    write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
    write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

!    do k = 1,ktot
       do j = 1,nout
          do i = 1,nin

             write(unit = ifich, iostat = iostatus) c12r(i,j)

          enddo
       enddo
!    enddo

    !---- 3. C13 ----!

    write(s_buffer, FMT = '(3A,I12)', iostat = iostatus) 'SCALARS ', 'debug3', ' double', 1
    write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
    write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

!    do k = 1,ktot
       do j = 1,nout
          do i = 1,nin

             write(unit = ifich, iostat = iostatus) c13r(i,j)

          enddo
       enddo
!    enddo

    !---- 4. C22 ----!

    write(s_buffer, FMT = '(A,A,A,I12)', iostat = iostatus) 'SCALARS ', 'debug4', ' double', 1
    write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
    write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

!    do k = 1,ktot
       do j = 1,nout
          do i = 1,nin

             write(unit = ifich, iostat = iostatus) c22r(i,j)

          enddo
       enddo
!    enddo

    !---- 5. C23 ----!

    write(s_buffer, FMT = '(A,A,A,I12)', iostat = iostatus) 'SCALARS ', 'debug5', ' double', 1
    write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
    write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

!    do k = 1,ktot
       do j = 1,nout
          do i = 1,nin

             write(unit = ifich, iostat = iostatus) c23r(i,j)

          enddo
       enddo
!    enddo

    !---- 6. C33 ----!

    write(s_buffer, FMT = '(A,A,A,I12)', iostat = iostatus) 'SCALARS ', 'debug6', ' double', 1
    write(unit = ifich, iostat = iostatus) trim(s_buffer) // newline
    write(unit = ifich, iostat = iostatus) 'LOOKUP_TABLE default' // newline

!    do k = 1,ktot
       do j = 1,nout
          do i = 1,nin

             write(unit = ifich, iostat = iostatus) c33r(i,j)

          enddo
       enddo
!    enddo

    write(unit = ifich, iostat = iostatus) newline
    close(ifich)


    return
  end subroutine todebug



end module mod_debug_write

