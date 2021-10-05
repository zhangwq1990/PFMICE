module mod_loadflds
use decomp_2d
use decomp_2d_io
use mod_cons
use mod_global
use mod_var
use mod_restart
implicit none
private
public loadflds
contains
!
subroutine loadflds(in,nr)
implicit none
integer :: in,nr
integer :: fh
integer(kind=MPI_OFFSET_KIND) :: filesize,disp
character(len=20) :: istepchar
real, dimension(3) :: fldinfo
real, dimension(imax,jmax,kmax) :: temp
integer:: i,j


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this is to read in data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (in.eq.0) then
  write(istepchar,'(I0)') nr

! then check if the restart file is there!
  call ck_restart(trim(istepchar))


  call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//'fld'//trim(istepchar), &
       MPI_MODE_RDONLY, MPI_INFO_NULL,fh, error)
  disp = 0_MPI_OFFSET_KIND
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! unew  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call decomp_2d_read_var(fh,disp,3,temp)
  unew(1:imax,1:jmax,1:kmax) = temp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! vnew  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call decomp_2d_read_var(fh,disp,3,temp)
  vnew(1:imax,1:jmax,1:kmax) = temp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! wnew  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call decomp_2d_read_var(fh,disp,3,temp)
  wnew(1:imax,1:jmax,1:kmax) = temp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! pnew  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call decomp_2d_read_var(fh,disp,3,temp)
  pnew(1:imax,1:jmax,1:kmax) = temp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! pold  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call decomp_2d_read_var(fh,disp,3,temp)
  pold(1:imax,1:jmax,1:kmax) = temp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Tnew  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call decomp_2d_read_var(fh,disp,3,temp)
  Tnew(1:imax,1:jmax,1:kmax) = temp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Phi  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call decomp_2d_read_var(fh,disp,3,temp)
  PFM_phi(1:imax,1:jmax,1:kmax) = temp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Phi boundary  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call decomp_2d_read_var(fh,disp,3,temp)
  do i=1,imax
   do j=1,jmax
    PFM_phi(i,j,0)  = temp(i,j,1)
    PFM_phi(i,j,k1) = temp(i,j,2)     
   enddo
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! C  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call decomp_2d_read_var(fh,disp,3,temp)
  PFM_c(1:imax,1:jmax,1:kmax) = temp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! C boundary  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call decomp_2d_read_var(fh,disp,3,temp)
  do i=1,imax
   do j=1,jmax
    PFM_c(i,j,0)  = temp(i,j,1)
    PFM_c(i,j,k1) = temp(i,j,2)     
   enddo
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Phi_c  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call decomp_2d_read_var(fh,disp,3,temp)
  Phi_c(1:imax,1:jmax,1:kmax) = temp

! the other information
  call decomp_2d_read_scalar(fh,disp,3,fldinfo)
  time = fldinfo(1)
  nr = int(fldinfo(2))
  dt = fldinfo(3)
  call MPI_FILE_CLOSE(fh,error)
endif
!


!this is to write files

if (in.eq.1) then
  !write(istepchar,'(i8.8)') nr
  write(istepchar,'(I0)') nr
  fldinfo = (/time,1.*nr,dt/)
  call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//'fld'//trim(istepchar), &
       MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
  filesize = 0_MPI_OFFSET_KIND
  call MPI_FILE_SET_SIZE(fh,filesize,error)  ! guarantee overwriting
  disp = 0_MPI_OFFSET_KIND
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! unew  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  temp = unew(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! vnew  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  temp = vnew(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! wnew  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  temp = wnew(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! pnew  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  temp = pnew(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! pold  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  temp = pold(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Tew  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  temp = Tnew(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Phi  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  temp = PFM_phi(1:imax,1:jmax,1:kmax) 
  call decomp_2d_write_var(fh,disp,3,temp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Phi bou  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  temp = 0.
  do i=1,imax
   do j=1,jmax
    temp(i,j,1) = PFM_phi(i,j,0)  
    temp(i,j,2) = PFM_phi(i,j,k1) 
   enddo
  enddo
  call decomp_2d_write_var(fh,disp,3,temp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! C  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  temp = PFM_c(1:imax,1:jmax,1:kmax) 
  call decomp_2d_write_var(fh,disp,3,temp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! C bou  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  temp = 0.
  do i=1,imax
   do j=1,jmax
    temp(i,j,1) = PFM_c(i,j,0)  
    temp(i,j,2) = PFM_c(i,j,k1) 
   enddo
  enddo
  call decomp_2d_write_var(fh,disp,3,temp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Phi_c  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  temp = Phi_c(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)


! write in other information
  call decomp_2d_write_scalar(fh,disp,3,fldinfo)
  call MPI_FILE_CLOSE(fh,error)
endif
!
return
end subroutine loadflds
!
end module mod_loadflds
