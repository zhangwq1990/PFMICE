program PFMICE
use mod_cons
use mod_var
use mod_global
use mod_parallel
use mod_init
use mod_old
use mod_logo
use mod_field
use mod_bou
use mod_loadflds
use mod_conserve
use mod_phi
use mod_phi_im
use mod_phi_RHS
use mod_revise
use mod_c
use mod_cim
use mod_c2
use mod_map
use mod_chknan
use mod_uvwstar
use mod_prefft
use mod_p_RHS
use mod_solver
use mod_correc
use mod_temp
use mod_energy
use mod_vtk_write
use mod_outlet
use mod_debug_write
use mod_chem
implicit none
integer :: begin,nstep,istep
real, dimension(0:i1,0:j1,0:k1) :: dPFMold
real, target, dimension(0:i1,0:j1,0:k1) :: p
real, pointer, dimension(:,:,:) :: ppp
real, target, dimension(0:i1,0:j1,0:k1) :: BigQ
real, pointer, dimension(:,:,:) :: qqq
integer :: i,j,k,save_vtk,output_screen,vtk_start,vtk_end,save_restart,line
real :: res0,res1,res2,res3,res4
real :: time1, time2
real:: mass,totalmass !, sum111,var_total,var_total1,var_total2,var_total3,var_total4,var_total5,var_total6
real:: phicmin,Tmax,Tmin,cmin,cmax,phimax,phimin
integer :: phi_scheme,c_scheme,T_scheme




! control parameters
phi_scheme=1   ! 1 is explicit
c_scheme  =1   ! 1 is explicit
T_scheme  =1   ! 1 is explicit
dt=2.0e-7
begin =0
!nstep =5000
nstep =10
save_vtk=500
!output_screen=500
output_screen=1
save_restart=10000
line=0  !this is to output a file with only 1 point, or more if you want
        ! if line=0, then no output; if it is not zero, then output with this frequency
kill=1  !this means the programme will go on
tec_view=0 !this means write out tecplot file
cishu=0 !tecplot file at first is 0
totalmass=0.

!!!!!!!!!!!the time step to start and end
vtk_start=1
vtk_end=nstep+begin
!!!!!!!!!!!!!!!!!!!!!!

!start program
call initmpi


call logo

if (myid .eq. 0) write(6,*) 'nr steps at beginning simulation = ',begin
if (myid .eq. 0) write(6,*) ' dt = ', dt


! start from 0 or restart from a solution
if(begin.eq.0) then
  call init
  time = 0.
else
  unew(:,:,:) = 0.
  vnew(:,:,:) = 0.
  wnew(:,:,:) = 0.
  pnew(:,:,:) = 0.
  PFM_phi(:,:,:) = 0.
  PFM_c(:,:,:) = 1.

variables(1)='PFM_phi'
variables(2)='PFM_c'
variables(3)='rhol'
variables(4)='Tnew'
variables(5)='phi_c'
variables(6)='residual'


  call loadflds(0,begin)
  call boundloadd(unew,vnew,wnew,pnew,PFM_phi)
endif




!init Poisson solver
call init_tri
call initsolver
ppp =>    p(1:imax,1:jmax,1:kmax) ! uncomment iff solver2d is used
qqq => BigQ(1:imax,1:jmax,1:kmax) ! uncomment iff solver2d is used

! init phi boundary and c boundary, calculate chemical potential
!this is needed at the 1st time step, but not needed at the restart
if (begin .ne. 0) then
!  call boundPFM(PFM_phi,dPFM_bound_old,dPFM_bound_old_top,unew,vnew,wnew)
  call boundChem(Phi_c)
endif

  call boundPFM(PFM_phi,dPFM_bound_old,dPFM_bound_old_top,unew,vnew,wnew)  !I think this is needed
  call boundChem(Phi_c) ! I think this is needed
  call ChemicalPotential
  call boundChem(chem_pot)



  call boundc(PFM_c)
  call boundT(Tnew)


! calculate the physical parameters

call updatefield(0)


call bounduvw(unew,vnew,wnew,PFM_phi,rhol,visl)
call boundp(pnew)


dPFM_bound_old(:,:)=0.0
dPFM_bound_old_top(:,:)=0.0
!alpha=0.
delta_u=0.0
delta_u_old=0.0
ct_change=0.
rholold=rhol
BigQ=0.


! output initial solution
C_Po(:,:,:,1) = PFM_phi(0:i1,0:j1,0:k1)
C_Po(:,:,:,2)=  PFM_c(0:i1,0:j1,0:k1) 
C_Po(:,:,:,3) = rhol(0:i1,0:j1,0:k1)
C_Po(:,:,:,4) = Tnew(0:i1,0:j1,0:k1)
C_Po(:,:,:,5) = pnew(:,:,:)
C_Po(:,:,:,6) = Phi_c(0:i1,0:j1,0:k1)
istep =0
!call post2dme(istep,1,itot/2,unew,vnew,wnew,pnew,C_Po)


!!debug output
!debug1=0.0
!debug2=0.
!debug3=0.
!debug4=0.
!debug5=0.
!debug6=0.

  use_c_map=0  !do not use map_c at the beginning
  use_phi_map=0  !do not use map_c at the beginning

print*, "myid=", myid, coords(1),coords(2)


!*********************************************************
 open(13, file="line", status='replace')
!**********************************************************


call CPU_time(time1)
! main loop below
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do istep = begin+1,begin+nstep  !here I changed the defination of nstep

! to provide a smaller time step at first
!if (istep > 3000) then
!   dt =5.0e-6
!endif

  time = time + dt




!! give old values
!rholold2=rholold
!rholold=rhol
!rholcp_old=rholcp
!PFM_phi_old=PFM_phi
!PFM_c_old=PFM_c
!kkold=kk
!vislold=visl
!Told=Tnew
!uold=unew
!vold=vnew
!wold=wnew
!pold=pnew
!chem_pot_old=chem_pot
!Phi_c_old=Phi_c
!delta_u_old=delta_u

call old



if (phi_scheme==1) then
!explicit phi
  call m_phi
  call updatePfm(dPFMold)


  call boundPFM(PFM_phi,dPFM_bound_old,dPFM_bound_old_top,unew,vnew,wnew)
  !here is the mapping process
  call check_map_phi
  
  call mpi_allreduce(MPI_IN_PLACE,use_phi_map,  1,mpi_integer,mpi_sum,comm_cart,error)   !calculate the demoninator and coefficient of equation50
  if (use_phi_map>0) then
     call mpi_allreduce(MPI_IN_PLACE,phi_wb,  1,mpi_real8,mpi_sum,comm_cart,error)   !calculate the demoninator and coefficient of equation50
     call mpi_allreduce(MPI_IN_PLACE,phi_diff,1,mpi_real8,mpi_sum,comm_cart,error)
     call mapphi   !this is to force phi back to [0,1]


      call phi_RHS(BigQ)
      call solver2d(qqq,tri1_l,tri1_m,tri1_r)
      call boundBigQ(BigQ)
      call revise(BigQ)
      call updthalosBig(m_y_phi,1)
      call updthalosBig(m_y_phi,2)
      call updthalosBig(m_x_phi,1)
      call updthalosBig(m_x_phi,2)
      call updthalosBig(m_z_phi,1)
      call updthalosBig(m_z_phi,2)

   do k=1,kmax
    do j=1,jmax
     do i=1,imax
     PFM_phi(i,j,k)=phi_b(i,j,k)
     enddo
    enddo
   enddo
     call boundPFM(PFM_phi,dPFM_bound_old,dPFM_bound_old_top,unew,vnew,wnew)
  endif
 !!end of the mapping process
  call ChemicalPotential
  call boundChem(chem_pot)


else   ! use implicit scheme

  !implicit
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  call updatephi_im
  call boundPFM(PFM_phi,dPFM_bound_old,dPFM_bound_old_top,unew,vnew,wnew)
  call ChemicalPotential_im
  call boundChem(chem_pot)
  call m_imphi

  call check_map_phi
  
  call mpi_allreduce(MPI_IN_PLACE,use_phi_map,  1,mpi_integer,mpi_sum,comm_cart,error)   !calculate the demoninator and coefficient of equation50
  if (use_phi_map>0) then
     !print*, "aaa"
     call mpi_allreduce(MPI_IN_PLACE,phi_wb,  1,mpi_real8,mpi_sum,comm_cart,error)   !calculate the demoninator and coefficient of equation50
     call mpi_allreduce(MPI_IN_PLACE,phi_diff,1,mpi_real8,mpi_sum,comm_cart,error)
     call mapphi   !this is to force phi back to [0,1]


      call phi_RHS(BigQ)
      call solver2d(qqq,tri1_l,tri1_m,tri1_r)
      call boundBigQ(BigQ)
      call revise(BigQ)
      call updthalosBig(m_y_phi,1)
      call updthalosBig(m_y_phi,2)
      call updthalosBig(m_x_phi,1)
      call updthalosBig(m_x_phi,2)
      call updthalosBig(m_z_phi,1)
      call updthalosBig(m_z_phi,2)


   do k=1,kmax
    do j=1,jmax
     do i=1,imax
        PFM_phi(i,j,k)=phi_b(i,j,k)
     enddo
    enddo
   enddo
   call boundPFM(PFM_phi,dPFM_bound_old,dPFM_bound_old_top,unew,vnew,wnew)
  endif


endif  ! end the ex or im scheme


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! then let us go to C
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call m_phi_c



!!calculate c value
if (c_scheme==1 ) then
!!explicit
  call updatephi_c
  call updatec
!!!implicit
else
  call updatecim
endif

! then let us map phi_c
  call check_map_c
  
  call mpi_allreduce(MPI_IN_PLACE,use_c_map,  1,mpi_integer,mpi_sum,comm_cart,error)   !calculate the demoninator and coefficient of equation50
  if (use_c_map>0) then
     !print*, "aaa"
     call mpi_allreduce(MPI_IN_PLACE,phi_c_wb,  1,mpi_real8,mpi_sum,comm_cart,error)   !calculate the demoninator and coefficient of equation50
     call mpi_allreduce(MPI_IN_PLACE,phi_c_diff,1,mpi_real8,mpi_sum,comm_cart,error)
     call mapc   !this is to force c back to [0,1]


      call c_RHS(BigQ)
      call solver2d(qqq,tri1_l,tri1_m,tri1_r)
      call boundBigQ(BigQ)
      call revisec(BigQ)
      call updthalosBig(m_y_phi_c,1)
      call updthalosBig(m_y_phi_c,2)
      call updthalosBig(m_x_phi_c,1)
      call updthalosBig(m_x_phi_c,2)
      call updthalosBig(m_z_phi_c,1)
      call updthalosBig(m_z_phi_c,2)

   do k=1,kmax
    do j=1,jmax
     do i=1,imax
     Phi_c(i,j,k)=Phi_c_b(i,j,k)
     enddo
    enddo
   enddo
endif
  call boundc(PFM_c)
  call boundChem(Phi_c)
  call updatec_t




call m_rho     !here it is still time step n, flux
call m_rho_cp



! update the physical parameters
! 1 is density, 2 is viscosity, 3 is kk, 4 is rho_cp, 0 is all
! here it is time step n+1, variable
call updatefield(4)    !update roh_cp
call updatefield(1)    !update rho






! then we solv temperature
if (T_scheme==1) then
  call temp
else
  call energy
endif
  call boundT(Tnew)

  call update_delta




! calculate intermediate velocity

  do k=0,k1
     do j=0,j1
        do i=0,i1
          residual1(i,j,k)=(rhol(i,j,k)-rholold(i,j,k))/dt+ & ! use ziyang consistency equation
                           (m_x_rho(i,j,k)-m_x_rho(i-1,j,k))/dx+&
                           (m_y_rho(i,j,k)-m_y_rho(i,j-1,k))/dy+&
                           (m_z_rho(i,j,k)-m_z_rho(i,j,k-1))/dz
        enddo
     enddo
  enddo

  call uvwstar
  call bounduvw(ustar,vstar,wstar,PFM_phi,rhol,visl)





call updatefield(2)   !viscosity
call updatefield(3)   ! kk




! calculate pressure and velocity for n+1
!if (istep.eq.begin+1) then
if (istep.eq.1) then  !I have restart
  pold2=pnew
  pold =pnew
else
  pold2=pold
  pold=pnew
endif

  do k=1,kmax
   do j=1,jmax
    do i=1,imax
      !if (istep.eq.begin+1) then
      if (istep.eq.1) then   !because I write pold into restart now
        phat(i,j,k) = pold(i,j,k)
      else
        phat(i,j,k) = 2.0*pold(i,j,k)-pold2(i,j,k)
      endif
    enddo
   enddo
  enddo
  call boundp(phat)
  call p_RHS(p)
  call solver2d(ppp,tri2_l,tri2_m,tri2_r)
  call boundp(p)
  call correc(p)


  call bounduvw(unew,vnew,wnew,PFM_phi,rhol,visl)
  pnew(:,:,:) = p(:,:,:)
  call boundp(pnew)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  !here, check the residuals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


if (mod(istep,output_screen).eq.0 .or. (istep==1)) then
  call chk(PFM_phi(0:i1,0:j1,0:k1))
  if (kill .eq. 0) then
     call CPU_time(time2)
     print*, "time used is ", time2-time1
     write(6,*) "NaN or Infinity appears in Phi, stop programme"
     call decomp_2d_finalize
     call MPI_FINALIZE(error)
  endif
  call chk(Tnew(0:i1,0:j1,0:k1))
  if (kill .eq. 0) then
     call CPU_time(time2)
     print*, "time used is ", time2-time1
     write(6,*) "NaN or Infinity appears in T, stop programme"
     call decomp_2d_finalize
     call MPI_FINALIZE(error)
  endif

    v_bulk =sum(vnew(1:imax,1:jmax,1:kmax))
    w_bulk =sum(wnew(1:imax,1:jmax,1:kmax))
    mass   =sum(rhol(1:imax,1:jmax,1:kmax))
    phicmin=minval(Phi_c(1:imax,1:jmax,1:kmax))
    Tmax   =maxval(Tnew(1:imax,1:jmax,1:kmax))
    Tmin   =minval(Tnew(1:imax,1:jmax,1:kmax))

    phimax   =maxval(PFM_phi(1:imax,1:jmax,1:kmax))
    phimin   =minval(PFM_phi(1:imax,1:jmax,1:kmax))
    cmax   =maxval(PFM_c(1:imax,1:jmax,1:kmax))
    cmin   =minval(PFM_c(1:imax,1:jmax,1:kmax))

    call mpi_allreduce(MPI_IN_PLACE,v_bulk,1,mpi_real8,mpi_sum,comm_cart,error)
    call mpi_allreduce(MPI_IN_PLACE,w_bulk,1,mpi_real8,mpi_sum,comm_cart,error)
    call mpi_allreduce(MPI_IN_PLACE,mass  ,1,mpi_real8,mpi_sum,comm_cart,error)
    call mpi_allreduce(MPI_IN_PLACE,phicmin,1,mpi_real8,mpi_min,comm_cart,error)
    call mpi_allreduce(MPI_IN_PLACE,Tmax  ,1,mpi_real8,mpi_max,comm_cart,error)
    call mpi_allreduce(MPI_IN_PLACE,Tmin  ,1,mpi_real8,mpi_min,comm_cart,error)

    call mpi_allreduce(MPI_IN_PLACE,phimax  ,1,mpi_real8,mpi_max,comm_cart,error)
    call mpi_allreduce(MPI_IN_PLACE,phimin  ,1,mpi_real8,mpi_min,comm_cart,error)
    call mpi_allreduce(MPI_IN_PLACE,cmax  ,1,mpi_real8,mpi_max,comm_cart,error)
    call mpi_allreduce(MPI_IN_PLACE,cmin  ,1,mpi_real8,mpi_min,comm_cart,error)

    v_bulk=v_bulk/(1.*itot*jtot*ktot)
    w_bulk=w_bulk/(1.*itot*jtot*ktot)

    !calculate the residual
    residual0=0.
    residual1=0.
    residual2=0.
    residual3=0.
    residual4=0.
    res0=0.
    res1=0.
    res2=0.
    res3=0.
    res4=0.
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
            residual0(i,j,k)=(unew(i,j,k)-unew(i-1,j,k))/dx+&
                             (vnew(i,j,k)-vnew(i,j-1,k))/dy+&
                             (wnew(i,j,k)-wnew(i,j,k-1))/dz
            residual1(i,j,k)=(rhol(i,j,k)-rholold(i,j,k))/dt+ & ! use ziyang consistency equation
                             (m_x_rho(i,j,k)-m_x_rho(i-1,j,k))/dx+&
                             (m_y_rho(i,j,k)-m_y_rho(i,j-1,k))/dy+&
                             (m_z_rho(i,j,k)-m_z_rho(i,j,k-1))/dz
            residual2(i,j,k)=(PFM_phi(i,j,k)-PFM_phi_old(i,j,k))/dt+&
                             (m_x_phi(i,j,k)-m_x_phi(i-1,j,k))/dx+&
                             (m_y_phi(i,j,k)-m_y_phi(i,j-1,k))/dy+&
                             (m_z_phi(i,j,k)-m_z_phi(i,j,k-1))/dz-&
                             PFM_phi_old(i,j,k)*((uold(i,j,k)-uold(i-1,j,k))/dx+&
                                                 (vold(i,j,k)-vold(i,j-1,k))/dy+&
                                                 (wold(i,j,k)-wold(i,j,k-1))/dz)
            residual3(i,j,k)=(Phi_c(i,j,k)-Phi_c_old(i,j,k))/dt+&
                             (m_x_phi_c(i,j,k)-m_x_phi_c(i-1,j,k))/dx+&
                             (m_y_phi_c(i,j,k)-m_y_phi_c(i,j-1,k))/dy+&
                             (m_z_phi_c(i,j,k)-m_z_phi_c(i,j,k-1))/dz-&
                             ct_change(i,j,k)
                             !ct_origin(i,j,k)
            !residual4(i,j,k)=rhol(i,j,k)*residual0(i,j,k)-delta_u(i,j,k)*(rho_2-rho_3)
            residual4(i,j,k)=rholold(i,j,k)*( (uold(i,j,k)-uold(i-1,j,k))/dx+&
                                              (vold(i,j,k)-vold(i,j-1,k))/dy+&
                                              (wold(i,j,k)-wold(i,j,k-1))/dz )-&
                                              delta_u_old(i,j,k)*(rho_2-rho_3)
                            
             res0=res0+residual1(i,j,k)
             res1=res1+abs(residual1(i,j,k))
             res2=res2+abs(residual2(i,j,k))
             res3=res3+abs(residual3(i,j,k))
             res4=res4+abs(residual4(i,j,k))
          enddo
       enddo
    enddo
    call boundBigQ(residual1)
    call boundBigQ(residual2)
    call boundBigQ(residual3)
    call mpi_allreduce(MPI_IN_PLACE,res1,1,mpi_real8,mpi_sum,comm_cart,error)
    call mpi_allreduce(MPI_IN_PLACE,res2,1,mpi_real8,mpi_sum,comm_cart,error)
    call mpi_allreduce(MPI_IN_PLACE,res3,1,mpi_real8,mpi_sum,comm_cart,error)
    if (myid.eq.0) write(6,*) 'time = ',time,'istep = ',istep,'mass=   ',mass, 'vbulk = ',v_bulk,'wbulk = ',w_bulk
    if (myid.eq.0) write(6,"(5es24.17)") res0,res1,res2,res3,res4,Tmax,Tmin
    if (myid.eq.0) write(6,"(2es24.17)") phimax,phimin,cmax,cmin
    if (myid.eq.0) write(6,*) "  "
    call CPU_time(time2)
    print*, "time used is ", time2-time1
    !if (myid.eq.0) then
    !    do k=1,kmax
    !       !print*, k,debug5(1,10,k),debug6(1,10,k),PFM_phi(1,10,k),Tnew(i,j,k)
    !       print*, k,debug2(1,10,k),debug3(1,10,k),debug4(1,10,k),debug1(1,10,k),debug5(1,10,k),PFM_c(1,10,k)
    !    enddo
    !endif
call outmass(totalmass)
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! output solution files
!if (istep>=vtk_start .and. istep<=vtk_end) then
!if(myid.eq.0) then
  if ((mod(istep,save_vtk).eq.0).or.(istep.eq.(begin+1))) then
     C_Po(:,:,:,1) = PFM_phi(0:i1,0:j1,0:k1)
     C_Po(:,:,:,2)=  PFM_c(0:i1,0:j1,0:k1) 
     C_Po(:,:,:,3) = rhol(0:i1,0:j1,0:k1)
     C_Po(:,:,:,4) = Tnew(0:i1,0:j1,0:k1)
     C_Po(:,:,:,5) = Phi_c(0:i1,0:j1,0:k1)
     C_Po(:,:,:,6) = residual1(0:i1,0:j1,0:k1)
     call post2dme(istep,1,itot/2,unew,vnew,wnew,pnew,C_Po)
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     call boundBigQ(debug1)
     call boundBigQ(debug2)
     call boundBigQ(debug3)
     call boundBigQ(debug4)
     call boundBigQ(debug5)
     call boundBigQ(debug6)
     debug(:,:,:,1) = debug1(0:i1,0:j1,0:k1)
     debug(:,:,:,2)=  debug2(0:i1,0:j1,0:k1) 
     debug(:,:,:,3) = debug3(0:i1,0:j1,0:k1)
     debug(:,:,:,4) = debug4(0:i1,0:j1,0:k1)
     debug(:,:,:,5) = debug5(0:i1,0:j1,0:k1)
     debug(:,:,:,6) = debug6(0:i1,0:j1,0:k1)
     call ddebug(istep,1,itot/2,debug)

  endif 
 
!endif
! output restart file

  if (mod(istep,save_restart).eq.0) then
    call loadflds(1,istep)
  endif



enddo

if(myid.eq.0) write(6,*) '*** Fim ***'
!
call decomp_2d_finalize
call MPI_FINALIZE(error)
!
close(13)
call CPU_time(time2)
print*, "time used is ", time2-time1
end program 
