module mod_var
use mod_cons
!
real ,dimension(0:i1,0:j1,0:k1) :: unew,vnew,wnew,pnew, &
                                   phat, &
                                   dudtold,dvdtold,dwdtold
real ,dimension(0:i1,0:j1,0:k1) :: delta_u,delta_u_old    !used for icing computation


real ,dimension(0:i1,0:j1,0:k1) :: uold,vold,wold,pold
real ,dimension(0:i1,0:j1,0:k1) :: ustar,vstar,wstar
real ,dimension(0:i1,0:j1,0:k1) :: uold2,vold2,wold2,pold2
real, dimension(0:i1,0:j1,0:k1) :: RHS_x,RHS_y,RHS_z,RHS_x_old,RHS_y_old,RHS_z_old



real ,dimension(0:i1,0:j1,0:k1) :: BigQold, BigQold2,BigQhat


real ,target ,dimension(0:i1,0:j1,0:k1) :: dudt,dvdt,dwdt
real ,dimension(0:i1,0:j1,0:k1,1:6) :: C_Po,debug
real(mytype) :: time,dt
real(mytype) ::  wi(itot+15), wj(jtot+15)
real, dimension(imax,jmax) :: xyrt
real, dimension(kmax) :: tri1_l,tri1_m,tri1_r
real, dimension(kmax) :: tri2_l,tri2_m,tri2_r
real, dimension(1:imax,1:jmax,1:kmax) :: aa,bb,cc    !so you are here, for Tommas
!real :: forcextot,forceytot,forceztot
real :: u_bulk,v_bulk,w_bulk
real :: kill   !used to stop programme

character (len=120) :: variables(6)

integer :: tec_view,cishu, ifile  !to control the tecplot file

!*********************** Phase Field Method*************
real ,dimension(0:i1,0:j1,0:k1) ::surf_tension_x, surf_tension_y,surf_tension_z,adduv
real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: visl,vislold,rhol,rholold,rholold2
real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: PFM_phi,PFM_phi_old,chem_pot,chem_pot_old
real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: alpha  !this is the coefficient used in the momentum equation, convection term
real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: Phi_c,Phi_c_old  !this is used in the icing
real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: PFM_c,PFM_c_old, Tnew,rholcp,rholcp_old,Told,kk,kkold,au
real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: Phi_c_bstar,phi_bstar,phi_b,Phi_c_b
real :: phi_c_wb,phi_c_diff,phi_wb,phi_diff
integer :: use_c_map,use_phi_map
real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: m_x_rho_cp,m_y_rho_cp,m_z_rho_cp
real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: m_x_rho,m_y_rho,m_z_rho
real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: m_x_phi,m_x_phi_c,m_y_phi,m_y_phi_c,m_z_phi,m_z_phi_c
real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: ct_origin,ct_change
!residual0 is delta dot U, 1 is mass conservation, 2 is phi conservation, 3 is phic conservation
!4 is delta dot U - mobility*(rho2-roh3)/rho*sigma
real, dimension(0:i1,0:j1,0:k1) :: residual0,residual1,residual2,residual3,residual4
real, dimension(0:i1,0:j1,0:k1) :: sm
real ,dimension(0:jtot)::PFM_phi_agg
real::PFM_lambda
real::vel_max
real ,dimension(0:i1,0:j1) :: dPFM_bound_old,dPFM_bound_old_top


!*************************  used for debug ***************
!real, dimension(-2:i1+2,-2:j1+2,-2:k1+2) :: debug1,debug2,debug3
real, dimension(0:i1,0:j1,0:k1) :: debug1,debug2,debug3,debug4,debug5,debug6
!******************************************************
end module mod_var
!
module mod_global
use mpi
use decomp_2d
implicit none
integer :: myid,xhalo,yhalo,restartp,rankcw,xhalo2,yhalo2,xhalo3,yhalo3
integer :: comm_cart
!
logical periods(3),reorder
integer error,request,status(MPI_STATUS_SIZE)
integer right,rightfront,front,leftfront,left,leftback,back,rightback
integer, dimension(0:8) :: neighbor
integer, dimension(1:2) :: coords
real :: boundleftmyid,boundfrontmyid
!
end module mod_global
