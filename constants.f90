module mod_cons
use decomp_2d
implicit none
!££££££££££££££££££££££££££
!file IO                  £
!££££££££££££££££££££££££££
 character(len=5), parameter :: datadir = 'data/'
 character(len=5), parameter :: partdir = 'part/'
!£££££££££££££££££££££££££
! domain constants       £
!£££££££££££££££££££££££££
integer, parameter :: ndims = 2
integer, dimension(ndims), parameter :: dims = (/1,4/)
integer,parameter:: number_of_processors=dims(1)*dims(2)
integer, parameter :: itot =2, jtot = 200, ktot = 200
integer, parameter :: it1 = itot+1, jt1 = jtot+1, kt1 = ktot+1
integer, parameter :: imax = itot/dims(1), jmax = jtot/dims(2),kmax = ktot
integer, parameter :: i1 = imax+1, j1 = jmax+1, k1 = kmax+1
real, parameter :: lx = 0.0008,ly = 0.008,lz = 0.008
real, parameter :: dxi = itot/lx, dyi = jtot/ly, dzi = ktot/lz
real, parameter :: dx = 1./dxi, dy = 1./dyi, dz = 1./dzi
real, parameter :: dx2 = dx*dx, dy2 = dy*dy, dz2 = dz*dz
!£££££££££££££££££££££££££
!math constants          £
!£££££££££££££££££££££££££
real, parameter :: pi = acos(-1.)
!real, parameter :: picon = acos(-1.)
!£££££££££££££££££££££££££
!case specific           £
!£££££££££££££££££££££££££
character(len=3), parameter :: iniu = 'wet' !
logical, parameter :: isfreeslip = .false.
!real, parameter :: bulk_v_sup = 1.
character(len=7), parameter :: PhaseField= 'Droplet' !'Droplet'! 'Couette'
character(len=7), parameter :: top = 'outlet' !'outlet'! 'wall'
!character(len=7), parameter :: top = 'wall' !'outlet'! 'wall'
real, parameter:: droplet_radius= lz/10. 
!££££££££££££££££££££££££
!physical propertities  £
!££££££££££££££££££££££££
real,parameter:: PFM_sigma_ref=0.0728
real,parameter :: u_ref = 1 
real,parameter::latent=334000
real,parameter::  PFM_l= droplet_radius*0.1
real,parameter::  PFM_l_c= PFM_l
real,parameter :: mobility = 2.5e-11
!real,parameter :: mobility = 2.5e-15
real,parameter :: mobility_c = 2.0e-3
real,parameter :: lamda_c = 3.0*sqrt(2.0)*PFM_l_c/2.0e5
real,parameter :: t_melt = 273   !all the temperature should be give as K
real,parameter:: rho_2= 1000  !liquid
real,parameter:: rho_1= 1
real,parameter:: rho_3= 900
real,parameter::vis_2=  1.0e-3
real,parameter::vis_1= 1.6e-5
real,parameter::vis_3= vis_2
real,parameter:: g_x = 0.0
real,parameter:: g_y = 0.0 
real,parameter:: g_z = -9.8
real,parameter:: kk_2=0.5918
real,parameter:: kk_1=0.0209
real,parameter:: kk_3=2.25
real,parameter:: cp_2=4200
real,parameter:: cp_1=1000
real,parameter:: cp_3=2018
real,parameter::PFM_thetta= 90*pi/180
real,parameter::PFM_mu_f= 1.28
real,parameter:: slip_length =  0.000025*droplet_radius
real,parameter:: vel_wall_top = 0.
real,parameter:: vel_wall_bot = -0.
logical, parameter :: Wetting = .true.
end module mod_cons
