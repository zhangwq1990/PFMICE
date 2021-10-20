module mod_old
use mod_var
use mod_var
implicit none
private
public old
contains
subroutine old
implicit none

! give old values
rholold2=rholold
rholold=rhol
rholcp_old=rholcp
PFM_phi_old=PFM_phi
PFM_c_old=PFM_c
kkold=kk
vislold=visl
Told=Tnew
uold=unew
vold=vnew
wold=wnew
pold=pnew
chem_pot_old=chem_pot
Phi_c_old=Phi_c
delta_u_old=delta_u



return
end subroutine old
!
end module mod_old
