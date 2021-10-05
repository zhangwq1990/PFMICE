module mod_chem
use mod_var
use mod_cons
use mod_global
implicit none
private
public ChemicalPotential, ChemicalPotential_im
contains

subroutine ChemicalPotential
implicit none
integer::i,j,k
real::Psi_prime,d2phidx2,d2phidy2,d2phidz2,lap_phi

do k=1,kmax
 do j=1,jmax
   do i=1,imax
     Psi_prime = 2.0*PFM_phi(i,j,k)*(PFM_phi(i,j,k)-1)*(2.0*PFM_phi(i,j,k)-1)
     d2phidx2 = (PFM_phi(i+1,j,k)-2.*PFM_phi(i,j,k)+PFM_phi(i-1,j,k))/dx2
     d2phidy2 = (PFM_phi(i,j+1,k)-2.*PFM_phi(i,j,k)+PFM_phi(i,j-1,k))/dy2
     d2phidz2 = (PFM_phi(i,j,k+1)-2.*PFM_phi(i,j,k)+PFM_phi(i,j,k-1))/dz2
     lap_phi= d2phidx2 + d2phidy2 + d2phidz2
     Chem_Pot(i,j,k) =3.0*sqrt(2.0)*(PFM_Sigma_ref/PFM_l*Psi_prime-PFM_Sigma_ref*PFM_l*lap_phi)
   enddo
  enddo
enddo
return

end subroutine ChemicalPotential



subroutine ChemicalPotential_im
implicit none
integer::i,j,k
real::Psi_prime,d2phidx2,d2phidy2,d2phidz2,lap_phi
real:: middle,SSSS

middle=3.0*sqrt(2.)*PFM_Sigma_ref*PFM_l  !this is the gamma1 in the paper
SSSS=1.0*PFM_l**2*sqrt(4.0*1.0/mobility/middle/dt)  !this is the S in the paper


do k=1,kmax
 do j=1,jmax
   do i=1,imax
     Psi_prime = 2.0*PFM_phi_old(i,j,k)*(PFM_phi_old(i,j,k)-1)*(2.0*PFM_phi_old(i,j,k)-1)
     d2phidx2 = (PFM_phi(i+1,j,k)-2.*PFM_phi(i,j,k)+PFM_phi(i-1,j,k))/dx2
     d2phidy2 = (PFM_phi(i,j+1,k)-2.*PFM_phi(i,j,k)+PFM_phi(i,j-1,k))/dy2
     d2phidz2 = (PFM_phi(i,j,k+1)-2.*PFM_phi(i,j,k)+PFM_phi(i,j,k-1))/dz2
     lap_phi= d2phidx2 + d2phidy2 + d2phidz2
     Chem_Pot(i,j,k) =3.0*sqrt(2.0)*(PFM_Sigma_ref/PFM_l*Psi_prime-PFM_Sigma_ref*PFM_l*lap_phi)+&
                      middle*SSSS/PFM_l**2*(PFM_phi(i,j,k)-PFM_phi_old(i,j,k))
   enddo
  enddo
enddo
return

end subroutine ChemicalPotential_im




end module mod_chem
