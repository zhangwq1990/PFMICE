module mod_logo
use mod_global
implicit none
private
public logo
contains
subroutine logo
implicit none


if (myid .eq. 0) then
write(6,*) "--------------------------------------"
write(6,*)"--██████╗ ███████╗██╗ ██████╗███████╗--"
write(6,*)"--██╔══██╗██╔════╝██║██╔════╝██╔════╝--"
write(6,*)"--██████╔╝█████╗  ██║██║     █████╗  --"
write(6,*)"--██╔═══╝ ██╔══╝  ██║██║     ██╔══╝  --"
write(6,*)"--██║     ██║     ██║╚██████╗███████╗--"
write(6,*)"--╚═╝     ╚═╝     ╚═╝ ╚═════╝╚══════╝--"
write(6,*) "-----------------------------------"
write(6,*) "--------Version 1.0----------------"
write(6,*) "-----Resease date: Jun 8 2021------"
write(6,*) "-----------------------------------"
write(6,*) "-----------------------------------"
endif
                                   

return
end subroutine logo
!
end module mod_logo
