module mod_thomas
implicit none
private
public thomas
contains
subroutine thomas(tz,aaa,bbb,ccc,kk)

!here default is z direction, but it can be other direction as well

implicit none
real, intent(inout), dimension(1:) :: tz
real, intent(in), dimension(1:) :: aaa
real, intent(in), dimension(1:) :: bbb
real, intent(in), dimension(1:) :: ccc
integer :: i,j,k,kk
real :: z,d(1:kk)
!!

!do j=1,jmax
!  do i=1,imax
    !z        = 1./bb(1)
    d(1) = ccc(1)/bbb(1)
    tz(1) = tz(1)/bbb(1)
!  enddo
!enddo
do k=2,kk-1
!   do j=1,jmax
!     do i=1,imax
       z        = 1./(bbb(k)-aaa(k)*d(k-1))
       d(k) = ccc(k)*z
       tz(k) = (tz(k)-aaa(k)*tz(k-1))*z
!     enddo
!  enddo
enddo
! backward substitution
!do j=1,jmax
!  do i=1,imax
    z        = bbb(kk)-aaa(kk)*d(kk-1)
    if(z.ne.0.) then
      tz(kk) = (tz(kk)-aaa(kk)*tz(kk-1))/z
    else
      tz(kk) =0.
    endif
!  enddo
!enddo

do k=kk-1,1,-1
!  do j=1,jmax
!    do i=1,imax
      tz(k) = tz(k)-d(k)*tz(k+1)
!    enddo
!  enddo
enddo

!!
return
end subroutine thomas
!
end module mod_thomas
