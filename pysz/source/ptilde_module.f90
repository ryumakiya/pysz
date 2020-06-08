!-----------------------------------------------------------------------------
! Read in, tabulate, and interpolate the \tilde{p} (Fourier transform of generalized NFW).
! The input file is in the following:
!
! [1st column] k (ell/l500) 
! [2nd column] tilde_p(x)
!-----------------------------------------------------------------------------
MODULE ptilde_module
  IMPLICIT none
  INTEGER :: ndata
  DOUBLE PRECISION, dimension(:), allocatable :: xa,ya,y2a
  DOUBLE PRECISION :: kmin, kmax
contains
  SUBROUTINE open_ptilde(n,ell_in,ptilde_in)
    IMPLICIT none
    INTEGER, intent(IN) :: n
    DOUBLE PRECISION, dimension(0:n-1), intent(IN) :: ell_in, ptilde_in
    kmin = ell_in(0)
    kmax = ell_in(n-1)
    ndata = n
    allocate(xa(ndata), ya(ndata), y2a(ndata))
    xa(1:ndata) = dlog(ell_in(0:n-1))
    ya(1:ndata) = dlog(ptilde_in(0:n-1))
    call spline(xa,ya,ndata,1.d30,1.d30,y2a)
    return
  END SUBROUTINE open_ptilde
  SUBROUTINE close_ptilde
    deallocate(xa,ya,y2a)
    return
  END SUBROUTINE close_ptilde
END MODULE ptilde_module
!------------------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION ptilde(ak)
  Use ptilde_module
  IMPLICIT none
  DOUBLE PRECISION :: a,b,h,x,y,ak
  INTEGER :: jlo

  ! cut off
  if(ak > kmax) then
    ptilde = 0d0
    return
  end if

  ! extrapolate
  if(ak <= kmin) then 
    ak = kmin
    jlo = 1
  else
    x  = dlog(ak)
    call hunt(xa,ndata,x,jlo)
  endif

  h=xa(jlo+1)-xa(jlo)
  a=(xa(jlo+1)-x)/h
  b=(x-xa(jlo))/h
  y=a*ya(jlo)+b*ya(jlo+1)+((a**3-a)*y2a(jlo)+(b**3-b)*y2a(jlo+1))*(h**2)/6.
  ptilde = dexp(y)

if (isnan(ptilde) .eqv. .True.) then
  print *, ak,x,jlo,xa(jlo+1),xa(jlo),h
endif

  return
END FUNCTION ptilde

