!-----------------------------------------------------------------------------
! Read in, tabulate, and interpolate the p_tilde (Fourier transform of generalized NFW).
! The input file is in the following:
!
! [1st column] k (ell/l500) 
! [2nd column] tilde_p(x)
!-----------------------------------------------------------------------------
MODULE mod_ptilde
  IMPLICIT none
  INTEGER :: ndata
  DOUBLE PRECISION, dimension(:), allocatable :: xa,ya,y2a
CONTAINS
!===================================================================================
  SUBROUTINE setup_ptilde()
    use input_ptilde
    IMPLICIT none
    CHARACTER(len=128) :: filename
    INTEGER :: i
    call ptilde_init
    ndata = ndata_ptilde
    ALLOCATE(xa(ndata),ya(ndata),y2a(ndata))
    xa=dlog(k_ptilde)
    ya=dlog(ptilde_in)
    DEALLOCATE(k_ptilde,ptilde_in)
    CALL spline(xa,ya,ndata,1.d30,1.d30,y2a)
    return
  END SUBROUTINE setup_ptilde
!===================================================================================
  SUBROUTINE close_ptilde
    DEALLOCATE(xa,ya,y2a)
    return
  END SUBROUTINE close_ptilde
!===================================================================================
END MODULE mod_ptilde
!===================================================================================
DOUBLE PRECISION FUNCTION ptilde(ak)
  Use mod_ptilde
  IMPLICIT none
  DOUBLE PRECISION :: a,b,h,x,y,ak
  INTEGER :: jlo
  ! cut off
  if(ak > 9.d3) then
    ptilde = 1d-30
    return
  end if

  if(ak <= 1.d-2) ak=1.d-2
  x  = dlog(ak)
  CALL hunt(xa,ndata,x,jlo)
  h=xa(jlo+1)-xa(jlo)
  a=(xa(jlo+1)-x)/h
  b=(x-xa(jlo))/h
  y=a*ya(jlo)+b*ya(jlo+1)+((a**3-a)*y2a(jlo)+(b**3-b)*y2a(jlo+1))*(h**2)/6.
  ptilde = dexp(y)
  return
END FUNCTION ptilde
!===================================================================================
