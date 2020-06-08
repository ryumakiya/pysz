PROGRAM calc_ptilde
  !$ use omp_lib
  IMPLICIT none
  double precision :: ptilde, k, lnk, xin=1d-5, xout=6d0, tol=1d-10
  integer :: i
  double precision :: ans, rombint1
  double precision :: k_arr(3002), ptilde_arr(3002)
  external integrand, rombint1

  k_arr(1) = 1d-20
  ptilde_arr(1) = rombint1(integrand,xin,xout,tol,k_arr(1))
  !$omp parallel
  !$omp do
  do i = 1, 3001
    lnk = -4d0+8d0/3000*(i-1)
    k_arr(i+1) = 10d0**(lnk)
    ptilde_arr(i+1) = rombint1(integrand,xin,xout,tol,k_arr(i+1))
  end do
  !$omp end do
  !$omp end parallel

  do i = 1, 3002
    print *, k_arr(i), ptilde_arr(i)
  end do

END PROGRAM calc_ptilde
!!$================================================================
  DOUBLE PRECISION FUNCTION integrand(x,y) ! x=r/r500 & y=l/l500
  ! Ref: Planck Collaboration, A&A, 550, A131 (2013)
    IMPLICIT none
    double precision :: x,y
    double precision :: c=1.81d0,g=0.31d0,a=1.33d0,b=4.13d0 ! Planck 2013
    double precision :: pgnfw, sinc
    external sinc

    pgnfw = 1.d0/(c*x)**g/(1d0+(c*x)**a)**((b-g)/a)
    integrand = pgnfw*x**2d0*sinc(y*x)

    return
  END FUNCTION integrand
!!$================================================================
  DOUBLE PRECISION FUNCTION sinc(x)
    double precision :: x
    if (abs(x)<=1d-2) then
      sinc = 1.-x*x/6.+x**4/120.
    else
      sinc = dsin(x)/x
    endif
    return
  END FUNCTION sinc
!!$================================================================
