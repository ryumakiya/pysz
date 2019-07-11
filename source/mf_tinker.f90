!-----------------------------------------------------------------------------
!
! Halo multiplicity function in double precision, mf(lnnu,z), 
! where lnnu=ln(nu) and nu is called the "threshold," 
! given by nu = [delta_c/sigma(R,z)]^2. Here, delta_c=1.6865, sigma(R,z) 
! is the r.m.s. mass fluctuation within a top-hat smoothing of scale R 
! at a redshift z. 
!
! The functional form of the mass function is from Tinker et al., ApJ, 
! 688, 709 (2008). Note that "mf(lnnu,z)" here corresponds to f(sigma)/2 
! in this reference, where sigma=sqrt(1.6865/nu).
!
!-----------------------------------------------------------------------------

DOUBLE PRECISION FUNCTION mf_T08_intp(lnnu,z,delta)
!!$ Halo multiplicity function per log(nu): 0.5*f(sigma), 
!!$ where sigma=deltac/sqrt(nu).
!!$ \int_{-infty}^infty mf(lnnu,z) dlnnu = 1 is NOT satisfied!
!!$ Ref: Eq.(3) of Tinker et al. (2008) 
!!$ parameters are interpolated by splint
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: lnnu,z,delta ! delta is in m, not c
  DOUBLE PRECISION :: nu,sigma
  double precision :: Ap,a,b,c,zp
  double precision :: Ap0, a0, b0, c0
  double precision :: Apz = -0.14d0, az = -0.06d0, bz, cz = 0d0

  if (z > 3d0) then
    zp = 3d0
  else
    zp = z
  endif

  call intp_mf_params(delta,Ap0,a0,b0,c0)
  bz = 10**(-(0.75/log10(delta/75.d0))**1.2)
  bz = -1*bz

  Ap = Ap0*(1d0+zp)**Apz
  a = a0*(1d0+zp)**az
  b = b0*(1d0+zp)**bz
  c = c0*(1d0+zp)**cz
  nu = exp(lnnu)
  sigma = 1.6865d0/dsqrt(nu)
  mf_T08_intp = 0.5d0*(Ap*((sigma/b)**(-a)+1d0)*dexp(-c/sigma**2d0))
  return
END FUNCTION mf_T08_intp
!-----------------------------------------------------------------------------
SUBROUTINE intp_mf_params(Delta,Ap,a,b,c)
  double precision, intent(IN) :: Delta
  double precision, intent(INOUT) :: Ap, a, b, c
  double precision :: Delta_arr(9), logDelta_arr(9), Ap_arr(9), a_arr(9), b_arr(9), c_arr(9)
  double precision :: d2_Ap_arr(9), d2_a_arr(9), d2_b_arr(9), d2_c_arr(9)
  integer :: i, nx = 9
  double precision :: splint_mf
  external splint_mf

  Delta_arr = (/200d0,300d0,400d0,600d0,800d0,1200d0,1600d0,2400d0,3200d0/)
  do i = 1, 9
    logDelta_arr(i) = log10(Delta_arr(i))
  end do

  ! f(sigma) parameters from Tinker et al. (2008)
  Ap_arr = (/0.186d0,0.200d0,0.212d0,0.218d0,0.248d0,0.255d0,0.260d0,0.260d0,0.260d0/)
  a_arr = (/1.47d0,1.52d0,1.56d0,1.61d0,1.87d0,2.13d0,2.30d0,2.53d0,2.66d0/)
  b_arr = (/2.57d0,2.25d0,2.05d0,1.87d0,1.59d0,1.51d0,1.46d0,1.44d0,1.41d0/)
  c_arr = (/1.19d0,1.27d0,1.34d0,1.45d0,1.58d0,1.80d0,1.97d0,2.24d0,2.44d0/)

  ! 2nd derivative of parameters
  d2_Ap_arr = (/0.00d0,0.50d0,-1.56d0,3.05d0,-2.95d0,1.07d0,-0.71d0,0.21d0,0.00d0/)
  d2_a_arr = (/0.00d0,1.19d0,-6.34d0,21.36d0,-10.95d0,2.59d0,-0.85d0,-2.07d0,0.00d0/)
  d2_b_arr = (/0.00d0,-1.08d0,12.61d0,-20.96d0,24.08d0,-6.64d0,3.84d0,-2.09d0,0.00d0/)
  d2_c_arr = (/0.00d0,0.94d0,-0.43d0,4.61d0,0.01d0,1.21d0,1.43d0,0.33d0,0.00d0/)

  Ap = splint_mf(log10(Delta),logDelta_arr,Ap_arr,d2_Ap_arr,nx)
  a = splint_mf(log10(Delta),logDelta_arr,a_arr,d2_a_arr,nx)
  b = splint_mf(log10(Delta),logDelta_arr,b_arr,d2_b_arr,nx)
  c = splint_mf(log10(Delta),logDelta_arr,c_arr,d2_c_arr,nx)

END SUBROUTINE intp_mf_params
!======================================================================================
DOUBLE PRECISION FUNCTION splint_mf(x,x_arr,y_arr,d2y_arr,nx)
  integer, intent(IN) :: nx
  double precision, intent(IN) :: x, x_arr(nx), y_arr(nx), d2y_arr(nx)
  integer :: k_low, k_hi, k
  double precision :: y

  if (x <= x_arr(nx)) then ! spline interpolation
    call hunt(x_arr,nx,x,k)
    h = x_arr(k+1)-x_arr(k)
    a = (x_arr(k+1)-x)/h
    b = (x-x_arr(k))/h
    y = a*y_arr(k)+b*y_arr(k+1)+((a**3-a)*d2y_arr(k)+(b**3-b)*d2y_arr(k+1))*(h**2)/6.
  else ! linear interpolation
    y = (y_arr(nx)-y_arr(nx-1))/(x_arr(nx)-x_arr(nx-1))*(x-x_arr(nx-1))+y_arr(nx-1)
  end if 

  splint_mf = y
  return
END FUNCTION splint_mf
!======================================================================================
