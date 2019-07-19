!!$================================================================
  DOUBLE PRECISION FUNCTION calc_uy(lnx, lnM, ell)
    USE global_var
    USE angular_distance
    USE cosmo
    USE mod_ptilde
    IMPLICIT none
    double precision, intent(IN) :: lnx, lnM, ell
    double precision :: m500, r500, l500
    double precision :: omega, rhoc
    double precision :: xin=1.d-5, xout=6d0, tol=1.d-4
    double precision :: da, rombint1
    double precision :: z
    double precision :: ptilde, P0=6.41d0

    m500 = dexp(lnM)
    z = dexp(lnx)-1d0
  
    ! compute omega and rhoc
    rhoc=2.775d11*Ez(z)**2d0 ! critical density in units of h^2 M_sun/Mpc^3
  
    m500 = m500/(mass_bias*(1.d0+z)**beta)
    r500=(3d0*m500/4d0/pi/500d0/rhoc)**(1d0/3d0) ! h^-1 Mpc
    l500=da(z)/r500
    calc_uy = 1.65d0*(h0/0.7d0)**2d0*Ez(z)**(8d0/3d0) &
            *(m500/3d14/0.7d0)**(2d0/3d0+alp_p) &
            *P0*(0.7d0/h0)**1.5d0*ptilde(ell/l500) &
            /0.5176d0*(4d0*3.14159d0)/l500**2d0*(r500/h0)
    calc_uy = 1.04d-4*(calc_uy/50d0)  

    return
  END FUNCTION calc_uy
!!$================================================================
  DOUBLE PRECISION FUNCTION integrand_by(lnM,z,ell)
    use cosmo
    use mf_module
    use global_var
    IMPLICIT none
    double precision, intent(IN) :: lnM,z,ell
    double precision :: calc_uy, lnx,omega
  
    omega = (om0_cb+onu)*(1d0+z)**3d0/Ez(z)**2d0

    lnx = dlog(z+1d0)
    integrand_by = calc_uy(lnx,lnM,ell) &
                  *bl_delta(lnnu_500c(lnM,z),500d0/omega) &
                  *dndlnMh_500c_T08(lnM,z)
 
  END FUNCTION integrand_by
!!$================================================================
!   FUNCTION pgnfw
!   ! Ref: Planck Collaboration, A&A, 550, A131 (2013)
!     USE cosmo
!     IMPLICIT none
!     double precision :: x,y
!     double precision :: c=1.81d0,g=0.31d0,a=1.33d0,b=4.13d0,P0=6.41d0 ! Planck 2013
!     double precision :: sinc
!     external sinc
!   
!     pgnfw=P0*(0.7d0/h0)**1.5d0/(c*x)**g/(1d0+(c*x)**a)**((b-g)/a) &
!          *x**2d0*sinc(y*x) !dsin(y*x)/(y*x)
!     return
!   END FUNCTION pgnfw
! !!$================================================================
!   DOUBLE PRECISION FUNCTION sinc(x)
!     double precision :: x
!     if (abs(x)<=1d-2) then
!       sinc = 1.-x*x/6.+x**4/120.
!     else
!       sinc = dsin(x)/x
!     endif
!     return
!   END FUNCTION sinc
! !!$================================================================
