MODULE cosmo
  ! cosmological parameters
  double precision :: h0, obh2, och2, om0_cb, ode0, w, mnu, onu
!!$================================================
CONTAINS
!!$================================================
DOUBLE PRECISION FUNCTION Ez(z)
  ! Ez = H(z)/H0 (dimensionless)
  ! x=a
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: z
  double precision :: x, ok0, or0=0d0 ! radiation density is ignored!
  ok0 = 1d0-om0_cb-ode0-onu
  x = 1d0/(1d0+z)
  Ez = dsqrt((om0_cb+onu)/x**3d0+ok0/x**2d0+or0/x**4d0+ode0/x**(3d0+3d0*w))
  return
END FUNCTION Ez
!!$================================================
END MODULE cosmo
