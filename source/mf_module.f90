MODULE mf_module
  USE cosmo
  USE sigma_z
  implicit none
CONTAINS
!======================================================================================
  DOUBLE PRECISION FUNCTION lnnu_500c(lnM500,z)
    USE global_var
    double precision, intent(IN) :: lnM500,z
    double precision :: m500, Rh
    double precision :: deltac=1.6865d0
    double precision :: c1, c2, dz
    real :: chebev,lnsigma2,lnRh
    integer :: iz

    m500 = dexp(lnM500)
 
    ! calc dndlnRh
    if (flag_nu == 0) then
      Rh = (3d0*m500/4d0/pi/(om0_cb+onu)/2.775d11)**(1d0/3d0) ! h^-1 Mpc
    else if (flag_nu == 1) then
      Rh = (3d0*m500/4d0/pi/om0_cb/2.775d11)**(1d0/3d0) ! h^-1 Mpc
    end if
    lnRh = real(dlog(Rh))

    dz = (z2-0d0)/(pk_nz-1)
    iz = int((z-1d-5)/dz)
    if (iz < pk_nz-1) then
      c1 = CHEBEV(lnR1,lnR2,c(:,iz+1),ndim,lnRh)
      c2 = CHEBEV(lnR1,lnR2,c(:,iz+2),ndim,lnRh)
      lnsigma2 = (c2-c1)/dz*(z-dz*iz)+c1 ! ln(sigma^2)
    else 
      lnsigma2 = CHEBEV(lnR1,lnR2,c(:,pk_nz),ndim,lnRh)          ! ln(sigma^2)
    end if

    lnnu_500c = 2d0*dlog(deltac)-dble(lnsigma2) ! ln(nu)
  END FUNCTION lnnu_500c
!======================================================================================
  ! CALCULATE MASS FUNCTION
  DOUBLE PRECISION FUNCTION dndlnMh_500c_T08(lnM500,z)
    USE global_var
    double precision, intent(IN) :: lnM500, z
    double precision :: m500,omz
    double precision :: Rh
    double precision :: deltac=1.6865d0,mf_500c,dndlnRh,lnnu,dlnnudlnRh
    real :: chebev,lnsigma2,lnRh,dlnsigma2dlnRh
    double precision :: mf_T08_intp
    double precision :: c1, c2, cder1, cder2, dz
    integer :: iz
    external mf_T08_intp

    m500 = dexp(lnM500)
    omz = (om0_cb+onu)*(1d0+z)**3d0/Ez(z)**2d0 ! Omega(z); E(z)=H(z)/H0

    ! calc dndlnRh
    if (flag_nu == 0) then
      Rh = (3d0*m500/4d0/pi/(om0_cb+onu)/2.775d11)**(1d0/3d0) ! h^-1 Mpc
    else if (flag_nu == 1) then
      Rh = (3d0*m500/4d0/pi/om0_cb/2.775d11)**(1d0/3d0) ! h^-1 Mpc
    end if
    lnRh = real(dlog(Rh))

    dz = (z2-0d0)/(pk_nz-1)
    iz = int((z-1d-5)/dz)
    if (iz < pk_nz-1) then
      c1 = CHEBEV(lnR1,lnR2,c(:,iz+1),ndim,lnRh)
      c2 = CHEBEV(lnR1,lnR2,c(:,iz+2),ndim,lnRh)
      lnsigma2 = (c2-c1)/dz*(z-dz*iz)+c1 ! ln(sigma^2)
      cder1 = CHEBEV(lnR1,lnR2,cder(:,iz+1),ndim,lnRh)
      cder2 = CHEBEV(lnR1,lnR2,cder(:,iz+2),ndim,lnRh)
      dlnsigma2dlnRh = (cder2-cder1)/dz*(z-dz*iz)+cder1 ! ln(sigma^2)
    else 
      lnsigma2 = CHEBEV(lnR1,lnR2,c(:,pk_nz),ndim,lnRh)          ! ln(sigma^2)
      dlnsigma2dlnRh = CHEBEV(lnR1,lnR2,cder(:,pk_nz),ndim,lnRh) ! dln(sigma^2)/dlnRh
    end if
    lnnu=2d0*dlog(deltac)-dble(lnsigma2) ! ln(nu
    dlnnudlnRh=-dble(dlnsigma2dlnRh)     ! dln(nu)/dlnRh
    dndlnRh = (3d0/4d0/pi)*dlnnudlnRh*mf_T08_intp(lnnu,z,500d0/omz)/Rh**3d0
    dndlnMh_500c_T08 = dndlnRh/3d0 ! in units of h^3 Mpc^-3
  END FUNCTION dndlnMh_500c_T08
!======================================================================================
DOUBLE PRECISION FUNCTION bl_delta(lnnu,delta) ! linear bias factor bl
  DOUBLE PRECISION, intent(IN) :: lnnu, delta
  DOUBLE PRECISION :: nu, deltac = 1.6865d0
  DOUBLE PRECISION :: y, A, as, B, bs, C, cs

  !! Cooray & Sheth 2002
  !nu = dexp(lnnu)
  !bl = 1.d0+(0.707*nu-1d0)/deltac+2d0*0.3/deltac/(1d0+(0.707*nu)**0.3)

  ! Tinker et al. 2010
  nu = dexp(lnnu/2.d0)
  y = log10(delta)
  A  = 1.0+0.24*y*exp(-(4/y)**4)
  as = 0.44*y-0.88
  B = 0.183
  bs = 1.5
  C = 0.019+0.107*y+0.19*exp(-(4/y)**4)
  cs = 2.4
  bl_delta = 1.d0-A*(nu**as)/(nu**as+deltac**as)+B*nu**bs+C*nu**cs

  return
END FUNCTION bl_delta
!======================================================================================
END MODULE mf_module
!======================================================================================
