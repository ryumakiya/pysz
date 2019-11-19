#define MAXSIZE 4096
SUBROUTINE calc_by(h0_in, obh2_in, och2_in, mnu_in, mass_bias_in,&
                   pk_nk_in, pk_nz_in, k_arr, pk_arr,&
                   nz_in, z_arr, by_arr, dydz_arr, flag_nu_in)
  !$ USE omp_lib
  use cosmo
  use global_var
  use linearpk_z
  use sigma_z
  use angular_distance
  use mod_ptilde
  IMPLICIT none
  integer, intent(IN) :: pk_nk_in, pk_nz_in, nz_in
  double precision, intent(IN) :: h0_in, obh2_in, och2_in, mnu_in
  double precision, intent(IN) :: mass_bias_in
  double precision, dimension(0:pk_nk_in-1,0:pk_nz_in-1), intent(IN) :: k_arr, pk_arr
  double precision, dimension(0:nz_in-1), intent(INOUT) :: z_arr
  double precision, dimension(0:nz_in-1), intent(INOUT) :: by_arr
  double precision, dimension(0:nz_in-1), intent(INOUT) :: dydz_arr
  integer, intent(IN) :: flag_nu_in

  ! flag for neutrino prescription
  flag_nu = flag_nu_in

  ! read in linear P(k,z)
  pk_nk = pk_nk_in
  pk_nz = pk_nz_in
  call open_linearpk(pk_nk,pk_nz,k_arr,pk_arr)

  ! input parameters
  !! cosmological parameters
  h0 = h0_in
  obh2 = obh2_in
  och2 = och2_in
  om0_cb = (obh2+och2)/h0**2
  mnu = mnu_in
  onu = mnu/93.14/h0**2
  ode0 = 1.d0-om0_cb-onu
  w = -1d0

  !! tSZ parameters
  mass_bias = mass_bias_in
  ! fixed parameters
  alp_p = 0.12d0
  beta = 0.d0

  ! preface
  !! ptilde
  call setup_ptilde

  !! fit sigma^2(R) to Chebyshev polynomials
  call compute_sigma2

  !! compute and tabulate da(z) 
  call setup_da

  ! calculate Cls
  !$OMP barrier
  call calc_by_dydz(z_arr,by_arr,dydz_arr)
  !$OMP barrier

  call close_linearpk
  call close_sigma
  call close_ptilde
  !$OMP barrier

!===============================================================
CONTAINS
!===============================================================
  SUBROUTINE calc_by_dydz(z_arr,by_arr,dydz_arr)
    !$ USE omp_lib
    USE global_var
    IMPLICIT none
    double precision, intent(INOUT) :: by_arr(0:nz_in-1)
    double precision, intent(INOUT) :: dydz_arr(0:nz_in-1)
    double precision, intent(IN) :: z_arr(0:nz_in-1)
    integer :: i, j, k, nm=101
    double precision :: lnx, dvdz, z
    double precision :: lnM, lnM1, lnM2, dlnM
    double precision :: a1,a2,da
    external da

    by_arr(:) = 0d0
    dydz_arr(:) = 0d0

    lnM1=dlog(Mmin); lnM2=dlog(Mmax)
    dlnM = (lnM2-lnM1)/nm

    do i = 0, nz_in-1
      z = z_arr(i)
      lnx = dlog(z+1d0)
      dvdz = (1d0+z)**2d0*da(z)**2d0*2998d0/Ez(z) ! h^-3 Mpc^3
      ! mass integration
      a1 = (integrand1(lnM1,lnx)+integrand1(lnM2,lnx))*(0.5d0*dlnM)
      a2 = (integrand2(lnM1,lnx)+integrand2(lnM2,lnx))*(0.5d0*dlnM)
      do j = 1, nm-1
        lnM = lnM1+dlnM*j
        a1 = a1+integrand1(lnM,lnx)*dlnM
        a2 = a2+integrand2(lnM,lnx)*dlnM
      end do
      by_arr(i) = a2/a1
      dydz_arr(i) = a1*dvdz
    end do
    
  END SUBROUTINE calc_by_dydz
!===============================================================
  DOUBLE PRECISION FUNCTION integrand1(lnM, lnx)
    USE global_var
    USE cosmo
    USE mf_module
    USE angular_distance
    IMPLICIT none
    double precision, intent(IN) :: lnx, lnM
    double precision :: uy, ell
    double precision :: z
    double precision :: calc_uy
    external calc_uy

    ell = 0d0
    uy = calc_uy(lnx,lnM,ell)
    z = dexp(lnx)-1d0
    integrand1 = uy*dndlnMh_500c_T08(lnM,z)

    return  
  END FUNCTION integrand1
!===============================================================
  DOUBLE PRECISION FUNCTION integrand2(lnM, lnx)
    USE global_var
    USE cosmo
    USE mf_module
    USE angular_distance
    IMPLICIT none
    double precision, intent(IN) :: lnx, lnM
    double precision :: uy, ell, bl
    double precision :: z, omega
    double precision :: calc_uy
    external calc_uy

    ell = 0d0
    uy = calc_uy(lnx,lnM,ell)
    z = dexp(lnx)-1d0
    omega = (om0_cb+onu)*(1d0+z)**3d0/Ez(z)**2d0
    bl = bl_delta(lnnu_500c(lnM,z),500d0/omega)
    integrand2 = uy*bl*dndlnMh_500c_T08(lnM,z)

    return  
  END FUNCTION integrand2
!===============================================================
END SUBROUTINE calc_by
