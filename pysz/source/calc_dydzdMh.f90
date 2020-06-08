#define MAXSIZE 4096
SUBROUTINE calc_dydzdMh(h0_in, obh2_in, och2_in, mnu_in, mass_bias_in,&
                   pk_nk_in, pk_nz_in, k_arr, pk_arr, &
                   nell_ptilde_in, ell_ptilde_arr, ptilde_arr, & 
                   z_in, Mh_in, dydzdMh, flag_nu_in)
  !$ USE omp_lib
  use cosmo
  use global_var
  use linearpk_z
  use sigma_z
  use angular_distance
  use ptilde_module
  IMPLICIT none
  integer, intent(IN) :: pk_nk_in, pk_nz_in, nell_ptilde_in
  double precision, intent(IN) :: h0_in, obh2_in, och2_in, mnu_in
  double precision, intent(IN) :: mass_bias_in
  double precision, dimension(0:pk_nk_in-1,0:pk_nz_in-1), intent(IN) :: k_arr, pk_arr
  double precision, dimension(0:nell_ptilde_in-1), intent(IN) :: ell_ptilde_arr, ptilde_arr
  double precision, intent(IN) :: z_in, Mh_in
  double precision, intent(INOUT) :: dydzdMh
  integer, intent(IN) :: flag_nu_in

  ! flag for neutrino prescription
  flag_nu = flag_nu_in

  ! read in linear P(k,z)
  allocate(pk_z_arr(0:0))
  pk_z_arr(:) = z_in
  pk_nk = pk_nk_in
  pk_nz = pk_nz_in
  call open_linearpk(pk_nk,pk_nz,k_arr,pk_arr)

  ! read ptilde
  call open_ptilde(nell_ptilde_in,ell_ptilde_arr,ptilde_arr)
 
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
  !! fit sigma^2(R) to Chebyshev polynomials
  call compute_sigma2

  !! compute and tabulate da(z) 
  call setup_da

  ! calculate Cls
  !$OMP barrier
  call get_dydzdMh(z_in,Mh_in,dydzdMh)
  !$OMP barrier

  call close_linearpk
  call close_sigma
  call close_ptilde
  deallocate(pk_z_arr)
  !$OMP barrier

!===============================================================
CONTAINS
!===============================================================
  SUBROUTINE get_dydzdMh(z,Mh,dydzdMh)
    !$ USE omp_lib
    USE global_var
    USE cosmo
    USE mf_module
    USE angular_distance
    IMPLICIT none
    double precision, intent(INOUT) :: dydzdMh
    double precision, intent(IN) :: z, Mh
    double precision :: lnx, lnMh, dvdz, uy
    double precision :: da, calc_uy, ell
    external da, calc_uy

    dydzdMh = 0.d0

    lnx = dlog(z+1d0)
    lnMh = dlog(Mh)
    dvdz = (1d0+z)**2d0*da(z)**2d0*2998d0/Ez(z) ! h^-3 Mpc^3

    ell = 0d0
    uy = calc_uy(lnx,lnMh,ell)
    dydzdMh = uy*dndlnMh_500c_T08(lnMh,z)*dvdz

    return
  END SUBROUTINE get_dydzdMh
!===============================================================
END SUBROUTINE calc_dydzdMh
