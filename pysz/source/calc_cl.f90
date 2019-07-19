#define MAXSIZE 4096
SUBROUTINE calc_cl(h0_in, obh2_in, och2_in, mnu_in, mass_bias_in,&
                   pk_nk_in, pk_nz_in, k_arr, pk_arr,&
                   nl_in, ell_arr, cl_yy, tll, flag_nu_in, flag_tll_in)
  !$ USE omp_lib
  use cosmo
  use global_var
  use linearpk_z
  use sigma_z
  use angular_distance
  use mod_ptilde
  IMPLICIT none
  integer, intent(IN) :: pk_nk_in, pk_nz_in, nl_in
  double precision, intent(IN) :: h0_in, obh2_in, och2_in, mnu_in
  double precision, intent(IN) :: mass_bias_in
  double precision, dimension(0:pk_nk_in-1,0:pk_nz_in-1), intent(IN) :: k_arr, pk_arr
  double precision, dimension(0:nl_in-1), intent(INOUT) :: ell_arr
  double precision, dimension(0:nl_in-1,0:1), intent(INOUT) :: cl_yy
  double precision, dimension(0:nl_in-1,0:nl_in-1), intent(INOUT) :: tll
  integer, intent(IN) :: flag_nu_in, flag_tll_in

  nl = nl_in

  ! flag for neutrino prescription
  flag_nu = flag_nu_in

  ! Calc Tll or not
  flag_tll = flag_tll_in

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
  call calc_cl_yy(ell_arr,cl_yy,tll)
  !$OMP barrier

  call close_linearpk
  call close_sigma
  call close_ptilde
  !$OMP barrier

!===============================================================
CONTAINS
!===============================================================
  SUBROUTINE calc_cl_yy(ell_arr,cl_yy,tll)
    !$ USE omp_lib
    USE global_var
    IMPLICIT none
    double precision, intent(INOUT) :: cl_yy(0:nl-1,0:1)
    double precision, intent(INOUT) :: tll(0:nl-1,0:nl-1)
    double precision, intent(IN) :: ell_arr(:)
    double precision, allocatable :: cls_th_1h(:,:), cls_th_2h(:,:)
    double precision :: ell
    integer :: i, j, k, nth, ith
    double precision :: lnx, lnx1, lnx2, dlnx
    double precision :: lnM, lnM1, lnM2, dlnM
    double precision :: fac1, fac2
    double precision :: intg_1h(nl), intg_2h(nl), intg_tll(nl,nl)

    cl_yy(:,:) = 0.d0
    intg_tll(:,:) = 0.d0

    lnx1=dlog(1d0+z1); lnx2 = dlog(1d0+z2)
    dlnx = (lnx2-lnx1)/nz
  
    lnM1=dlog(Mmin); lnM2=dlog(Mmax)
    
    ! 1-halo term
    !$OMP parallel private(j,lnx,fac1,lnM,intg_1h,intg_2h,ith), shared(cls_th_1h,cls_th_2h,nth,dlnx)
    nth = omp_get_num_threads()
    ith = omp_get_thread_num()
    !$OMP single
    allocate(cls_th_1h(nl,nth),cls_th_2h(nl,nth))
    cls_th_1h(:,:) = 0.d0
    cls_th_2h(:,:) = 0.d0
    !$OMP end single
    !$OMP do
    do j = 1, nz+1
      lnx = lnx1+dlnx*(j-1)
      fac1 = 1.d0
      if (j == 1 .or. j == nz+1) fac1 = 0.5d0
      ! mass integration
      call qgaus2_n20_arr(integrand_1h,lnM1,lnM2,intg_1h,lnx,ell_arr)
      cls_th_1h(:,ith+1) = cls_th_1h(:,ith+1)+intg_1h*dlnx*fac1
    end do 
    !$OMP end do
    !$OMP barrier

    !$OMP barrier
    ! 2-halo term
    dlnx = (lnx2-lnx1)/nz_2h
    !$OMP do
    do j = 1, nz_2h+1
      lnx = lnx1+dlnx*(j-1)
      fac1 = 1.d0
      if (j == 1 .or. j == nz_2h+1) fac1 = 0.5d0
      call calc_integrand_2h(lnx,ell_arr,intg_2h)
      cls_th_2h(:,ith+1) = cls_th_2h(:,ith+1)+intg_2h*dlnx*fac1
    end do
    !$OMP end do
    !$OMP barrier
    !$OMP end parallel

    do i = 1, nl
      cl_yy(i-1,0) = sum(cls_th_1h(i,1:nth))
      cl_yy(i-1,1) = sum(cls_th_2h(i,1:nth))
    end do

    !$OMP single
    deallocate(cls_th_1h,cls_th_2h)
    !$OMP end single

    ! Tll
    if (flag_tll == 1) then
      dlnM = (lnM2-lnM1)/nm_tll
      dlnx = (lnx2-lnx1)/nz_tll
      do j = 1, nz_tll+1
        lnx = lnx1+dlnx*(j-1)
        fac1 = 1.d0
        if (j == 1 .or. j == nz+1) fac1 = 0.5d0
        do k =1, nm_tll+1
          lnM = lnM1+dlnM*(k-1)
          fac2 = 1.d0
          if (k == 1 .or. k == nm_tll+1) fac2 = 0.5d0
          intg_tll = integrand_tll(lnM,lnx,ell_arr)
          tll(0:nl-1,0:nl-1) = tll(0:nl-1,0:nl-1)+intg_tll*dlnM*dlnx*fac1*fac2
        end do
      end do
      tll = tll/4/pi
    end if

  END SUBROUTINE calc_cl_yy
!===============================================================
  FUNCTION integrand_1h(lnM, lnx, ell_arr)
    USE global_var
    USE cosmo
    USE mf_module
    USE angular_distance
    IMPLICIT none
    double precision :: integrand_1h(nl)
    double precision, intent(IN) :: lnx, lnM, ell_arr(nl)
    double precision :: uy, ell
    double precision :: z, dvdz, fac, da
    double precision :: calc_uy
    integer :: i
    external calc_uy, da

    integrand_1h(:) = 0d0

    do i = 1, nl
      ell = ell_arr(i)
      uy = calc_uy(lnx,lnM,ell)
      integrand_1h(i) = uy*uy
    end do

    z = dexp(lnx)-1d0
    dvdz = (1d0+z)**2d0*da(z)**2d0*2998d0/Ez(z) ! h^-3 Mpc^3
    fac = dvdz*(1+z)*dndlnMh_500c_T08(lnM,z)

    integrand_1h = integrand_1h*fac

    return  
  END FUNCTION integrand_1h
!===============================================================
  SUBROUTINE calc_integrand_2h(lnx, ell_arr, intg_2h)
    USE global_var
    USE cosmo
    USE mf_module
    USE angular_distance
    USE linearpk_z
    IMPLICIT none
    double precision, intent(IN) :: lnx, ell_arr(nl)
    double precision, intent(INOUT) :: intg_2h(nl)
    double precision :: ell
    double precision :: z, dvdz, fac, da, chi
    double precision :: by, lnM1, lnM2, lnM, dlnM
    double precision :: integrand_by
    double precision :: linear_pk
    double precision :: pk1, pk2, pk_z, dz
    integer :: il, iz
    external linear_pk, integrand_by

    intg_2h = 0d0
    z = dexp(lnx)-1d0
    dz = (z2-0d0)/(pk_nz-1)
    iz = int((z-1d-5)/dz)

    do il = 1, nl
      ell = ell_arr(il)

      ! SZ term
      lnM1 = dlog(Mmin); lnM2 = dlog(Mmax)
      call qgaus2(integrand_by,lnM1,lnM2,by,z,ell) ! Gaussian quadrature, SZ bias
  
      dvdz=(1d0+z)**2d0*da(z)**2d0*2998d0/Ez(z)
      chi = da(z)*(1.d0+z)
      if (iz < pk_nz-1) then
        pk1 = linear_pk(ell/chi,iz+1)
        pk2 = linear_pk(ell/chi,iz+2)
        pk_z = (pk2-pk1)/dz*(z-dz*iz)+pk1
      else
        pk_z = linear_pk(ell/chi,pk_nz) 
      end if

      fac = pk_z*dvdz*(1d0+z)
  
      intg_2h(il) = by*by*fac
    end do

    return  
  END SUBROUTINE calc_integrand_2h
!===============================================================
 FUNCTION integrand_tll(lnM, lnx, ell_arr)
    USE global_var
    USE cosmo
    USE mf_module
    USE angular_distance
    IMPLICIT none
    double precision :: integrand_tll(nl,nl), integrand_1h(nl,1)
    double precision, intent(IN) :: lnx, lnM, ell_arr(nl)
    double precision :: uy, ell
    double precision :: z, dvdz, fac, da
    double precision :: calc_uy
    integer :: i
    external calc_uy, da
    double precision :: fsky_yy = 0.494d0

    integrand_tll(:,:) = 0d0

    do i = 1, nl
      ell = ell_arr(i)
      uy = calc_uy(lnx,lnM,ell)
      integrand_1h(i,1) = uy*uy
    end do

    integrand_1h = integrand_1h/sqrt(fsky_yy)
    integrand_tll = matmul(integrand_1h,transpose(integrand_1h))

    z = dexp(lnx)-1d0
    dvdz = (1d0+z)**2d0*da(z)**2d0*2998d0/Ez(z) ! h^-3 Mpc^3
    fac = dvdz*(1+z)*dndlnMh_500c_T08(lnM,z)

    integrand_tll = integrand_tll*fac

    return  
  END FUNCTION integrand_tll
!===============================================================
END SUBROUTINE calc_cl
