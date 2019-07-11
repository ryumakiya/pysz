MODULE global_var
  integer :: nl
  double precision :: mass_bias, alp_p, beta
  double precision :: z1 = 1.d-5, z2 = 4d0
  double precision :: Mmin = 1d11, Mmax = 5d15
  integer :: nz = 2001, nz_2h = 100
  integer :: flag_nu = 0 ! 0:'tot', 1:'cb'
  integer :: pk_nz, pk_nk
  double precision :: pi = 3.14159265359d0
END MODULE global_var
