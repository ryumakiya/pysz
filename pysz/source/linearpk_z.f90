!-----------------------------------------------------------------------------
!
! Read in, tabulate, and interpolate the linear matter power spectrum.
! The input file is in the following format (same as CAMB format):
!
! [1st column] k/h (i.e., wavenumber in units of h Mpc^-1) 
! [2nd column] h^3 P(k) (i.e., P(k) in units of h^-3 Mpc^3)
!
! << sample code >>
! 
! USE linearpk
! double precision :: linear_pk, k_ov_h
! character(len=128) :: filename
! integer :: n
! external linear_pk
! filename='wmap5baosn_max_likelihood_matterpower.dat'
! n=896 ! # of lines in the file
! CALL open_linearpk(filename,n)
! k_ov_h=1d0
! print*,'P(k) at k=',k_ov_h,' h Mpc^-1 is',linear_pk(k_ov_h),' h^3 Mpc^-3'
! end
!
! August 23, 2008
! E. Komatsu
!
!-----------------------------------------------------------------------------
MODULE LinearPk_z
  IMPLICIT none
  INTEGER :: ndata,nz
  DOUBLE PRECISION, dimension(:,:), allocatable :: xa,ya,y2a
  DOUBLE PRECISION :: kmin_in, kmax_in
contains
  SUBROUTINE Open_LinearPk(n,nz,aklin,pk_lin)
    IMPLICIT none
    INTEGER, intent(IN) :: n
    INTEGER, intent(IN) :: nz
    DOUBLE PRECISION, dimension(0:n-1,0:nz-1), intent(IN) :: aklin, pk_lin
    INTEGER :: iz
    kmin_in=aklin(0,0)
    kmax_in=aklin(n-1,0)
    ndata=n
    ALLOCATE(xa(ndata,nz),ya(ndata,nz),y2a(ndata,nz))
    do iz = 1, nz
      xa(:,iz) = dlog(aklin(0:n-1,iz-1))
      ya(:,iz) = dlog(pk_lin(0:n-1,iz-1))
      CALL spline(xa(:,iz),ya(:,iz),ndata,1.d30,1.d30,y2a(:,iz))
    end do
    return
  END SUBROUTINE Open_LinearPk
  SUBROUTINE Close_LinearPk
    DEALLOCATE(xa,ya,y2a)
    return
  END SUBROUTINE Close_LinearPk
END MODULE LinearPk_z
DOUBLE PRECISION FUNCTION Linear_Pk(ak,iz)
  Use LinearPk_z
  IMPLICIT none
  DOUBLE PRECISION :: a,b,h,x,y,ak
  INTEGER :: jlo,iz

  if((ak > kmax_in) .or. (ak < kmin_in)) then
    Linear_Pk = 0d0
    return
  endif

  x  = dlog(ak)
  CALL hunt(xa(:,iz),ndata,x,jlo)
  h=xa(jlo+1,iz)-xa(jlo,iz)
  a=(xa(jlo+1,iz)-x)/h
  b=(x-xa(jlo,iz))/h
  y=a*ya(jlo,iz)+b*ya(jlo+1,iz)+((a**3-a)*y2a(jlo,iz)+(b**3-b)*y2a(jlo+1,iz))*(h**2)/6.
  Linear_Pk = dexp(y)
  return
END FUNCTION Linear_Pk
DOUBLE PRECISION FUNCTION dlnPkdlnk(ak,iz)
  Use LinearPk_z
  IMPLICIT none
  DOUBLE PRECISION :: a,b,h,x,y,ak
  INTEGER :: jlo,iz
  x  = dlog(ak)
  CALL hunt(xa(:,iz),ndata,x,jlo)
  h=xa(jlo+1,iz)-xa(jlo,iz)
  a=(xa(jlo+1,iz)-x)/h
  b=(x-xa(jlo,iz))/h
  y=(ya(jlo+1,iz)-ya(jlo,iz))/h+(-(3.*a**2-1.)*y2a(jlo,iz)+(3.*b**2-1.)*y2a(jlo+1,iz))*h/6.
  dlnPkdlnk = y
  return
END FUNCTION DlnPkdlnk
