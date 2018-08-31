module vec
  use iso_c_binding
  implicit none

#ifndef VEC_NCELL
# define VEC_NCELL 128
#endif
#ifndef VEC_DP
# define VEC_DP 1
#endif
#if VEC_DP
# define MYREAL 8
# define c_real c_double
#else
# define MYREAL 4
# define c_real c_float
#endif

  integer, parameter :: &
       ncell = VEC_NCELL, &
       myreal = MYREAL

  real(kind=myreal), parameter :: &
       L = 2.0_myreal, &
       xl = 0.0_myreal, &
       xr = xl + L, &
       u_max = 0.1_myreal, &
       dx = L / ncell, &
       rho_ref = 1.0_myreal

contains
  function get_x_ctr(i) result(xc)
    integer, intent(in) :: i
    real(kind=myreal) :: xc

    xc = xl + ((i - 0.5_myreal)/ncell)*L;
  end function get_x_ctr

  function map_x_to_n11(x) result(r)
    real(kind=myreal), intent(in) :: x
    real(kind=myreal) :: r

    r = 2.0_myreal*(x - xl)/L - 1.0_myreal
  end function map_x_to_n11

  function get_src(rho) result(src)
    real(kind=myreal), intent(in) :: rho
    integer, parameter :: ratesz = 5
    real(kind=myreal) :: src, trhodiffmax, rhodiff
    real(kind=myreal), parameter :: rate(ratesz) = &
         (/ 0.2_myreal, 0.1_myreal, 0.05_myreal, 0.025_myreal, 0.0125_myreal /)
    integer :: tsize, idx

    tsize = ratesz - 1.0_myreal;
    trhodiffmax = 1.0_myreal;
    rhodiff = abs(rho - rho_ref)
    if (rhodiff >= trhodiffmax) then
       src = 0.0_myreal
    else
       idx = min(tsize, int(tsize*(rhodiff/trhodiffmax)) + 1)
       src = rhodiff*rate(idx);
    end if
  end function get_src

  function get_u(x) result(u)
    real(kind=myreal), intent(in) :: x
    real(kind=myreal) :: u

    u = u_max*cos(0.2_myreal*map_x_to_n11(x) + 0.25_myreal);
  end function get_u

  function get_ic() result(ic)
    real(kind=myreal) :: ic

    ic = rho_ref
  end function get_ic

  subroutine calc_numerical_flux(rho, flux_bdy, flux_int)
    real(kind=myreal), intent(in) :: rho(ncell)
    real(kind=myreal), intent(out) :: flux_bdy, flux_int(ncell)
    integer :: i

    flux_bdy = rho_ref * get_u(xl - 0.5_myreal*dx);
    do i = 1, ncell
       flux_int(i) = rho(i) * get_u(get_x_ctr(i))
    end do
  end subroutine calc_numerical_flux

  subroutine step1(dt, rho, flux_int)
    real(kind=myreal), intent(in) :: dt
    real(kind=myreal), intent(inout) :: rho(ncell), flux_int(ncell)
    real(kind=myreal) :: rdx, flux_bdy, neg_flux_div, src
    integer :: i

    call calc_numerical_flux(rho, flux_bdy, flux_int)
    rdx = 1.0_myreal/dx
    neg_flux_div = rdx*(flux_bdy - flux_int(1));
    src = get_src(rho(1));
    rho(1) = rho(1) + dt*(src + neg_flux_div);
    do i = 2, ncell
       neg_flux_div = rdx*(flux_int(i-1) - flux_int(i));
       src = get_src(rho(i));
       rho(i) = rho(i) + dt*(src + neg_flux_div);
    end do
  end subroutine step1

  subroutine f90_step(ncol, nstep, dt, rho, work) bind(c)
    integer(kind=c_int), intent(in) :: ncol, nstep
    real(kind=c_real), intent(in) :: dt
    real(kind=c_real), intent(inout) :: rho(ncell,ncol), work(ncell,ncol)
    integer :: c, i

    !$omp parallel do
    do c = 1, ncol
       do i = 1, nstep
          call step1(dt, rho(:,c), work(:,c))
       end do
    end do
  end subroutine f90_step
end module vec
