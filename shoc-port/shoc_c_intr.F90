module shoc_c_intr

  use iso_c_binding

contains

  subroutine hi(shcol, nlev, nlevi, nqtracers, dtime, arr) bind(c)
    integer, value, intent(in) :: shcol, nlev, nlevi, nqtracers
    real(kind=c_double), value, intent(in) :: dtime
    real(kind=c_double), dimension(shcol,nlev), intent(inout) :: arr
    integer :: i, j

    print *, 'hi', shcol, nlev, nlevi, nqtracers, dtime
    do j = 1, nlev
       do i = 1, shcol
          arr(i,j) = 10*i + j
       end do
    end do
  end subroutine hi

  subroutine shoc_c_init(nlev, gravit, rair, rh2o, cpair, zvir, latvap, latice, &
       karman, pref_mid) bind(c)
    use shoc, only: shoc_init
    implicit none

    integer(kind=c_int), value, intent(in) :: nlev
    real(kind=c_double), value, intent(in) :: gravit, rair, rh2o, cpair, zvir, latvap, &
         latice, karman
    real(kind=c_double), intent(in) :: pref_mid(nlev)

    call shoc_init(nlev, gravit, rair, rh2o, cpair, zvir, latvap, latice, karman, &
         pref_mid, nlev, 1)
  end subroutine shoc_c_init

  subroutine shoc_c_main( &
     shcol, nlev, nlevi, dtime, nadv, &   ! Input
     host_dx, host_dy,thv, &              ! Input
     zt_grid,zi_grid,pres,presi,pdel,&    ! Input
     wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, & ! Input
     wtracer_sfc,num_qtracers,w_field, &  ! Input
     exner, phis, &                       ! Input
     host_temp, tke, thetal, qw, &        ! Input/Output
     u_wind, v_wind,qtracers,&            ! Input/Output
     wthv_sec,tkh,tk,shoc_ql, &           ! Input/Output
     shoc_cldfrac,pblh,&                  ! Output
     shoc_mix, isotropy,&                 ! Output (diagnostic)
     w_sec, thl_sec, qw_sec, qwthl_sec,&  ! Output (diagnostic)
     wthl_sec, wqw_sec, wtke_sec,&        ! Output (diagnostic)
     uw_sec, vw_sec, w3,&                 ! Output (diagnostic)    
     wqls_sec, brunt &                    ! Output (diagnostic)
     ) bind(c)

    use shoc, only: shoc_main
    implicit none

    integer, value, intent(in) :: shcol, nlev, nlevi, num_qtracers, nadv
    real(kind=c_double), value, intent(in) :: dtime
    real(kind=c_double), dimension(shcol), intent(in) :: &
         host_dx, host_dy, wthl_sfc, wqw_sfc, phis
    real(kind=c_double), dimension(shcol), intent(inout) :: uw_sfc, vw_sfc
    real(kind=c_double), dimension(shcol,nlev), intent(in) :: &
         zt_grid, pres, pdel, thv, w_field, exner
    real(kind=c_double), dimension(shcol,nlevi), intent(in) :: &
         zi_grid, presi
    real(kind=c_double), dimension(shcol,num_qtracers), intent(in) :: &
         wtracer_sfc

    real(kind=c_double), dimension(shcol,nlev), intent(inout) :: &
         tke, thetal, qw, u_wind, v_wind, wthv_sec, tk, tkh, host_temp
    real(kind=c_double), dimension(shcol,nlev,num_qtracers), intent(inout) :: &
         qtracers

    real(kind=c_double), dimension(shcol,nlev), intent(out) :: &
         shoc_cldfrac, shoc_ql, shoc_mix, w_sec, wqls_sec, brunt, isotropy
    real(kind=c_double), dimension(shcol,nlevi), intent(out) :: &
         thl_sec, qw_sec, qwthl_sec, wthl_sec, wqw_sec, wtke_sec, uw_sec, vw_sec, w3
    real(kind=c_double), dimension(shcol), intent(out) :: pblh

    integer :: i, c
    real(8), parameter :: ustar2 = 0.28d0**2
    real(8) :: speed

    call shoc_main( &
         shcol, nlev, nlevi, dtime, nadv, &   ! Input
         host_dx, host_dy,thv, &              ! Input
         zt_grid, zi_grid, pres,presi,pdel, & ! Input
         wthl_sfc, wqw_sfc, uw_sfc, vw_sfc,&  ! Input
         wtracer_sfc, num_qtracers, w_field,& ! Input     
         exner,phis, &                        ! Input
         host_temp, tke, thetal, qw, &        ! Input/Output
         u_wind, v_wind, qtracers, &          ! Input/Output
         wthv_sec,tkh,tk, shoc_ql, &          ! Input/Output
         shoc_cldfrac, pblh, &                ! Output
         shoc_mix, isotropy, &                ! Output (diagnostic)
         w_sec, thl_sec, qw_sec, qwthl_sec, & ! Output (diagnostic)
         wthl_sec, wqw_sec, wtke_sec, &       ! Output (diagnostic)
         uw_sec, vw_sec, w3, &                ! Output (diagnostic)    
         wqls_sec, brunt)
  end subroutine shoc_c_main

end module shoc_c_intr
