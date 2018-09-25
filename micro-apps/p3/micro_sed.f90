#include "types.f.h"

module micro_sed_mod

  implicit none

  !
  ! Constants
  !
  real, parameter ::                   &
       INV_RHOW = 1.e-3,               &
       RHOW     = 997.,                & ! density liquid water [kg/m^3]
       THRD     = 1./3.,               &
       SXTH     = 1./6.,               &
       PI       = 3.14159265,          &
       PIOV6    = PI*SXTH,             &
       CONS1    = PIOV6*RHOW,          &
       QSMALL   = 1.e-14,              &
       NSMALL   = 1.e-16,              &
       RD       = 287.15,              & ! specific gas constant of dry air [J/kg/K]
       RHOSUR   = 100000./(RD*273.15), & ! density air at ~ sea level: (p [Pa])/(R_d (T [K])) [kg/m^3]
       CP       = 1005.,               & ! specific heat of air, const pressure, at 300 K [J/kg/K]
       INV_CP   = 1./CP

  ! Quantities:

  !   qr (q_r) is mass mixing ratio for rain [kg/kg], i.e., dimensionless.

  !   nr (N_r) is total number mixing ratio for rain [dimensionless]. It is not,
  ! e.g., the number concentration [1/m^3].

  !   Fall terminal velocities V_* have units [m/s]. the expressions for V are
  ! weighted by [(air reference density)/(air actual density)]^0.54 (see
  ! doi.org/10.1175/JAS3890.1 eq 5).

  !
  ! Globals
  !
  real, dimension(300,10), target :: VN_TABLE,VM_TABLE  ! lookup table values for rain number- and mass-weighted fallspeeds
  real, dimension(150), target    :: MU_R_TABLE         ! lookup table values for rain shape parameter mu_r

contains

  function c_get_vn_table() bind(c)
    use iso_c_binding
    type(c_ptr) :: c_get_vn_table
    c_get_vn_table = c_loc(VN_TABLE)
  end function c_get_vn_table

  function c_get_vm_table() bind(c)
    use iso_c_binding
    type(c_ptr) :: c_get_vm_table
    c_get_vm_table = c_loc(VM_TABLE)
  end function c_get_vm_table

  function c_get_mu_r_table() bind(c)
    use iso_c_binding
    type(c_ptr) :: c_get_mu_r_table
    c_get_mu_r_table = c_loc(MU_R_TABLE)
  end function c_get_mu_r_table

  subroutine append_precision(string, prefix)
    use iso_c_binding

    character(kind=c_char, len=128), intent(inout) :: string
    character(*), intent(in) :: prefix
    real(kind=c_real) :: s

    write (string, '(a,i1,a1)') prefix, sizeof(s), C_NULL_CHAR
  end subroutine append_precision

#ifdef TRACE
  !=============================================================================!
  subroutine trace_loop_impl(name, b, e)
  !=============================================================================!
    character(*), intent(in) :: name
    integer, intent(in) :: b, e


    print '(A," LOOP ",I0," -> ",I0)', name, b, e
  end subroutine trace_loop_impl

  !=============================================================================!
  subroutine trace_data_impl(name, i, k, item)
  !=============================================================================!
    character(*), intent(in) :: name
    integer, intent(in) :: i, k
    real, intent(in) :: item

    print '(A,"[",I0,"][",I0,"] = ",F20.12)', name, i, k, item
  end subroutine trace_data_impl

#define trace_loop(n, b, e) call trace_loop_impl(n, b, e)
#define trace_data(n, i, k, v) call trace_data_impl(n, i, k, v)
#else
#define trace_loop(n, b, e)
#define trace_data(n, i, k, v)
#endif

  !=============================================================================!
  subroutine p3_init() bind(c)
  !=============================================================================!
    ! Generate lookup table for rain fallspeed and ventilation parameters
    ! the lookup table is two dimensional as a function of number-weighted mean size
    ! proportional to qr/Nr and shape parameter mu_r

    use cpp_bridge
    use iso_c_binding

    integer                      :: i,ii,jj,kk
    real                         :: lamr,lamold,mu_r,dum,dm,dum1,dum2,dum3,dum4,dum5,dd,amg,vt,dia,initlamr
    logical :: ok
    character(kind=c_char, len=128) :: mu_r_filename, vn_filename, vm_filename

    ! Generate lookup table for rain shape parameter mu_r
    ! this is very fast so it can be generated at the start of each run
    ! make a 150x1 1D lookup table, this is done in parameter
    ! space of a scaled mean size proportional qr/Nr -- initlamr

    !print*, '   Generating rain lookup-table ...'

    call append_precision(mu_r_filename, c_char_"mu_r_table.dat")
    call append_precision(vn_filename, c_char_"vn_table.dat")
    call append_precision(vm_filename, c_char_"vm_table.dat")
    if ( array_io_file_exists(mu_r_filename) .and. &
         array_io_file_exists(vn_filename) .and. &
         array_io_file_exists(vm_filename)) then
       ok = array_io_read(mu_r_filename, c_loc(MU_R_TABLE), size(MU_R_TABLE)) .and. &
            array_io_read(vn_filename, c_loc(VN_TABLE), size(VN_TABLE)) .and. &
            array_io_read(vm_filename, c_loc(VM_TABLE), size(VM_TABLE))
       if (.not. ok) then
          print *, 'p3_init: One or more table files exists but gave a read error; computing from scratch.'
       else
          return
       end if
    end if

    do i = 1,150              ! loop over lookup table values
       initlamr = 1./((real(i)*2.)*1.e-6 + 250.e-6)

       ! iterate to get mu_r
       ! mu_r-lambda relationship is from Cao et al. (2008), eq. (7)

       ! start with first guess, mu_r = 0

       mu_r = 0.

       do ii=1,50
          lamr = initlamr*((mu_r+3.)*(mu_r+2.)*(mu_r+1.)/6.)**thrd

          ! new estimate for mu_r based on lambda
          ! set max lambda in formula for mu_r to 20 mm-1, so Cao et al.
          ! formula is not extrapolated beyond Cao et al. data range
          dum  = min(20.,lamr*1.e-3)
          mu_r = max(0.,-0.0201*dum**2+0.902*dum-1.718)

          ! if lambda is converged within 0.1%, then exit loop
          if (ii.ge.2) then
             if (abs((lamold-lamr)/lamr).lt.0.001) goto 111
          end if

          lamold = lamr

       enddo

111    continue

       ! assign lookup table values
       MU_R_TABLE(i) = mu_r

    enddo

    mu_r_loop: do ii = 1,10   !** change 10 to 9, since range of mu_r is 0-8  CONFIRM
       !mu_r_loop: do ii = 1,9   !** change 10 to 9, since range of mu_r is 0-8

       mu_r = real(ii-1)  ! values of mu

       ! loop over number-weighted mean size
       meansize_loop: do jj = 1,300

          if (jj.le.20) then
             dm = (real(jj)*10.-5.)*1.e-6      ! mean size [m]
          elseif (jj.gt.20) then
             dm = (real(jj-20)*30.+195.)*1.e-6 ! mean size [m]
          endif

          lamr = (mu_r+1)/dm

          ! do numerical integration over PSD

          dum1 = 0. ! numerator,   number-weighted fallspeed
          dum2 = 0. ! denominator, number-weighted fallspeed
          dum3 = 0. ! numerator,   mass-weighted fallspeed
          dum4 = 0. ! denominator, mass-weighted fallspeed
          dum5 = 0. ! term for ventilation factor in evap
          dd   = 2.

          ! loop over PSD to numerically integrate number and mass-weighted mean fallspeeds
          do kk = 1,10000

             dia = (real(kk)*dd-dd*0.5)*1.e-6  ! size bin [m]
             amg = PIOV6*997.*dia**3           ! mass [kg]
             amg = amg*1000.                   ! convert [kg] to [g]
             !get fallspeed as a function of size [m s-1]
             if (dia*1.e+6.le.134.43)      then
                vt = 4.5795e+3*amg**(2.*THRD)
             elseif (dia*1.e+6.lt.1511.64) then
                vt = 4.962e+1*amg**THRD
             elseif (dia*1.e+6.lt.3477.84) then
                vt = 1.732e+1*amg**SXTH
             else
                vt = 9.17
             endif

             !note: factor of 4.*mu_r is non-answer changing and only needed to
             !      prevent underflow/overflow errors, same with 3.*mu_r for dum5
             dum1 = dum1 + vt*10.**(mu_r*alog10(dia)+4.*mu_r)*exp(-lamr*dia)*dd*1.e-6
             dum2 = dum2 + 10.**(mu_r*alog10(dia)+4.*mu_r)*exp(-lamr*dia)*dd*1.e-6
             dum3 = dum3 + vt*10.**((mu_r+3.)*alog10(dia)+4.*mu_r)*exp(-lamr*dia)*dd*1.e-6
             dum4 = dum4 + 10.**((mu_r+3.)*alog10(dia)+4.*mu_r)*exp(-lamr*dia)*dd*1.e-6
             dum5 = dum5 + (vt*dia)**0.5*10.**((mu_r+1.)*alog10(dia)+3.*mu_r)*exp(-lamr*dia)*dd*1.e-6

          enddo ! kk-loop (over PSD)

          dum2 = max(dum2, 1.e-30)  !to prevent divide-by-zero below
          dum4 = max(dum4, 1.e-30)  !to prevent divide-by-zero below
          dum5 = max(dum5, 1.e-30)  !to prevent log10-of-zero below

          VN_TABLE(jj,ii)    = dum1/dum2
          VM_TABLE(jj,ii)    = dum3/dum4

       enddo meansize_loop

    enddo mu_r_loop

    ok = array_io_write(mu_r_filename, c_loc(MU_R_TABLE), size(MU_R_TABLE))
    ok = array_io_write(vn_filename, c_loc(VN_TABLE), size(VN_TABLE))
    ok = array_io_write(vm_filename, c_loc(VM_TABLE), size(VM_TABLE))

  end subroutine p3_init

  !=============================================================================!
  subroutine populate_input(ni, nk, qr, nr, th, dzq, pres, kdir)
  !=============================================================================!
    use cpp_bridge
    use iso_c_binding

    implicit none

    integer, intent(in) :: ni, nk, kdir

    real, dimension(1:ni,1:nk), intent(inout), target :: qr, nr, th, dzq, pres

    call fully_populate_input_data(ni, nk, kdir, c_loc(qr), c_loc(nr), c_loc(th), c_loc(dzq), c_loc(pres))

  end subroutine populate_input

  !=============================================================================!
  subroutine micro_sed_func_wrap(ni, nk, dt, ts, kdir, repeat) bind(c)
  !=============================================================================!
    use cpp_bridge
    use iso_c_binding
    use omp_lib

    implicit none

    integer, intent(in) :: ni, nk, ts, kdir, repeat
    real, intent(in) :: dt

    integer, parameter :: chunksize = CHUNKSIZE
    real, dimension(chunksize,nk) :: cqr, cnr, cth, cdzq, cpres
    real, dimension(:,:), allocatable, target :: qr, nr, th, dzq, pres
    real, dimension(chunksize) :: cprt_liq
    real, dimension(ni), target :: prt_liq
    real(8) :: start, finish
    integer :: ti, ci, nchunk, cni, cnk, ws, i
    logical :: ok

    character (kind=c_char, len=*), parameter :: filename = c_char_"fortran"//char(0)

    allocate(qr(ni,nk))
    allocate(nr(ni,nk))
    allocate(th(ni,nk))
    allocate(dzq(ni,nk))
    allocate(pres(ni,nk))

    call dump_arch_f90()
    print '("Running with ni=",I0," nk=",I0," dt=",F6.2," ts=",I0)', ni, nk, dt, ts

    call populate_input(ni, nk, qr, nr, th, dzq, pres, kdir)
    print *, 'chunksize',chunksize

    start = omp_get_wtime()

    prt_liq(:) = 0

#if CHUNKSIZE > 0
    nchunk = (ni + chunksize - 1) / chunksize
    do ti = 1, ts
!$OMP PARALLEL DO DEFAULT(SHARED)
       do ci = 1, nchunk
          call micro_sed_func_chunk(ni, nk, kdir, dt, qr, nr, th, dzq, pres, prt_liq, ci)
       end do
!$OMP END PARALLEL DO
    end do
#else
    do ti = 1, ts
!$OMP PARALLEL DO DEFAULT(SHARED)
       do i = 1, ni
          call micro_sed_func(1, nk, kdir, 1, 1, dt, &
               qr(i,:), nr(i,:), th(i,:), dzq(i,:), pres(i,:), prt_liq(i))
       end do
!$OMP END PARALLEL DO
    end do
#endif

    finish = omp_get_wtime()

    print '("Time = ",E20.3," seconds.")', finish - start

    ok = dump_all(filename, c_loc(qr), c_loc(nr), c_loc(th), c_loc(dzq), c_loc(pres), c_loc(prt_liq), ni, nk, dt, ts)

    deallocate(qr, nr, th, dzq, pres)

  end subroutine micro_sed_func_wrap

  !=============================================================================!
  subroutine micro_sed_func_c(kts, kte, kdir, its, ite, dt, qr, nr, th, dzq, pres, prt_liq) bind(c)
  !=============================================================================!
    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: kts, kte, kdir, its, ite
    real(kind=c_real), value, intent(in) :: dt
    real(kind=c_real), dimension(its:ite,kts:kte), intent(inout) :: qr, nr
    real(kind=c_real), intent(in), dimension(its:ite,kts:kte) :: th, dzq, pres
    real(kind=c_real), dimension(its:ite), intent(out) :: prt_liq

    call micro_sed_func(kts, kte, kdir, its, ite, dt, qr, nr, th, dzq, pres, prt_liq)
  end subroutine micro_sed_func_c

  !=============================================================================!
  subroutine micro_sed_func_chunk(ni, nk, kdir, dt, qr, nr, th, dzq, pres, prt_liq, ci)
  !=============================================================================!
    implicit none

    integer, intent(in) :: ni, nk, kdir, ci
    real, intent(in) :: dt

    real, dimension(ni,nk), intent(inout) :: qr, nr, th, dzq, pres

    real, dimension(ni), intent(inout) :: prt_liq

    integer, parameter :: chunksize = CHUNKSIZE
    integer :: cni, cnk, ws
    real, dimension(chunksize,nk) :: cqr, cnr, cth, cdzq, cpres
    real, dimension(chunksize) :: cprt_liq

    cni = chunksize*(ci - 1) + 1
    cnk = min(ni, cni + chunksize - 1)
    ws = cnk - cni + 1
    cqr(1:ws,:) = qr(cni:cnk,:)
    cnr(1:ws,:) = nr(cni:cnk,:)
    cth(1:ws,:) = th(cni:cnk,:)
    cdzq(1:ws,:) = dzq(cni:cnk,:)
    cpres(1:ws,:) = pres(cni:cnk,:)
    cprt_liq(1:ws) = prt_liq(1:ws)

    call micro_sed_func(1, nk, kdir, 1, ws, dt, &
         cqr(1:ws,:), cnr(1:ws,:), cth(1:ws,:), cdzq(1:ws,:), cpres(1:ws,:), cprt_liq(1:ws))

    qr(cni:cnk,:) = cqr(1:ws,:)
    nr(cni:cnk,:) = cnr(1:ws,:)
    th(cni:cnk,:) = cth(1:ws,:)
    dzq(cni:cnk,:) = cdzq(1:ws,:)
    pres(cni:cnk,:) = cpres(1:ws,:)
    prt_liq(cni:cnk) = cprt_liq(1:ws)

  end subroutine micro_sed_func_chunk

  !=============================================================================!
  subroutine micro_sed_func(kts, kte, kdir, its, ite, dt, qr, nr, th, dzq, pres, prt_liq)
  !=============================================================================!
    implicit none

    !
    ! Arg explanation
    !
    ! kts: vertical array bound
    ! kte: vertical array bound; kte > kts
    ! kdir: 1 if k=1 is the surface, -1 if k=1 is the top
    ! its: horizontal array bound
    ! ite: horizontal array bound
    ! dt: time step
    ! qr: rain, mass mixing ratio  (in/out)
    ! nr: rain, number mixing ratio (in/out)
    ! th: potential temperature                    K
    ! dzq: vertical grid spacing                   m
    ! pres: pressure                               Pa
    ! prt_liq: precipitation rate, total liquid    m s-1  (output)

    integer, intent(in) :: kts, kte, kdir, its, ite
    real, intent(in) :: dt

    real, dimension(its:ite,kts:kte), intent(inout) :: qr, nr

    real, intent(in),    dimension(its:ite,kts:kte) :: th, dzq, pres

    real, dimension(its:ite), intent(out) :: prt_liq

    integer :: i,k,ktop,kbot,k_qxbot,k_qxtop,k_temp,tmpint1,dumii,dumjj

    logical :: log_qxpresent

    real :: Co_max, dt_left, dt_sub, dum1, dum2, rdumii, rdumjj, &
            fluxdiv_nx, fluxdiv_qx, inv_dum3, odt, prt_accum, tmp1, tmp2

    real, dimension(kts:kte) :: V_qr, V_nr, flux_qx, flux_nx
    real, dimension(its:ite,kts:kte) :: mu_r, lamr, rhofacr, inv_dzq, rho, inv_rho, t, tmparr1

    inv_dzq    = 1./dzq  ! inverse of thickness of layers

    dumii = 0

    ! constants
    odt      = 1./dt   ! inverse model time step

    tmparr1 = (pres*1.e-5)**(RD*INV_CP)
    t       = th    *tmparr1    !compute temperature from theta (value at beginning of microphysics step)

    ! direction of vertical leveling:
    if (kdir == -1) then
       ktop = kts        !k of top level
       kbot = kte        !k of bottom level
    else
       ktop = kte        !k of top level
       kbot = kts        !k of bottom level
    endif

    ! Rain sedimentation:  (adaptive substepping)
    trace_loop("i_loop_main", its, ite)
    i_loop_main: do i = its,ite

       trace_loop("  k_loop_1", kbot, ktop)
       k_loop_1: do k = kbot,ktop,kdir
          rho(i,k)     = pres(i,k)/(RD*t(i,k))
          inv_rho(i,k) = 1./rho(i,k)
          rhofacr(i,k) = (RHOSUR*inv_rho(i,k))**0.54
          trace_data("    rhofacr", i, k, rhofacr(i,k))
       end do k_loop_1

       ! Note, we are skipping supersaturation checks

       log_qxpresent = .false.
       k_qxtop       = kbot

       !find top, determine qxpresent
       do k = ktop,kbot,-kdir
          if (qr(i,k).ge.QSMALL) then
             log_qxpresent = .true.
             k_qxtop = k
             exit
          endif !
       enddo

       qr_present: if (log_qxpresent) then

          dt_left   = dt  !time remaining for sedi over full model (mp) time step
          prt_accum = 0.  !precip rate for individual category

          !find bottom
          do k = kbot,k_qxtop,kdir
             if (qr(i,k).ge.QSMALL) then
                k_qxbot = k
                exit
             endif
          enddo

          substep_sedi_r: do while (dt_left.gt.1.e-4)

             Co_max = 0.
             V_qr = 0.
             V_nr = 0.

             trace_loop("  k_loop_sedi_r1", k_qxtop, k_qxbot)
             kloop_sedi_r1: do k = k_qxtop,k_qxbot,-kdir

                qr_notsmall_r1: if (qr(i,k)>QSMALL) then

                   !Compute Vq, Vn:
                   nr(i,k)  = max(nr(i,k),NSMALL)
                   trace_data("    nr", i, k, nr(i, k))
                   call get_rain_dsd2(qr(i,k),nr(i,k),mu_r(i,k),rdumii,dumii,lamr(i,k),     &
                        mu_r_table,tmp1,tmp2)
                   call find_lookupTable_indices_3(dumii,dumjj,dum1,rdumii,rdumjj,inv_dum3, &
                        mu_r(i,k),lamr(i,k))
                   !mass-weighted fall speed:
                   dum1 = VM_TABLE(dumii,dumjj)+(rdumii-real(dumii))*inv_dum3*              &
                        (VM_TABLE(dumii+1,dumjj)-VM_TABLE(dumii,dumjj))         !at mu_r
                   dum2 = VM_TABLE(dumii,dumjj+1)+(rdumii-real(dumii))*inv_dum3*            &
                        (VM_TABLE(dumii+1,dumjj+1)-VM_TABLE(dumii,dumjj+1))   !at mu_r+1
                   V_qr(k) = dum1 + (rdumjj-real(dumjj))*(dum2-dum1)         !interpolated
                   V_qr(k) = V_qr(k)*rhofacr(i,k)               !corrected for air density
                   trace_data("    V_qr", 1, k, V_qr(k))

                   ! number-weighted fall speed:
                   dum1 = VN_TABLE(dumii,dumjj)+(rdumii-real(dumii))*inv_dum3*              &
                        (VN_TABLE(dumii+1,dumjj)-VN_TABLE(dumii,dumjj)         ) !at mu_r

                   dum2 = VN_TABLE(dumii,dumjj+1)+(rdumii-real(dumii))*inv_dum3*            &
                        (VN_TABLE(dumii+1,dumjj+1)-VN_TABLE(dumii,dumjj+1))    !at mu_r+1
                   V_nr(k) = dum1+(rdumjj-real(dumjj))*(dum2-dum1)            !interpolated
                   V_nr(k) = V_nr(k)*rhofacr(i,k)                !corrected for air density
                   trace_data("    V_nr", 1, k, V_nr(k))

                endif qr_notsmall_r1

                Co_max = max(Co_max, V_qr(k)*dt_left*inv_dzq(i,k))
                !            Co_max = max(Co_max, max(V_nr(k),V_qr(k))*dt_left*inv_dzq(i,k))

                trace_data("  Co_max", 1, 1, Co_max)

             enddo kloop_sedi_r1

             !-- compute dt_sub
             tmpint1 = int(Co_max+1.)    !number of substeps remaining if dt_sub were constant
             dt_sub  = min(dt_left, dt_left/float(tmpint1))

             if (k_qxbot.eq.kbot) then
                k_temp = k_qxbot
             else
                k_temp = k_qxbot-kdir
             endif

             !-- calculate fluxes
             trace_loop("  k_flux_loop", k_temp, k_qxtop)
             do k = k_temp,k_qxtop,kdir
                flux_qx(k) = V_qr(k)*qr(i,k)*rho(i,k)
                trace_data("    flux_qx", 1, k, flux_qx(k))
                flux_nx(k) = V_nr(k)*nr(i,k)*rho(i,k)
                trace_data("    flux_nx", 1, k, flux_nx(k))
             enddo

             !accumulated precip during time step
             if (k_qxbot.eq.kbot) prt_accum = prt_accum + flux_qx(kbot)*dt_sub
             !or, optimized: prt_accum = prt_accum - (k_qxbot.eq.kbot)*dt_sub

             !--- for top level only (since flux is 0 above)
             k = k_qxtop
             !- compute flux divergence
             fluxdiv_qx = -flux_qx(k)*inv_dzq(i,k)
             fluxdiv_nx = -flux_nx(k)*inv_dzq(i,k)
             !- update prognostic variables
             qr(i,k) = qr(i,k) + fluxdiv_qx*dt_sub*inv_rho(i,k)
             trace_data("  qr", i, k, qr(i,k))
             nr(i,k) = nr(i,k) + fluxdiv_nx*dt_sub*inv_rho(i,k)
             trace_data("  nr", i, k, nr(i,k))

             trace_loop("  k_flux_div_loop", k_qxtop-kdir, k_temp)
             do k = k_qxtop-kdir,k_temp,-kdir
                !-- compute flux divergence
                fluxdiv_qx = (flux_qx(k+kdir) - flux_qx(k))*inv_dzq(i,k)
                fluxdiv_nx = (flux_nx(k+kdir) - flux_nx(k))*inv_dzq(i,k)
                !-- update prognostic variables
                qr(i,k) = qr(i,k) + fluxdiv_qx*dt_sub*inv_rho(i,k)
                trace_data("    qr", i, k, qr(i,k))
                nr(i,k) = nr(i,k) + fluxdiv_nx*dt_sub*inv_rho(i,k)
                trace_data("    nr", i, k, nr(i,k))
             enddo

             dt_left = dt_left - dt_sub  !update time remaining for sedimentation
             if (k_qxbot.ne.kbot) k_qxbot = k_qxbot - kdir
             !or, optimzed: k_qxbot = k_qxbot +(k_qxbot.eq.kbot)*kdir

          enddo substep_sedi_r

          trace_data("  prt_liq", i, 1, prt_liq(i))
          trace_data("  prt_accum", 1, 1, prt_accum)
          prt_liq(i) = prt_liq(i) + prt_accum*INV_RHOW*odt
          trace_data("  prt_liq", i, 1, prt_liq(i))

       endif qr_present

    enddo i_loop_main

  end subroutine micro_sed_func

  !=============================================================================!
  subroutine find_lookupTable_indices_3(dumii,dumjj,dum1,rdumii,rdumjj,inv_dum3,mu_r,lamr)
  !=============================================================================!

    !------------------------------------------------------------------------------------------!
    ! Finds indices in rain lookup table (3)
    !------------------------------------------------------------------------------------------!

    implicit none

    ! arguments:
    integer, intent(out) :: dumii,dumjj
    real,    intent(out) :: dum1,rdumii,rdumjj,inv_dum3
    real,    intent(in)  :: mu_r,lamr

    !------------------------------------------------------------------------------------------!

    ! find location in scaled mean size space
    dum1 = (mu_r+1.)/lamr
    if (dum1.le.195.e-6) then
       inv_dum3  = 0.1
       rdumii = (dum1*1.e6+5.)*inv_dum3
       rdumii = max(rdumii, 1.)
       rdumii = min(rdumii,20.)
       dumii  = int(rdumii)
       dumii  = max(dumii, 1)
       dumii  = min(dumii,20)
    elseif (dum1.gt.195.e-6) then
       inv_dum3  = THRD*0.1            !i.e. 1/30
       rdumii = (dum1*1.e+6-195.)*inv_dum3 + 20.
       rdumii = max(rdumii, 20.)
       rdumii = min(rdumii,300.)
       dumii  = int(rdumii)
       dumii  = max(dumii, 20)
       dumii  = min(dumii,299)
    endif

    ! find location in mu_r space
    rdumjj = mu_r+1.
    rdumjj = max(rdumjj,1.)
    rdumjj = min(rdumjj,10.)
    dumjj  = int(rdumjj)
    dumjj  = max(dumjj,1)
    dumjj  = min(dumjj,9)

  end subroutine find_lookupTable_indices_3

  !=============================================================================!
  subroutine get_rain_dsd2(qr,nr,mu_r,rdumii,dumii,lamr,mu_r_table,cdistr,logn0r)
  !=============================================================================!

    ! Computes and returns rain size distribution parameters

    implicit none

    !arguments:
    real, dimension(:), intent(in)  :: mu_r_table
    real,     intent(in)            :: qr
    real,     intent(inout)         :: nr
    real,     intent(out)           :: rdumii,lamr,mu_r,cdistr,logn0r
    integer,  intent(out)           :: dumii

    !local variables:
    real                            :: inv_dum,lammax,lammin

    !--------------------------------------------------------------------------

    if (qr.ge.QSMALL) then

       ! use lookup table to get mu
       ! mu-lambda relationship is from Cao et al. (2008), eq. (7)

       ! find spot in lookup table
       ! (scaled N/q for lookup table parameter space_
       nr      = max(nr,NSMALL)
       inv_dum = (qr/(CONS1*nr*6.))**THRD

       if (inv_dum.lt.282.e-6) then
          mu_r = 8.282
       elseif (inv_dum.ge.282.e-6 .and. inv_dum.lt.502.e-6) then
          ! interpolate
          rdumii = (inv_dum-250.e-6)*1.e+6*0.5
          rdumii = max(rdumii,1.)
          rdumii = min(rdumii,150.)
          dumii  = int(rdumii)
          dumii  = min(149,dumii)
          mu_r   = mu_r_table(dumii)+(mu_r_table(dumii+1)-mu_r_table(dumii))*(rdumii-  &
               real(dumii))
       elseif (inv_dum.ge.502.e-6) then
          mu_r = 0.
       endif

       lamr   = (CONS1*nr*(mu_r+3.)*(mu_r+2)*(mu_r+1.)/(qr))**THRD  ! recalculate slope based on mu_r
       lammax = (mu_r+1.)*1.e+5   ! check for slope
       lammin = (mu_r+1.)*1250.   ! set to small value since breakup is explicitly included (mean size 0.8 mm)

       ! apply lambda limiters for rain
       if (lamr.lt.lammin) then
          lamr = lammin
          nr   = exp(3.*log(lamr)+log(qr)+log(gamma(mu_r+1.))-log(gamma(mu_r+4.)))/(CONS1)
       elseif (lamr.gt.lammax) then
          lamr = lammax
          nr   = exp(3.*log(lamr)+log(qr)+log(gamma(mu_r+1.))-log(gamma(mu_r+4.)))/(CONS1)
       endif

       cdistr  = 0. !nr/gamma(mu_r+1.)
       logn0r  = 0. !alog10(nr)+(mu_r+1.)*alog10(lamr)-alog10(gamma(mu_r+1)) !note: logn0r is calculated as log10(n0r)

    else

       lamr   = 0.
       cdistr = 0.
       logn0r = 0.

    endif

  end subroutine get_rain_dsd2

end module micro_sed_mod
