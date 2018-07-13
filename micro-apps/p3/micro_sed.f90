program micro_sed

  implicit none

  integer, parameter :: kts=0, kte=42, ni=0, nk=42, its=0, ite=42


  integer :: i,k,ktop,kbot,kdir,i_strt,k_strt,k_qxbot,k_qxtop,k_temp,tmpint1

  logical :: log_qxpresent

  real :: Co_max, dt, dt_left, dt_sub, dum1, dum2, dumii, dumjj, rdumii, rdumjj, &
          fluxdiv_nx, fluxdix_qx, fluxdiv_qx, inv_dum3, odt, prt_accum, tmp1, tmp2, &
          cons1, rhow, inv_rhow, sxth, pi, piov6, nsmall, qsmall

  real, dimension(150) :: mu_r_table
  real, dimension(kts:kte) :: V_qr, V_nr
  real, dimension(ni,nk)  :: qr
  real, dimension(ni,nk)  :: nr
  real, dimension(its:ite,kts:kte) :: mu_r
  real, dimension(its:ite,kts:kte) :: lamr

  ! constants
  inv_rhow = 1.e-3   ! inverse of (max.) density of liquid water
  odt      = 1./dt   ! inverse model time step
  rhow     = 997.
  sxth     = 1./6.
  pi       = 3.14159265
  piov6    = pi*sxth
  cons1    = piov6*rhow
  qsmall   = 1.e-14
  nsmall   = 1.e-16

  ! Rain sedimentation:  (adaptivive substepping)
  i_loop_main: do i = its,ite

     log_qxpresent = .false.
     k_qxtop       = kbot

     !find top, determine qxpresent
     do k = ktop,kbot,-kdir
        if (qr(i,k).ge.qsmall) then
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
           if (qr(i,k).ge.qsmall) then
              k_qxbot = k
              exit
           endif
        enddo

        substep_sedi_r: do while (dt_left.gt.1.e-4)

           Co_max = 0.
           V_qr = 0.
           V_nr = 0.

           kloop_sedi_r1: do k = k_qxtop,k_qxbot,-kdir

              qr_notsmall_r1: if (qr(i,k)>qsmall) then

                 !Compute Vq, Vn:
                 nr(i,k)  = max(nr(i,k),nsmall)
                 call get_rain_dsd2(qr(i,k),nr(i,k),mu_r(i,k),rdumii,dumii,lamr(i,k),     &
                      mu_r_table,tmp1,tmp2)
                 call find_lookupTable_indices_3(dumii,dumjj,dum1,rdumii,rdumjj,inv_dum3, &
                      mu_r(i,k),lamr(i,k))
                 !mass-weighted fall speed:
                 dum1 = vm_table(dumii,dumjj)+(rdumii-real(dumii))*inv_dum3*              &
                      (vm_table(dumii+1,dumjj)-vm_table(dumii,dumjj))         !at mu_r
                 dum2 = vm_table(dumii,dumjj+1)+(rdumii-real(dumii))*inv_dum3*            &
                      (vm_table(dumii+1,dumjj+1)-vm_table(dumii,dumjj+1))   !at mu_r+1
                 V_qr(k) = dum1 + (rdumjj-real(dumjj))*(dum2-dum1)         !interpolated
                 V_qr(k) = V_qr(k)*rhofacr(i,k)               !corrected for air density

                 ! number-weighted fall speed:
                 dum1 = vn_table(dumii,dumjj)+(rdumii-real(dumii))*inv_dum3*              &
                      (vn_table(dumii+1,dumjj)-vn_table(dumii,dumjj)         ) !at mu_r

                 dum2 = vn_table(dumii,dumjj+1)+(rdumii-real(dumii))*inv_dum3*            &
                      (vn_table(dumii+1,dumjj+1)-vn_table(dumii,dumjj+1))    !at mu_r+1
                 V_nr(k) = dum1+(rdumjj-real(dumjj))*(dum2-dum1)            !interpolated
                 V_nr(k) = V_nr(k)*rhofacr(i,k)                !corrected for air density

              endif qr_notsmall_r1

              Co_max = max(Co_max, V_qr(k)*dt_left*inv_dzq(i,k))
              !            Co_max = max(Co_max, max(V_nr(k),V_qr(k))*dt_left*inv_dzq(i,k))

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
           do k = k_temp,k_qxtop,kdir
              flux_qx(k) = V_qr(k)*qr(i,k)*rho(i,k)
              flux_nx(k) = V_nr(k)*nr(i,k)*rho(i,k)
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
           nr(i,k) = nr(i,k) + fluxdiv_nx*dt_sub*inv_rho(i,k)

           do k = k_qxtop-kdir,k_temp,-kdir
              !-- compute flux divergence
              fluxdiv_qx = (flux_qx(k+kdir) - flux_qx(k))*inv_dzq(i,k)
              fluxdiv_nx = (flux_nx(k+kdir) - flux_nx(k))*inv_dzq(i,k)
              !-- update prognostic variables
              qr(i,k) = qr(i,k) + fluxdiv_qx*dt_sub*inv_rho(i,k)
              nr(i,k) = nr(i,k) + fluxdiv_nx*dt_sub*inv_rho(i,k)
           enddo

           dt_left = dt_left - dt_sub  !update time remaining for sedimentation
           if (k_qxbot.ne.kbot) k_qxbot = k_qxbot - kdir
           !or, optimzed: k_qxbot = k_qxbot +(k_qxbot.eq.kbot)*kdir

        enddo substep_sedi_r

        prt_liq(i) = prt_liq(i) + prt_accum*inv_rhow*odt

     endif qr_present

  enddo i_loop_main

end program micro_sed

real function gamma(X)
  !----------------------------------------------------------------------
  ! THIS ROUTINE CALCULATES THE gamma FUNCTION FOR A REAL ARGUMENT X.
  !   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
  !   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE gamma
  !   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
  !   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
  !   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
  !   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
  !   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
  !   MACHINE-DEPENDENT CONSTANTS.
  !----------------------------------------------------------------------
  !
  ! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
  !
  ! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
  ! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
  ! XBIG   - THE LARGEST ARGUMENT FOR WHICH gamma(X) IS REPRESENTABLE
  !          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
  !                  gamma(XBIG) = BETA**MAXEXP
  ! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
  !          APPROXIMATELY BETA**MAXEXP
  ! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
  !          1.0+EPS .GT. 1.0
  ! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
  !          1/XMININ IS MACHINE REPRESENTABLE
  !
  !     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
  !
  !                            BETA       MAXEXP        XBIG
  !
  ! CRAY-1         (S.P.)        2         8191        966.961
  ! CYBER 180/855
  !   UNDER NOS    (S.P.)        2         1070        177.803
  ! IEEE (IBM/XT,
  !   SUN, ETC.)   (S.P.)        2          128        35.040
  ! IEEE (IBM/XT,
  !   SUN, ETC.)   (D.P.)        2         1024        171.624
  ! IBM 3033       (D.P.)       16           63        57.574
  ! VAX D-FORMAT   (D.P.)        2          127        34.844
  ! VAX G-FORMAT   (D.P.)        2         1023        171.489
  !
  !                            XINF         EPS        XMININ
  !
  ! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
  ! CYBER 180/855
  !   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
  ! IEEE (IBM/XT,
  !   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
  ! IEEE (IBM/XT,
  !   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
  ! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
  ! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
  ! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
  !
  !----------------------------------------------------------------------
  !
  ! ERROR RETURNS
  !
  !  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
  !     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
  !     TO BE FREE OF UNDERFLOW AND OVERFLOW.
  !
  !
  !  INTRINSIC FUNCTIONS REQUIRED ARE:
  !
  !     INT, DBLE, EXP, log, REAL, SIN
  !
  !
  ! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
  !              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
  !              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
  !              (ED.), SPRINGER VERLAG, BERLIN, 1976.
  !
  !              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
  !              SONS, NEW YORK, 1968.
  !
  !  LATEST MODIFICATION: OCTOBER 12, 1989
  !
  !  AUTHORS: W. J. CODY AND L. STOLTZ
  !           APPLIED MATHEMATICS DIVISION
  !           ARGONNE NATIONAL LABORATORY
  !           ARGONNE, IL 60439
  !
  !----------------------------------------------------------------------
  implicit none
  integer :: I,N
  logical :: l_parity
  real ::                                                       &
       CONV,EPS,FACT,HALF,ONE,res,sum,TWELVE,                    &
       TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
  real, dimension(7) :: C
  real, dimension(8) :: P
  real, dimension(8) :: Q
  real, parameter    :: constant1 = 0.9189385332046727417803297

  !----------------------------------------------------------------------
  !  MATHEMATICAL CONSTANTS
  !----------------------------------------------------------------------
  data ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/
  !----------------------------------------------------------------------
  !  MACHINE DEPENDENT PARAMETERS
  !----------------------------------------------------------------------
  data XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,XINF/3.4E38/
  !----------------------------------------------------------------------
  !  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
  !     APPROXIMATION OVER (1,2).
  !----------------------------------------------------------------------
  data P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,  &
       -3.79804256470945635097577E+2,6.29331155312818442661052E+2,  &
       8.66966202790413211295064E+2,-3.14512729688483675254357E+4,  &
       -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
  data Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,  &
       -1.01515636749021914166146E+3,-3.10777167157231109440444E+3, &
       2.25381184209801510330112E+4,4.75584627752788110767815E+3,  &
       -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
  !----------------------------------------------------------------------
  !  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
  !----------------------------------------------------------------------
  data C/-1.910444077728E-03,8.4171387781295E-04,                      &
       -5.952379913043012E-04,7.93650793500350248E-04,                 &
       -2.777777777777681622553E-03,8.333333333333333331554247E-02,    &
       5.7083835261E-03/
  !----------------------------------------------------------------------
  !  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
  !----------------------------------------------------------------------
  CONV(I) = REAL(I)
  l_parity=.FALSE.
  FACT=ONE
  N=0
  Y=X
  if (Y.LE.ZERO) then
     !----------------------------------------------------------------------
     !  ARGUMENT IS NEGATIVE
     !----------------------------------------------------------------------
     Y=-X
     Y1=AINT(Y)
     res=Y-Y1
     if (res.NE.ZERO) then
        if(Y1.NE.AINT(Y1*HALF)*TWO)l_parity=.TRUE.
        FACT=-PI/SIN(PI*res)
        Y=Y+ONE
     else
        res=XINF
        goto 900
     endif
  endif
  !----------------------------------------------------------------------
  !  ARGUMENT IS POSITIVE
  !----------------------------------------------------------------------
  if (Y.LT.EPS) then
     !----------------------------------------------------------------------
     !  ARGUMENT .LT. EPS
     !----------------------------------------------------------------------
     if (Y.GE.XMININ) then
        res=ONE/Y
     else
        res=XINF
        goto 900
     endif
  elseif (Y.LT.TWELVE) then
     Y1=Y
     if (Y.LT.ONE) then
        !----------------------------------------------------------------------
        !  0.0 .LT. ARGUMENT .LT. 1.0
        !----------------------------------------------------------------------
        Z=Y
        Y=Y+ONE
     else
        !----------------------------------------------------------------------
        !  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
        !----------------------------------------------------------------------
        N=INT(Y)-1
        Y=Y-CONV(N)
        Z=Y-ONE
     endif
     !----------------------------------------------------------------------
     !  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
     !----------------------------------------------------------------------
     XNUM=ZERO
     XDEN=ONE
     do I=1,8
        XNUM=(XNUM+P(I))*Z
        XDEN=XDEN*Z+Q(I)
     enddo
     res=XNUM/XDEN+ONE
     if (Y1.LT.Y) then
        !----------------------------------------------------------------------
        !  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
        !----------------------------------------------------------------------
        res=res/Y1
     elseif (Y1.GT.Y) then
        !----------------------------------------------------------------------
        !  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
        !----------------------------------------------------------------------
        do I=1,N
           res=res*Y
           Y=Y+ONE
        enddo
     endif
  else
     !----------------------------------------------------------------------
     !  EVALUATE FOR ARGUMENT .GE. 12.0,
     !----------------------------------------------------------------------
     if (Y.LE.XBIG) then
        YSQ=Y*Y
        sum=C(7)
        do I=1,6
           sum=sum/YSQ+C(I)
        enddo
        sum=sum/Y-Y+constant1
        sum=sum+(Y-HALF)*log(Y)
        res=exp(sum)
     else
        res=XINF
        goto 900
     endif
  endif
  !----------------------------------------------------------------------
  !  FINAL ADJUSTMENTS AND RETURN
  !----------------------------------------------------------------------
  if (l_parity)res=-res
  if (FACT.NE.ONE)res=FACT/res
900 gamma=res
  return
  ! ---------- LAST LINE OF gamma ----------

end function gamma

!======================================================================================!
subroutine find_lookupTable_indices_3(dumii,dumjj,dum1,rdumii,rdumjj,inv_dum3,mu_r,lamr)

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
     inv_dum3  = thrd*0.1            !i.e. 1/30
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

!===========================================================================================
subroutine get_rain_dsd2(qr,nr,mu_r,rdumii,dumii,lamr,mu_r_table,cdistr,logn0r)

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
  real ::           cons1, rhow, inv_rhow, sxth, pi, piov6, nsmall, qsmall, thrd

    ! constants
  inv_rhow = 1.e-3   ! inverse of (max.) density of liquid water
  rhow     = 997.
  sxth     = 1./6.
  thrd     = 1./3.
  pi       = 3.14159265
  piov6    = pi*sxth
  cons1    = piov6*rhow
  nsmall   = 1.e-16
  qsmall   = 1.e-14

  !--------------------------------------------------------------------------

  if (qr.ge.qsmall) then

     ! use lookup table to get mu
     ! mu-lambda relationship is from Cao et al. (2008), eq. (7)

     ! find spot in lookup table
     ! (scaled N/q for lookup table parameter space_
     nr      = max(nr,nsmall)
     inv_dum = (qr/(cons1*nr*6.))**thrd

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

     lamr   = (cons1*nr*(mu_r+3.)*(mu_r+2)*(mu_r+1.)/(qr))**thrd  ! recalculate slope based on mu_r
     lammax = (mu_r+1.)*1.e+5   ! check for slope
     lammin = (mu_r+1.)*1250.   ! set to small value since breakup is explicitly included (mean size 0.8 mm)

     ! apply lambda limiters for rain
     if (lamr.lt.lammin) then
        lamr = lammin
        nr   = exp(3.*log(lamr)+log(qr)+log(gamma(mu_r+1.))-log(gamma(mu_r+4.)))/(cons1)
     elseif (lamr.gt.lammax) then
        lamr = lammax
        nr   = exp(3.*log(lamr)+log(qr)+log(gamma(mu_r+1.))-log(gamma(mu_r+4.)))/(cons1)
     endif

     cdistr  = nr/gamma(mu_r+1.)
     logn0r  = alog10(nr)+(mu_r+1.)*alog10(lamr)-alog10(gamma(mu_r+1)) !note: logn0r is calculated as log10(n0r)

  else

     lamr   = 0.
     cdistr = 0.
     logn0r = 0.

  endif

end subroutine get_rain_dsd2
